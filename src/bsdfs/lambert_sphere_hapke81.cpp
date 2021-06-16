/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

#define Power powf
#define Dot dot
#define Abs fabsf
#define Sqrt sqrtf
#define Pi M_PI

MTS_NAMESPACE_BEGIN

// This implements Hapke's 1981 approximation to the Lambert-sphere BRDF for a half space comprised of spherical Lambertian particles [d'Eon 2021 - EGSR].

// Albedo mapping from diffuse color kd to single-scattering albedo c [d'Eon 2021, Eq.(46)]
float c_from_kd( const float kd )
{
   return (1 - Power(1 - kd, 2.73556))/(1 - 0.184096*Power(1 - kd, 2.48423));
}

Spectrum reflectanceMap( const Spectrum kd )
{
    return Color3(
        c_from_kd( kd[0] ),
        c_from_kd( kd[1] ),
        c_from_kd( kd[2] )
    );
}

float lambertSpherePhaseFunction( const float u )
{
    return (2*(sqrtf(1 - u*u) - u*acosf(u)))/(3.*M_PI * M_PI);
}

Spectrum non_negative_color3( const Spectrum& c )
{
    return Color3( std::max( 0.0f, c[0] ), std::max( 0.0f, c[1] ), std::max( 0.0f, c[2] ) );
}

float safesqrt( const float x )
{
    if( x > 0.0 )
    {
        return sqrt(x);
    }
    return 0.0;
}

float safeacos( const float x )
{
    if( x > -1.0 && x < 1.0  )
    {
        return acos(x);
    }
    return 0.0;
}

// H-function for the two-stream approximation in 3D with iso scatter [Hapke 1981]
float H(const float u, const float c)
{
    return (1 + 2 * u) / (1 + 2 * sqrt(1-c) * u);
}

// [Hapke 1981 approximation]
float brdfGray( const float c, const float ui, const float uo, const float graySS, const float kd, const float iodot )
{
    const float SS = c * graySS;
    return SS + c * (H(ui, c) * H(uo, c) - 1) / (ui + uo) / (4 * Pi);
}

Spectrum lambertSphereBRDF( Spectrum c, const float ui, const float uo, const Vector3& wi, const Vector3& wo, Spectrum kd )
{
    const float iodot = dot(wi,wo);
    const float graySS = lambertSpherePhaseFunction(-iodot) / ( ui + uo );
    return Color3( brdfGray(c[0], ui, uo, graySS, kd[0], iodot ), 
                   brdfGray(c[1], ui, uo, graySS, kd[1], iodot ), 
                   brdfGray(c[2], ui, uo, graySS, kd[2], iodot ) );
}

class LambertSphereHapke81 : public BSDF {
public:
    LambertSphereHapke81(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'reflectance' and 'diffuseReflectance' as parameter names */
        m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("reflectance") ? "reflectance"
                : "diffuseReflectance", Spectrum(.5f)));
    }

    LambertSphereHapke81(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_reflectance = static_cast<Texture *>(manager->getInstance(stream));

        configure();
    }

    void configure() {
        /* Verify the input parameter and fix them if necessary */
        m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

        m_components.clear();
        if (m_reflectance->getMaximum().max() > 0)
            m_components.push_back(EDiffuseReflection | EFrontSide
                | (m_reflectance->isConstant() ? 0 : ESpatiallyVarying));
            m_usesRayDifferentials = m_reflectance->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return m_reflectance->eval(its);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);

        const Float ui = Frame::cosTheta(bRec.wi);
        const Float uo = Frame::cosTheta(bRec.wo);

        Spectrum kd = m_reflectance->eval(bRec.its);
        Spectrum c = reflectanceMap( kd );

        return (lambertSphereBRDF(c,ui,uo,bRec.wi,bRec.wo,kd) * Frame::cosTheta(bRec.wo));
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        return warp::squareToCosineHemispherePdf(bRec.wo);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;

        const Float ui = Frame::cosTheta(bRec.wi);
        const Float uo = Frame::cosTheta(bRec.wo);
        Spectrum kd = m_reflectance->eval(bRec.its);
        Spectrum c = reflectanceMap( kd );

        return lambertSphereBRDF(c,ui,uo,bRec.wi,bRec.wo,kd) / Spectrum( INV_PI );
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        pdf = warp::squareToCosineHemispherePdf(bRec.wo);

        const Float ui = Frame::cosTheta(bRec.wi);
        const Float uo = Frame::cosTheta(bRec.wo);
        Spectrum kd = m_reflectance->eval(bRec.its);
        Spectrum c = reflectanceMap( kd );

        return lambertSphereBRDF(c,ui,uo,bRec.wi,bRec.wo,kd) / Spectrum( INV_PI );
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
                && (name == "reflectance" || name == "diffuseReflectance")) {
            m_reflectance = static_cast<Texture *>(child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_reflectance.get());
    }

    Float getRoughness(const Intersection &its, int component) const {
        return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "LambertSphereHapke81[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  reflectance = " << indent(m_reflectance->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
};

// ================ Hardware shader implementation ================

class LambertSphereHapke81Shader : public Shader {
public:
    LambertSphereHapke81Shader(Renderer *renderer, const Texture *reflectance)
        : Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
        m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
    }

    bool isComplete() const {
        return m_reflectanceShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_reflectance.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_reflectanceShader.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << evalName << "(uv, wi, wo);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_reflectance;
    ref<Shader> m_reflectanceShader;
};

Shader *LambertSphereHapke81::createShader(Renderer *renderer) const {
    return new LambertSphereHapke81Shader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(LambertSphereHapke81Shader, false, Shader)
MTS_IMPLEMENT_CLASS_S(LambertSphereHapke81, false, BSDF)
MTS_EXPORT_PLUGIN(LambertSphereHapke81, "Lambert-sphere BRDF (Hapke approximation)")
MTS_NAMESPACE_END
