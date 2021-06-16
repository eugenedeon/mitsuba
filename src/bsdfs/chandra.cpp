/*  

    THIS VERSION: takes single-scattering albedo c - but has importance sampling

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

mitsuba::Spectrum powf( const mitsuba::Spectrum a, const int b )
{
    return a.pow(double(b));
}

MTS_NAMESPACE_BEGIN

// Albedo inversion: map spherical/bond albedo kd to single-scattering albedo c [Hitchhiker's guide approximation 3]
float c_from_kd( const float kd )
{
   return (-1.0893564101764108 - 2.3011443820100275*kd - 0.8247913150511265*Power(kd,2) + 
     0.9234846594011388*(0.39417098289205227 + kd)*Sqrt(0.5667087129883742 + kd)*
      Sqrt(15.803432476716953 + kd))/Power(0.5461928294878655 + kd,2);
}

Spectrum reflectanceMap( const Spectrum kd )
{
    return Color3(c_from_kd( kd[0] ), c_from_kd( kd[1] ), c_from_kd( kd[2] ));
}

// Hapke's improved approximation for the 3D iso scatter H function [Hapke 2002]
mitsuba::Spectrum H2( const mitsuba::Spectrum& c, const float u){
    const mitsuba::Spectrum y = (mitsuba::Spectrum(1) - c).sqrt();
    const mitsuba::Spectrum n = (mitsuba::Spectrum(1) - y) / (mitsuba::Spectrum(1) + y);
    return mitsuba::Spectrum(1)/(mitsuba::Spectrum(1) - u*(mitsuba::Spectrum(1) - y)*(n + (mitsuba::Spectrum(1) - 0.5*n - n*u)*logf((1 + u)/u)));
}

mitsuba::Spectrum diffuse_transport_BRDF( mitsuba::Spectrum kd, const float ui, const float uo )
{
    Spectrum c = reflectanceMap(kd);
    return INV_PI * c * (1.0 / 4.0 ) * ( H2(c, ui) * H2(c, uo )) / ( ui + uo );
}

class Chandra : public BSDF {
public:
    Chandra(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'reflectance' and 'diffuseReflectance' as parameter names */
        m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("reflectance") ? "reflectance"
                : "diffuseReflectance", Spectrum(.5f)));
    }

    Chandra(Stream *stream, InstanceManager *manager)
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

        return (diffuse_transport_BRDF(kd,ui,uo) * Frame::cosTheta(bRec.wo));
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        const Float ui = Frame::cosTheta(bRec.wi);
        const Float u = Frame::cosTheta(bRec.wo);

        return (powf(u,Sqrt(ui))*(1 + sqrtf(ui)))/(2.*M_PI);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        const Float ui = Frame::cosTheta(bRec.wi);
        const Float uo = powf( sample.y, 1 / ( 1 + sqrtf(ui) ) );
        const Float sinAlpha = sqrtf( 1.0 - uo * uo );
        Float phi = (2.0f * M_PI) * sample.x;
        Vector localDir = Vector(
            sinAlpha * std::cos(phi),
            sinAlpha * std::sin(phi),
            uo
        );

        
        bRec.wo = localDir;
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;

        Spectrum kd = m_reflectance->eval(bRec.its);

        return diffuse_transport_BRDF(kd,ui,uo) * uo / ( (powf(uo,Sqrt(ui))*(1 + sqrtf(ui)))/(2.*M_PI) );
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        const Float ui = Frame::cosTheta(bRec.wi);
        const Float uo = powf( sample.y, 1 / ( 1 + sqrtf(ui) ) );
        const Float sinAlpha = sqrtf( 1.0 - uo * uo );
        Float phi = (2.0f * M_PI) * sample.x;
        Vector localDir = Vector(
            sinAlpha * std::cos(phi),
            sinAlpha * std::sin(phi),
            uo
        );

        
        bRec.wo = localDir;
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        pdf = (powf(uo,Sqrt(ui))*(1 + sqrtf(ui)))/(2.*M_PI);

        Spectrum kd = m_reflectance->eval(bRec.its);

        return diffuse_transport_BRDF(kd,ui,uo) * uo / pdf;
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
        oss << "Chandra[" << endl
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

class ChandraShader : public Shader {
public:
    ChandraShader(Renderer *renderer, const Texture *reflectance)
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

Shader *Chandra::createShader(Renderer *renderer) const {
    return new ChandraShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(ChandraShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Chandra, false, BSDF)
MTS_EXPORT_PLUGIN(Chandra, "Chandrasekhar's H-function BRDF")
MTS_NAMESPACE_END
