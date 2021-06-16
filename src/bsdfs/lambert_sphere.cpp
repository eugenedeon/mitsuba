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

// This implements the Lambert-sphere BRDF for a half space comprised of spherical Lambertian particles [d'Eon 2021 - EGSR].

// Albedo mapping from diffuse color kd to single-scattering albedo c [d'Eon 2021, Eq.(46)]
float c_from_kd( const float kd )
{
   return (1 - Power(1 - kd, 2.73556))/(1 - 0.184096*Power(1 - kd, 2.48423));
}

Spectrum reflectanceMap( const Spectrum kd )
{
    return Color3(c_from_kd( kd[0] ), c_from_kd( kd[1] ), c_from_kd( kd[2] ));
}

// far-field phase function for a spherical Lambertian particle [Lambert 1760,Blinn 1982]
float lambertSpherePhaseFunction( const float u )
{
    return (2*(sqrtf(1 - u*u) - u*acosf(u)))/(3.*M_PI * M_PI);
}

// zeroth-order H function [d'Eon 2021, Eq.(37)]
float H0(const float u, const float c)
{
    const float C = sqrtf( 1.0f - c );
    const float a = (8.216443463470172 + 1.501116252342486*Power(C,6.054351911707695))/(4.175932557866179 - 1.2122216865199813*C);
    const float d = (7.773097842804312 - 0.5658108102188075*Power(C,0.9615460393324836))/
   (8.659120811349371 - 0.15997430541238322*Power(C,7));
   return (1 + a*Power(u,d))/(1 + a*Sqrt(1 - 2*((89*c)/288. + (59*c*c)/288. - Power(c,3)/72.))*Power(u,d));
}

// zeroth-order H function [d'Eon 2021, Eq.(28)]
float H1(const float u, const float c)
{
    return exp((-0.14483935877260054*c + 0.024285125615255733*c*c)*
    Power(u,0.45944184025347456 - 1.0787917226917565*u + 1.8572804924335546*Power(u,2) - 
      1.1283119660147671*Power(u,3)));
}

// zeroth-order component of the BRDF
float f0(const float ui, const float uo, const float c)
{
    const float A = (69*c)/128.;
    const float B = (-0.08446297239646791 + 0.5153566145541554*Sqrt(1 - c) - 0.77757371002123*(1 - c) + 
     0.34668869623791543*Power(1 - c,1.5))/
   (0.9648932738041012 - 0.6655015709030936*Sqrt(1 - c) + 0.1826019462608555*(1 - c));
    const float C = (682.8479477533338 - 2567.7368047535556*Sqrt(1 - c) + 7487.987105705168*(1 - c) - 
     5602.448801045478*Power(1 - c,1.5))/
   (5850.602606063662 - 4008.3309624647227*Sqrt(1 - c) + 1480.250743805733*(1 - c));
    const float D = (0.2855294320307508 + 160.39651500649123*Sqrt(1 - c) - 327.42799697993706*(1 - c) + 
     166.88327184107732*Power(1 - c,1.5))/
   (674.1908010450103 - 412.9837444306491*Sqrt(1 - c) + 596.4232294419696*(1 - c));
    const float E = (15*(1 - c)*c*(3 + (4*c)/3.))/128.;
    const float F = (-1.9208967199948512 - 242.16001167844007*Sqrt(1 - c) - 21.914139454773085*(1 - c) + 
     266.06342182761813*Power(1 - c,1.5))/
   (1499.904420175135 + 457.4200839912641*Sqrt(1 - c) + 215.77342164754094*(1 - c));
    return 0.5 * INV_PI * ( H0( ui, c ) ) * ( H0( uo, c ) ) / ( ui + uo ) * ( A + C*ui*uo + E*Power(ui,2)*Power(uo,2) + B*(ui + uo) + D*ui*uo*(ui + uo) + 
   F*(Power(ui,2) + Power(uo,2)) );
}

// 1st-order component of the BRDF divided by sin(theta_i)*sin(theta_o)
float f1m(const float ui, const float uo, const float c)
{
    const float l = -0.05890107167021953*c - 0.004740543453386809*Power(c,2);

    return 2.0 * INV_PI * (-(1*
      ((c*(64 + 45*ui*uo))/48. - (15*(1 + 0.44*c)*c*ui*uo*H1(ui,c)*H1(uo,c))/16. - 
        (4*c*(1 + l*ui)*(1 + l*uo)*H1(ui,c)*H1(uo,c))/3.))/(8.*(ui + uo)));
}

Spectrum non_negative_color3( const Spectrum& c )
{
    return Color3( std::max( 0.0f, c[0] ), std::max( 0.0f, c[1] ), std::max( 0.0f, c[2] ) );
}

Spectrum lambertSphereBRDF( Spectrum c, const float ui, const float uo, const Vector3& wi, const Vector3& wo )
{
    const float iodot = dot(wi,wo);
    Spectrum f_single = c * lambertSpherePhaseFunction(-iodot) / ( ui + uo );
    const float cosphisinisino = ( iodot - ui * uo ); // cos(phi) * sin(theta_i) * sin(theta_o)
    
    const Spectrum F0single = INV_PI * c * (207 + 256*ui*uo - 45*Power(uo,2) + 45*Power(ui,2)*(-1 + 3*Power(uo,2)))/(768.*(ui + uo));

    Spectrum F0 = Color3( f0( ui, uo, c[0] ), f0( ui, uo, c[1] ), f0( ui, uo, c[2] ) );
    Spectrum F1m = Color3( f1m( ui, uo, c[0] ), f1m( ui, uo, c[1] ), f1m( ui, uo, c[2] ) );

    return non_negative_color3( f_single + F0 - F0single + F1m * cosphisinisino );
}

class LambertSphere : public BSDF {
public:
    LambertSphere(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'reflectance' and 'diffuseReflectance' as parameter names */
        m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("reflectance") ? "reflectance"
                : "diffuseReflectance", Spectrum(.5f)));
    }

    LambertSphere(Stream *stream, InstanceManager *manager)
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

        Spectrum c = reflectanceMap( m_reflectance->eval(bRec.its) );

        return lambertSphereBRDF(c,ui,uo,bRec.wi,bRec.wo) * Frame::cosTheta(bRec.wo);
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
        Spectrum c = reflectanceMap( m_reflectance->eval(bRec.its) );

        return lambertSphereBRDF(c,ui,uo,bRec.wi,bRec.wo) / Spectrum( INV_PI );
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
        Spectrum c = reflectanceMap( m_reflectance->eval(bRec.its) );

        return lambertSphereBRDF(c,ui,uo,bRec.wi,bRec.wo) / Spectrum( INV_PI );
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
        oss << "LambertSphere[" << endl
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

class LambertSphereShader : public Shader {
public:
    LambertSphereShader(Renderer *renderer, const Texture *reflectance)
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

Shader *LambertSphere::createShader(Renderer *renderer) const {
    return new LambertSphereShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(LambertSphereShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(LambertSphere, false, BSDF)
MTS_EXPORT_PLUGIN(LambertSphere, "Lambert-Sphere BRDF")
MTS_NAMESPACE_END
