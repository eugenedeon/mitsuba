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

mitsuba::Spectrum powf( const mitsuba::Spectrum a, const int b )
{
    return a.pow(double(b));
}

// This implements the "Diffusion-transport" BRDF for a half space with isotropic scattering
// and Gamma/Erlang-2 steps between collisions [d'Eon 2019, d'Eon and Krivanek 2020]

mitsuba::Spectrum diffuse_transport_BRDF( mitsuba::Spectrum c, const float ui, const float uo )
{
    const mitsuba::Spectrum sqrt1c = ( mitsuba::Spectrum(1.0f) - c ).sqrt();

    return (c*(2*Power(uo,2)*(mitsuba::Spectrum(1 + 2*uo) + sqrt1c*Power(uo,2)) + 
       Power(ui,4)*(2*sqrt1c + (mitsuba::Spectrum(1.0f) + 3*sqrt1c - c)*uo - 
          2*(mitsuba::Spectrum(-1.0f) + c)*Power(uo,2)) + 
       ui*uo*(mitsuba::Spectrum(6) + 4*(mitsuba::Spectrum(3) + sqrt1c)*uo + (mitsuba::Spectrum(6) + 8*sqrt1c)*Power(uo,2) + 
          (mitsuba::Spectrum(1) + 3*sqrt1c - c)*Power(uo,3)) + 
       Power(ui,3)*(mitsuba::Spectrum(4) + (mitsuba::Spectrum(6) + 8*sqrt1c)*uo + 
          (mitsuba::Spectrum(3) + 13*sqrt1c - 3*c)*Power(uo,2) - 6*(mitsuba::Spectrum(-1.0) + c)*Power(uo,3)) + 
       Power(ui,2)*(mitsuba::Spectrum(2) + 4*(mitsuba::Spectrum(3) + sqrt1c)*uo + 
          2*(mitsuba::Spectrum(8) + 6*sqrt1c - c)*Power(uo,2) + 
          (mitsuba::Spectrum(3) + 13*sqrt1c - 3*c)*Power(uo,3) - 2*(mitsuba::Spectrum(-1) + c)*Power(uo,4))))/
   (8.*Pi*Power(mitsuba::Spectrum(1) + sqrt1c*ui,2)*Power(ui + uo,3)*Power(mitsuba::Spectrum(1) + sqrt1c*uo,2));
}

mitsuba::Spectrum cFromR( const mitsuba::Spectrum & in_R )
{
    return (4*in_R)/((mitsuba::Spectrum(1) + in_R)*(mitsuba::Spectrum(1) + in_R));
}

MTS_NAMESPACE_BEGIN

class DiffusionTransport : public BSDF {
public:
    DiffusionTransport(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'reflectance' and 'diffuseReflectance' as parameter names */
        m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("reflectance") ? "reflectance"
                : "diffuseReflectance", Spectrum(.5f)));
    }

    DiffusionTransport(Stream *stream, InstanceManager *manager)
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

        Spectrum c = cFromR( m_reflectance->eval(bRec.its) );

        return (diffuse_transport_BRDF(c,ui,uo) * Frame::cosTheta(bRec.wo));
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

        Spectrum c = cFromR( m_reflectance->eval(bRec.its) );

        return diffuse_transport_BRDF(c,ui,uo) * uo / ( (powf(uo,Sqrt(ui))*(1 + sqrtf(ui)))/(2.*M_PI) );
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

        Spectrum c = cFromR( m_reflectance->eval(bRec.its) );

        return diffuse_transport_BRDF(c,ui,uo) * uo / pdf;
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
        oss << "DiffusionTransport[" << endl
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

class DiffusionTransportShader : public Shader {
public:
    DiffusionTransportShader(Renderer *renderer, const Texture *reflectance)
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

Shader *DiffusionTransport::createShader(Renderer *renderer) const {
    return new DiffusionTransportShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(DiffusionTransportShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(DiffusionTransport, false, BSDF)
MTS_EXPORT_PLUGIN(DiffusionTransport, "Diffusion Transport BRDF")
MTS_NAMESPACE_END
