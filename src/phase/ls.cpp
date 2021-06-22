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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

// Far-field scattering from a Lambertian sphere [Lambert 1760, Blinn 1982, d'Eon 2021]
class LSPhaseFunction : public PhaseFunction {
public:
    LSPhaseFunction(const Properties &props)
        : PhaseFunction(props) {
    }

    LSPhaseFunction(Stream *stream, InstanceManager *manager)
        : PhaseFunction(stream, manager) {
        configure();
    }

    virtual ~LSPhaseFunction() { }

    void configure() {
        PhaseFunction::configure();
        m_type = EAngleDependence;
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        PhaseFunction::serialize(stream, manager);
    }

    Float sample(PhaseFunctionSamplingRecord &pRec,
            Sampler *sampler) const {
        Point2 sample(sampler->next2D());

        // approximate CDF inverse: [d'Eon 2021 - EGSR]
        Float cosTheta = 1.0f - 2.0f * powf(1.0f - powf(sample.x, 1.01938f + 0.0401885 * sample.x), 0.397225f);

        Float sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
              sinPhi, cosPhi;

        math::sincos(2*M_PI*sample.y, &sinPhi, &cosPhi);

        pRec.wo = Frame(-pRec.wi).toWorld(Vector(
            sinTheta * cosPhi,
            sinTheta * sinPhi,
            cosTheta
        ));

        return 1.0f;
    }

    Float sample(PhaseFunctionSamplingRecord &pRec,
            Float &pdf, Sampler *sampler) const {
        LSPhaseFunction::sample(pRec, sampler);
        pdf = LSPhaseFunction::eval(pRec);
        return 1.0f;
    }

    Float eval(const PhaseFunctionSamplingRecord &pRec) const {
        const Float u = -dot(pRec.wi, pRec.wo);
        return (2.0f*(sqrt(1.0f - u * u) - u * acos(u)))/(3.0f * M_PI * M_PI);
    }

    Float getMeanCosine() const {
        return -0.4444444f;
    }

    std::string toString() const {
        return "LSPhaseFunction[]";
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(LSPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(LSPhaseFunction, "Lambert-sphere phase function");
MTS_NAMESPACE_END
