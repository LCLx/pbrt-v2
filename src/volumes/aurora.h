/*
Tao Du
taodu@stanford.edu
Jun 4, 2014
*/

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_VOLUMES_AURORA_H
#define PBRT_VOLUMES_AURORA_H

// volumes/aurora.h*
#include "volume.h"

// AuroraDensity Declarations
class AuroraDensity : public DensityRegion {
public:
    // AuroraDensity Public Methods
    AuroraDensity(const Spectrum &sa, const Spectrum &ss,
                       float gg, const Spectrum &emit, const BBox &e,
                       const Transform &v2w)
        : DensityRegion(sa, ss, gg, emit, v2w), extent(e) 
	{
		//	TODO
    }
    BBox WorldBound() const { return Inverse(WorldToVolume)(extent); }
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
    float Density(const Point &Pobj) const {
        if (!extent.Inside(Pobj)) return 0;
		//	TODO: return the density in Pobj
		/*
		float height = Dot(Pobj - extent.pMin, upDir);
        return a * expf(-b * height);
		*/
		return 0.f;
    }
private:
    // AuroraDensity Private Data
    BBox extent;
	//	TODO
};

AuroraDensity *CreateAuroraVolumeRegion(const Transform &volume2world,
        const ParamSet &params);

#endif // PBRT_VOLUMES_AURORA_H
