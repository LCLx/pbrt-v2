/*
Tao Du
taodu@stanford.edu
Jun 4, 2014
*/

// volumes/aurora.cpp*
#include "stdafx.h"
#include "volumes/aurora.h"
#include "paramset.h"

// AuroraDensity Method Definitions
AuroraDensity *CreateAuroraVolumeRegion(const Transform &volume2world,
        const ParamSet &params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.);
    Spectrum Le = params.FindOneSpectrum("Le", 0.);
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
	//	TODO: what is our parameter for AuroraDensity?
	return new AuroraDensity(sigma_a, sigma_s, g, Le, BBox(p0, p1),
        volume2world);
}