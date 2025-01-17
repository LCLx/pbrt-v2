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
#include "../core/splines.h"
#include "../core/perlin2d.h"

#include <fstream>

struct AuroraPhoton
{
	AuroraPhoton(const Point& pp, float rr, float gg, float bb)
		: p(pp), r(rr), g(gg), b(bb)
	{
	}
	Point p;
	float r, g, b;
};

struct AuroraVoxel
{
	AuroraVoxel()
	{
		photons.clear();
	}
	~AuroraVoxel()
	{
		photons.clear();
	}
	int PhotonNum() {return (int)photons.size();}
	vector<AuroraPhoton> photons;
};

class AuroraGrid
{
public:
	AuroraGrid(const BBox &e, int x, int y, int z, float r);
	~AuroraGrid(){delete []grids;}

	//	add a new photon into the grid
	void AddPhoton(const AuroraPhoton &photon);
	//	given a position p and a search radius, find the nearby photon
	//	then use gaussian kernel to average them
	void SearchInGrid(const Point &p, float &r, float &g, float &b) const;

	float LoadFactor();

private:
	BBox extent;
	float step;
	float radius;
	int nx, ny, nz;
	AuroraVoxel *grids;

    Perlin2D noise2d;
};


// AuroraDensity Declarations
class AuroraDensity : public VolumeRegion {
public:
	// AuroraRegion Public Methods
    AuroraDensity(const Spectrum &sa, const Spectrum &ss, float gg, const BBox &e,
                  const Transform &VolumeToWorld, float aa, float bb,
                       const Vector &up, int x, int y, int z, float radius, const float *d,
					   float cr, float cg, float cb, float is,
					   const string rcolor, const string gcolor, const string bcolor, const string intensity,
					   const Vector &b, int bn, float a, float l, float dd, float dA, float ee,
                       float _persistence, float _frequency, int _level, int _seed, float _magnitude, float _perturb)
        : sig_a(sa), sig_s(ss), g(gg), extent(e),
          WorldToVolume(Inverse(VolumeToWorld)), a(aa), b(bb), 
		  nx(x), ny(y), nz(z), grid(e, nx, ny, nz, radius),
		  B(b), beamNum(bn), alphaD(a), L(l), dt(dd), deltaAlpha(dA), eleThreshold(ee),
          persistence(_persistence), frequency(_frequency), level(_level), seed(_seed), magnitude(_magnitude), perturb(_perturb)
	{
		upDir = Normalize(up);
		eleDensity = new float[nx*ny*nz];
        memcpy(eleDensity, d, nx * ny * nz * sizeof(float));
		//	color and intensity
		//	read aurora color information
		std::ifstream fin;
		float height, value;
		//	r color
		fin.open(rcolor.c_str());
		while (fin >> height)
		{
			fin >> value;
			value *= cr;
			auroraColor[0].Add_Sample(height, value);
		}
		auroraColor[0].Build();
		fin.close();
		//	g color
		fin.open(gcolor.c_str());
		while (fin >> height)
		{
			fin >> value;
			value *= cg;
			auroraColor[1].Add_Sample(height, value);
		}
		auroraColor[1].Build();
		fin.close();
		//	b color
		fin.open(bcolor.c_str());
		while (fin >> height)
		{
			fin >> value;
			value *= cb;
			auroraColor[2].Add_Sample(height, value);
		}
		auroraColor[2].Build();
		fin.close();
		//	intensity
		fin.open(intensity.c_str());
		while (fin >> height)
		{
			fin >> value;
			value *= is;
			auroraIntensity.Add_Sample(height, value);
		}
		auroraIntensity.Build();
		fin.close();
		//	generate photons
		GeneratePhotons();
	}
	~AuroraDensity();
	BBox WorldBound() const { return Inverse(WorldToVolume)(extent); }
	bool IntersectP(const Ray &r, float *t0, float *t1) const 
	{
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
	float Density(const Point &Pobj) const;
    Spectrum sigma_a(const Point &p, const Vector &, float) const 
	{
        return Density(WorldToVolume(p)) * sig_a;
    }
    Spectrum sigma_s(const Point &p, const Vector &, float) const 
	{
        return Density(WorldToVolume(p)) * sig_s;
    }
    Spectrum sigma_t(const Point &p, const Vector &, float) const 
	{
        return Density(WorldToVolume(p)) * (sig_a + sig_s);
    }
    Spectrum Lve(const Point &p, const Vector &, float) const;	//	TODO
    float p(const Point &p, const Vector &w, const Vector &wp, float) const 
	{
        return PhaseHG(w, wp, g);
    }
    Spectrum tau(const Ray &r, float stepSize, float offset) const;

	float EleDensity(const Point &Pobj) const;

	float D(int x, int y, int z) const {
        x = Clamp(x, 0, nx-1);
        y = Clamp(y, 0, ny-1);
        z = Clamp(z, 0, nz-1);
        return eleDensity[z*nx*ny + y*nx + x];
    }

	void GeneratePhotons();
private:
	Spectrum sig_a, sig_s;
    float g;
    Transform WorldToVolume;
	//BBox is used to calculate the max and min point of 2 points.
	//(0,0,2) and (1,1,1): extent.pMax = (1,1,2), extent.pMin = (0,0,1)
	BBox extent;
	
	//	use a, b and upDir to compute density
	float a, b;
    Vector upDir;

	//	eleDensity represents whether aurora 'exists' or not
	//	use eleDensity to generate the photon in aurora
	float *eleDensity;
	float eleThreshold;
	//	the number of bins in x, y, z in eDensity
	int nx, ny, nz;

	//	aurora grid
	AuroraGrid grid;

	//	the direction of the geomagnetic field
	Vector B;
	int beamNum;	//	number of beams to generate
	float alphaD;	//	alpha angle
	float L;		//	the length of the beam
	float dt;			
	float deltaAlpha;
	//Catmull_Rom is a class used to interpolate, defined in core/splines.h
	Catmull_Rom auroraColor[3];
	Catmull_Rom auroraIntensity;
                       
    float persistence;
    float frequency;
    int level;
    int seed;
    float magnitude;
    float perturb;
};


AuroraDensity *CreateAuroraVolumeRegion(const Transform &volume2world,
        const ParamSet &params);

#endif // PBRT_VOLUMES_AURORA_H
