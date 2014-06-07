/*
Tao Du
taodu@stanford.edu
Jun 4, 2014
*/

// volumes/aurora.cpp*
#include "stdafx.h"
#include "volumes/aurora.h"
#include "paramset.h"

void AuroraGrid::AddPhoton(const AuroraPhoton &photon)
{
	Point p = photon.p;
	if (!extent.Inside(p))
		return;
	Vector vox = photon.p - extent.pMin;
    int vx = int(vox.x / step);
	int vy = int(vox.y / step);
	int vz = int(vox.z / step);
	if (vx == nx || vy == ny || vz == nz)
		return;
	grids[vz*nx*ny + vy*nx + vx].photons.push_back(photon);
}

void AuroraGrid::SearchInGrid(const Point &p, float &r, float &g, float &b) const
{
	if (!extent.Inside(p))
		return;
	Vector vox = p - extent.pMin;
	float fx = vox.x / step;
	float fy = vox.y / step;
	float fz = vox.z / step;
	int vx = int(fx); int vy = int(fy); int vz = int(fz);
	int xmin = fx - vx > 0.5f ? vx : vx - 1;
	int ymin = fy - vy > 0.5f ? vy : vy - 1;
	int zmin = fz - vz > 0.5f ? vz : vz - 1;
	xmin = Clamp(xmin, 0, nx - 1);
	ymin = Clamp(ymin, 0, ny - 1);
	zmin = Clamp(zmin, 0, nz - 1);
	int xmax = xmin + 1; int ymax = ymin + 1; int zmax = zmin + 1;
	xmax = Clamp(xmax, 0, nx - 1);
	ymax = Clamp(ymax, 0, ny - 1);
	zmax = Clamp(zmax, 0, nz - 1);

	//	scan all the grids within the range above
	r = 0.f;
	g = 0.f;
	b = 0.f;
	float radiusSq = radius * radius;
	//	sigam in gaussian kernel
	float sigma = radius;
	float halfInvSigmaSq = .5f / sigma / sigma;
	float sum = 0.f;
	for (int i = xmin; i <= xmax; i++)
		for (int j = ymin; j <= ymax; j++)
			for (int k = zmin; k <= zmax; k++)
			{
				AuroraVoxel voxel = grids[k*nx*ny + j*nx + i];
				int photonNum = voxel.PhotonNum();
				for (int s = 0; s < photonNum; s++)
				{
					//	get the position of the photon
					AuroraPhoton photon = voxel.photons[s];
					Vector v = photon.p - p;
					float lenSq = v.LengthSquared();
					if (lenSq >= radiusSq)
						continue;
					//	compute the gaussian weight
					float weight = expf(-lenSq * halfInvSigmaSq);
					r += (weight * photon.r);
					g += (weight * photon.g);
					b += (weight * photon.b);
					sum += weight;
				}
			}
	r /= sum;
	g /= sum;
	b /= sum;
}

AuroraDensity::~AuroraDensity()
{
	delete []eleDensity;
}

float AuroraDensity::EleDensity(const Point &Pobj) const
{
	if (!extent.Inside(Pobj)) return 0;
    // Compute voxel coordinates and offsets for _Pobj_
    Vector vox = extent.Offset(Pobj);
    vox.x *= nx;
    vox.y *= ny;
    vox.z *= nz;
    int vx = (int)vox.x;
	int vy = (int)vox.y;
	int vz = (int)vox.z;
	if (vx == nx || vy == ny || vz == nz)
		return 0;
    float dx = vox.x - vx, dy = vox.y - vy, dz = vox.z - vz;
    // Trilinearly interpolate density values to compute local density
    float d00 = Lerp(dx, D(vx, vy, vz),     D(vx+1, vy, vz));
    float d10 = Lerp(dx, D(vx, vy+1, vz),   D(vx+1, vy+1, vz));
    float d01 = Lerp(dx, D(vx, vy, vz+1),   D(vx+1, vy, vz+1));
    float d11 = Lerp(dx, D(vx, vy+1, vz+1), D(vx+1, vy+1, vz+1));
    float d0 = Lerp(dy, d00, d10);
    float d1 = Lerp(dy, d01, d11);
    return Lerp(dz, d0, d1);
}

float AuroraDensity::Density(const Point &Pobj) const
{
    if (!extent.Inside(Pobj)) return 0;
    float height = Dot(Pobj - extent.pMin, upDir);
    return a * expf(-b * height);
}

Spectrum AuroraDensity::Lve(const Point &p, const Vector &, float) const
{
	const Point Pobj = WorldToVolume(p);
	float rgb[3];
	grid.SearchInGrid(Pobj, rgb[0], rgb[1], rgb[2]);
	rgb[0] = 0.5f;
	rgb[1] = 0.1f;
	rgb[2] = 0.2f;
	return Spectrum::FromRGB(rgb);
}

Spectrum AuroraDensity::tau(const Ray &r, float stepSize, float u) const 
{
    float t0, t1;
    float length = r.d.Length();
    if (length == 0.f) return 0.f;
    Ray rn(r.o, r.d / length, r.mint * length, r.maxt * length, r.time);
    if (!IntersectP(rn, &t0, &t1)) return 0.;
    Spectrum tau(0.);
    t0 += u * stepSize;
    while (t0 < t1) {
        tau += sigma_t(rn(t0), -rn.d, r.time);
        t0 += stepSize;
    }
    return tau * stepSize;
}

// AuroraDensity Method Definitions
AuroraDensity *CreateAuroraVolumeRegion(const Transform &volume2world,
        const ParamSet &params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.);
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
	float a = params.FindOneFloat("a", 1.);
    float b = params.FindOneFloat("b", 1.);
    Vector up = params.FindOneVector("updir", Vector(0,1,0));
	float radius = params.FindOneFloat("radius", 1.);
	int nitems;
    const float *data = params.FindFloat("density", &nitems);
    if (!data) {
        Error("No \"density\" values provided for volume grid?");
        return NULL;
    }
    int nx = params.FindOneInt("nx", 1);
    int ny = params.FindOneInt("ny", 1);
    int nz = params.FindOneInt("nz", 1);
    if (nitems != nx*ny*nz) {
        Error("VolumeGridDensity has %d density values but nx*ny*nz = %d",
            nitems, nx*ny*nz);
        return NULL;
    }
	string rcolor = params.FindOneFilename("aurora_r", "aurora_r.txt.even");
	string gcolor = params.FindOneFilename("aurora_g", "aurora_g.txt.even");
	string bcolor = params.FindOneFilename("aurora_b", "aurora_b.txt.even");
	string intensity = params.FindOneFilename("aurora_intensity", "rays.txt.even");
	//	TODO: what is our parameter for AuroraDensity?
	return new AuroraDensity(sigma_a, sigma_s, g, BBox(p0, p1),
        volume2world, a, b, up, nx, ny, nz, radius, data, rcolor, gcolor, bcolor, intensity);
}