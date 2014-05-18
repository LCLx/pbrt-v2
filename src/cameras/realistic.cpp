// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"

#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // Realistic camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   string autofocusfile = params.FindOneString("af_zones", "");
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, autofocusfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter,
                                 const string &specfile,
								 const string &autofocusfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   film(f)
{

	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.

	//	compute the film size
	filmDiag = filmdiag;
	float xRes = float(f->xResolution);
	float yRes = float(f->yResolution);
	float diagRes = sqrt(xRes * xRes + yRes * yRes);
	filmLengthInX = filmdiag / diagRes * xRes;
	filmLengthInY = filmdiag / diagRes * yRes;
	//	store the film position
	filmDist = filmdistance;
	//	store the aperture
	aperture = aperture_diameter;
	//	parse the specfile
	ParseLens(specfile);

	// If 'autofocusfile' is the empty string, then you should do
	// nothing in any subsequent call to AutoFocus()
	autofocus = false;

	if (autofocusfile.compare("") != 0)  {
		ParseAfZones(autofocusfile);
		autofocus = true;
	}
}

void RealisticCamera::ParseLens(const string& filename)
{
	ifstream specfile(filename.c_str());
	if (!specfile)
	{
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		exit(-1);
	}
	char line[512];
	filmPos = 0.0;
	while (!specfile.eof()) 
	{
		specfile.getline(line, 512);
		if (line[0] != '\0' && line[0] != '#' && 
		line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
		{
			lenses.resize(lenses.size() + 1);
			Lens& lens = lenses[lenses.size() - 1];
			sscanf(line, "%f %f %f %f\n", &lens.radius, &lens.thickness, &lens.refraction, &lens.aperture);
			//	if we find the aperture stop
			if (lens.radius == 0.0f)
			{
				lens.aperture = min(lens.aperture, aperture);
				//	aperture stop should be in air
				lens.refraction = 1.f;			
			}
			lens.zPos = filmPos;
			filmPos -= lens.thickness;
			//printf("zPos = %f aperture = %f\n", lens.zPos, lens.aperture);
		}
	}
	filmPos -= filmDist;
	//printf("filmPos = %f\n", filmPos);
	printf("Read in %zu lens from %s\n", lenses.size(), filename.c_str());

}

// parses the AF zone file
void RealisticCamera::ParseAfZones(const string& filename)
{
  ifstream specfile(filename.c_str());
   if (!specfile) {
      fprintf(stderr, "Cannot open file %s\n", filename.c_str());
      exit (-1);
   }

   char line[512];

   while (!specfile.eof()) {
      specfile.getline(line, 512);
      if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
      {
		afZones.resize(afZones.size()+1);
		AfZone& zone = afZones[afZones.size()-1];
		sscanf(line, "%f %f %f %f\n", &zone.left, &zone.right, &zone.top, &zone.bottom);
      }
   }

	printf("Read in %zu AF zones from %s\n", afZones.size(), filename.c_str());
}

RealisticCamera::~RealisticCamera()
{
	
}

bool RealisticCamera::RefractFromLens(int lensId, Ray Rin, Ray &Rout) const
{
	//	compute the intersection of lenses[lensId] and Rin
	//	if they intersect
	//	return true and build a refracted ray in Rout
	//	otherwise return false
	//	compute the center of the circle

	//	if we encounter the aperture stop
	if (lenses[lensId].radius == 0.0f)
	{
		//	the aperture shot
		//	Rin.o + t * Rin.d
		float zPos = lenses[lensId].zPos;
		float t = (zPos - Rin.o.z) / Rin.d.z;
		Point Pinter = Rin(t);
		//	decide whether Pinter is inside the aperture
		float aper = lenses[lensId].aperture;
		if (Pinter.x * Pinter.x + Pinter.y * Pinter.y > aper * aper / 4)
			return false;
		else
		{
			Rout = Rin;
			return true;
		}
	}
	//	else
	Lens lens = lenses[lensId];
	float lensAper = lens.aperture;
	float radius = lens.radius;
	float radiusSq = radius * radius;
	float zCenter = lens.zPos - radius;
	//	compute the distance between center and the ray
	Point center(0.0, 0.0, zCenter);
	//	(Rin.o + t * Rin.d - center) * Rin.d = 0
	//	(Rin.o - center) * Rin.d + t * Rin.d * Rin.d = 0;
	float dSq = Rin.d.LengthSquared();
	float t = Dot(Rin.o - center, Rin.d) / -dSq;
	Point midPoint = Rin(t);
	float midSq = (midPoint - center).LengthSquared();
	if (midSq > radiusSq)	//	no intersection
		return false;
	//	find the intersection point
	float dLen = Rin.d.Length();
	float deltaT = sqrt((radiusSq - midSq)) / dLen;
	float t0 = t - deltaT;
	float t1 = t + deltaT;
	float tInt;
	if (t0 > 0.0)
		tInt = t0;
	else if (t1 > 0.0)
		tInt = t1;
	else
		return false;		//	no intersection with t > 0.0
	Point Pinter = Rin(tInt);
	float apSq = Pinter.x * Pinter.x + Pinter.y * Pinter.y;
	if (apSq > lensAper * lensAper / 4)
		return false;		//	no intersection within the aperture
	Rout.o = Pinter;
	//	now decide the direction of the refracted light
	Vector normal = Normalize(Pinter - center);
	if (radius < 0.0f)
		normal = -normal;
	float n2 = lens.refraction;
	float n1 = 1.0f;
	if (lensId > 0)
		n1 = lenses[lensId - 1].refraction;
	float cos2 = Dot(normal, Rin.d) / dLen;
	float sin2 = sqrt(1 - cos2 * cos2);
	float sin1 = n2 * sin2 / n1;
	if (sin1 > 1.0f)
		return false;
	//	now let's build the dir vector
	float cos1 = sqrt(1 - sin1 * sin1);
	float tan1 = sin1 / cos1;
	Vector v1 = Cross(Rin.d, normal);
	Vector v2 = Normalize(Cross(normal, v1));
	Rout.d = normal + tan1 * v2;
	return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{

  // YOUR CODE HERE -- make that ray!

  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens
	
	//	the range of imageX: 0 to xResolution
	//	the range of imageY: 0 to yResolution
	float xCamera, yCamera;
	xCamera = -(sample.imageX * 1.0 / film->xResolution - 0.5) * filmLengthInX;
	yCamera = (sample.imageY * 1.0 / film->yResolution - 0.5) * filmLengthInY;
	Point Pcamera(xCamera, yCamera, filmPos);

	//	generate lens samples
	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	//	the last lens
	int nLens = (int)lenses.size();
	double LensRadius = lenses[nLens - 1].aperture / 2;
	lensU *= LensRadius;
	lensV *= LensRadius;
	//	now we have a ray from Pcamera, and it points towards 
	Point Phit(lensU, lensV, lenses[nLens - 1].zPos);
	Ray Rin(Pcamera, Phit - Pcamera, 0.f, INFINITY);
	//	start to iterate all the lenses
	Ray Rout;
	for (int i = nLens - 1; i >= 0; i--)
	{
		bool succeed = RefractFromLens(i, Rin, Rout);
		//	if Rin and lens won't intersect
		//	the function will return false		
		if (!succeed)
			return 0.0f;
		//	update Rin
		Rin = Rout;
	}
	//	return Rin
	*ray = Ray(Rin.o, Normalize(Vector(Rin.d)), 0.f, INFINITY);
	ray->time = sample.time;
	CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);
	//	compute the weight
	Vector XX = Phit - Pcamera;
	float XXlenSq = XX.LengthSquared();
	float w = XX.z * LensRadius / XXlenSq;
	return M_PI * w * w;
}

void  RealisticCamera::AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample) {
	// YOUR CODE HERE:
	// The current code shows how to create a new Sampler, and Film cropped to the size of the auto focus zone.
	// It then renders the film, producing rgb values.  You need to:
	//
	// 1. Modify this code so that it can adjust film plane of the camera
	// 2. Use the results of raytracing to evaluate whether the image is in focus
	// 3. Search over the space of film planes to find the best-focused plane.

	if(!autofocus)
		return;

	for (size_t i=0; i<afZones.size(); i++) {

		AfZone & zone = afZones[i];

		RNG rng;
		MemoryArena arena;
		Filter * filter = new BoxFilter(.5f,.5f);
		const float crop[] = {zone.left,zone.right,zone.top,zone.bottom};
		ImageFilm sensor(film->xResolution, film->yResolution, filter, crop,"foo.exr",false);
		int xstart,xend,ystart,yend;
		sensor.GetSampleExtent(&xstart,&xend,&ystart,&yend);

		StratifiedSampler sampler(xstart, xend, ystart, yend,
		                          16, 16, true, ShutterOpen, ShutterClose);

		// Allocate space for samples and intersections
		int maxSamples = sampler.MaximumSampleCount();
		Sample *samples = origSample->Duplicate(maxSamples);
		RayDifferential *rays = new RayDifferential[maxSamples];
		Spectrum *Ls = new Spectrum[maxSamples];
		Spectrum *Ts = new Spectrum[maxSamples];
		Intersection *isects = new Intersection[maxSamples];

		// Get samples from _Sampler_ and update image
		int sampleCount;
		while ((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
			// Generate camera rays and compute radiance along rays
			for (int i = 0; i < sampleCount; ++i) {
				// Find camera ray for _sample[i]_

				float rayWeight = this->GenerateRayDifferential(samples[i], &rays[i]);
				rays[i].ScaleDifferentials(1.f / sqrtf(sampler.samplesPerPixel));


				// Evaluate radiance along camera ray

				if (rayWeight > 0.f)
					Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
													 arena, &isects[i], &Ts[i]);
				else {
					Ls[i] = 0.f;
					Ts[i] = 1.f;
				}

				// Issue warning if unexpected radiance value returned
				if (Ls[i].HasNaNs()) {
					Error("Not-a-number radiance value returned "
						  "for image sample.  Setting to black.");
					Ls[i] = Spectrum(0.f);
				}
				else if (Ls[i].y() < -1e-5) {
					Error("Negative luminance value, %f, returned"
						  "for image sample.  Setting to black.", Ls[i].y());
					Ls[i] = Spectrum(0.f);
				}
				else if (isinf(Ls[i].y())) {
					Error("Infinite luminance value returned"
						  "for image sample.  Setting to black.");
					Ls[i] = Spectrum(0.f);
				}

			}

			// Report sample results to _Sampler_, add contributions to image
			if (sampler.ReportResults(samples, rays, Ls, isects, sampleCount))
			{
				for (int i = 0; i < sampleCount; ++i)
				{

					sensor.AddSample(samples[i], Ls[i]);

				}
			}

			// Free _MemoryArena_ memory from computing image sample values
			arena.FreeAll();
		}

		float * rgb;
		int width;
		int height;
		sensor.WriteRGB(&rgb,&width,&height,1.f);
		// YOUR CODE HERE! The rbg contents of the image for this zone
		// are now stored in the array 'rgb'.  You can now do whatever
		// processing you wish


		//you own rgb  now so make sure to delete it:
		delete [] rgb;
		//if you want to see the output rendered from your sensor, uncomment this line (it will write a file called foo.exr)
		//sensor.WriteImage(1.f);


		delete[] samples;
		delete[] rays;
		delete[] Ls;
		delete[] Ts;
		delete[] isects;
	}
}
