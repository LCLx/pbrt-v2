
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// accelerators/bvh.cpp*
#include "stdafx.h"
#include "accelerators/bvh.h"
#include "probes.h"
#include "paramset.h"
#include "core/timer.h"

//  LBVH Morton Code data structure
struct LBVHPrimitiveInfo {
    LBVHPrimitiveInfo() {}
    LBVHPrimitiveInfo(uint32_t pn, const BBox &b)
        : primitiveNumber(pn), bounds(b)
    {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }

    uint32_t primitiveNumber;
    Point centroid;
    BBox bounds;
    uint32_t mortonCode;
};

uint32_t part1by2(uint32_t x)
{
    x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    return x;
}

uint32_t compMortonCode(uint32_t x, uint32_t y, uint32_t z)
{
    return (part1by2(z) << 2) + (part1by2(y) << 1) + (part1by2(x));
}

// BVHAccel Local Declarations
struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() { }
    BVHPrimitiveInfo(int pn, const BBox &b)
        : primitiveNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int primitiveNumber;
    Point centroid;
    BBox bounds;
};


void countSortMortonCode(vector<LBVHPrimitiveInfo> &buildData,
        vector<uint32_t> &index, int byteId)
{
    //  current order is given by index
    int count[256];
    for (int i = 0; i < 256; i++)
        count[i] = 0;
    uint32_t length = buildData.size();
    vector<uint32_t> keys;
    keys.reserve(length);
    for (uint32_t i = 0; i < length; i++)
    {
        uint32_t key = (buildData[index[i]].mortonCode >> (byteId * 8)) & 0x000000ff;
        keys.push_back(key);
        count[key]++;
    }

    int total = 0;
    for (int i = 0; i < 256; i++)
    {
        int oldCount = count[i];
        count[i] = total;
        total += oldCount;
    }
    vector<uint32_t> nIndex = index;
    for (uint32_t i = 0; i < length; i++)
    {
        uint32_t key = keys[i];
        index[count[key]] = nIndex[i];
        count[key]++;
    }
    nIndex.clear();
}

void radixSortMortonCode(vector<LBVHPrimitiveInfo> &buildData)
{
    //  sort the buildData.mortonCode by using radix sort + counting sort
    //  init an index vector
    vector<uint32_t> index;
    index.reserve(buildData.size());
	
    for (uint32_t i = 0; i < buildData.size(); i++)
        index.push_back(i);
    for (int round = 0; round < 4; round++)
    {
        countSortMortonCode(buildData, index, round);
    }
    vector<LBVHPrimitiveInfo> nBuildData;
    nBuildData.reserve(buildData.size());
    for (uint32_t i = 0; i < buildData.size(); i++)
        nBuildData.push_back(buildData[index[i]]);
    buildData.swap(nBuildData);
}

struct BVHBuildNode {
    // BVHBuildNode Public Methods
    BVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(uint32_t axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    BVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};


struct CompareToMid {
    CompareToMid(int d, float m) { dim = d; mid = m; }
    int dim;
    float mid;
    bool operator()(const BVHPrimitiveInfo &a) const {
        return a.centroid[dim] < mid;
    }
};


struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHPrimitiveInfo &a,
                    const BVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};


struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const BVHPrimitiveInfo &p) const;

    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};


bool CompareToBucket::operator()(const BVHPrimitiveInfo &p) const {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    Assert(b >= 0 && b < nBuckets);
    return b <= splitBucket;
}


struct LinearBVHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };

    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};


static inline bool IntersectP(const BBox &bounds, const Ray &ray,
        const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}



// BVHAccel Method Definitions
BVHAccel::BVHAccel(const vector<Reference<Primitive> > &p,
                   uint32_t mp, const string &sm) {
    maxPrimsInNode = min(255u, mp);
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    if (sm == "sah")         splitMethod = SPLIT_SAH;
    else if (sm == "middle") splitMethod = SPLIT_MIDDLE;
    else if (sm == "equal")  splitMethod = SPLIT_EQUAL_COUNTS;
    else if (sm == "lbvh")   splitMethod = SPLIT_LBVH;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                sm.c_str());
        splitMethod = SPLIT_SAH;
    }

    if (primitives.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build BVH from _primitives_
    PBRT_BVH_STARTED_CONSTRUCTION(this, primitives.size());

    // Recursively build BVH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    vector<Reference<Primitive> > orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root = NULL;
    if (splitMethod == SPLIT_LBVH)
    {
        vector<LBVHPrimitiveInfo> buildData;
        buildData.reserve(primitives.size());
        for (uint32_t i = 0; i < primitives.size(); i++)
        {
            BBox bbox = primitives[i]->WorldBound();
            buildData.push_back(LBVHPrimitiveInfo(i, bbox));
        }
		Timer lbvhTimer;
		lbvhTimer.Start();
        root = buildLBVH(buildArena, buildData, &totalNodes, orderedPrims);
		//printf("Elapsed time to build LBVH: %.2f seconds\n", lbvhTimer.Time());
    }
    else
    {
		//	sah
        vector<BVHPrimitiveInfo> buildData;
        buildData.reserve(primitives.size());
        for (uint32_t i = 0; i < primitives.size(); ++i) {
            BBox bbox = primitives[i]->WorldBound();
            buildData.push_back(BVHPrimitiveInfo(i, bbox));
        }
		Timer sahTimer;
		sahTimer.Start();
        root = recursiveBuild(buildArena, buildData, 0,
                                        primitives.size(), &totalNodes,
                                        orderedPrims);
		//printf("Elapsed time to build BVH: %.2f seconds\n", sahTimer.Time());
    }
    primitives.swap(orderedPrims);
	orderedPrims.clear();
    Info("BVH created with %d nodes for %d primitives (%.2f MB)", totalNodes,
             (int)primitives.size(), float(totalNodes * sizeof(LinearBVHNode))/(1024.f*1024.f));

    // Compute representation of depth-first traversal of BVH tree
    nodes = AllocAligned<LinearBVHNode>(totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearBVHNode;
    uint32_t offset = 0;
    flattenBVHTree(root, &offset);
	Assert(offset == totalNodes);
    PBRT_BVH_FINISHED_CONSTRUCTION(this);
}

BBox BVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : BBox();
}

BVHBuildNode *BVHAccel::buildLBVH(MemoryArena &buildArena,
        vector<LBVHPrimitiveInfo> &buildData, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims)
{
    //  compute the bounding box of the centroid points
    BBox bbox;
    for (uint32_t i = 0; i < primitives.size(); i++)
        bbox = Union(bbox, buildData[i].bounds);
    //  quantize them to integers
	//	morton code
	//Timer mortonCodeTimer;
	//mortonCodeTimer.Start();
    float invXStep, invYStep, invZStep;
    Vector v = bbox.pMax - bbox.pMin;
	uint32_t num = 1024;    
	invXStep = num / v.x;
    invYStep = num / v.y;
    invZStep = num / v.z;
    for (uint32_t i = 0; i < primitives.size(); i++)
    {
        //  fill the member in buildData.mortonCodeX, mortonCodeY, mortonCodeZ
        //  and mortonCode
        Vector p = buildData[i].centroid - bbox.pMin;
        uint32_t x = uint32_t(p.x * invXStep);
        x = (x < num) ? x : num - 1;
        uint32_t y = uint32_t(p.y * invYStep);
        y = (y < num) ? y : num - 1;
        uint32_t z = uint32_t(p.z * invZStep);
        z = (z < num) ? z : num - 1;
        //  compute the morton code
        buildData[i].mortonCode = compMortonCode(x, y, z);
    }
	//printf("Elapsed time to compute morton code: %.2f seconds\n", mortonCodeTimer.Time());
    //  use radix sort to sort them
    radixSortMortonCode(buildData);
	//printf("Elapsed time for radix sort: %.2f seconds\n", myTimer.Time());
    //  top-down approach to build lbvh
    BVHBuildNode *node = recursiveBuildLBVH(buildArena, buildData, 0, primitives.size(), totalNodes, orderedPrims, 29);
	return node;
}

BVHBuildNode* BVHAccel::recursiveBuildLBVH(MemoryArena &buildArena,
        vector<LBVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,
        uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims, int bitLevel)
{
	assert(start != end);
    (*totalNodes)++;
    BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
	//	if nPrimitives is too small
	//	or it's the last bit and nPrimitives is no more than 255
	//	we can build a leaf node
	//	else if nPrimitives is not so small and it is not the last bit
	//	or it is the last bit but nPrimitive is very large  
    if (nPrimitives <= 4 
		|| (bitLevel < 0 && nPrimitives <= 255)) {
        // Create leaf _BVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
	else if (bitLevel < 0)
	{
		uint32_t mid = (start + end) / 2;
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();
        node->InitInterior(dim,
                        recursiveBuildLBVH(buildArena, buildData, start, mid,
                            totalNodes, orderedPrims, -1),
                        recursiveBuildLBVH(buildArena, buildData, mid, end,
                            totalNodes, orderedPrims, -1));
	}
    else
    {
        //  compute dimension
        //  0 for x, 1 for y and 2 for z
        //  find mid here!
        //  the mid satisfies that
        //  the bitLevel-th bit from start to mid - 1 is 0
        //  but from mid to end-1 is 1
        //  note that buildData has been sorted
        //  if mid - start == 0 || end - mid == 0
        //  then we should go to next level until bitLevel == -1
        //  if bitLevel == -1, then we should init a leaf node
        uint32_t mid = 0;
		bool found = false;
        while (bitLevel >= 0 && !found)
        {
            //  binary search
            uint32_t value = ((buildData[start].mortonCode >> (bitLevel+1)) << (bitLevel + 1)) + (1 << bitLevel);
			//  search from mid where buildData[mid].mortonCode >= value and buildData[mid - 1].mortoncode < value

            if (buildData[start].mortonCode >= value || buildData[end - 1].mortonCode < value)
            {
                bitLevel--;
                continue;
            }
            uint32_t p = start, q = end - 1;
            while (true)
            {
				mid = (p + q) / 2;
				if (buildData[mid].mortonCode < value)
                    p = mid + 1;
                else
                    q = mid;
                if (p == q)
                {
                    mid = p;
					found = true;
                    break;
                }
            }
        }
        if (bitLevel < 0)
        {
        	//  build a leaf node
			if (nPrimitives <= 255)
			{
        		uint32_t firstPrimOffset = orderedPrims.size();
        		for (uint32_t i = start; i < end; ++i) {
            		uint32_t primNum = buildData[i].primitiveNumber;
            		orderedPrims.push_back(primitives[primNum]);
        		}
        		node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
        	}
			else
			{
				//	split the node
				uint32_t mid = (start + end) / 2;
        		// Compute bound of primitive centroids, choose split dimension _dim_
        		BBox centroidBounds;
        		for (uint32_t i = start; i < end; ++i)
            		centroidBounds = Union(centroidBounds, buildData[i].centroid);
        		int dim = centroidBounds.MaximumExtent();
        		node->InitInterior(dim,
                        recursiveBuildLBVH(buildArena, buildData, start, mid,
                            totalNodes, orderedPrims, -1),
                        recursiveBuildLBVH(buildArena, buildData, mid, end,
                            totalNodes, orderedPrims, -1));
			}
		}
        else
        {
            node->InitInterior(bitLevel % 3,
                        recursiveBuildLBVH(buildArena, buildData, start, mid,
                            totalNodes, orderedPrims, bitLevel - 1),
                        recursiveBuildLBVH(buildArena, buildData, mid, end,
                            totalNodes, orderedPrims, bitLevel - 1));
        }
    }
    return node;
}

BVHBuildNode *BVHAccel::recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start,
        uint32_t end, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims) {
    Assert(start != end);
    (*totalNodes)++;
    BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
    // Compute bounds of all primitives in BVH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // Create leaf _BVHBuildNode_
            uint32_t firstPrimOffset = orderedPrims.size();
            for (uint32_t i = start; i < end; ++i) {
                uint32_t primNum = buildData[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
            return node;
        }

        // Partition primitives based on _splitMethod_
        switch (splitMethod) {
        case SPLIT_MIDDLE: {
            // Partition primitives through node's midpoint
            float pmid = .5f * (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]);
            BVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
                                                      &buildData[end-1]+1,
                                                      CompareToMid(dim, pmid));
            mid = midPtr - &buildData[0];
            if (mid != start && mid != end)
                // for lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall through
                // to SPLIT_EQUAL_COUNTS
                break;
        }
        case SPLIT_EQUAL_COUNTS: {
            // Partition primitives into equally-sized subsets
            mid = (start + end) / 2;
            std::nth_element(&buildData[start], &buildData[mid],
                             &buildData[end-1]+1, ComparePoints(dim));
            break;
        }
        case SPLIT_SAH: default: {
            // Partition primitives using approximate SAH
            if (nPrimitives <= 4) {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                                 &buildData[end-1]+1, ComparePoints(dim));
            }
            else {
                // Allocate _BucketInfo_ for SAH partition buckets
                const int nBuckets = 12;
                struct BucketInfo {
                    BucketInfo() { count = 0; }
                    int count;
                    BBox bounds;
                };
                BucketInfo buckets[nBuckets];

                // Initialize _BucketInfo_ for SAH partition buckets
                for (uint32_t i = start; i < end; ++i) {
                    int b = nBuckets *
                        ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                         (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                    if (b == nBuckets) b = nBuckets-1;
                    Assert(b >= 0 && b < nBuckets);
                    buckets[b].count++;
                    buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
                }

                // Compute costs for splitting after each bucket
                float cost[nBuckets-1];
                for (int i = 0; i < nBuckets-1; ++i) {
                    BBox b0, b1;
                    int count0 = 0, count1 = 0;
                    for (int j = 0; j <= i; ++j) {
                        b0 = Union(b0, buckets[j].bounds);
                        count0 += buckets[j].count;
                    }
                    for (int j = i+1; j < nBuckets; ++j) {
                        b1 = Union(b1, buckets[j].bounds);
                        count1 += buckets[j].count;
                    }
                    cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
                              bbox.SurfaceArea();
                }

                // Find bucket to split at that minimizes SAH metric
                float minCost = cost[0];
                uint32_t minCostSplit = 0;
                for (int i = 1; i < nBuckets-1; ++i) {
                    if (cost[i] < minCost) {
                        minCost = cost[i];
                        minCostSplit = i;
                    }
                }

                // Either create leaf or split primitives at selected SAH bucket
                if (nPrimitives > maxPrimsInNode ||
                    minCost < nPrimitives) {
                    BVHPrimitiveInfo *pmid = std::partition(&buildData[start],
                        &buildData[end-1]+1,
                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                    mid = pmid - &buildData[0];
                }

                else {
                    // Create leaf _BVHBuildNode_
                    uint32_t firstPrimOffset = orderedPrims.size();
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t primNum = buildData[i].primitiveNumber;
                        orderedPrims.push_back(primitives[primNum]);
                    }
                    node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                    return node;
                }
            }
            break;
        }
        }
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
    return node;
}


uint32_t BVHAccel::flattenBVHTree(BVHBuildNode *node, uint32_t *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    uint32_t myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1],
                                                       offset);
    }
    return myOffset;
}


BVHAccel::~BVHAccel() {
    FreeAligned(nodes);
}


bool BVHAccel::Intersect(const Ray &ray, Intersection *isect) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTION_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    // Follow ray through BVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        // Check ray against BVH node
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                for (uint32_t i = 0; i < node->nPrimitives; ++i)
                {
                    PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
                    {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        hit = true;
                    }
                    else {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                   }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                // Put far BVH node on _todo_ stack, advance to near node
                PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}


bool BVHAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTIONP_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    uint32_t todo[64];
    uint32_t todoOffset = 0, nodeNum = 0;
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                  for (uint32_t i = 0; i < node->nPrimitives; ++i) {
                    PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        return true;
                    }
                else {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   /// second child first
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTIONP_FINISHED();
    return false;
}


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    string splitMethod = ps.FindOneString("splitmethod", "sah");
    uint32_t maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return new BVHAccel(prims, maxPrimsInNode, splitMethod);
}


