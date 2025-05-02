#include <vector>
#include <string>
#include "vector.h"
#include "bounding.h"

#pragma once

class BVHNode {
    public:
        BoundingBox bbox;
        BVHNode* childLeft = nullptr; 
        BVHNode* childRight = nullptr;
        int start;
        int end; 
    
        BVHNode() = default;
        ~BVHNode() {
            delete childLeft;
            delete childRight;
        }
    
        bool isLeaf() const {
            return childLeft == nullptr && childRight == nullptr;
        }
    };

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false);
    int vtxi, vtxj, vtxk;
    int uvi, uvj, uvk; 
    int ni, nj, nk;
    int group;
};



class TriangleMesh : public Geometry {
    public:
        TriangleMesh(const Vector& albedo, MaterialType material = NORMAL) : Geometry(albedo, material) {}
        ~TriangleMesh();
    
        bool intersect(const Ray& ray, double& t, Vector& P, Vector& N) const override;

        void readOBJ(const char* obj);
        void computeBoundingBox(); 

        void buildBVH();
        void buildBVHRecursive(BVHNode* node, int start, int end);
        bool intersectBVH(const BVHNode* node, const Ray& ray, double& t, Vector& P, Vector& N) const; 

        std::vector<TriangleIndices> indices;
        std::vector<Vector> vertices;
        std::vector<Vector> normals;
        std::vector<Vector> uvs;
        std::vector<Vector> vertexcolors;
    
        BoundingBox box; 
        BVHNode* bvhRoot = nullptr;
    };

    