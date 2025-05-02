
#include <cstdio>
#include <algorithm>
#include <limits>
#include <list>
#include <cstring>

#include "bounding.h"
#include "vector.h"
#include "trianglemesh.h"



TriangleIndices::TriangleIndices(int vtxi, int vtxj, int vtxk, int ni, int nj, int nk, int uvi, int uvj, int uvk, int group, bool added)
    : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {}

TriangleMesh::~TriangleMesh() {}


void TriangleMesh::buildBVH() {
    bvhRoot = new BVHNode();

    buildBVHRecursive(bvhRoot, 0, indices.size());
}


void TriangleMesh::buildBVHRecursive(BVHNode* node, int start, int end) {
    if (start >= end) {
        return;
    }
    
    node->start = start;
    node->end = end;

    node->bbox = BoundingBox();
    for (int i = start; i < end; i++) {
        const auto& tri = indices[i];
        if (tri.vtxi >= 0 && tri.vtxi < vertices.size() &&
            tri.vtxj >= 0 && tri.vtxj < vertices.size() &&
            tri.vtxk >= 0 && tri.vtxk < vertices.size()) {
            node->bbox.next(vertices[tri.vtxi]);
            node->bbox.next(vertices[tri.vtxj]);
            node->bbox.next(vertices[tri.vtxk]);
        }
    }

    if (end - start <= 4) {
        return;
    }

    Vector diag = node->bbox.Bmax - node->bbox.Bmin;
    int axis = 0;
    if (diag[1] > diag[0]) axis = 1;
    if (diag[2] > diag[axis]) axis = 2;

    double mid = 0.5 * (node->bbox.Bmin[axis] + node->bbox.Bmax[axis]);

    int pivot = std::partition(indices.begin() + start, indices.begin() + end,
        [&](const TriangleIndices& tri) {
            if (tri.vtxi >= 0 && tri.vtxi < vertices.size() &&
                tri.vtxj >= 0 && tri.vtxj < vertices.size() &&
                tri.vtxk >= 0 && tri.vtxk < vertices.size()) {
                Vector centroid = (vertices[tri.vtxi] + vertices[tri.vtxj] + vertices[tri.vtxk]) / 3.0;
                return centroid[axis] < mid;
            }
            return false;
        }) - indices.begin();

    if (pivot == start || pivot == end) {
        return;
    }

    node->childLeft = new BVHNode();
    node->childRight = new BVHNode();
    buildBVHRecursive(node->childLeft, start, pivot);
    buildBVHRecursive(node->childRight, pivot, end);
}


bool TriangleMesh::intersect(const Ray& ray, double& t, Vector& P, Vector& N) const {
    return intersectBVH(bvhRoot, ray, t, P, N);
}


bool TriangleMesh::intersectBVH(const BVHNode* node, const Ray& ray, double& t, Vector& P, Vector& N) const {
    if (!node) {
        return false;
    }

    double tmin, tmax;
    if (!node->bbox.intersect(ray, tmin, tmax)) {
        return false;
    }

    if (node->isLeaf()) {
        bool hit = false;
        double tMin = 1e30;

        int validStart = std::max(0, node->start);
        int validEnd = std::min((int)indices.size(), node->end);
        
        for (int i = validStart; i < validEnd; i++) {
            const auto& tri = indices[i];

            if (tri.vtxi < 0 || tri.vtxi >= vertices.size() ||
                tri.vtxj < 0 || tri.vtxj >= vertices.size() ||
                tri.vtxk < 0 || tri.vtxk >= vertices.size()) {
                continue;
            }
            
            const Vector& A = vertices[tri.vtxi];
            const Vector& B = vertices[tri.vtxj];
            const Vector& C = vertices[tri.vtxk];

            Vector e1 = B - A;
            Vector e2 = C - A;
            Vector pvec = cross(ray.direction, e2);
            double det = dot(e1, pvec);

            if (fabs(det) < 1e-8) continue;

            double invDet = 1.0 / det;
            Vector tvec = ray.origin - A;
            double u = dot(tvec, pvec) * invDet;
            if (u < 0 || u > 1) continue;

            Vector qvec = cross(tvec, e1);
            double v = dot(ray.direction, qvec) * invDet;
            if (v < 0 || u + v > 1) continue;

            double tTemp = dot(e2, qvec) * invDet;
            if (tTemp < 0 || tTemp >= tMin) continue;

            hit = true;
            tMin = tTemp;
            P = ray.origin + tTemp * ray.direction;
            N = cross(e1, e2);
            N.normalize();

            if (dot(N, ray.direction) > 0) {
                N = -N;
            }
        }

        t = tMin;
        return hit;
    }

    double t_left = 1e30;
    Vector P_left, N_left;
    bool hit_left = node->childLeft && intersectBVH(node->childLeft, ray, t_left, P_left, N_left);
    
    double t_right = 1e30;
    Vector P_right, N_right;
    bool hit_right = node->childRight && intersectBVH(node->childRight, ray, t_right, P_right, N_right);
    
    if (hit_left && hit_right) {
        if (t_left < t_right) {
            t = t_left;
            P = P_left;
            N = N_left;
        } else {
            t = t_right;
            P = P_right;
            N = N_right;
        }
        return true;
    } else if (hit_left) {
        t = t_left;
        P = P_left;
        N = N_left;
        return true;
    } else if (hit_right) {
        t = t_right;
        P = P_right;
        N = N_right;
        return true;
    }
    
    return false;
}

void TriangleMesh::readOBJ(const char* obj) {

    char matfile[255];
    char grp[255];

    FILE* f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);

            } else {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char* consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                indices.push_back(t);
            } else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true) {
                if (consumedline[0] == '\n') break;
                if (consumedline[0] == '\0') break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                    if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                    if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                    if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                    if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                    if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                    if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                } else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    } else {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            } else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }

        }

    }
    fclose(f);

}
