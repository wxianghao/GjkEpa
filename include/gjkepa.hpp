#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#define MAX_GJK_ITERS 30
#define MAX_EPA_VERTS 30

using Vec3 = Eigen::Vector3f;
using Real = Eigen::Vector3f::Scalar;
// using Vec3 = Eigen::Vector3d;
// using Real = Eigen::Vector3d::Scalar;

struct Polytope
{
    std::vector<Vec3> verts;
    Polytope() = default;
    Polytope(std::initializer_list<Vec3> verts_list)
        : verts{verts_list}
    {
    }

    Polytope &operator=(std::initializer_list<Vec3> verts_list)
    {
        verts.clear();
        verts.assign(verts_list.begin(), verts_list.end());
        return *this;
    }
    ~Polytope() = default;
};

struct Simplex
{
    std::vector<Vec3> verts;
    Simplex() { verts.reserve(4); }
    Simplex(std::initializer_list<Vec3> verts_list)
        : verts{verts_list}
    {
    }
    ~Simplex() = default;
    Simplex &operator=(std::initializer_list<Vec3> verts_list)
    {
        verts.clear();
        verts.assign(verts_list.begin(), verts_list.end());
        return *this;
    }
};

struct EPATriangle
{
    std::array<int, 3> verts;     // Indices of vertices
    std::array<int, 3> adjacents; // Indices of adjacent faces
    Vec3               normal;    // Unit normal of the face pointing the outside
    Real               distance;  // Signed distance from the face to the origin

    void rebuild(const std::vector<Vec3> &vertices)
    {
        const Vec3 &p0  = vertices[verts[0]];
        const Vec3 &p1  = vertices[verts[1]];
        const Vec3 &p2  = vertices[verts[2]];
        Vec3        e01 = p1 - p0;
        Vec3        e02 = p2 - p0;
        Vec3        e12 = p2 - p1;

        if (e01.dot(e01) < e02.dot(e02)) {
            this->normal = e01.cross(e12);
        }
        else {
            this->normal = e02.cross(e12);
        }
        this->normal.normalize();
        this->distance = this->normal.dot(p0);
    }
};


void gjk(const Polytope &A, const Polytope &B, Simplex &simplex, Real &distance, const Vec3 *initDirection = nullptr);

#ifdef BUILD_TESTS
void handleSimplexLine(Simplex &simplex, Vec3 &closest);
void handleSimplexTri(Simplex &simplex, Vec3 &closest);
void handleSimplexTetra(Simplex &simplex, Vec3 &closest);
void epaStart4(const std::vector<Vec3> &vertices, std::vector<EPATriangle> &triangles);
bool replaceTriangle(int apexIdx, int faceIdx, const std::vector<Vec3> &vertices, std::vector<EPATriangle> &triangles);
#endif