#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#define MAX_GJK_ITERS 30

using Vec3 = Eigen::Vector3f;
using Real = Eigen::Vector3f::Scalar;

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

#ifdef BUILD_TESTS
void handleSimplexLine(Simplex &simplex, Vec3 &closest);
void handleSimplexTri(Simplex &simplex, Vec3 &closest);
void handleSimplexTetra(Simplex &simplex, Vec3 &closest);
void gjk(const Polytope &A, const Polytope &B, Simplex &simplex, Real &distance, const Vec3 *initDirection = nullptr);
#endif