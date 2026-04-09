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

    Polytope(std::initializer_list<Vec3> verts_list)
        : verts{verts_list}
    {
        if (verts_list.size() < 4) {
            verts.clear();
            throw std::invalid_argument("Polytope requires at least 4 vertices, got "
                                        + std::to_string(verts_list.size()));
        }
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
        if (verts_list.size() > 4) {
            verts.clear();
            throw std::invalid_argument("Simplex requires at most 4 vertices, got "
                                        + std::to_string(verts_list.size()));
        }
    }
    ~Simplex() = default;
    Simplex &operator=(std::initializer_list<Vec3> verts_list)
    {
        if (verts_list.size() > 4) {
            throw std::invalid_argument("Simplex requires at most 4 vertices, got "
                                        + std::to_string(verts_list.size()));
        }
        verts.clear();
        verts.assign(verts_list.begin(), verts_list.end());
        return *this;
    }
};

#ifdef BUILD_TESTS
void handleSimplexLine(Simplex &simplex, Vec3 &closest);
void handleSimplexTri(Simplex &simplex, Vec3 &closest);
void handleSimplexTetra(Simplex &simplex, Vec3 &closest);
#endif