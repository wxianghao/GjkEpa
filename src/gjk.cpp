#include <cassert>
#include <cstdio>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>

#include "common.hpp"
#include "gjkepa.hpp"

Vec3 support_gjk(const Polytope &A, const Vec3 &direction)
{
    int  idx = -1;
    Real len = Eigen::NumTraits<Real>::lowest();

    for (int i = 0; i < A.verts.size(); ++i) {
        Real t = A.verts[i].dot(direction);
        if (t > len) {
            idx = i;
            len = t;
        }
    }

    return A.verts[idx];
}

void handleSimplexLine(Simplex &simplex, Vec3 &closest)
{
    const Vec3 &A = simplex.verts[0];
    const Vec3 &B = simplex.verts[1];

    Vec3 AB = B - A;
    Real pB = B.dot(AB);

    // O is in the B region
    if (pB < 0) {
        closest = B;
        simplex = {B};
        return;
    }

    closest = B - AB * pB / AB.dot(AB);
}

void handleSimplexTri(Simplex &simplex, Vec3 &closest)
{
    const Vec3 &A = simplex.verts[0];
    const Vec3 &B = simplex.verts[1];
    const Vec3 &C = simplex.verts[2];

    Vec3 AC = C - A;
    Vec3 BC = C - B;
    Vec3 AB = B - A;

    // Choose the smaller edge to compute the plane ABC's normal
    Real dAC     = AC.dot(AC);
    Real dBC     = BC.dot(BC);
    Vec3 shorter = dAC < dBC ? AC : BC;
    Vec3 normal  = AB.cross(shorter);

    // Origin's position relative to the edges
    bool outsideAC = normal.cross(AC).dot(C) <= 0;
    bool outsideBC = BC.cross(normal).dot(C) <= 0;

    // Case 1: origin's projection falls inside ABC
    if (!outsideAC && !outsideBC) {
        Real x  = normal.dot(C);
        closest = normal * (x / normal.dot(normal));
        // Ensure the normal always points the origin
        if (x < 0) {
            // Equiv. to:
            // simplex = {A, B, C};
        }
        else {
            simplex = {B, A, C};
        }
        return;
    }

    // Case 2: origin projects outside ABC, but falls on region of either BC or AC
    Real pBC = BC.dot(C);
    Real pAC = AC.dot(C);
    if (outsideBC && pBC > 0) {
        closest = C - BC * (pBC / dBC);
        simplex = {B, C};
        return;
    }
    if (outsideAC && pAC > 0) {
        closest = C - AC * (pAC / dAC);
        simplex = {A, C};
        return;
    }

    // Case 3: origin projects on the apex region
    closest = C;
    simplex = {C};
}

void handleSimplexTetra(Simplex &simplex, Vec3 &closest)
{
    const Vec3 &s0   = simplex.verts[0];
    const Vec3 &s1   = simplex.verts[1];
    const Vec3 &s2   = simplex.verts[2];
    const Vec3 &apex = simplex.verts[3];

    // Edges
    Vec3 edge0 = apex - s0;
    Vec3 edge1 = apex - s1;
    Vec3 edge2 = apex - s2;

    // Squared norm of edges pointing the apex
    Real l0 = edge0.dot(edge0);
    Real l1 = edge1.dot(edge1);
    Real l2 = edge2.dot(edge2);

    // Normals of each face
    Vec3 shorter0 = l0 < l1 ? edge0 : edge1;
    Vec3 shorter1 = l1 < l2 ? edge1 : edge2;
    Vec3 shorter2 = l2 < l0 ? edge2 : edge0;
    Vec3 normal0  = (s1 - s0).cross(shorter0);
    Vec3 normal1  = (s2 - s1).cross(shorter1);
    Vec3 normal2  = (s0 - s2).cross(shorter2);

    // Origin's SDF from each face
    Real a0 = normal0.dot(apex);
    Real a1 = normal1.dot(apex);
    Real a2 = normal2.dot(apex);

    // Case 1: origin is inside the tetrahedron
    if (a0 > 0 && a1 > 0 && a2 > 0) {
        closest = Vec3{0, 0, 0};
        return;
    }

    // Projection onto edges
    Real p0 = edge0.dot(apex);
    Real p1 = edge1.dot(apex);
    Real p2 = edge2.dot(apex);

    // Case 2: origin is in the apex region
    if (p0 <= 0 && p1 <= 0 && p2 <= 0) {
        closest = apex;
        simplex = {apex};
        return;
    }


    // Case 3: origin is in the edge region
    Real u0 = normal0.cross(edge0).dot(apex);
    Real u1 = normal1.cross(edge1).dot(apex);
    Real u2 = normal2.cross(edge2).dot(apex);
    Real w0 = normal0.cross(edge1).dot(apex);
    Real w1 = normal1.cross(edge2).dot(apex);
    Real w2 = normal2.cross(edge0).dot(apex);
    if (u0 <= 0 && w2 > 0 && p0 >= 0) {
        closest = apex - edge0 * (p0 / l0);
        simplex = {apex, s0};
        return;
    }
    if (u1 <= 0 && w0 > 0 && p1 >= 0) {
        closest = apex - edge1 * (p1 / l1);
        simplex = {apex, s1};
        return;
    }
    if (u2 <= 0 && w1 > 0 && p2 >= 0) {
        closest = apex - edge2 * (p2 / l2);
        simplex = {apex, s2};
        return;
    }

    // Case 4: origin is in the face region
    if (u0 > 0 && w0 <= 0 && a0 <= 0) {
        closest = normal0 * (a0 / normal0.dot(normal0));
        simplex = {apex, s0, s1};
        return;
    }

    if (u1 > 1 && w1 <= 0 && a1 <= 0) {
        closest = normal1 * (a1 / normal1.dot(normal1));
        simplex = {apex, s1, s2};
        return;
    }

    if (u2 > 2 && w2 <= 0 && a2 <= 0) {
        closest = normal2 * (a2 / normal2.dot(normal2));
        simplex = {apex, s2, s0};
        return;
    }

    assert(false);
}

void handleSimplex(Simplex &simplex, Vec3 &closest)
{
    switch (simplex.verts.size()) {
    case 1:
        closest = -closest;
        return;
    case 2:
        handleSimplexLine(simplex, closest);
        return;
    case 3:
        handleSimplexTri(simplex, closest);
    case 4:
        handleSimplexTetra(simplex, closest);
    default:
        throw std::invalid_argument("Illegal simplex");
    }
}


void gjk(const Polytope &A, const Polytope &B, Simplex &simplex, Real &distance, const Vec3 *initDirection = nullptr)
{
    Vec3 closest;
    int  k       = 0;
    Real normMax = 0.f;

    // Initialize direction
    if (initDirection == nullptr) {
        closest = A.verts[0] - B.verts[0];
    }
    else {
        closest = support_gjk(A, *initDirection) - support_gjk(B, -*initDirection);
    }

    simplex.verts.clear();
    simplex.verts.push_back(closest);

    do {
        // Exit condition: the nearest point approaches the origin
        if (closest.dot(closest) < abstol * abstol) {
            break;
        }

        // Calculate support point
        Vec3 support = support_gjk(A, -closest) - support_gjk(B, closest);

        // Exit condition: the new simplex vertex cannot move further
        if (closest.dot(closest) - support.dot(closest) < abstol) {
            break;
        }

        // Add new vertex to the simplex
        simplex.verts.push_back(support);

        // Determine the new closest point
        handleSimplex(simplex, closest);

        // Exit condition: the nearest point approaches the origin in terms of very large simplex shape
        for (const auto &v : simplex.verts) {
            normMax = Eigen::numext::maxi(normMax, v.dot(v));
        }
        if (closest.dot(closest) <= abstol * abstol * normMax) {
            break;
        }

        ++k;
    } while (simplex.verts.size() == 4 || k < MAX_GJK_ITERS);

    distance = closest.norm();
}