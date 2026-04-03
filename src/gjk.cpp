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

void handleSimplex(Simplex &simplex, Vec3 &closest)
{
    switch (simplex.verts.size()) {
    case 1:
        closest = -closest;
        return;
    case 2:
        handleSimplexLine(simplex, closest);
        return;
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
    } while (k < MAX_GJK_ITERS);

    distance = closest.norm();
}