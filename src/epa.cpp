#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common.hpp"
#include "gjkepa.hpp"


bool replaceTriangle(int apexIdx, int faceIdx, const std::vector<Vec3> &vertices, std::vector<EPATriangle> &triangles)
{
    // Collect faces which can be seen by the apex
    const Vec3             &apex = vertices[apexIdx];
    std::unordered_set<int> seens;
    for (int i = 0; i < triangles.size(); ++i) {
        const EPATriangle &tri = triangles[i];
        if (apex.dot(tri.normal) - tri.distance >= epsRel) {
            seens.insert(i);
        }
    }

    //! Skip degenerate check for fast implementation

    // Build the silhouette
    std::vector<int> next(vertices.size() - 1);
    std::vector<int> faces(vertices.size() - 1);
    std::vector<int> silVerts;
    int              loopStart = -1;
    // Build linked list of silhouette nodes
    for (int i = 0; i < triangles.size(); ++i) {
        if (seens.count(i)) {
            continue;
        }
        const EPATriangle &tri  = triangles[i];
        const auto        &v    = tri.verts;
        const auto        &adjs = tri.adjacents;

        if (seens.count(adjs[0])) {
            next[v[2]]  = v[1];
            faces[v[2]] = i;
            loopStart   = v[2];
        }

        if (seens.count(adjs[1])) {
            next[v[0]]  = v[2];
            faces[v[0]] = i;
            loopStart   = v[0];
        }

        if (seens.count(adjs[2])) {
            next[v[1]]  = v[0];
            faces[v[1]] = i;
            loopStart   = v[1];
        }
    }
    // Extract the list into an array
    assert(loopStart != -1);
    silVerts.push_back(loopStart);
    for (int vIdx = next[loopStart]; vIdx != loopStart; vIdx = next[vIdx]) {
        silVerts.push_back(vIdx);
    }

    // Build new triangle faces
    int newStartIdx = triangles.size();
    for (int i = 0; i < silVerts.size(); ++i) {
        int from = silVerts[i], to = silVerts[(i + 1) % silVerts.size()];
        // Build triangle
        EPATriangle tri;
        tri.verts     = {apexIdx, from, to};
        tri.adjacents = {faces[from], -1, -1};
        // Compute new normal and distance
        tri.rebuild(vertices);

        triangles.push_back(tri);
        faces[from] = newStartIdx + i;
    }
    // Find adjacents of new triangles
    for (int i = 0; i < silVerts.size(); ++i) {
        EPATriangle &tri = triangles[newStartIdx + i];
        tri.adjacents[1] = newStartIdx + (i + 1) % silVerts.size();
        tri.adjacents[2] = newStartIdx + (i - 1 + silVerts.size()) % silVerts.size();
    }

    // Find adjacents of old triangles
    for (int i = 0; i < newStartIdx; ++i) {
        if (seens.count(i)) {
            continue;
        }
        EPATriangle &tri = triangles[i];
        for (int adj = 0; adj < 3; ++adj) {
            if (seens.count(tri.adjacents[adj])) {
                tri.adjacents[adj] = faces[tri.verts[(adj + 2) % 3]];
            }
        }
    }

    // Eliminate seen triangles
    std::unordered_map<int, int> newIdxMap;
    for (int slow = 0, fast = 0; fast < triangles.size(); ++fast) {
        if (seens.count(fast)) {
            continue;
        }
        triangles[slow] = triangles[fast];
        newIdxMap[fast] = slow;
        ++slow;
    }
    triangles.resize(triangles.size() - seens.size());
    for (auto &tri : triangles) {
        for (int i = 0; i < 3; ++i) {
            int adj          = tri.adjacents[i];
            tri.adjacents[i] = newIdxMap.at(adj);
        }
    }

    return true;
}

int supportIndex(const Polytope &poly, const Vec3 &dir)
{
    int  bestIdx = -1;
    Real bestVal = -inf;

    for (int i = 0; i < poly.verts.size(); ++i) {
        const auto &v   = poly.verts[i];
        Real        val = v.dot(dir);
        if (bestVal < val) {
            bestVal = val;
            bestIdx = i;
        }
    }

    return bestIdx;
}

void epaStart4(const std::vector<Vec3> &vertices, std::vector<EPATriangle> &triangles)
{
    triangles.clear();

    constexpr int triVerts[4][3] = {{0, 1, 2}, {3, 1, 0}, {3, 2, 1}, {3, 0, 2}};
    constexpr int triAdjs[4][3]  = {{2, 3, 1}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

    triangles.resize(4);
    for (int i = 0; i < 4; ++i) {
        auto &tri     = triangles[i];
        tri.verts     = {triVerts[i][0], triVerts[i][1], triVerts[i][2]};
        tri.adjacents = {triAdjs[i][0], triAdjs[i][1], triAdjs[i][2]};
        tri.rebuild(vertices);

        // Check whether to flip
        if (tri.distance < 0) {
            std::swap(tri.verts[0], tri.verts[1]);
            std::swap(tri.adjacents[0], tri.adjacents[1]);
            tri.rebuild(vertices);
        }
    }
}

Vec3 newPoint(const Polytope &A, const Polytope &B, std::vector<Vec3> &vertices, const Vec3 &dir, const Vec3 &ref)
{
    Vec3 vp = A.verts[supportIndex(A, dir)] - B.verts[supportIndex(B, -dir)];
    Vec3 vq = A.verts[supportIndex(A, -dir)] - B.verts[supportIndex(B, dir)];
    if (abs((vp - ref).dot(dir)) > abs((vq - ref).dot(dir))) {
        return vp;
    }
    else {
        return vq;
    }
}

// void epaStart123(const Polytope &A, const Polytope &B, std::vector<Vec3> &vertices, std::vector<EPATriangle>
// &triangles)
// {
//     if (vertices.size() == 1) {
//         Vec3 v = {1, 0, 0};
//         vertices.push_back(newPoint(A, B, vertices, v, vertices[0]));
//     }

//     if (vertices.size() == 2) {
//         const Vec3 &p0 = vertices[0];
//         const Vec3 &p1 = vertices[1];
//         Vec3        v  = p1 - p0;

//         if (abs(v.x()) < abs(v.y())) {
//             v = Vec3{1, 0, 0}.cross(v);
//         }
//         else {
//             v = Vec3{0, 1, 0}.cross(v);
//         }

//         vertices.push_back(newPoint(A, B, vertices, v, vertices[1]));
//     }

//     const Vec3 &p0  = vertices[0];
//     const Vec3 &p1  = vertices[1];
//     const Vec3 &p2  = vertices[2];
//     Vec3        e01 = p1 - p0;
//     Vec3        e02 = p2 - p0;
//     Vec3        e12 = p2 - p1;

//     Vec3 normal;
//     if (e01.dot(e01) < e02.dot(e02)) {
//         normal = e01.cross(e12);
//     }
//     else {
//         normal = e02.cross(e12);
//     }

//     Real distance = normal.dot(p0);
// }


// //!: should use binary operations in GPU implementation
// bool replaceTriangle(Vec3 &apex, int apexIdx, std::vector<EpaTriangle> &triangles, int closestIdx)
// {

//     // Calculate whether apex is on the outside of each triangle face
//     int               nThreadtens = 0;
//     std::vector<bool> threatened(triangles.size());
//     for (int i = 0; i < triangles.size(); ++i) {
//         threatened[i] = triangles[i].normal.dot(apex) - triangles[i].distance >= epsRel;
//         nThreadtens += threatened[i];
//     }

//     // Degenerate if illegal cases happen
//     if (nThreadtens == 0 || nThreadtens == triangles.size() || !threatened[closestIdx]) {
//         return false;
//     }

//     // Flood fill the threatened faces. The set of triangles to kill = triangles
//     // reachable from closestIdx by walking adjacency through threatened faces.
//     std::vector<bool> flooded(triangles.size(), false);
//     std::queue<int>   floodQ;
//     floodQ.push(closestIdx);
//     flooded[closestIdx] = true;
//     while (!floodQ.empty()) {
//         int idx = floodQ.front();
//         floodQ.pop();
//         for (int adj : triangles[idx].adj) {
//             if (!flooded[adj] && threatened[adj]) {
//                 floodQ.push(adj);
//                 flooded[adj] = true;
//             }
//         }
//     }

//     int loopStart = -1;

//     // Build a list of silhouette
//     std::vector<int> next(MAX_EPA_VERTS, -1);
//     std::vector<int> edges(MAX_EPA_VERTS, -1);
//     int              nSilhouette = 0;
//     for (int i = 0; i < (int)triangles.size(); ++i) {
//         if (flooded[i])
//             continue;

//         const auto &v   = triangles[i].v;
//         const auto &adj = triangles[i].adj;

//         bool d0 = flooded[adj[0]];
//         bool d1 = flooded[adj[1]];
//         bool d2 = flooded[adj[2]];

//         if (d0) {
//             next[v[2]]  = v[1];
//             edges[v[2]] = i;
//             ++nSilhouette;
//             loopStart = v[2];
//         }
//         if (d1) {
//             next[v[0]]  = v[2];
//             edges[v[0]] = i;
//             ++nSilhouette;
//             loopStart = v[0];
//         }
//         if (d2) {
//             next[v[1]]  = v[0];
//             edges[v[1]] = i;
//             ++nSilhouette;
//             loopStart = v[1];
//         }
//     }

//     if (nSilhouette == 0) {
//         return false;
//     }

//     // Build new triangle faces with those silhouette edges and the apex
//     std::vector<int> loop;
//     loop.push_back(loopStart);
//     for (int vertex = next[loopStart]; vertex != loopStart; vertex = next[vertex]) {
//         loop.push_back(vertex);
//     }

//     if (loop.size() != nSilhouette) {
//         return false;
//     }

//     for (int i = 0; i < loop.size(); ++i) {
//         EpaTriangle t;
//         int         v1 = loop[i], v2 = loop[(i + 1) % loop.size()];
//         t.v = {apexIdx, v1, v2};

//         int oldAdj = edges[v1];
//         t.adj[0]   = oldAdj;
//     }


//     return true;
// }


// // Initialize epa triangles from simplex of 4 vertices
// static void epaStart4(const std::vector<Vec3> &verts, std::vector<EpaTriangle> &triangles)
// {
//     // Build faces
//     constexpr int faceVerts[4][3] = {
//         {0, 1, 2},
//         {3, 1, 0},
//         {3, 2, 1},
//         {3, 0, 2},
//     };

//     // Build adjacent faces
//     constexpr int faceAdj[4][3] = {
//         {2, 3, 1},
//         {0, 3, 2},
//         {0, 1, 3},
//         {0, 2, 1},
//     };

//     // All normals should point outside
//     Vec3 e01  = verts[1] - verts[0];
//     Vec3 e02  = verts[2] - verts[0];
//     Vec3 e03  = verts[3] - verts[0];
//     bool flip = e01.cross(e02).dot(e03) > 0;

//     triangles.resize(4);
//     for (int i = 0; i < 4; i++) {
//         auto &tri = triangles[i];
//         tri.v     = {faceVerts[i][0], faceVerts[i][1], faceVerts[i][2]};
//         tri.adj   = {faceAdj[i][0], faceAdj[i][1], faceAdj[i][2]};

//         if (flip) {
//             std::swap(tri.v[1], tri.v[2]);
//             std::swap(tri.adj[1], tri.adj[2]);
//         }

//         Vec3 p0 = verts[tri.v[0]], p1 = verts[tri.v[1]], p2 = verts[tri.v[2]];
//         Vec3 edge01 = p1 - p0, edge02 = p2 - p0, edge12 = p2 - p1;

//         // Compute each face's normal
//         Vec3 n;
//         if (edge01.squaredNorm() < edge02.squaredNorm())
//             n = edge01.cross(edge12);
//         else
//             n = edge02.cross(edge12);

//         tri.normal   = n.normalized();
//         tri.distance = tri.normal.dot(p0);
//     }
// }

void epa(const Polytope &A, const Polytope &B, Simplex &simplex)
{
    std::vector<EPATriangle> triangles;
    std::vector<Vec3>        vertices(simplex.verts.begin(), simplex.verts.end());
    if (simplex.verts.size() == 4) {
        epaStart4(vertices, triangles);
    }
    else {
        throw std::runtime_error("Not implemented!");
    }

    Real bestDist      = inf;
    int  bestFace      = -1;
    bool nondegenerate = true;
    while (nondegenerate && vertices.size() < MAX_EPA_VERTS) {
        // Search for the closest face to the origin
        int  closestFace = 0;
        Real distance    = triangles[0].distance;
        for (size_t i = 1; i < triangles.size(); ++i) {
            if (triangles[i].distance < distance) {
                distance    = triangles[i].distance;
                closestFace = i;
            }
        }

        // Search for a new support point
        Vec3 v    = triangles[closestFace].normal;
        int  idxA = supportIndex(A, v), idxB = supportIndex(B, -v);
        Vec3 apex = A.verts[idxA] - B.verts[idxB];

        Real apexDist = apex.dot(v);
        if (abs(apexDist - distance) <= epsRel) {
            bestFace = closestFace;
            bestDist = distance;
            break;
        }

        // if (distance > bestDist) {
        //     break;
        // }

        if (apexDist < bestDist)
        {
            
        }
    }


    // Real lowerBound = -inf;
    // Real upperBound = inf;
    // int  nvertices;

    // std::vector<EPATriangle> triangles;
    // int                      bestFaceIdx = -1;

    // if (simplex.verts.size() == 4) {
    //     epaStart4(simplex.verts, triangles);
    //     nvertices = 4;
    // }
    // else {
    //     throw std::runtime_error("Not implemented!");
    // }

    // bool nondegenerate = true;
    // while (nondegenerate && nvertices < MAX_EPA_VERTS) {
    //     /* Find closest face */
    //     int  closestIdx = 0;
    //     Real distance   = triangles[0].distance;
    //     for (int i = 1; i < (int)triangles.size(); i++) {
    //         if (triangles[i].distance < distance) {
    //             distance   = triangles[i].distance;
    //             closestIdx = i;
    //         }
    //     }

    //     /* Search for a new support point */
    //     Vec3 v       = triangles[closestIdx].normal;
    //     int  idxA    = supportIndex(A, v);
    //     int  idxB    = supportIndex(B, -v);
    //     Vec3 support = A.verts[idxA] - B.verts[idxB];

    //     /* Check convergence */
    //     Real supportDistance = support.dot(v);

    //     // Check termination
    //     if (abs(supportDistance - distance) <= epsRel) {
    //         upperBound  = supportDistance;
    //         lowerBound  = distance;
    //         bestFaceIdx = closestIdx;
    //         break;
    //     }

    //     if (distance > upperBound) {
    //         break;
    //     }

    //     // Check update
    //     if (supportDistance < upperBound) {
    //         upperBound  = supportDistance;
    //         lowerBound  = distance;
    //         bestFaceIdx = closestIdx;
    //     }

    //     /* Add new vertex via reconstructing faces */
    //     bool nondegenerate = true;
    //     Vec3 apex          = support;
    //     if (supportDistance - distance < epsRel) {
    //         nondegenerate = false;
    //     }

    //     // Reconstruct faces
    //     // replaceTriangle(, int faceIdx, const std::vector<Vec3> &vertices, std::vector<EPATriangle> &triangles)
    // }
}