#include <cstdio>
#include <vector>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "gjkepa.hpp"

TEST_SUITE("EPA helper functions")
{
    TEST_CASE("replaceTriangle")
    {
        std::vector<Vec3> vertices = {{1.32013, -3.59795, 0},
                                      {-3.23657, 4.408, 12},
                                      {2.10064, 5.69651, -3.44659},
                                      {-7.99261, 3.09085, 0},
                                      {-5.8034, 9.29242, 7.09037}};

        std::vector<EPATriangle> triangles;
        epaStart4(vertices, triangles);
        replaceTriangle(4, 2, vertices, triangles);

        for (auto &tri : triangles) {
            printf("Triangle verts: %d %d %d\n", tri.verts[0], tri.verts[1], tri.verts[2]);
            printf("Triangle distance:  %f\n", tri.distance);
        }
    }
}