#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <sstream>

#include "doctest.h"
#include "gjkepa.hpp"

std::string vecstr(const Vec3 &v)
{
    std::ostringstream oss;
    oss << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return oss.str();
}

void transformSimplex(const Eigen::Matrix3f &trans, Simplex &simplex)
{
    for (auto &v : simplex.verts) {
        v = trans * v;
    }
}

TEST_SUITE("GJK helper functions")
{
    TEST_CASE("handleSimplexLine()")
    {
        Simplex simplex;
        Vec3    closest;
        Vec3    correct;

        SUBCASE("Apex region 1")
        {
            simplex = {{1.5695644479381963, 2.888135673370931, 2.426578170072121},
                       {-0.3227685063381752, 1.7944596712810046, 1.3019346106146925}};
            correct = {-0.3227685063381752, 1.7944596712810046, 1.3019346106146925};
        }

        SUBCASE("Apex region 2")
        {
            simplex = {{1.5695644479381963, 2.888135673370931, -5.565083853908408},
                       {0.8173082025418976, 0.19334958423338233, -4.2489385907111155}};
            correct = {0.8173082025418976, 0.19334958423338233, -4.2489385907111155};
        }

        SUBCASE("Edge region 1")
        {
            simplex = {{5.657664790366462, 3.822531595273899, -2.7466729184312673},
                       {-4.609189575824734, 2.515301916674902, 3.1790994850367618}};
            correct = {-0.07081339329407843, 3.0931517455365354, 0.5596618871336689};
        }

        SUBCASE("Edge region 2")
        {
            simplex = {{1.5695644479381963, 2.888135673370931, 2.426578170072121},
                       {-2.965762802691635, -5.220990483557489, 0.9349958216279619}};
            correct = {-0.17990095565067743, -0.23989343696807275, 1.8512124751612373};
        }

        INFO("simplex A: ", vecstr(simplex.verts[0]));
        INFO("simplex B: ", vecstr(simplex.verts[1]));

        handleSimplexLine(simplex, closest);
        INFO("closest: ", vecstr(closest));
        INFO("correct: ", vecstr(correct));
        REQUIRE(closest.isApprox(correct));
    }

    TEST_CASE("handleSimplexTri()")
    {
        Simplex simplex;
        Vec3    closest;
        Vec3    correct;

        SUBCASE("Inside triangle")
        {
            simplex = {{5.75, 1.04, 3.14}, {-1.34, -6.2, 3.14}, {-5.13, 4.43, 3.14}};
            correct = {0, 0, 3.14};
            SUBCASE("Plane") {}
            SUBCASE("Rotated")
            {
                auto randomQ = Eigen::Quaterniond::UnitRandom();
                auto rotate  = randomQ.toRotationMatrix().cast<Real>();
                correct      = rotate * correct;
                transformSimplex(rotate, simplex);
            }

            // Check correct normal direction
            auto A      = simplex.verts[0];
            auto B      = simplex.verts[1];
            auto C      = simplex.verts[2];
            auto normal = (B - A).cross(C - A);
            REQUIRE(normal.dot(C) <= 0);
        }

        SUBCASE("Edge AC")
        {
            simplex = {
                {3.56, -0.34, -1.2},
                {4.57, 2.42, -1.2},
                {-2.02, 4, -1.2},
            };
            correct = {1.17707692, 1.51338462, -1.2};
            SUBCASE("Plane") {}
            SUBCASE("Rotated")
            {
                auto randomQ = Eigen::Quaterniond::UnitRandom();
                auto rotate  = randomQ.toRotationMatrix().cast<Real>();
                correct      = rotate * correct;
                transformSimplex(rotate, simplex);
            }
        }

        SUBCASE("Edge BC")
        {
            simplex = {{10, 10, 0.05}, {3, 0, 0.05}, {-5, 4, 0.05}};
            correct = {0.6, 1.2, 0.05};
            SUBCASE("Plane") {}
            SUBCASE("Rotated")
            {
                auto randomQ = Eigen::Quaterniond::UnitRandom();
                auto rotate  = randomQ.toRotationMatrix().cast<Real>();
                correct      = rotate * correct;
                transformSimplex(rotate, simplex);
            }
        }

        SUBCASE("Apex")
        {
            simplex = {{5, 10, 19}, {15, 3, 19}, {2, 1, 19}};
            correct = {2, 1, 19};
            SUBCASE("Plane") {}
            SUBCASE("Rotated")
            {
                auto randomQ = Eigen::Quaterniond::UnitRandom();
                auto rotate  = randomQ.toRotationMatrix().cast<Real>();
                correct      = rotate * correct;
                transformSimplex(rotate, simplex);
            }
        }

        SUBCASE("Exact apex")
        {
            simplex = {{2, 6.6, 0}, {5.2, 2.1, 0}, {0, 0, 0}};
            correct = {0, 0, 0};
        }


        INFO("simplex A: ", vecstr(simplex.verts[0]));
        INFO("simplex B: ", vecstr(simplex.verts[1]));
        INFO("simplex C: ", vecstr(simplex.verts[2]));
        handleSimplexTri(simplex, closest);

        // Ensure correct closest point
        INFO("closest: ", vecstr(closest));
        INFO("correct: ", vecstr(correct));
        REQUIRE(closest.isApprox(correct));
    }

    TEST_CASE("handleSimplexTetra()")
    {
        Simplex simplex;
        Vec3    closest;
        Vec3    correct;

        handleSimplexTetra(simplex, closest);

        SUBCASE("Inside tetra")
        {
            simplex = {
                {-2.92593, -2.25338, 0}, {3.21313, -0.99584, -1}, {-3, 2, -1.76271}, {-1.34325, 0.9676, 3.71282}};
            correct = {0, 0, 0};
        }

        SUBCASE("Edge")
        {
            simplex = {{0.73248, 0.11043, -1.68946},
                       {5.83595, -1.21708, -1},
                       {3.29811, 2.54955, -1.76271},
                       {1.18249, 0.04732, 3.4885}};
            correct = {0.87281, 0.09075, -0.07475};
        }

        SUBCASE("Apex")
        {
            simplex = {{-2.49432, -0.44896, -1.68946},
                       {5.83595, -1.21708, -1},
                       {3.29811, 2.54955, -1.76271},
                       {0, 0, -0.36963}};
            correct = {0, 0, -0.36963};
        }

        SUBCASE("Face")
        {
            simplex = {{-0.44258, 1.15203, -1.68946},
                       {5.83595, -1.21708, -1},
                       {3.29811, 2.54955, -1.76271},
                       {2.67042, 1.21494, 4.25063}};
            correct = {0.44728, 1.1137, -0.2462};
        }

        INFO("simplex A: ", vecstr(simplex.verts[0]));
        INFO("simplex B: ", vecstr(simplex.verts[1]));
        INFO("simplex C: ", vecstr(simplex.verts[2]));
        INFO("simplex D: ", vecstr(simplex.verts[3]));

        handleSimplexTetra(simplex, closest);

        // Ensure correct closest point
        INFO("closest: ", vecstr(closest));
        INFO("correct: ", vecstr(correct));
        REQUIRE(closest.isApprox(correct));
    }
}