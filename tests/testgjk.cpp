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

        handleSimplexLine(simplex, closest);
        INFO("simplex A: ", vecstr(simplex.verts[0]));
        INFO("simplex B: ", vecstr(simplex.verts[1]));
        INFO("closest: ", vecstr(closest));
        INFO("correct: ", vecstr(correct));
        REQUIRE(closest.isApprox(correct));
    }

    TEST_CASE("handleSimplexTri()") {}
}