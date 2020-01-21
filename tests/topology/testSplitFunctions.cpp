/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file testSplitFunctions.cpp
 * @ingroup Tests
 * @author Pablo Hernandez-Cerdan (\c pablo.hernandez.cerdan@outlook.com)
 *
 * @date 2019/01/08
 *
 * Testing class for SplitFunctions
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/SplitFunctions.h"
#include "DGtalCatch.h"
#include "DGtal/topology/KhalimskyCellHashFunctions.h"
#include "DGtal/topology/VoxelComplex.h"
#include "DGtal/topology/VoxelComplexFunctions.h"
// #include <DGtal/io/viewers/Viewer3D.h>
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;

TEST_CASE("computeSplitsEasy", "[computeSplits, getSplit]") {
    using Point = DGtal::Z3i::Point;
    Point lowerBound = {0,0,0};
    Point upperBound = {4,6,8};
    size_t requested_number_of_splits = 2;
    std::vector<unsigned int> splits(lowerBound.dimension);
    trace.beginBlock("computeSplits Easy");
    trace.info() << "lowerBound: " << lowerBound << std::endl;
    trace.info() << "upperBound: " << upperBound << std::endl;
    auto number_of_splits = functions::computeSplits(
            requested_number_of_splits,
            lowerBound, upperBound,
            splits.data());
    trace.info() << "number_of_splits: " << number_of_splits << std::endl;
    trace.info() << "splits: " << splits[ 0 ] << ", " << splits[ 1 ] << ", "
                 << splits[ 2 ] << std::endl;
    CHECK( number_of_splits == requested_number_of_splits );
    CHECK( splits[0] == 1 );
    CHECK( splits[1] == 1 );
    CHECK( splits[2] == 2 );
    trace.endBlock();
    trace.beginBlock("getSplit Easy");
    for (size_t split_index = 0; split_index < number_of_splits; split_index++) {
        auto outputBounds =
            DGtal::functions::getSplit(split_index, requested_number_of_splits,
                    lowerBound, upperBound);
        const auto & outputLowerBound = outputBounds[0];
        const auto & outputUpperBound = outputBounds[1];
        trace.info() << "split_index: " << split_index << std::endl;
        trace.info() << "outputLowerBound: " << outputLowerBound << std::endl;
        trace.info() << "outputUpperBound: " << outputUpperBound << std::endl;
    }
    size_t split_index = 1;
    auto outputBounds =
        DGtal::functions::getSplit(split_index, requested_number_of_splits,
                lowerBound, upperBound);
    const auto & outputLowerBound = outputBounds[0];
    const auto & outputUpperBound = outputBounds[1];
    CHECK( outputLowerBound[2] == 4 );
    CHECK( outputUpperBound[2] == 8 );
    trace.endBlock();
}
////////////////////// splitComplex ///////////////////////////
// Fixture for a X
struct Fixture_X {
    ///////////////////////////////////////////////////////////
    // type aliases
    ///////////////////////////////////////////////////////////
    using Point = DGtal::Z3i::Point;
    using Domain = DGtal::Z3i::Domain;
    using KSpace = DGtal::Z3i::KSpace;

    using FixtureDigitalTopology = DGtal::Z3i::DT26_6;
    using FixtureDigitalSet = DGtal::DigitalSetByAssociativeContainer<
        DGtal::Z3i::Domain,
        std::unordered_set<typename DGtal::Z3i::Domain::Point>>;
    using FixtureMap = std::unordered_map<KSpace::Cell, CubicalCellData>;
    using FixtureComplex = DGtal::VoxelComplex<KSpace, FixtureMap>;

    ///////////////////////////////////////////////////////////
    // fixture data
    FixtureComplex complex_fixture;
    FixtureDigitalSet set_fixture;
    KSpace ks_fixture; // needed because ConstAlias in CC constructor.
    ///////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////
    // Constructor
    ///////////////////////////////////////////////////////////
    Fixture_X() : complex_fixture(ks_fixture), set_fixture(create_set()) {
        create_complex_from_set(set_fixture);
    }

    ///////////////////////////////////////////////////////////
    // Function members
    ///////////////////////////////////////////////////////////
    FixtureDigitalSet create_set() {
        using namespace DGtal;

        Point p1(-16, -16, -16);
        Point p2(16, 16, 16);
        Domain domain(p1, p2);

        FixtureDigitalSet a_set(domain);
        std::vector<Point> center_set;
        center_set.reserve(9);

        Point c00(0, 0, 0);
        center_set.push_back(c00);
        Point c01x(-1, 0, 0);
        center_set.push_back(c01x);
        Point c10x(1, 0, 0);
        center_set.push_back(c10x);
        Point c02x(-2, 0, 0);
        center_set.push_back(c02x);
        Point c20x(2, 0, 0);
        center_set.push_back(c20x);

        Point c01y(0, -1, 0);
        center_set.push_back(c01y);
        Point c10y(0, 1, 0);
        center_set.push_back(c10y);
        Point c02y(0, -2, 0);
        center_set.push_back(c02y);
        Point c20y(0, 2, 0);
        center_set.push_back(c20y);

        Point z_pos(0, 0, 1);
        int branch_length(4);
        std::vector<Point> diagonals;
        diagonals.reserve(6);
        for (const auto &p : center_set) {
            diagonals.clear();
            for (int l = 0; l <= branch_length; ++l) {
                diagonals.push_back({l, l, 0});
                diagonals.push_back({l, -l, 0});
                diagonals.push_back({-l, l, 0});
                diagonals.push_back({-l, -l, 0});
                for (int z = -1; z <= 1; ++z)
                    for (const auto &d : diagonals)
                        a_set.insert(p + d + (z * z_pos));
            }
        }

		return a_set;
    }

    FixtureComplex &create_complex_from_set(FixtureDigitalSet &input_set) {

        ks_fixture.init(input_set.domain().lowerBound(),
                        input_set.domain().upperBound(), true);
        complex_fixture = FixtureComplex(ks_fixture);
        complex_fixture.construct(input_set);
        return complex_fixture;
    }
};

TEST_CASE_METHOD(Fixture_X, "splitComplex", "[parallel]") {
    using namespace DGtal::functions;
    auto & vc = complex_fixture;
    // CHECK(vc.nbCells(0) == 528);
    // CHECK(vc.nbCells(1) == 1276);
    // CHECK(vc.nbCells(2) == 1016);
    // CHECK(vc.nbCells(3) == 267);
    // Modify data of a cell to check if data is copied to sub complexes
    Point p{2, 0, 0};
    const auto modified_cell = vc.space().uSpel(p);
    const auto modified_cell_dim = vc.space().uDim(modified_cell);
    auto it_cell= vc.findCell(modified_cell_dim, modified_cell);
    CHECK(it_cell != vc.end(modified_cell_dim));
    DGtal::uint32_t magic_number = 15;
    it_cell->second.data = magic_number;
    size_t requested_number_of_splits = 2;
    SECTION("SplitComplex")  {
        trace.beginBlock("SplitComplex");
        trace.info() << "lowerBound" <<  vc.space().lowerBound() << std::endl;
        trace.info() << "upperBound" <<  vc.space().upperBound() << std::endl;
        auto out = splitComplex(vc, requested_number_of_splits);
        auto & sub_complexes = out.sub_complexes;
        CHECK(out.number_of_splits == requested_number_of_splits);
        CHECK(out.splits_domain_with_ghost_layers.empty() == true);

        CHECK(sub_complexes.size() == requested_number_of_splits);
        {
            size_t sub_index = 0;
            auto & sc = sub_complexes[sub_index];
            CHECK(sc.nbCells(0) != 0);
            CHECK(sc.nbCells(1) != 0);
            CHECK(sc.nbCells(2) != 0);
            CHECK(sc.nbCells(3) != 0);
            trace.info() << "lowerBound S0" <<  sc.space().lowerBound() << std::endl;
            trace.info() << "upperBound S0" <<  sc.space().upperBound() << std::endl;
            CHECK(sc.space().lowerBound() == typename KSpace::Point(-16,-16,-16));
            CHECK(sc.space().upperBound() == typename KSpace::Point(-1,16,16));
            CHECK(sc.space().lowerBound() == out.splits_domain[sub_index][0]);
            CHECK(sc.space().upperBound() == out.splits_domain[sub_index][1]);
        }
        {
            size_t sub_index = 1;
            auto & sc = sub_complexes[sub_index];
            CHECK(sc.nbCells(0) != 0);
            CHECK(sc.nbCells(1) != 0);
            CHECK(sc.nbCells(2) != 0);
            CHECK(sc.nbCells(3) != 0);
            trace.info() << "lowerBound S1" <<  sc.space().lowerBound() << std::endl;
            trace.info() << "upperBound S1" <<  sc.space().upperBound() << std::endl;
            CHECK(sc.space().lowerBound() == typename KSpace::Point(0,-16,-16));
            CHECK(sc.space().upperBound() == typename KSpace::Point(16,16,16));
            CHECK(sc.space().lowerBound() == out.splits_domain[sub_index][0]);
            CHECK(sc.space().upperBound() == out.splits_domain[sub_index][1]);
            // Check cell exist in sub_complex
            auto sc_it_cell = sc.findCell(it_cell->first);
            CHECK(sc_it_cell != sc.end(modified_cell_dim));
            // Check data is copied into sub_complexes
            CHECK(sc[modified_cell].data == magic_number);
        }
        const size_t row_voxels_sc0 = 7;
        CHECK(sub_complexes[0].nbCells(0) + sub_complexes[1].nbCells(0) ==
                vc.nbCells(0) + (row_voxels_sc0 + 1) * 4);
        CHECK(sub_complexes[0].nbCells(1) + sub_complexes[1].nbCells(1) ==
                vc.nbCells(1) + ( (row_voxels_sc0 + 1) * 3 + row_voxels_sc0 * 4));
        CHECK(sub_complexes[0].nbCells(2) + sub_complexes[1].nbCells(2) ==
                vc.nbCells(2) + row_voxels_sc0 * 3);
        CHECK(sub_complexes[0].nbCells(3) + sub_complexes[1].nbCells(3) ==
                vc.nbCells(3));
        trace.endBlock();
    }
    SECTION("SplitComplex with Ghost Layers")  {
        trace.beginBlock("SplitComplex with Ghost Layers");
        size_t wide_of_ghost_layer = 1;
        auto out = splitComplex(vc, requested_number_of_splits, wide_of_ghost_layer);
        auto & sub_complexes = out.sub_complexes;
        {
            size_t sub_index = 0;
            auto & sc = sub_complexes[sub_index];
            trace.info() << "lowerBound S0" <<  sc.space().lowerBound() << std::endl;
            trace.info() << "upperBound S0" <<  sc.space().upperBound() << std::endl;
            const auto sc_expected_lowerBound = typename KSpace::Point(-16,-16,-16);
            const auto sc_expected_upperBound = typename KSpace::Point(0,16,16);
            CHECK(sc.space().lowerBound() == sc_expected_lowerBound);
            CHECK(sc.space().upperBound() == sc_expected_upperBound);
            CHECK(out.splits_domain_with_ghost_layers[sub_index][0] == sc_expected_lowerBound);
            CHECK(out.splits_domain_with_ghost_layers[sub_index][1] == sc_expected_upperBound);
            CHECK(out.splits_domain[sub_index][0] == typename KSpace::Point(-16,-16,-16));
            CHECK(out.splits_domain[sub_index][1] == typename KSpace::Point(-1,16,16));
        }
        {
            size_t sub_index = 1;
            auto & sc = sub_complexes[sub_index];
            trace.info() << "lowerBound S1" <<  sc.space().lowerBound() << std::endl;
            trace.info() << "upperBound S1" <<  sc.space().upperBound() << std::endl;
            const auto sc_expected_lowerBound = typename KSpace::Point(-1,-16,-16);
            const auto sc_expected_upperBound = typename KSpace::Point(16, 16, 16);
            CHECK(sc.space().lowerBound() == sc_expected_lowerBound);
            CHECK(sc.space().upperBound() == sc_expected_upperBound);
            CHECK(out.splits_domain_with_ghost_layers[sub_index][0] == sc_expected_lowerBound);
            CHECK(out.splits_domain_with_ghost_layers[sub_index][1] == sc_expected_upperBound);
            CHECK(out.splits_domain[sub_index][0] == typename KSpace::Point(0,-16,-16));
            CHECK(out.splits_domain[sub_index][1] == typename KSpace::Point(16,16,16));
        }
        trace.endBlock();
    }
    SECTION("setBorderVoxelsAsFixed")  {
        trace.beginBlock("setBorderVoxelsAsFixed");
        auto out = splitComplex(vc, requested_number_of_splits);
        auto & sub_complexes = out.sub_complexes;
        const auto lowerBound = vc.space().lowerBound();
        const auto upperBound = vc.space().upperBound();
        const auto fixed_data = FixtureComplex::Parent::FIXED;
        {
            size_t sub_index = 1;
            auto & sc = sub_complexes[sub_index];
            const auto sub_lowerBound = sc.space().lowerBound();
            const auto sub_upperBound = sc.space().upperBound();
            trace.info() << "lowerBound S0" <<  sub_lowerBound << std::endl;
            trace.info() << "upperBound S0" <<  sub_upperBound << std::endl;
            const auto nullptr_point = static_cast<Point *>(nullptr);
            auto border_iterators =
                getBorderVoxels(sc, sub_lowerBound, sub_upperBound, nullptr_point,
                        &lowerBound, &upperBound);
            CHECK(border_iterators.size() == 15);
            const auto sc_expected_fixed = typename KSpace::Point(0,0,0);
            auto it_expected = sc.findCell(3, sc.space().uSpel(sc_expected_fixed));
            REQUIRE(it_expected != sc.end(3));
            CHECK(it_expected->second.data == 0);
            setBorderData(sc, border_iterators, fixed_data);
            CHECK(it_expected->second.data == fixed_data);
            {
            auto inside_voxel = vc.space().uSpel(typename KSpace::Point(0,0,0));
            auto faces = vc.cellBoundary(inside_voxel);
            auto coFaces = vc.cellCoBoundary(inside_voxel);
            std::cout << faces.size() << std::endl;
            std::cout << coFaces.size() << std::endl;
            }

            // SECTION( "visualize the split" ){
            //     int argc(1);
            //     char** argv(nullptr);
            //     QApplication app(argc, argv);
            //     Viewer3D<> viewer(vc.space());
            //     viewer.show();
            //
            //     viewer.setFillColor(Color(200, 0, 0, 100));
            //     auto & sc0 = sub_complexes[0];
            //     for ( auto it = sc0.begin(3); it!= sc0.end(3); ++it )
            //         viewer << it->first;
            //     viewer.setFillColor(Color(0, 0, 100, 100));
            //     auto & sc1 = sub_complexes[1];
            //     for ( auto it = sc1.begin(3); it!= sc1.end(3); ++it )
            //         viewer << it->first;
            //
            //     // All kspace voxels
            //     viewer.setFillColor(Color(40, 40, 40, 10));
            //     for ( auto it = vc.begin(3); it!= vc.end(3); ++it )
            //         viewer << it->first;
            //
            //     viewer << Viewer3D<>::updateDisplay;
            //     app.exec();
            // }
        }
        trace.endBlock();
    }
}
////////////////////// end splitComplex ////////////////////////