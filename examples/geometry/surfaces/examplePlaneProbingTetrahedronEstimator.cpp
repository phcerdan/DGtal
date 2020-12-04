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
 * @file
 * @ingroup Examples
 * @author Jocelyn Meyron (\c jocelyn.meyron@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systemes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2020/12/04
 *
 * An example file named examplePlaneProbingTetrahedronEstimator.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "ConfigExamples.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/geometry/surfaces/DigitalPlanePredicate.h"
#include "DGtal/geometry/surfaces/estimation/PlaneProbingHNeighborhood.h"
#include "DGtal/geometry/surfaces/estimation/PlaneProbingTetrahedronEstimator.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

using DigitalPlane = DigitalPlanePredicate<Z3i::Space>;
using Point = DigitalPlane::Vector;
using Vector = DigitalPlane::Point;
using Neighborhood = PlaneProbingHNeighborhood<DigitalPlane>;
using Estimator = PlaneProbingTetrahedronEstimator<DigitalPlane, ProbingMode::H>;

///////////////////////////////////////////////////////////////////////////////

int main(void)
{
  Vector n(2, 6, 15);
  DigitalPlane plane(n, 0, n.norm1());
  Point o(0, 0, 0);
  detail::Triplet<Point> m = { Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1) };
  Estimator estimator(o, m, plane);

  int it = 0;
  while (estimator.advance()) {
      it++;

      std::clog << "it = " << it << " "
          << estimator.m(0) << " " << estimator.m(1) << " " << estimator.m(2) << " "
          << estimator.getNormal() << std::endl;
  }

  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
