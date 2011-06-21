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
 * @file   dgtalCairo-4-modes.cpp
 * @author Martial Tola <http://liris.cnrs.fr/martial.tola/>
 * @date   mercredi 25 mai 2011
 * 
 * @brief
 *
 * Simple example of class DGtalCairo.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/io/CairoViewers/DGtalCairo.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/helpers/Shapes.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  DGtalCairo viewer;

  Point p1( -1, -1, -2 );
  Point p2( 2, 2, 3 );
  Domain domain( p1, p2 );

  Point p3( 1, 1, 1 );
  Point p4( 2, -1, 3 );
  Point p5( -1, 2, 3 );
  Point p6( 0, 0, 0 );
  Point p0( 0, 2, 1 );

  //viewer << SetMode3DCairo( p1.styleName(), "Grid" );

  viewer << p1 << p2 << p3 << p4 << p5 << p6 << p0;

  //viewer << SetMode3DCairo(domain.styleName(), "PavingGrids");
  viewer << domain;
  
  viewer << Cairo3dCameraPosition(0.500000, 0.500000, 11.274194)
	<< Cairo3dCameraDirection(0.000000, 0.000000, -1.000000)
	<< Cairo3dCameraUpVector(0.000000, 1.000000, 0.000000);
  
  //viewer << DGtalCairo::Cairo3dWireFrame::yes;
  viewer.saveCairo("dgtalCairo-4-modes.png", DGtalCairo::CairoPNG, 600, 400);
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////




