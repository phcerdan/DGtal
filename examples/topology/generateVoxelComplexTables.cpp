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
 * @file generateVoxelComplexTables.cpp
 * @ingroup Examples
 * @author Pablo Hernandez-Cerdan (\c pablo.hernandez.cerdan@outlook.com)
 * Institute of Fundamental Sciences. Massey University.
 * Palmerston North, New Zealand
 *
 * @date 2016/03/20
 *
 * Creates precomputed tables for all possible configurations of the
 * neighborhood of a voxel in 3D with a topology of 26 neighbors.
 * The table is a map between all this configurations and the result of applying
 * an input selected predicate function to a center voxel.
 * The options for the functions are the Skel function from
 * @see VoxelComplexFunctions.h.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <bitset>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/VoxelComplexFunctions.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////

/**
   Generates a table mapping the number of configuration of a 26 topology voxel
   neighborhood, and the boolean result of a predicate function applied
   to the central point for each configuration.
   The configuration is determined by a sequence of bits, the
   first bit for the point in the neighborhood, the second bit for the second
   point, etc. When set to one, the point is in the neighborhood.

   @tparam TVoxelComplex the type of the VoxelComplex whose
   property we wish to precompute.

   @tparam TMap the type used to store the mapping configuration -> bool.

   @param map (modified) the mapping configuration -> bool.
   @param skelFunction a predicate function related to the property we want to check.
*/
template <typename TVoxelComplex, typename TMap>
void
generateVoxelComplexTable(
	     TMap & map,
	     std::function< bool(
	       const TVoxelComplex & ,
	       const typename TVoxelComplex::Cell & )
	     > skelFunction
	     )
{
  using Object = typename TVoxelComplex::Object;
  using DigitalSet = typename Object::DigitalSet ;
  using Point = typename Object::Point ;
  using Domain = typename DigitalSet::Domain ;
  using DomainConstIterator = typename Domain::ConstIterator ;
  using KSpace = typename TVoxelComplex::KSpace;
  using DigitalTopology = typename Object::DigitalTopology;
  using ForegroundAdjacency = typename Object::ForegroundAdjacency;
  using BackgroundAdjacency = typename Object::BackgroundAdjacency;
  ForegroundAdjacency adjF;
  BackgroundAdjacency adjB;
  DigitalTopology dt( adjF, adjB,
      DigitalTopologyProperties::JORDAN_DT);

  Point p1 = Point::diagonal( -1 );
  Point p2 = Point::diagonal(  1 );
  Point c = Point::diagonal( 0 );
  Domain domain( p1, p2 );
  DigitalSet shapeSet( domain );
  Object shape( dt, shapeSet );
  unsigned int k = 0;
  for ( DomainConstIterator it = domain.begin(); it != domain.end(); ++it )
    if ( *it != c ) ++k;
  ASSERT( ( k < 32 )
	  && "[generateVoxelComplexTable] number of configurations is too high." );
  unsigned int nbCfg = 1 << k;

  KSpace ks;
  // Pad KSpace domain.
  ks.init(shape.domain().lowerBound() + Point::diagonal( -1 ) ,
          shape.domain().upperBound() + Point::diagonal( 1 ),
	  true);
  TVoxelComplex vc(ks);
  vc.construct(shape);
  for ( unsigned int cfg = 0; cfg < nbCfg; ++cfg ){
    if ( ( cfg % 1000 ) == 0 )
      trace.progressBar( (double) cfg, (double) nbCfg );
    vc.clear();
    vc.insertVoxelPoint(c);
    unsigned int mask = 1;
    for ( DomainConstIterator it = domain.begin(); it != domain.end(); ++it ){
      if ( *it != c ) {
	if ( cfg & mask ) vc.insertVoxelPoint( *it );
	mask <<= 1;
      }
    }
    const auto &kcell = vc.space().uSpel(c);
    bool predicate_output = skelFunction(vc, kcell);
    map[ cfg ] = predicate_output;
  }
}


int main( int argc, char** argv )
{
  typedef std::bitset<67108864> ConfigMap; // 2^26

  using namespace Z3i;
  using DigitalSet = DigitalSetByAssociativeContainer<
	    Domain, std::unordered_set< typename Domain::Point> >;
  using Object = Object<DT26_6, DigitalSet>;
  using VoxelComplex = VoxelComplex<KSpace, Object>;

  std::function< bool(
		 const VoxelComplex & ,
		 const typename VoxelComplex::Cell & )
	       > skelFunction;
  string error_message(
      "Provide one of the following arguments for select function:\n"
      "- skelIsthmus \n "
      "- oneIsthmus \n "
      "- twoIsthmus \n ");
  if (argc != 2 ){
    cout << error_message << std::endl;
    return 1;
  }
  std::string input_str = std::string(argv[1]);
  if (input_str == "skelIsthmus")
    skelFunction = functions::skelIsthmus<VoxelComplex>;
  else if (input_str == "oneIsthmus")
    skelFunction = functions::oneIsthmus<VoxelComplex>;
  else if (input_str == "twoIsthmus")
    skelFunction = functions::twoIsthmus<VoxelComplex>;
  else{
    cout << error_message << endl;
    return 1;
  }

  trace.beginBlock ( "Generate " + input_str + " table for 26_6 topology" );
  // Too big for stack. Use heap instead.
  auto table26_6 = make_shared<ConfigMap>();
  generateVoxelComplexTable< VoxelComplex >(
      *table26_6,
      skelFunction );

  string filename = input_str + "_table26_6.txt";
  trace.info() << "Save to file... " + filename << std::endl;
  ofstream file26_6( filename );
  file26_6 << *table26_6;
  file26_6.close();

  trace.endBlock();

  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
