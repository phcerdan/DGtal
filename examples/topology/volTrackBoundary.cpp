/**
 * @file topology/volTrackBoundary.cpp
 * @ingroup Examples
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/02/06
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

/**
 * In many circumstances, it is better to use the presented graph
 * structure of digital surfaces. For instance it may be used to find
 * the surface just by searching it by adjacencies. This process is
 * called \b tracking. This is done for you by
 * static method Surfaces::trackBoundary.
 * 
 * @see \ref dgtal_digsurf_sec2_2
 * 
 * On the lobser.vol volume, volTrackBoundary.cpp extracts 148364
 * surfels in 351ms.
 * 
 * @verbatim
 * # Commands
 * $ ./examples/topology/volTrackBoundary ../examples/samples/lobster.vol 50 255
 * @endverbatim
 * 
 * @image html volTrackBoundary-lobster.png "Digital surface that is the boundary of a (6,18)-connected component in image lobst* er.vol, extracted by tracking from an initial surfel in 351ms."
 * 
 * \example topology/volTrackBoundary.cpp
 */


///////////////////////////////////////////////////////////////////////////////
//! [volTrackBoundary-basicIncludes]
#include <iostream>
#include <queue>
#include "DGtal/base/Common.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/Color.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/IntervalForegroundPredicate.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/CodedKhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"
//! [volTrackBoundary-basicIncludes]

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////

void usage( int /*argc*/, char** argv )
{
  std::cerr << "Usage: " << argv[ 0 ] << " <fileName.vol> <minT> <maxT>" << std::endl;
  std::cerr << "\t - displays the boundary of the shape stored in vol file <fileName.vol>." << std::endl;
  std::cerr << "\t - voxel v belongs to the shape iff its value I(v) follows minT <= I(v) <= maxT." << std::endl;
}

int main( int argc, char** argv )
{
  if ( argc < 4 )
    {
      usage( argc, argv );
      return 1;
    }
  std::string inputFilename = argv[ 1 ];
  unsigned int minThreshold = atoi( argv[ 2 ] );
  unsigned int maxThreshold = atoi( argv[ 3 ] );

  //! [volTrackBoundary-typedefs]
  typedef SpaceND< 3, int32_t >                         Space;
  typedef HyperRectDomain< Space >                      Domain;
  typedef DigitalSetBySTLSet<Domain>                    DigitalSet;
  typedef CodedKhalimskySpaceND< 3, int32_t, uint64_t > KSpace;
  typedef KSpace::SCell                                 SCell;
  typedef KSpace::SCellSet                              SCellSet;
  typedef SurfelAdjacency<KSpace::dimension>            MySurfelAdjacency;
  typedef KhalimskySpaceND< 3, int32_t >                DisplayKSpace;  // used for display.
  typedef SignedKhalimskyCell< 3, int32_t >             DisplaySCell;
  //! [volTrackBoundary-typedefs]

  //! [volTrackBoundary-readVol]
  trace.beginBlock( "Reading vol file into an image." );
  typedef ImageSelector < Domain, int>::Type Image;
  Image image = VolReader<Image>::importVol(inputFilename);
  typedef IntervalForegroundPredicate<Image> ThresholdedImage;
  ThresholdedImage thresholdedImage( image, minThreshold, maxThreshold );
  trace.endBlock();
  //! [volTrackBoundary-readVol]
  
  
  //! [volTrackBoundary-KSpace]
  trace.beginBlock( "Construct the Khalimsky space from the image domain." );
  KSpace ks;
  bool space_ok = ks.init( image.domain().lowerBound(), 
                           image.domain().upperBound(), true );
  if (!space_ok)
    {
      trace.error() << "Error in the Khamisky space construction."<<std::endl;
      return 2;
    }
  trace.endBlock();
  //! [volTrackBoundary-KSpace]

  //! [volTrackBoundary-SurfelAdjacency]
  MySurfelAdjacency surfAdj( true ); // interior in all directions.
  //! [volTrackBoundary-SurfelAdjacency]

  //! [volTrackBoundary-ExtractingSurface]
  trace.beginBlock( "Extracting boundary by tracking from an initial bel." );
  SCellSet boundary;
  SCell bel = Surfaces<KSpace>::findABel( ks, thresholdedImage, 100000 );
  Surfaces<KSpace>::trackBoundary( boundary, ks, 
                                   surfAdj,
                                   thresholdedImage, bel );
  trace.endBlock();
  //! [volTrackBoundary-ExtractingSurface]

  //! [volTrackBoundary-DisplayingSurface]
  trace.beginBlock( "Displaying surface in Viewer3D." );
  QApplication application(argc,argv);
  DisplayKSpace dks;
  space_ok = dks.init( image.domain().lowerBound(), image.domain().upperBound(), true );
  Viewer3D<> viewer( dks );
  viewer.show(); 
  DisplaySCell scell;
  viewer << SetMode3D( scell.className(), "Basic" );
  viewer << CustomColors3D(Color(250, 0, 0 ), Color( 128, 128, 128 ) );
  unsigned long nbSurfels = 0;
  for ( SCellSet::const_iterator it = boundary.begin(),
          it_end = boundary.end(); it != it_end; ++it, ++nbSurfels )
    {
      scell = dks.sCell( ks.sKCoords( *it ), ks.sSign( *it ) );
      viewer << scell;
    }
  viewer << Viewer3D<>::updateDisplay;
  trace.info() << "nb surfels = " << nbSurfels << std::endl;
  trace.endBlock();
  return application.exec();
  //! [volTrackBoundary-DisplayingSurface]
}

