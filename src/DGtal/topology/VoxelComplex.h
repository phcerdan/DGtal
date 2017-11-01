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

#pragma once

/**
* @file VoxelComplex.h
* @author Pablo Hernandez-Cerdan (\c pablo.hernandez.cerdan@outlook.com)
* Insitute of Fundamental Sciences, Massey University, New Zealand.
*
* @date 2016/01/21
*
* Header file for module VoxelComplex.cpp
*
* This file is part of the DGtal library.
*/
#if defined(VoxelComplex_RECURSES)
#error Recursive header files inclusion detected in VoxelComplex.h
#else // defined(VoxelComplex_RECURSES)
/** Prevents recursive inclusion of headers. */
#define VoxelComplex_RECURSES

#if !defined VoxelComplex_h
/** Prevents repeated inclusion of headers. */
#define VoxelComplex_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <DGtal/topology/CubicalComplex.h>
#include <DGtal/topology/Object.h>
#include <DGtal/kernel/sets/DigitalSetBySTLSet.h>
#include <DGtal/topology/DigitalTopology.h>
#include "boost/dynamic_bitset.hpp"

namespace DGtal
{
  /**
  * Any cell is stored within a cubical complex with an associated
  * data. You may define your own associated data but the type must
  * derive from this class CubicalCellData. Its basic usage is to
  * store flags associated to the cells, but it may store other
  * values.
  *
  * @note Predefined flags are CubicalComplex::REMOVED,
  * CubicalComplex::COLLAPSIBLE, CubicalComplex::FIXED,
  * CubicalComplex::USER1. Other bits can be used to associate an
  * integer to the cell. The corresponding mask is
  * CubicalComplex::VALUE.
  *
  * @note Such data is notably used in collapse operation
  * (CubicalComplex::collapse).
  */
    // BOOST_STATIC_CONSTANT( uint32_t, NOT_DEATH     = 0x10000000 );
  struct VoxelComplexCellData: public CubicalCellData {
    inline VoxelComplexCellData() :
      CubicalCellData(), birth_date( 0 ), mature( 0 ) {}
    VoxelComplexCellData( uint32_t d ) :
      CubicalCellData( d ), birth_date( 0 ), mature( 0 ) {}
    VoxelComplexCellData( uint32_t d, uint32_t birth, uint32_t lifetime ) :
      CubicalCellData( d ), birth_date( birth ), mature( lifetime ) {}
    uint32_t birth_date;
    uint32_t mature;
  };
  /////////////////////////////////////////////////////////////////////////////
  // template class VoxelComplex
  /**
  * Description of template class 'VoxelComplex' <p> \brief
  * This class represents a voxel complex living in some
  * Khalimsky space. Voxel complexes are derived from @sa cubical complexes,
  * with specialized methods to deal with spels.
  *
  * The aim is to implement critical kernels, ie, cliques of spels, as shown
  * by M.Couprie and G.Bertrand (2014).
  *
  * Implemented using resources from \a CubicalComplex and
  * (Digital) \a Object.
  *
  * @tparam TKSpace any model of concepts::CCellularGridSpaceND, i.e. a type
  * that models a Khalimsky space.
  *
  * @tparam TCellContainer any model of associative container, mapping
  * a KSpace::Cell to a CubicalCellData or any type deriving from
  * it. It could be for instance a std::map or a
  * std::unordered_map. Note that unfortunately, unordered_map are
  * (strangely) not models of boost::AssociativeContainer, hence we
  * cannot check concepts here.
  *
  */

template < typename TKSpace, typename TObject, typename TCellContainer =
               typename TKSpace::template CellMap< VoxelComplexCellData >::Type >
class VoxelComplex :
  public CubicalComplex< TKSpace, TCellContainer >
{
  // The TObject::DigitalSet::Container must be associative.
  BOOST_CONCEPT_ASSERT(( concepts::CSTLAssociativeContainer< typename TObject::DigitalSet::Container > ));
    // ----------------------- associated types ------------------------------
public:
  using Parent          = CubicalComplex< TKSpace, TCellContainer > ; ///< Type of the parent class CubicalComplex.
  using Self            = VoxelComplex< TKSpace, TObject, TCellContainer > ; ///< Type of this instance of VoxelComplex.
  using KSpace          = TKSpace ;  ///< Type of the cellular grid space.
  using CellContainer   = TCellContainer ; ///< Type for storing cells, an associative container Cell -> Data
  using Data            = typename CellContainer::mapped_type ; ///< Type of data associated to each cell.
  using Object          = TObject ; ///< Type for input Object storing the digital set of spels with a topology.
  using DigitalTopology = typename TObject::DigitalTopology ; ///< Type for the type of the digital topology.

  /// The dimension of the embedding space.
  static const Dimension dimension = KSpace::dimension;
  using Integer              = typename KSpace::Integer ; ///< Type for integers in the space.
  using Cell                 = typename KSpace::Cell ;    ///< Type for a cell in the space.
  using Cells                = typename KSpace::Cells ;   ///< Type for a sequence of cells in the space.
  using Space                = typename KSpace::Space ;   ///< Type of the digital space
  using Size                 = typename KSpace::Size ;    ///< Type for a number of elements
  using Point                = typename KSpace::Point ;   ///< Type for a point in the digital space
  using DirIterator          = typename KSpace::DirIterator ; ///< Type for iterating over cell directions
  using CellMap              = CellContainer ; ///< Type for storing cells, an associative container Cell -> Data
  using CellMapConstIterator = typename CellMap::const_iterator ; ///< Const iterator for visiting type CellMap
  using CellMapIterator      = typename CellMap::iterator ; ///< Iterator for visiting type CellMap

  //Clique alias
  using Clique = Parent;
  using CliqueContainer = std::vector<Clique>;

  //Tables
  using ConfigMap = boost::dynamic_bitset<>;
  using PointToMaskMap = std::unordered_map<Point, unsigned int>;
public:
  /// Inherit all constructors from parent CubicalComplex.
  using CubicalComplex<KSpace, CellContainer>::CubicalComplex;
  /**
   * Assignment.
   * @param other the object to copy.
   * @return a reference on 'this'.
   */
  Self & operator= ( const Self & other );

  /**
   * Construct the VoxelComplex with target \a obj
   *
   * @param obj input object to copy to the complex.
   *
   */
  void construct ( const TObject& obj );

  /**
   * Construct from digital set and precomputed look up table for simplicity.
   * Object and methods involving object will be empty/invalid.
   * But they might not me needed thanks to the loaded table.
   *
   * @note if you require object operation but also wants the speed up of the table
   * call construct from object, and then loadTable.
   *
   * @param input_set points to construct the CubicalComplex.
   * Calls @ref CubicalComplex::construct
   * @param input_table input table[conf]->bool
   *
   * @see LookUpTableFunctions.h
   */
  void construct (const typename Object::DigitalSet & input_set,
                  const Alias<ConfigMap> input_table );

  /**
   * Set precomputed look up table for simplicity.
   *
   * @param input_table input table[conf]->bool
   *
   * @see LookUpTableFunctions.h
   */
  void setSimplicityTable(const Alias<ConfigMap> input_table);

  /**
   * Get const reference to table[conf]->bool for simplicity.
   *
   * @return table[conf]->bool for simplicity.
   */
  const ConfigMap & table() const;
  /**
   * Get const reference to isTableLoaded bool member.
   *
   * @return true if table for simplicity has been loaded.
   */
  const bool & isTableLoaded() const;
  void voxelClose( const Cell & kcell );
  void cellsClose( Dimension k, const Cells & cells );
  void insertVoxelCell( const Cell & kcell, const bool & close_it = true, const Data& data = Data() );
  void insertVoxelPoint( const Point & dig_point, const bool & close_it = true, const Data& data = Data());
    /**
    * Clears the voxel complex, which becomes empty.
    * This includes the khalmisky cells and also the object points.
    */
  void clear();
  /**
   * Get reference to a point in the DigitalSet of myObject
   * corresponding to input voxel
   *
   * @param voxel input voxel
   *
   * @return A point in myObject
   */
  const typename Object::Point & objPointFromVoxel(
      const Cell & voxel) const;
  //------ Spels ------//
  /**
   * Get pointels that are Faces of input_cell.
   * @note If input_cell is a spel, then it will be emplaced only if it belongs
   * to the complex.
   *
   * @param pointels_out in/out pointels.
   * @param input_cell input cell from which get the surrounding pointels.
   *
   * @see KhalimskySpaceND::uFaces
   */
  void pointelsFromCell(std::set<Cell>& pointels_out,
                        const Cell & input_cell) const;
  /**
   * Get spels that are coFaces of input_cell.
   * @note If input_cell is a spel, then it will be emplaced only if
   * it belongs to the complex.
   *
   * @param spels_out in/out spels.
   * @param input_cell input cell from which get the surrounding spels.
   *
   * @see KhalimskySpaceND::uCoFaces
   */
  void spelsFromCell(std::set<Cell>& spels_out, const Cell & input_cell) const;
  Clique Kneighborhood(const Cell & input_cell) const ;

  /**
   * Populate this complex object from the spels belonging to the Khalimsky space.
   *
   * @return const reference to object().
   */
  const Object & populateObjectFromCells();
  /**
   * Check if the input_spel from khalimsky space is simple using object properties.
   *
   * @param input_spel khalimsky space spel.
   *
   * @return true if input_spel is simple.
   * @note It uses isSimple from Object.
   * There are no guarantees than objectSet and khalimsky space are synchronized
   * so user must take care of the sync of space and kspace,
   * using insertion methods such as @see insertVoxelCell, insertVoxelPoint
   * @see Object::isSimple
   */
  bool isSimple(const Cell& input_spel) const ;

  /**
   * @return true if object is connected, false if disconnected.
   *
   * @note connectedness::unkwown is not possible.
   *
   * @see Object::computeConnectedness.
   */
  bool isConnected() const ;

  /**
   * @return Object representing the spels.
   */
  Object & object() ;
  /**
   * @return Object representing the spels, read only.
   */
  const Object & object() const ;

  /**
   * @return digitalSet of Object representing the spels.
   */
  typename Object::DigitalSet & objectSet() ;

  /**
   * @return digitalSet of Object representing the spels, read only.
   */
  const typename Object::DigitalSet & objectSet() const ;


  //------ Cliques ------//
  // Cliques, union of adjacent spels.
  // The intersection of all spels of the clique define the type.
  // 0-clique, 1-clique, 2-clique. 3-clique are isolated spels.
public:
  std::pair<bool, Clique> criticalCliquePair(
      const Dimension d, const CellMapConstIterator & cellMapIterator) const;

  std::array<CliqueContainer, dimension + 1> criticalCliques(
      bool verbose = false) const;

  std::array<CliqueContainer, dimension + 1> criticalCliques(
      const Parent & cubical, bool verbose = false) const;

  CliqueContainer criticalCliquesForD(
      const Dimension d, const Parent & cubical, bool verbose = false) const;
  /**
   * Compute the criticality of the surfel between A,B voxels and returns the associated 2-clique.
   *
   * @param A Spel (voxel).
   * @param B Spel (voxel).
   * @param verbose flag for verbose output
   *
   * @return <is_critical, 2-clique>
   * @note A,B can be swaped with same result.
   */
  std::pair<bool, Clique> K_2(
      const Cell & A,
      const Cell & B,
      bool verbose = false) const;
  /**
   * K_2 from two DigitalSet Points (uCoords)
   * @see VoxelComplex::K_2
   *
   * @param A voxel point
   * @param B voxel point
   * @param verbose flag for verbose output
   *
   * @return <is_critical, 2-clique>
   */
  std::pair<bool, Clique> K_2(
      const typename Object::Point & A,
      const typename  Object::Point & B,
      bool verbose = false ) const;

  /**
   * K_2 from a 2-face (surfel)
   * @see VoxelComplex::K_2
   *
   * @param face2 surfel between 2 spels
   * @param verbose flag for verbose output
   *
   * @return <is_critical, 2-clique>
   */
  std::pair<bool, Clique> K_2( const Cell & face2,
                               bool verbose = false ) const;

  /**
   * Compute the criticality of the linel and the associated 1-clique.
   *
   * @param face1 linel between two pointels
   * @param verbose flag for verbose output
   *
   * @return <is_critical, 1-clique>
   */
  std::pair<bool, Clique> K_1( const Cell &face1,
                               bool verbose = false ) const;

  /**
   * Compute the criticality of the pointel and the associated 0-clique.
   *
   * @param face0 a pointel cell.
   * @param verbose flag for verbose output
   *
   * @return <is_critical, 0-clique>
   */
  std::pair<bool, Clique> K_0( const Cell &face0,
                               bool verbose = false ) const;

  /**
   * Compute the criticality of the spel and the associated 3-clique.
   * It uses isSimple to check criticality.
   *
   * @param input_spel input spel
   * @param verbose bool flag
   *
   * @return <is_critical, 3-clique (=spel)>
   */
  std::pair<bool, Clique> K_3( const Cell &input_spel,
                               bool verbose = false ) const;
  /**
   * True if input cell is a cell with max dimension.
   *
   * @param c input Cell
   *
   * @return true if c is a Spel
   */

  bool isSpel( const Cell & c ) const;
  /**
   * Surfel between two adjacent spels.
   *
   * @param A input spel cell
   * @param B input spel cell
   *
   * @return Surfel between adjacent spels A,B
   */
  Cell surfelBetweenAdjacentSpels( const Cell & A, const Cell & B) const;

 /*------------- Data --------------*/
protected:
  Object myObject; ///< Object with a topology representing spels.
  CountedPtrOrPtr<ConfigMap> myTablePtr;
  CountedPtrOrPtr<PointToMaskMap> myPointToMaskPtr;
  bool myIsTableLoaded{false};

 /*------------- Internal Methods --------------*/
  /**
   * pointToMask map, used internally in @ref VoxelComplex::isSimple for
   * @ref LookUpTableFunctions.h::getSpelNeighborhoodConfigurationOccupancy
   *
   * @return reference to pointToMaskMap member.
   *
   * @see LookUpTableFunctions.h::pointToBitMaskMap()
   */
  const PointToMaskMap & pointToMask() const;

  /**
   * Populate myObject member with an empty set with valid domain and topology.
   * Used in @ref VoxelComplex::construct with digital sets.
   *
   * @param dig_set input digital set.
   */
  void instantiateEmptyObject(const typename Object::DigitalSet & dig_set);

    // ----------------------- Interface --------------------------------------
public:

  /**
   * Writes/Displays the object on an output stream.
   * @param out the output stream where the object is written.
   */
  void selfDisplay ( std::ostream & out ) const;

  /**
   * Checks the validity/consistency of the object.
   * @return 'true' if the object is valid, 'false' otherwise.
   */
  bool isValid() const;

  // --------------- CDrawableWithBoard2D realization ------------------
public:
  /**
   * @return the style name used for drawing this object.
   */
  std::string className() const;



}; // end of class CubicalComplex

  /**
  * Overloads 'operator<<' for displaying objects of class 'VoxelComplex'.
  * @param out the output stream where the object is written.
  * @param object the object of class 'VoxelComplex' to write.
  * @return the output stream after the writing.
  */
  template <typename TKSpace, typename TObject, typename TCellContainer>
  std::ostream&
  operator<< ( std::ostream & out,
               const VoxelComplex<TKSpace, TObject, TCellContainer> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/topology/VoxelComplex.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined VoxelComplex_h

#undef VoxelComplex_RECURSES
#endif // else defined(VoxelComplex_RECURSES)