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
 * @file DigitalSetBySTLVector.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @author Sebastien Fourey (\c Sebastien.Fourey@greyc.ensicaen.fr )
 * Groupe de Recherche en Informatique, Image, Automatique et
 * Instrumentation de Caen - GREYC (CNRS, UMR 6072), ENSICAEN, France
 *
 * @date 2010/07/01
 *
 * Header file for module DigitalSetBySTLVector.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DigitalSetBySTLVector_RECURSES)
#error Recursive header files inclusion detected in DigitalSetBySTLVector.h
#else // defined(DigitalSetBySTLVector_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DigitalSetBySTLVector_RECURSES

#if !defined DigitalSetBySTLVector_h
/** Prevents repeated inclusion of headers. */
#define DigitalSetBySTLVector_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <string>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////
/*#ifdef WITH_VISU3D_QGLVIEWER
#include "DGtal/io/Display3D.h"
#endif*/


namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class DigitalSetBySTLVector
  /**
   * Description of template class 'DigitalSetBySTLVector' <p> \brief
   * Aim: Realizes the concept CDigitalSet by using the STL container
   * std::vector.
   *
   * It thus describes a modifiable set of points within the given
   * domain [Domain].
   *
   * @tparam Domain a realization of the concept CDomain.
   * @see CDigitalSet,CDomain
   */
  template <typename TDomain>
  class DigitalSetBySTLVector
  {
  public:
    typedef TDomain Domain;
    typedef typename Domain::Point Point;
    typedef typename Domain::Size Size;
    typedef typename std::vector<Point>::iterator Iterator;
    typedef typename std::vector<Point>::const_iterator ConstIterator;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~DigitalSetBySTLVector();

    /**
     * Constructor.
     * Creates the empty set in the domain [d].
     *
     * @param d any domain.
     */
    DigitalSetBySTLVector( const Domain & d );

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    DigitalSetBySTLVector ( const DigitalSetBySTLVector & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    DigitalSetBySTLVector & operator= ( const DigitalSetBySTLVector & other );

    /**
     * @return the embedding domain.
     */
    const Domain & domain() const;


    // ----------------------- Standard Set services --------------------------
  public:

    /**
     * @return the number of elements in the set.
     */
    Size size() const;

    /**
     * @return 'true' iff the set is empty (no element).
     */
    bool empty() const;
     
    /**
     * Adds point [p] to this set.
     *
     * @param p any digital point.
     * @pre p should belong to the associated domain.
     */
    void insert( const Point & p );

    /**
     * Adds the collection of points specified by the two iterators to
     * this set.
     *
     * @param first the start point in the collection of Point.
     * @param last the last point in the collection of Point.
     * @pre all points should belong to the associated domain.
     */
    template <typename PointInputIterator>
    void insert( PointInputIterator first, PointInputIterator last );

    /**
     * Adds point [p] to this set if the point is not already in the
     * set. There is no defined behavior if the point is already in
     * the set (for instance, may be present twice).
     *
     * @param p any digital point.
     *
     * @pre p should belong to the associated domain.
     * @pre p should not belong to this.
     */
    void insertNew( const Point & p );

    /**
     * Adds the collection of points specified by the two iterators to
     * this set. The collection should contain distinct points. Each
     * of these points should also not belong already to the set.
     * set. There is no defined behavior if the preceding requisites
     * are not satisfied (for instance, points may be present several
     * times in the set).
     *
     * @param first the start point in the collection of Point.
     * @param last the last point in the collection of Point.
     *
     * @pre all points should belong to the associated domain.
     * @pre each point should not belong to this.
     */
    template <typename PointInputIterator>
    void insertNew( PointInputIterator first, PointInputIterator last );

    /**
     * Removes point [p] from the set.
     * 
     * @param p the point to remove.
     * @return the number of removed elements (0 or 1).
     */
    Size erase( const Point & p );

    /**
     * Removes the point pointed by [it] from the set.
     * 
     * @param it an iterator on this set.
     * @pre it should point on a valid element ( it != end() ).
     * Note: generally faster than giving just the point.
     */
    void erase( Iterator it );

    /**
     * Removes the collection of points specified by the two iterators from
     * this set.
     *
     * @param first the start point in this set.
     * @param last the last point in this set.
     */
    void erase( Iterator first, Iterator last );

    /**
     * Clears the set.
     * @post this set is empty.
     */
    void clear();

    /**
     * @param p any digital point.
     * @return a const iterator pointing on [p] if found, otherwise end().
     */
    ConstIterator find( const Point & p ) const;

    /**
     * @param p any digital point.
     * @return an iterator pointing on [p] if found, otherwise end().
     */
    Iterator find( const Point & p );

    /**
     * @return a const iterator on the first element in this set.
     */
    ConstIterator begin() const;

    /**
     * @return a const iterator on the element after the last in this set.
     */
    ConstIterator end() const;

    /**
     * @return an iterator on the first element in this set.
     */
    Iterator begin();

    /**
     * @return a iterator on the element after the last in this set.
     */
    Iterator end();

    /**
     * set union to left.
     * @param aSet any other set.
     */
    DigitalSetBySTLVector<Domain> & operator+=
    ( const DigitalSetBySTLVector<Domain> & aSet );

    // ----------------------- Other Set services -----------------------------
  public:
    
    /**
     * @return the complement of this set in the domain.
     *
     * NB: be aware of the overhead cost when returning the object.
     */
    DigitalSetBySTLVector<Domain> computeComplement() const; 

    /**
     * Builds the complement in the domain of the set [other_set] in
     * this.
     *
     * @param other_set defines the set whose complement is assigned to 'this'.
     */
    void assignFromComplement( const DigitalSetBySTLVector<Domain> & other_set ); 
    
    /**
     * Computes the bounding box of this set.
     *
     * @param lower the first point of the bounding box (lowest in all
     * directions).
     * @param upper the last point of the bounding box (highest in all
     * directions).
     */
    void computeBoundingBox( Point & lower, Point & upper ) const;

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

    // ------------------------- Protected Datas ------------------------------
  protected:
    
    /**
     * The associated domain.
     */
    const Domain & myDomain;

    /**
     * The container storing the points of the set.
     */
    std::vector<Point> myVector;

    // ------------------------- Private Datas --------------------------------
  private:

    /** 
     * Default Style Functor for selfDraw methods
     * 
     * @param aBoard 
     */
    struct SelfDrawStyle
    {
      SelfDrawStyle(Board2D & aBoard) 
      {
#if(0)
  aBoard.setFillColorRGBi(160,160,160);
  aBoard.setPenColorRGBi(80,80,80);
#endif
      }
    };

  public:
    /** 
     * Default style.
     */
    struct DefaultDrawStyle : public DrawableWithBoard2D
    {
      virtual void selfDraw(Board2D & aBoard) const
      {
#if(0)
  aBoard.setFillColorRGBi(160,160,160);
  aBoard.setPenColorRGBi(80,80,80);
#endif
      }
    };

    // --------------- CDrawableWithBoard2D realization --------------------
  public:

    /**
     * Default drawing style object.
     * @return the dyn. alloc. default style for this object. 
     */
    DrawableWithBoard2D* defaultStyle( std::string mode = "" ) const;

    /**
     * @return the style name used for drawing this object.
     */
    std::string styleName() const;

#if(0)
    /**
     * Draw the object on a Board2D board.
     * @param board the output board where the object is drawn.
     */
    void selfDraw(Board2D & board ) const;
#endif



#if(0)
         /** 
     * Default style.
     */
    struct DefaultDrawStyleDisplay3D : public  DrawableWithDisplay3D 
    {
       virtual void selfDrawDisplay3D(Display3D & display) const
        {
    display.myModes[ "DigitalSetBySTLVector" ] = "";
  }

    };

    /**
     * Default drawing style object.
     * @return the dyn. alloc. default style for this object.
     */
  DrawableWithDisplay3D* defaultStyleDisplay3D( std::string mode = "" ) const;

    /**
     * Draw the object on a Board2D board.
     * @param board the output board where the object is drawn.
     */
    void selfDrawDisplay3D(  Display3D & display ) const;
    void selfDrawAsGridDisplay3D( Display3D & display  ) const;
    void selfDrawAsPavingDisplay3D( Display3D & display ) const;
    void selfDrawAsPavingTransparentDisplay3D( Display3D & display ) const;
#endif




  public:
    
    /**
     * Draw the object on a LiBoard board
     * @param board the output board where the object is drawn.
     * @tparam Functor a Functor to specialize the Board style
     */
    template<typename Functor>
    void selfDraw(Board2D & board ) const;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Default Constructor.
     * Forbidden since a Domain is necessary for defining a set.
     */
    DigitalSetBySTLVector();

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class DigitalSetBySTLVector


  /**
   * Overloads 'operator<<' for displaying objects of class 'DigitalSetBySTLVector'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'DigitalSetBySTLVector' to write.
   * @return the output stream after the writing.
   */
  template <typename Domain>
  std::ostream&
  operator<< ( std::ostream & out, 
         const DigitalSetBySTLVector<Domain> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/kernel/sets/DigitalSetBySTLVector.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DigitalSetBySTLVector_h

#undef DigitalSetBySTLVector_RECURSES
#endif // else defined(DigitalSetBySTLVector_RECURSES)
