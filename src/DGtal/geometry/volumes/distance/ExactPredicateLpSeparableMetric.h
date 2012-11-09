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
 * @file ExactPredicateLpSeparableMetric.h
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/11/02
 *
 * Header file for module ExactPredicateLpSeparableMetric.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ExactPredicateLpSeparableMetric_RECURSES)
#error Recursive header files inclusion detected in ExactPredicateLpSeparableMetric.h
#else // defined(ExactPredicateLpSeparableMetric_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ExactPredicateLpSeparableMetric_RECURSES

#if !defined ExactPredicateLpSeparableMetric_h
/** Prevents repeated inclusion of headers. */
#define ExactPredicateLpSeparableMetric_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cmath>
#include "DGtal/base/Common.h"
#include "DGtal/math/BasicMathFunctions.h"
#include "DGtal/base/CQuantity.h"
#include "DGtal/kernel/CSpace.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{
  
/////////////////////////////////////////////////////////////////////////////
// template class ExactPredicateLpSeparableMetric
/**
 * Description of template class 'ExactPredicateLpSeparableMetric' <p>
   * \brief Aim:
   */
  template <typename TSpace, DGtal::uint32_t p,  
            typename TValue=double, typename TPromoted=DGtal::int64_t>
  class ExactPredicateLpSeparableMetric
  {
    // ----------------------- Standard services ------------------------------
  public:
    

    ///Copy the space type
    typedef TSpace Space;
    BOOST_CONCEPT_ASSERT(( CSpace<TSpace> ));

    ///Type for points
    typedef typename Space::Point Point;
    ///Type for points
    typedef typename Point::Coordinate Abscissa;
    ///Type for vectors
    typedef typename Space::Vector Vector;
    
    ///Type for internal distance values
    typedef TPromoted Promoted;
    BOOST_CONCEPT_ASSERT(( CSignedInteger<Promoted> ));
    
    ///Type for distance values
    typedef TValue Value;
    BOOST_CONCEPT_ASSERT(( CQuantity<TValue> ));

    /**
     * Constructor.
     */
    ExactPredicateLpSeparableMetric();

   
    /**
     * Destructor.
     */
    ~ExactPredicateLpSeparableMetric();

  /**
     * Copy constructor.
     * @param other the object to clone.
     */
    ExactPredicateLpSeparableMetric ( const ExactPredicateLpSeparableMetric & other ) {}
    
    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    ExactPredicateLpSeparableMetric & operator= ( const ExactPredicateLpSeparableMetric & other ) { return *this;}

    // ----------------------- Interface --------------------------------------
  public:

    // ----------------------- CLocalMetric --------------------------------------
    /** 
     * Compute the local distance between @a aP and its displacement
     * along the direction @a aDir. 
     * 
     * @param aP a point.
     * @param aDir a direction.
     * 
     * @return the distance between @a aP and @a aP+@a aDir. 
     */
    Value localDistance(const Point & aP, const Vector &aDir) const;

   
    // ----------------------- CMetric --------------------------------------
    /** 
     * Compute the distance between @a aP and @a aQ.
     * 
     * @param aP a first point.
     * @param aQ a second point.
     * 
     * @return the distance between aP and aQ.
     */
    Value distance(const Point & aP, const Point &aQ) const;
    
    /** 
     * Given an origin and two points, this method decides which one
     * is closest to the origin. This method should be faster than
     * comparing distance values.
     * 
     * @param origin the origin
     * @param first  the first point
     * @param second the second point
     * 
     * @return a Closest enum: FIRST, SECOND or BOTH.
     */  
    Closest closest(const Point &origin, 
		    const Point &first,
		    const Point &second) const;
    
      // ----------------------- CSeparableMetric --------------------------------------
    /** 
     * Given three sites (a,b,c) and a straight segment
     * [startingPoint,endPoint] along dimension dim, we detect if the
     * voronoi cells of a and c @e hide the voronoi cell of c on the
     * straight line.
     *
     * @pre both voronoi cells associated with @a a and @a b must
     * intersect the straight line. 
     * 
     * @param u a site
     * @param v a site
     * @param w a site
     * @param startingPoint starting point of the segment
     * @param endPoint end point of the segment
     * @param dim direction of the straight line
     * 
     * @return true if (a,c) hides b.
     */ 
    bool hiddenBy(const Point &u, 
                  const Point &v,
                  const Point &w, 
                  const Point &startingPoint,
                  const Point &endPoint,
                  const typename Point::UnsignedComponent dim) const;


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
  private:

    /** 
     * Compute the Lp distance without the computation of the power
     * 1/p. I.e. only @f$ \sum |p_i- q_i|^p@f$ is given.
     * 
     * @param aP a first point
     * @param aQ a second point
     * 
     * @return the power p of the l_p distance between aP and aQ.
     */    
    Promoted exactDistanceRepresentation(const Point &aP, const Point &aQ) const;

     /** 
     * Perform a binary search on the interval [lower,upper] to
     * detect the mid-point between u and v according to the l_p
     * distance.
     * 
     * @param udim coordinate of u along dimension dim
     * @param vdim coordinate of v along dimension dim
     * @param nu  partial distance of u (sum of |xj-x_i|^p) discarding
     * the term along the dimension dim
     * @param nv partial distance of v (sum of |xj-x_i|^p) discarding
     * the term along the dimension dim
     * @param lower interval lower bound 
     * @param upper interval upper bound
     * 
     * @return the Voronoi boundary point coordinates along dimension dim.
     */
    Abscissa binarySearchHidden(const Abscissa &udim, 
                                const Abscissa &vdim,
                                const Promoted &nu,
                                const Promoted &nv,
                                const Abscissa &lower,
                                const Abscissa &upper) const;

      
    // ------------------------- Private Datas --------------------------------
  private:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ExactPredicateLpSeparableMetric


///////////////////////////////////////////////////////////////////////////////
//            L_2 specialization                                             //
///////////////////////////////////////////////////////////////////////////////


 template <typename TSpace,  
           typename TValue, typename TPromoted>
 class ExactPredicateLpSeparableMetric<TSpace, 2, TValue, TPromoted>
  {
    // ----------------------- Standard services ------------------------------
  public:
    
    ///Copy the space type
    typedef TSpace Space;
    ///Type for points
    typedef typename Space::Point Point;
    ///Type for points
    typedef typename Point::Coordinate Abscissa;
    ///Type for vectors
    typedef typename Space::Vector Vector;
    
    ///Type for internal distance values
    typedef TPromoted Promoted;
    
    ///Type for distance values
    typedef TValue Value;
    BOOST_CONCEPT_ASSERT(( CQuantity<TValue> ));
    
    
    /**
     * Constructor.
     */
    ExactPredicateLpSeparableMetric();

   
    /**
     * Destructor.
     */
    ~ExactPredicateLpSeparableMetric();

  
    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    ExactPredicateLpSeparableMetric ( const ExactPredicateLpSeparableMetric & other ) {}
    
    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    ExactPredicateLpSeparableMetric & operator= ( const ExactPredicateLpSeparableMetric & other ) { return *this;}

   
    // ----------------------- Interface --------------------------------------
  public:

    // ----------------------- CLocalMetric --------------------------------------
    /** 
     * Compute the local distance between @a aP and its displacement
     * along the direction @a aDir. 
     * 
     * @param aP a point.
     * @param aDir a direction.
     * 
     * @return the distance between @a aP and @a aP+@a aDir. 
     */
    Value localDistance(const Point & aP, const Vector &aDir) const;

   
    // ----------------------- CMetric --------------------------------------
    /** 
     * Compute the distance between @a aP and @a aQ.
     * 
     * @param aP a first point.
     * @param aQ a second point.
     * 
     * @return the distance between aP and aQ.
     */
    Value distance(const Point & aP, const Point &aQ) const;
    
    /** 
     * Given an origin and two points, this method decides which one
     * is closest to the origin. This method should be faster than
     * comparing distance values.
     * 
     * @param origin the origin
     * @param first  the first point
     * @param second the second point
     * 
     * @return a Closest enum: FIRST, SECOND or BOTH.
     */  
    Closest closest(const Point &origin, 
		    const Point &first,
		    const Point &second) const;
    

    /** 
     * Compute an exact representation of  Lp distances. We just
     *  return the sum of powers without computing the last 1/p power.
     *  I.e. only @f$ \sum |p_i- q_i|^p@f$ is given.
     * 
     * @param aP a first point
     * @param aQ a second point
     * 
     * @return the power p of the l_p distance between aP and aQ.
     */    
    Promoted exactDistanceRepresentation(const Point &aP, const Point &aQ) const;
    
    // ----------------------- CSeparableMetric --------------------------------------
    /** 
     * Given three sites (a,b,c) and a straight segment
     * [startingPoint,endPoint] along dimension dim, we detect if the
     * voronoi cells of a and c @e hide the voronoi cell of c on the
     * straight line.
     *
     * @pre both voronoi cells associated with @a a and @a b must
     * intersect the straight line. 
     * 
     * @param u a site
     * @param v a site
     * @param w a site
     * @param startingPoint starting point of the segment
     * @param endPoint end point of the segment
     * @param dim direction of the straight line
     * 
     * @return true if (a,c) hides b.
     */ 
    bool hiddenBy(const Point &u, 
                  const Point &v,
                  const Point &w, 
                  const Point &startingPoint,
                  const Point &endPoint,
                  const typename Point::UnsignedComponent dim) const;

   // ----------------------- Other services --------------------------------------
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
  private:

   
     /** 
     * Perform a binary search on the interval [lower,upper] to
     * detect the mid-point between u and v according to the l_p
     * distance.
     * 
     * @param udim coordinate of u along dimension dim
     * @param vdim coordinate of v along dimension dim
     * @param nu  partial distance of u (sum of |xj-x_i|^p) discarding
     * the term along the dimension dim
     * @param nv partial distance of v (sum of |xj-x_i|^p) discarding
     * the term along the dimension dim
     * @param lower interval lower bound 
     * @param upper interval upper bound
     * 
     * @return the Voronoi boundary point coordinates along dimension dim.
     */
    Abscissa binarySearchHidden(const Abscissa &udim, 
                                const Abscissa &vdim,
                                const Promoted &nu,
                                const Promoted &nv,
                                const Abscissa &lower,
                                const Abscissa &upper) const;

      
    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ExactPredicateLpSeparableMetric

  /**
   * Overloads 'operator<<' for displaying objects of class 'ExactPredicateLpSeparableMetric'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ExactPredicateLpSeparableMetric' to write.
   * @return the output stream after the writing.
   */
  template <typename T, DGtal::uint32_t p, typename D, typename P>
  std::ostream&
  operator<< ( std::ostream & out, const ExactPredicateLpSeparableMetric<T,p,D,P> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ExactPredicateLpSeparableMetric_h

#undef ExactPredicateLpSeparableMetric_RECURSES
#endif // else defined(ExactPredicateLpSeparableMetric_RECURSES)
