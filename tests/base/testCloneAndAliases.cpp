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
 * @file testCloneAndAliases.cpp
 * @ingroup Tests
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 *
 * @date 2012/07/02
 *
 * This file is part of the DGtal library
 */

//#define TRACE_BITS

#include <cstdio>
#include <cmath>
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/Clone.h"
#include "DGtal/base/Alias.h"
#include "DGtal/base/ConstAlias.h"

using namespace DGtal;
using namespace std;

// Dummy class to test clones and aliases.
struct A1
{
private:
  A1();
public:
  ~A1()
  {
    std::cout << "  ~A1() " << std::endl;
    ++nbDeleted;
  }
  A1( int i ) : data( i )
  {
    std::cout << "  A1( int i ) " << std::endl;
    ++nbCreated;
  }
  A1( const A1 & a ) : data( a.data )
  {
    std::cout << "  A1( const A1 & a ) " << std::endl;
    ++nbCreated;
  }
  A1& operator=( const A1 & a ) 
  {
    data = a.data;
    std::cout << "  A1::operator=( const A1 & a ) " << std::endl;
    return *this;
  }
  static void reset() {
    nbCreated = 0; 
    nbDeleted = 0;
  }

  int data;

  static int nbCreated;
  static int nbDeleted;
};

int A1::nbCreated = 0;
int A1::nbDeleted = 0;

// This class uses standard by-value parameter passing 
// The data member is an instance of A1.
struct DByValue {
  DByValue( A1 a )
    : myA1( a )
  {
    std::cout << "  DByValue( A1 a) " << std::endl;
  }

  int value() const
  {
    return myA1.data;
  }

  A1 myA1;
};

// This class uses the provided explicit by-value parameter passing (with Clone).
// The data member is an instance of A1.
struct DByClone {
  DByClone( Clone<A1> a )
    : myA1( a )
  {
    std::cout << "  DByClone( Clone<A1> a) " << std::endl;
  }

  int value() const
  {
    return myA1.data;
  }

  A1 myA1;
};

// This class uses the provided explicit by-reference parameter passing (with Alias).
// The data member is a pointer to the given instance.
struct EByAlias {
  EByAlias( Alias<A1> a )
    : myA1( a )
  {
    std::cout << "  EByAlias( A1lias<A1> a ) " << myA1 << std::endl;
  }

  int value() const
  {
    return myA1->data;
  }
  
  A1* myA1;
  // ou
  // A1lias<A1> myA1;
};

// This class uses the provided explicit by-reference parameter passing (with Alias).
// The data member is a const pointer to the given instance.
struct EByConstAlias {
  EByConstAlias( ConstAlias<A1> a )
    : myA1( a )
  {
    std::cout << "  EByConstAlias( A1lias<A1> a ) " << myA1 << std::endl;
  }

  int value() const
  {
    return myA1->data;
  }
  
  const A1* myA1;
  // ou
  // A1lias<A1> myA1;
};

// This class uses the provided explicit by-reference parameter passing (with Alias).
// The data member is an Alias of the given instance (should be same behaviour).
struct F {
  F( Alias<A1> a )
    : myA1( a )
  {
    std::cout << "  F( A1lias<A1> a ) " << myA1 << std::endl;
  }

  int value() const
  {
    A1 & a = myA1; // in fact, necessary
    return a.data; 
  }
  
  Alias<A1> myA1;
};


int main()
{
  unsigned int nb = 0;
  unsigned int nbok = 0;
  trace.beginBlock ( "Number of A1 instances with standard by-value parameter passing." );
  A1 a1( 10 ); // +1/ 0
  DByValue d1( a1 ); //  +2/+1
  trace.info() << "D: d1.value() = " << d1.value() << std::endl;
  ++nb, nbok += A1::nbCreated==3 ? 1 : 0;
  ++nb, nbok += A1::nbDeleted==1 ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ")"
               << " nbCreated=" << A1::nbCreated 
               << " nbDeleted=" << A1::nbDeleted << std::endl; 
  trace.endBlock();
  trace.beginBlock ( "Number of A1 instances with explicit by-value parameter passing (Clone)." );
  DByClone dd1( a1 ); // +1/0
  ++nb, nbok += A1::nbCreated==4 ? 1 : 0;
  ++nb, nbok += A1::nbDeleted==1 ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ")"
               << " nbCreated=" << A1::nbCreated 
               << " nbDeleted=" << A1::nbDeleted << std::endl; 
  trace.endBlock();
  trace.beginBlock ( "Number of A1 instances with explicit by-reference parameter passing (Alias)." );
  EByAlias e1( a1 ); // +0/0
  ++nb, nbok += A1::nbCreated==4 ? 1 : 0;
  ++nb, nbok += A1::nbDeleted==1 ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ")"
               << " nbCreated=" << A1::nbCreated 
               << " nbDeleted=" << A1::nbDeleted << std::endl; 
  trace.endBlock();
  trace.beginBlock ( "Number of A1 instances with explicit by-const reference parameter passing (Alias)." );
  EByConstAlias ee1( a1 ); // +0/0
  ++nb, nbok += A1::nbCreated==4 ? 1 : 0;
  ++nb, nbok += A1::nbDeleted==1 ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ")"
               << " nbCreated=" << A1::nbCreated 
               << " nbDeleted=" << A1::nbDeleted << std::endl; 
  trace.endBlock();


  // These two lines should not compile.
  // Clone<A1> clone1( a1 );
  // Clone<A1> clone2( clone2 );  

  return ( nb == nbok ) ? 0 : 1;
}
/** @ingroup Tests **/
