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
 * @file StdDefsExtern.h
 * @author Pablo Hernandez-Cerdan (pablo.hernandez.cerdan@outlook.com)
 *
 * @date 2020/03/06
 *
 * Declare all the types used in StdDef as extern template (c++11)
 * To avoid recompilation of these templates in every TU.
 * It requires to link with the TU generated from StdDefs.cpp
 *
 * This file is part of the DGtal library.
 */

#include "DGtal/helpers/StdDefs.h"

//////////////////
/// Z2i Explicit instanciation
//////////////////
extern template class DGtal::SpaceND<2, DGtal::Z2i::Integer>;
extern template class DGtal::PointVector<2, DGtal::Z2i::Integer>;
extern template class DGtal::KhalimskySpaceND< 2, DGtal::Z2i::Integer > ;
extern template class DGtal::HyperRectDomain< DGtal::Z2i::Space > ; 
extern template class DGtal::MetricAdjacency< DGtal::Z2i::Space, 1> ;
extern template class DGtal::MetricAdjacency< DGtal::Z2i::Space, 2> ;
extern template class DGtal::DigitalTopology< DGtal::Z2i::Adj4, DGtal::Z2i::Adj8 > ;
extern template class DGtal::DigitalTopology< DGtal::Z2i::Adj8, DGtal::Z2i::Adj4 > ;
extern template class DGtal::DigitalSetByAssociativeContainer<DGtal::Z2i::Domain,
                                                       std::unordered_set< typename DGtal::Z2i::Domain::Point> >;
extern template class DGtal::GridCurve<DGtal::Z2i::K2> ;
extern template class DGtal::ExactPredicateLpSeparableMetric<DGtal::Z2i::Space,2> ;
extern template class DGtal::ExactPredicateLpSeparableMetric<DGtal::Z2i::Space,1> ;
extern template class DGtal::ExactPredicateLpPowerSeparableMetric<DGtal::Z2i::Space,2> ;
extern template class DGtal::ExactPredicateLpPowerSeparableMetric<DGtal::Z2i::Space,1> ;


//////////////////
/// Z3i Explicit instanciation
//////////////////
extern template class DGtal::SpaceND<3, DGtal::Z3i::Integer>;
extern template class DGtal::PointVector<3, DGtal::Z3i::Integer>;
extern template class DGtal::KhalimskySpaceND< 3, DGtal::Z3i::Integer > ;
extern template class DGtal::HyperRectDomain< DGtal::Z3i::Space > ; 
extern template class DGtal::MetricAdjacency< DGtal::Z3i::Space, 1> ;
extern template class DGtal::MetricAdjacency< DGtal::Z3i::Space, 2> ;
extern template class DGtal::MetricAdjacency< DGtal::Z3i::Space, 3> ;
extern template class DGtal::DigitalTopology< DGtal::Z3i::Adj6, DGtal::Z3i::Adj26 > ;
extern template class DGtal::DigitalTopology< DGtal::Z3i::Adj26, DGtal::Z3i::Adj6 > ;
extern template class DGtal::DigitalTopology< DGtal::Z3i::Adj18, DGtal::Z3i::Adj26 > ;
extern template class DGtal::DigitalTopology< DGtal::Z3i::Adj26, DGtal::Z3i::Adj18 > ;
extern template class DGtal::DigitalSetByAssociativeContainer<DGtal::Z3i::Domain,
                                                       std::unordered_set< typename DGtal::Z3i::Domain::Point> >;
extern template class DGtal::ExactPredicateLpSeparableMetric<DGtal::Z3i::Space,2> ;
extern template class DGtal::ExactPredicateLpSeparableMetric<DGtal::Z3i::Space,1> ;
extern template class DGtal::ExactPredicateLpPowerSeparableMetric<DGtal::Z3i::Space,2> ;
extern template class DGtal::ExactPredicateLpPowerSeparableMetric<DGtal::Z3i::Space,1> ;
