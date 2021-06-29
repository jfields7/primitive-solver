#ifndef PS_TYPES_HPP
#define PS_TYPES_HPP

//! \file ps_types.hpp
//  \brief contains some basic type definitions consistent with Athena++.
// 
//  Ideally this file shouldn't be required when the code is dropped into Athena.
//  Therefore, all type definitions should be consistent with Athena.
//

#define NHYDRO 6

using Real = double;
// FIXME: The location of IYF seems bad. Think about it some more and fix it later.
enum ConsIndex {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4, IYD=5};
enum PrimIndex {IVX=1, IVY=2, IVZ=3, IPR=4, ITM=5, IBY=(NHYDRO), IBZ=((NHYDRO)+1), IYF=((NHYDRO)+2)};
enum MagneticIndex {IB1=0, IB2=1, IB3=2};
enum SpatialMetricIndex{S11=0, S12=1, S13=2, S22=3, S23=4, S33=5, NSPMETRIC=6};
enum MetricIndex{I00=0, I01=1, I02=2, I03=3, I11=4, I12=5, I13=6, I22=7, I23=8, I33=9,
                 NMETRIC=10};

#endif
