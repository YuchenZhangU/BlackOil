#include "Point.hpp"
#include <cmath>

double distance( const Point &_v, const Point &_w )
{
   return std::sqrt( std::pow(_v.p[0] - _w.p[0],2)+
		     std::pow(_v.p[1] - _w.p[1],2)+
		     std::pow(_v.p[2] - _w.p[2],2) );
};

double dot_product( const Point &_v, const Point &_w )
{
   return _v.p[0]*_w.p[0] + _v.p[1]*_w.p[1] + _v.p[2]*_w.p[2];
};

Point operator+ ( const Point &_v, const Point &_w )
{
   Point nrvo;
   nrvo.p[0] = _v.p[0]+_w.p[0];
   nrvo.p[1] = _v.p[1]+_w.p[1];
   nrvo.p[2] = _v.p[2]+_w.p[2];
   return nrvo;
};

Point operator- ( const Point &_v, const Point &_w )
{
   Point nrvo;
   nrvo.p[0]=_v.p[0]-_w.p[0];
   nrvo.p[1]=_v.p[1]-_w.p[1];
   nrvo.p[2]=_v.p[2]-_w.p[2]; 
   return nrvo;
};

double norm( const Point &_v )
{
   return std::sqrt(_v.p[0]*_v.p[0] + _v.p[1]*_v.p[1] + _v.p[2]*_v.p[2] );
};

