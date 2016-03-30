#pragma once
#include <cmath>

typedef struct
{
	double p[3];
} Point;

static inline
double distance(const Point &_v, const Point &_w)
{
	return std::sqrt(std::pow(_v.p[0] - _w.p[0], 2) +
		std::pow(_v.p[1] - _w.p[1], 2) +
		std::pow(_v.p[2] - _w.p[2], 2));
};

static inline
double dot_product(const Point &_v, const Point &_w)
{
	return _v.p[0] * _w.p[0] + _v.p[1] * _w.p[1] + _v.p[2] * _w.p[2];
};

static inline
Point operator+ (const Point &_v, const Point &_w)
{
	Point nrvo;
	nrvo.p[0] = _v.p[0] + _w.p[0];
	nrvo.p[1] = _v.p[1] + _w.p[1];
	nrvo.p[2] = _v.p[2] + _w.p[2];
	return nrvo;
};

static inline
Point operator* (double _a, const Point &_w)
{
	Point nrvo;
	nrvo.p[0] = _a*_w.p[0];
	nrvo.p[1] = _a*_w.p[1];
	nrvo.p[2] = _a*_w.p[2];
	return nrvo;
};

static inline
Point operator- (const Point &_v, const Point &_w)
{
	Point nrvo;
	nrvo.p[0] = _v.p[0] - _w.p[0];
	nrvo.p[1] = _v.p[1] - _w.p[1];
	nrvo.p[2] = _v.p[2] - _w.p[2];
	return nrvo;
};

static inline
double norm(const Point &_v)
{
	return std::sqrt(_v.p[0] * _v.p[0] + _v.p[1] * _v.p[1] + _v.p[2] * _v.p[2]);
};

