#pragma once
#ifndef smoothing_spline_h
#define smoothing_spline_h

#include "point/point.h"
#include <vector>

using namespace std;

class SmoothingSpline {
    double smooth;

    vector<Point> points;

    vector<double> alpha;

    void transitionToMasterElement(int segNum, const double &X, double &ksi) const;

    double basisFunction(int number, const double &ksi) const;

    double derBasisFunction(int number) const;

public:
    SmoothingSpline(const double &smooth);

    void updateSpline(const vector<Point> &points, const vector<double> &fValue);

    void getValue(const Point &point, double *res) const;
};

#endif
