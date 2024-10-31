#pragma once
#ifndef Smoothing_Spline_1D_h
#define Smoothing_Spline_1D_h

#include "point.h"
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
