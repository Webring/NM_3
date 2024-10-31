#pragma once
#ifndef point_h
#define point_h

class Point {
    double X, Y, Z;

public:
    Point(double x = 0, double y = 0, double z = 0) : X(x), Y(y), Z(z) {
    }

    double x() const;

    double y() const;

    double z() const;
};


#endif
