#include <functional>
#include "smoothing_spline.h"
#include <iostream>
#include "cmath"
#include "vector"


double calculateMean(const vector<double> &arr) {
    int n = arr.size();
    if (n == 0) return 0.0;

    double sum = 0.0;
    for (double num: arr) {
        sum += num;
    }

    return sum / n;
}

double calculateSigma(const vector<double> &arr) {
    int n = arr.size();
    if (n == 0) return 0.0;

    double mean = calculateMean(arr), variance = 0.0;

    for (double num: arr) {
        variance += pow(num - mean, 2);
    }
    variance /= n;

    return sqrt(variance);
}


SmoothingSpline::SmoothingSpline(const double &smooth) {
    this->smooth = smooth;
}


void SmoothingSpline::transitionToMasterElement(int segmentIndex, const double &X, double &ksi) const {
    ksi = 2.0 * (X - points[segmentIndex].x()) / (points[segmentIndex + 1].x() - points[segmentIndex].x()) - 1.0;
}

double SmoothingSpline::basisFunction(int number, const double &ksi) const {
    if (number == 1) return 0.5 * (1 - ksi);
    if (number == 2) return 0.5 * (1 + ksi);
    throw runtime_error("Error in the basis function number...");
}


double SmoothingSpline::derBasisFunction(int number) const {
        if (number == 1) return -0.5;
        if (number == 2) return 0.5;
        throw runtime_error("Error in the basis function derivative number...");
}

void SmoothingSpline::updateSpline(const vector<Point> &points, const vector<double> &values) {
    double sigma = calculateSigma(values);
    double mean = calculateMean(values);
    cout << "mean " << mean << endl;
    cout << "sigma " << sigma << endl;


    this->points.clear();
    for (auto &elem: points) this->points.push_back(elem);

    int numberOfSegments = points.size() - 1;

    alpha.resize(numberOfSegments + 1);

    vector<double> a, b, c;
    a.resize(numberOfSegments + 1);
    b.resize(numberOfSegments + 1);
    c.resize(numberOfSegments + 1);


    function<void(int Num_Segment, const Point &point, const double &values, const double &w)>
            Assembling = [&](int i, const Point &point, const double &values, const double &w) {
                double X = point.x(), ksi;
                transitionToMasterElement(i, X, ksi);
                double f1 = basisFunction(1, ksi);
                double f2 = basisFunction(2, ksi);

                b[i] += (1.0 - smooth) * w * f1 * f1;
                b[i + 1] += (1.0 - smooth) * w * f2 * f2;
                a[i + 1] += (1.0 - smooth) * w * f1 * f2;
                c[i] += (1.0 - smooth) * w * f2 * f1;
                alpha[i] += (1.0 - smooth) * w * f1 * values;
                alpha[i + 1] += (1.0 - smooth) * w * f2 * values;
            };

    for (int i = 0; i < numberOfSegments; i++) {
        double W = 1.0;

        // if (577 <= i and i <= 1154) {
        //     W = 0.1;
        // } else if (i > 1154) {
        //     W = 0.01;
        // }

        // if (fabs(mean - values[i]) > fabs(3 * sigma)) {
        //     W = 0.01;
        // }

        Assembling(i, this->points[i], values[i], W);
        Assembling(i, this->points[i + 1], values[i + 1], W);

        double h = points[i + 1].x() - points[i].x();
        b[i] += 1.0 / h * smooth;
        b[i + 1] += 1.0 / h * smooth;
        a[i + 1] -= 1.0 / h * smooth;
        c[i] -= 1.0 / h * smooth;
    }

    for (int j = 1; j < numberOfSegments + 1; j++) {
        b[j] -= a[j] / b[j - 1] * c[j - 1];
        alpha[j] -= a[j] / b[j - 1] * alpha[j - 1];
    }

    alpha[numberOfSegments] /= b[numberOfSegments];
    for (int j = numberOfSegments - 1; j >= 0; j--)
        alpha[j] = (alpha[j] - alpha[j + 1] * c[j]) / b[j];
}

void SmoothingSpline::getValue(const Point &point, double *res) const {
    double eps = 1e-7;
    int numberOfSegments = points.size() - 1;
    double X = point.x();

    for (int i = 0; i < numberOfSegments; i++) {
        if (X > points[i].x() && X < points[i + 1].x() ||
            fabs(X - points[i].x()) < eps ||
            fabs(X - points[i + 1].x()) < eps) {
            double h = points[i + 1].x() - points[i].x();
            double ksi;
            transitionToMasterElement(i, X, ksi);
            res[0] = alpha[i] * basisFunction(1, ksi) +
                     alpha[i + 1] * basisFunction(2, ksi);
            res[1] = (alpha[i] * derBasisFunction(1) +
                      alpha[i + 1] * derBasisFunction(2)) * 2.0 / h;
            res[2] = 0.0;
            return;
        }
    }
    throw runtime_error("The point is not found in the segments...");
}
