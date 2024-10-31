#include <assert.h>
#include "point.h"
#include "smoothing_spline.h"
#include "fstream"

#define OUTPUT_FILE_PATH "../output.txt"
#define VALUES_FILE_PATH "../values.txt"
#define NUMBER_OF_TESTS 1731
#define SMOOTH_COEF 0.99

using namespace std;

void readValues(vector<double> &values) {
    ifstream file(VALUES_FILE_PATH);

    assert(file.is_open());

    for (int i = 0; i < NUMBER_OF_TESTS; i++) {
        double value;
        file >> value;
        values.push_back(value);
    }
    file.close();
}

void writeResults(const SmoothingSpline &spline) {
    ofstream file(OUTPUT_FILE_PATH);
    assert(file.is_open());

    double res[3];

    for (int i = 1; i <= NUMBER_OF_TESTS; i++) {
        spline.getValue(Point(i, 0.0, 0.0), res);
        file << res[0] << endl;
    }

    file.close();
}

int main() {
    vector<Point> Mesh;
    for (int i = 1; i <= NUMBER_OF_TESTS; i++) {
        Mesh.push_back(Point(i, 0.0, 0.0));
    }

    vector<double> values;
    readValues(values);

    SmoothingSpline spline(SMOOTH_COEF);

    spline.updateSpline(Mesh, values);

    writeResults(spline);
}
