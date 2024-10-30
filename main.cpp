#include <assert.h>
#include "Point.h"
#include "Smoothing_Spline_1D.h"
#include <iomanip>
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

void writeResults(const Com_Methods::Smoothing_Spline_1D &spline) {
    ofstream file(OUTPUT_FILE_PATH);
    assert(file.is_open());

    double res[3];

    for (int i = 1; i <= NUMBER_OF_TESTS; i++) {
        spline.Get_Value(Com_Methods::Point(i, 0.0, 0.0), res);
        file << res[0] << endl;
    }

    file.close();
}

int main() {
    vector<Com_Methods::Point> Mesh;
    for (int i = 1; i <= NUMBER_OF_TESTS; i++) {
        Mesh.push_back(Com_Methods::Point(i, 0.0, 0.0));
    }

    vector<double> values;
    readValues(values);

    Com_Methods::Smoothing_Spline_1D spline(SMOOTH_COEF);

    spline.Update_Spline(Mesh, values);

    writeResults(spline);
}
