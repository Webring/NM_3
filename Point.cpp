#include "Point.h"

namespace Com_Methods
{
	//конструктор
	Point::Point(double x = 0, double y = 0, double z = 0) : X(x), Y(y), Z(z) {}
	//доступ к полям
	double Point::x()const { return X; }
	double Point::y()const { return Y; }
	double Point::z()const { return Z; }
}
