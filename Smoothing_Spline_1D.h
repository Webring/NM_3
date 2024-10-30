#pragma once
#ifndef Smoothing_Spline_1D_h
#define Smoothing_Spline_1D_h

#include "Spline.h"

namespace Com_Methods
{
	//сглаживающий сплайн
	class Smoothing_Spline_1D : public Spline
	{
	private:
		//параметр сглаживания
		double SMOOTH;
		//точки сетки
		std::vector<Point> Points;
		//коэффициенты разложения по базису
		std::vector<double> alpha;
		//переход на мастер элемент [-1, 1]: 
		//Seg_Num - номер сегмента, Х - абсцисса, Ksi - координата на мастер-элементе
		void Transition_To_Master_Element(int Seg_Num, const double &X, double &Ksi) const;
		//базисные функции на [-1, 1]:
		//Number - номер функции, Ksi - координата на мастер-элементе
		double Basis_Function(int Number, const double &Ksi) const;
		//производные базисных функций на [-1, 1]:
		//Number - номер функции, Ksi - координата на мастер-элементе
		double Der_Basis_Function(int Number, const double & Ksi) const;
	public:
		//конструктор: SMOOTH - параметр сглаживания
		Smoothing_Spline_1D(const double &SMOOTH);
		//обновить сплайн
		void Update_Spline(const std::vector<Point> &Points, const std::vector<double> &F_Value) override;
		//вычислить значение сплайна в точке P
		void Get_Value(const Point &P, double * Res)const override;
	};
}

#endif