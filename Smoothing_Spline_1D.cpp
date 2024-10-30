#include <functional>
#include "Smoothing_Spline_1D.h"
#include <iostream>
#include "cmath"
#include "vector"

double calculateMean(const std::vector<double> &arr) {
    int n = arr.size();
    if (n == 0) return 0.0;

    double sum = 0.0;
    for (double num: arr) {
        sum += num;
    }

    return sum / n;
}

double calculateSigma(const std::vector<double> &arr) {
    int n = arr.size();
    if (n == 0) return 0.0;

    double mean = calculateMean(arr), variance = 0.0;

    // Считаем дисперсию
    for (double num: arr) {
        variance += pow(num - mean, 2);
    }
    variance /= n;

    // Стандартное отклонение
    return sqrt(variance);
}

namespace Com_Methods {
    //конструктор: SMOOTH - параметр сглаживания
    Smoothing_Spline_1D::Smoothing_Spline_1D(const double &SMOOTH) {
        this->SMOOTH = SMOOTH;
    }

    //переход на мастер элемент [-1, 1]:
    //Seg_Num - номер сегмента, Х - абсцисса, Ksi - координата на мастер-элементе
    void Smoothing_Spline_1D::Transition_To_Master_Element(int Seg_Num, const double &X, double &Ksi) const {
        Ksi = 2.0 * (X - Points[Seg_Num].x()) / (Points[Seg_Num + 1].x() - Points[Seg_Num].x()) - 1.0;
    }

    //базисные функции на [-1, 1]:
    //Number - номер функции, Ksi - координата на мастер-элементе
    double Smoothing_Spline_1D::Basis_Function(int Number, const double &Ksi) const {
        switch (Number) {
            case 1: {
                return 0.5 * (1 - Ksi);
                break;
            }
            case 2: {
                return 0.5 * (1 + Ksi);
                break;
            }
            default: {
                throw std::runtime_error("Error in the basis function number...");
                break;
            }
        }
    }

    //производные базисных функций на [-1, 1]:
    //Number - номер функции, Ksi - координата на мастер-элементе
    double Smoothing_Spline_1D::Der_Basis_Function(int Number, const double &Ksi) const {
        switch (Number) {
            case 1: {
                return -0.5;
                break;
            }
            case 2: {
                return 0.5;
                break;
            }
            default: {
                throw std::runtime_error("Error in the basis function derivative number...");
                break;
            }
        }
    }

    //обновить сплайн
    void Smoothing_Spline_1D::Update_Spline(const std::vector<Point> &Points,
                                            const std::vector<double> &F_Value) {
        double sigma = calculateSigma(F_Value);
        double mean = calculateMean(F_Value);
        std::cout << "mean " << mean << std::endl;
        std::cout << "sigma " << sigma << std::endl;


        //обновление списка точек сплайна
        this->Points.clear();
        for (auto &elem: Points) this->Points.push_back(elem);

        //число отрезков разбиения
        int Num_Segments = Points.size() - 1;

        //коэффициенты разложения по базису
        alpha.resize(Num_Segments + 1);

        //диагонали матрицы
        std::vector<double> a, b, c;
        a.resize(Num_Segments + 1);
        b.resize(Num_Segments + 1);
        c.resize(Num_Segments + 1);

        //процедура для ассемблирования СЛАУ:
        //Num_Segment - номер отрезка, P - точка, F_Val - значение табличной функции в точке, w - вес
        std::function<void(int Num_Segment, const Point &P, const double &F_Val, const double &w)>
                Assembling = [&](int i, const Point &P, const double &F_Val, const double &w) {
                    double X = P.x(), Ksi;
                    //переход на мастер-элемент
                    Transition_To_Master_Element(i, X, Ksi);
                    //вычисление значений базисных функций на мастер-элементе
                    double f1 = Basis_Function(1, Ksi);
                    double f2 = Basis_Function(2, Ksi);

                    //внесение вкладов в матрицу СЛАУ
                    b[i] += (1.0 - SMOOTH) * w * f1 * f1;
                    b[i + 1] += (1.0 - SMOOTH) * w * f2 * f2;
                    a[i + 1] += (1.0 - SMOOTH) * w * f1 * f2;
                    c[i] += (1.0 - SMOOTH) * w * f2 * f1;
                    //внесение вкладов в вектор правой части СЛАУ
                    alpha[i] += (1.0 - SMOOTH) * w * f1 * F_Val;
                    alpha[i + 1] += (1.0 - SMOOTH) * w * f2 * F_Val;
                };

        //сборка СЛАУ по сетке: сумма вкладов от каждого сегмента разбиения
        for (int i = 0; i < Num_Segments; i++) {
            //добавление узла сетки в СЛАУ
            double W = 1.0;

            // if (577 <= i and i <= 1154) {
            //     W = 0.1;
            // } else if (i > 1154) {
            //     W = 0.01;
            // }

            // if (fabs(mean - F_Value[i]) > fabs(3 * sigma)) {
            //     W = 0.01;
            // }

            Assembling(i, this->Points[i], F_Value[i], W);
            Assembling(i, this->Points[i + 1], F_Value[i + 1], W);

            //вклад от сглаживания по первой производной
            double h = Points[i + 1].x() - Points[i].x();
            b[i] += 1.0 / h * SMOOTH;
            b[i + 1] += 1.0 / h * SMOOTH;
            a[i + 1] -= 1.0 / h * SMOOTH;
            c[i] -= 1.0 / h * SMOOTH;
        }

        //метод прогонки: прямой ход
        for (int j = 1; j < Num_Segments + 1; j++) {
            //диагональ
            b[j] -= a[j] / b[j - 1] * c[j - 1];
            //правая часть
            alpha[j] -= a[j] / b[j - 1] * alpha[j - 1];
        }

        //метод прогонки: обратный ход
        alpha[Num_Segments] /= b[Num_Segments];
        for (int j = Num_Segments - 1; j >= 0; j--)
            alpha[j] = (alpha[j] - alpha[j + 1] * c[j]) / b[j];
    }

    //вычислить значение сплайна в точке P
    void Smoothing_Spline_1D::Get_Value(const Point &P, double *Res) const {
        //машинный ноль
        double eps = 1e-7;
        //число отрезков
        int Num_Segments = Points.size() - 1;
        //координата точки
        double X = P.x();

        //поиск отрезка, которому принадлежит точка
        for (int i = 0; i < Num_Segments; i++) {
            if (X > Points[i].x() && X < Points[i + 1].x() ||
                fabs(X - Points[i].x()) < eps ||
                fabs(X - Points[i + 1].x()) < eps) {
                //длина отрезка
                double h = Points[i + 1].x() - Points[i].x();
                //переход на местер-элемент, Ksi - координата на мастер-элементе
                double Ksi;
                Transition_To_Master_Element(i, X, Ksi);
                //вычисляем значение сплайна и производной по базисным функциям
                Res[0] = alpha[i] * Basis_Function(1, Ksi) +
                         alpha[i + 1] * Basis_Function(2, Ksi);
                Res[1] = (alpha[i] * Der_Basis_Function(1, Ksi) +
                          alpha[i + 1] * Der_Basis_Function(2, Ksi)) * 2.0 / h;
                Res[2] = 0.0;
                return;
            }
        }
        throw std::runtime_error("The point is not found in the segments...");
    }
}
