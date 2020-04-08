#include "pch.h"
#include <iostream>
#include <fstream>
using namespace std;
// Создание матрицы

double **Create(int n, int m) {
	double **M = new double *[n];
	for (int i = 0; i < n; ++i) {
		M[i] = new double[m];
	}
	return M;
}

// Удаление матрицы

void Free(double **M, int n) {
	for (int i = 0; i < n; ++i) {
		delete[] M[i];
	}
	delete[] M;
}

//ввод матрицы

void Input(double **M, int n, int m) {

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < m; j++) {
			double c = 1.0 + j;

			if (i == j) {
				M[i][j] = 10;
			}
			else if (j > i + 1) {
				M[i][j] = 0;
			}
			else {
				M[i][j] = 1 / c;
			}

		}
	}
}

//вывод матрицы (я его отключил)

void Print(double ** M, int n, int m, double *f) {

	for (int i = 0; i < n; ++i) {

		for (int j = 0; j < m; ++j) {
			if (M[i][j] != 0) {
				//cout << M[i][j] << "*x" << j + 1 << "  +  ";
			}
		}
		f[i] = i + 1;
		//cout << "=" << f[i] << endl;
	}
}


//Метод гауса

double *gauss(double **M, double *f, int n)
{
	double *x, max;
	int k, index;
	const double eps = 0.0001;  // точность
	x = new double[n];
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(M[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(M[i][k]) > max)
			{
				max = abs(M[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			return 0;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = M[k][j];
			M[k][j] = M[index][j];
			M[index][j] = temp;
		}
		double temp = f[k];
		f[k] = f[index];
		f[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = M[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				M[i][j] = M[i][j] / temp;
			f[i] = f[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				M[i][j] = M[i][j] - M[k][j];
			f[i] = f[i] - f[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = f[k];
		for (int i = 0; i < k; i++)
			f[i] = f[i] - M[i][k] * x[k];
	}

	return x;
}


// Метод Якоби

void Jacobi(int n, double** M)
{
	const double eps = 0.0001;
	double* TempX = new double[n];
	double norm; //норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций.


	double* F = new double[n];
	double* X = new double[n];

	for (int i = 0;i < n;i++) {
		F[i] = i + 1;  //Значение функции
		X[i] = 1;  //Начальное приближение
	}


	do {
		for (int i = 0; i < n; i++) {
			TempX[i] = F[i]/10;
			for (int g = 0; g < n; g++) {
				if (i != g)
					TempX[i] -= M[i][g] * X[g];
			}
			TempX[i] /= M[i][i];
		}
		norm = fabs(X[0] - TempX[0]);
		for (int h = 0; h < n; h++) {
			if (fabs(X[h] - TempX[h]) > norm)
				norm = fabs(X[h] - TempX[h]);
			X[h] = TempX[h];
		}
	} while (norm > eps);


	cout << endl << "Metod Jacobi" << endl;
	ofstream out("ans2.dat");
	
	for (int i = 0;i < 100;i++) {
		cout << TempX[i] << endl;
		out	 << TempX[i] << endl;
	
	}
	cout << "norma vectora nevyazki";
		for (int i = 0;i < 100;i++) {
			cout << i + 1 - TempX[i] * 10 << endl;
		}
	out.close();
	delete[] TempX;
}

int main()
{
	double *x;

	//слишком много придется менять если изменить на инт, у меня несколько функций завязано на этих значениях
	int n;
	int m;

	n = 100;
	m = 100;


	double **A = Create(n, m);
	double *f;
	f = new double[n];

	Input(A, n, m);
	Print(A, n, m, f);


	x = gauss(A, f, n);

	cout << endl << "Metod Gaussa" << endl;
	ofstream out("ans1.dat");
		for (int i = 0; i < n; i++) {

			cout << x[i] << endl;
			out << x[i] << endl;
		}
		cout << "norma vectora nevyazki";
		for (int i = 0;i < 100;i++) {
			cout << i + 1 - x[i] * 10 << endl;
		}
	out.close();

	Jacobi(n, A);

	cin.get();


	Free(A, n);


	system("pause");
	return 0;
}
