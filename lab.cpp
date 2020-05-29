//#include "pch.h"
#include <iostream>
#include <fstream>
#include <cmath>
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

			double c = i + 1;
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

void NVN(int n, int z, double *x) {
	
	double **A = Create(n, n);
	Input(A, n, n);

	double* S = new double[n];
	for (int i = 0; i < n; i++) {
		S[i] = 0;
		for (int j = 0; j < n; j++) {

			if (A[i][j] != 0) {
				S[i]+= A[i][j]*x[j];
			}
		}		
		
	}
	for (int i = 0; i < n; i++) {
		 S[i]=i + 1 - S[i];
		 S[i] = S[i] ;
		 //cout << S[i] << endl;;
	}

	double SS=0;
	for (int i = 0; i < n; i++) {
		SS += S[i] * S[i];
	}
	cout << sqrt(SS) << endl;
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
double *Jacobi(int N, double** A, double* F, double* X)
{
	double S = 0;
	double* TempX = new double[N];
	double norm; // норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций.
	const double eps = 0.0001;
	do {
		for (int i = 0; i < N; i++) {
			TempX[i] = F[i];
			for (int g = 0; g < N; g++) {
				if (i != g)
					TempX[i] -= A[i][g] * X[g];
			}
			TempX[i] /= A[i][i];
		}
		norm = fabs(X[0] - TempX[0]);
		for (int h = 0; h < N; h++) {
			if (fabs(X[h] - TempX[h]) > norm)
				norm = fabs(X[h] - TempX[h]);
			X[h] = TempX[h];

		}
		

	} while (norm > eps);
	cout << norm;
	
	delete[] TempX;
	return X;
}


int main()
{
	double *x;
	int n, m;
	n = 100;
	m = 100;


	double **A = Create(n, m);
	double *f;
	f = new double[n];
	double F1 = 0;
	Input(A, n, m);
	Print(A, n, m, f);

	/////////////Metod Gaussa
	x = gauss(A, f, n);

	cout << endl << "Metod Gaussa" << endl;
	ofstream out("ans1.dat");
	for (int i = 0; i < n; i++) {

		//cout << x[i] << endl;
		out << x[i] << endl;
	

	}
	int z = 1;
	NVN(n, z,x);
	out.close();

	





	F1 = 0;
	/////////////Metod Jacobi
	double *x2;
	x2 = new double[n];
	cout << endl << "Metod Jacobi" << endl;
	x2 = Jacobi(n, A, f, x2);
	ofstream out1("ans2.dat");
	for (int i = 0; i < n; i++) {

		//cout << x2[i] << endl;
		out1 << x2[i] << endl;

	}	
	//z = 100000;
	//NVN(n, z,x2);

	out1.close();












	cin.get();
	Free(A, n);
	system("pause");
	return 0;
}
