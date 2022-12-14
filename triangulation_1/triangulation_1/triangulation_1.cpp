#include <ctime>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
int const n = 10, m = 50, p = (n + 1) * (m + 1), N = 2 * n * m;
double X[p], Y[p];
int T[N][3];

int savetostl()
{
	int k = 0, ia, ib, ic;
	fstream TRstl;
	TRstl.open("triangulation.stl", ios::out | ios::app);
	TRstl.close();
	TRstl.open("triangulation.stl", ios::out | ios::in);
	TRstl << "solid <Triangulation>\n";
	for (k = 0; k < N; k++)
	{
		ia = T[k][0];
		ib = T[k][1];
		ic = T[k][2];
		TRstl << "facet normal " << 0.0 << " " << 0.0 << " " << 1.0 << "\n";
		TRstl << "outer loop\n";
		TRstl << "vertex ";
		TRstl << X[ia] << " " << Y[ia] << " " << 0.0 << "\n";
		TRstl << "vertex ";
		TRstl << X[ib] << " " << Y[ib] << " " << 0.0 << "\n";
		TRstl << "vertex ";
		TRstl << X[ic] << " " << Y[ic] << " " << 0.0 << "\n";
		TRstl << "endloop\n";
		TRstl << "endfacet\n";
	}
	TRstl << "endsolid";
	TRstl.close();
	return 0;
}

/*int main()
{
double a = 0.0, b = 4.0, c = 1.0, d = 3.0;
int k = 0;
for (int i = 0; i <= n; i++)
{
for (int j = 0; j <= m; j++)
{
X[k] = a + i * (b - a) / n;
Y[k] = c + j * (d - c) / m;
k++;
}
}
k = 0;
for (int i = 0; i < n; i++)
{
for (int j = 0; j < m; j++)
{
T[k][0] = (m + 1) * i + j;
T[k][1] = (m + 1) * (i + 1) + j;
T[k][2] = (m + 1) * i + j + 1;
k++;

T[k][0] = (m + 1) * i + j + 1;
T[k][1] = (m + 1) * (i + 1) + j + 1;
T[k][2] = (m + 1) * (i + 1) + j;
k++;
}
}
double S = 0;
for (k = 0; k < N; k++)
{
S = S + (0.5 * abs((X[T[k][1]] - X[T[k][0]]) * (Y[T[k][2]] - Y[T[k][0]]) - (X[T[k][2]] - X[T[k][0]]) * (Y[T[k][1]] - Y[T[k][0]])));
}
cout « S « endl;
savetostl();
}*/

int main()
{
	double a = 0.0, b = 3.0, c, d;
	int k = 0, h = 0;
	for (double i = 0; i <= n; i++) // Нахождение необходимых значений
	{

		for (double j = 0; j <= m; j++)
		{
			X[k] = a + (b - a) * (i / n);
			d = ((3 * X[k]) - (X[k] * X[k]));
			c = ((X[k] * X[k]) - (3 * X[k]) - 2);
			Y[k] = c + (d - c) * (j / m);
			k++;
		}

	}
	k = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			T[k][0] = (m + 1) * i + j;
			T[k][1] = (m + 1) * (i + 1) + j;
			T[k][2] = (m + 1) * i + j + 1;
			k++;

			T[k][0] = (m + 1) * i + j + 1;
			T[k][1] = (m + 1) * (i + 1) + j + 1;
			T[k][2] = (m + 1) * (i + 1) + j;
			k++;
		}
	}
	double S = 0;
	for (k = 0; k < N; k++) // Функция вычисления площади
	{
		S = S + (0.5 * abs((X[T[k][1]] - X[T[k][0]]) * (Y[T[k][2]] - Y[T[k][0]]) - (X[T[k][2]] - X[T[k][0]]) * (Y[T[k][1]] - Y[T[k][0]])));
	}
	cout << S << endl;
	savetostl();
}