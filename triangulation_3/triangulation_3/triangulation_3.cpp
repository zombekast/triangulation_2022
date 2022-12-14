#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;
int const n = 5, m = 5, p = 1 + n * m * (m + 1) / 2, N = n * (m*m);
double X[p], Y[p];
int T[N][3];
fstream TRstl;
#define PI 3.14159265358979323846
double pio(double);
double pio(double x)
{
	return cos(2*x) + 3;
}


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

int NumbPoint(int i, int j)
{
	int Nij;
	if (i == 0 && j == 0)
	{
		Nij = 0;
	}
	else
	{
		Nij = n * j * (j - 1) / 2 + i + 1;
	}
	//cout << Nij << endl;
	return Nij;
}

int main()
{
	double alpha = 0, beta = 2 * PI, s, tau, rho, theta, fi; // asin(1)
	int k = 0;
	X[k] = 0;
	Y[k] = 0;
	k++;
	for (double j = 1; j <= m; j++)
	{
		for (double i = 0; i <= n*j - 1; i++)
		{
			tau = j / m;
			fi = (beta * i) / (j * n);
			X[k] = (tau * pio(fi) * cos(fi));
			Y[k] = (tau * pio(fi) * sin(fi));
			k++;
		}
	}
	k = 0;
	for (int i = 1; i < n; i++)
	{
		T[k][0] = NumbPoint(0, 0);
		T[k][1] = NumbPoint(i-1, 1);
		T[k][2] = NumbPoint(i, 1);
		k++;
	}
	T[k][0] = NumbPoint(0, 0);
	T[k][1] = NumbPoint(n - 1, 1);
	T[k][2] = NumbPoint(0, 1);
	k++;
	for (int s = 1; s < n; s++)
	{
		for (int j = 1; j < m; j++)
		{
			for (int i = (s-1)*j; i <= (s*j) ; i++)
			{
				T[k][0] = NumbPoint(i, j);
				T[k][1] = NumbPoint((i + s - 1), (j + 1));
				T[k][2] = NumbPoint((i + s), (j + 1));
				k++;
				//T[k][0] = NumbPoint(i, j);
				//T[k][1] = NumbPoint((i + s), (j + 1));
				//T[k][2] = NumbPoint((i + 1), (j + 1));
				//k++;
			}
		}
	}
	for (int s = n; s <= n; s++)
	{
		for (int j = 1; j < m; j++)
		{
			for (int i = (s - 1) * j; i < (s * j); i++)
			{
				T[k][0] = NumbPoint(i, j);
				T[k][1] = NumbPoint((i + s - 1), (j + 1));
				T[k][2] = NumbPoint((i + s), (j + 1));
				k++;
			}
		}
	}
	for (int s = n; s <= n; s++)
	{
		for (int j = 1; j < m; j++)
		{
			for (int i = (s * j); i <= (s * j); i++)
			{
				T[k][0] = NumbPoint(0, j);
				T[k][1] = NumbPoint(s * (j + 1) - 1, (j + 1));
				T[k][2] = NumbPoint(0, (j + 1));
				k++;
			}
		}
	}
	//for (int s = n; s <= n; s++)
	//{
	//	for (int j = 1; j < m; j++)
	//	{
	//		for (int i = (s * j) - 1; i <= (s * j) - 1; i++)
	//		{
	//			T[k][0] = NumbPoint(s * j - 1, j);
	//			T[k][1] = NumbPoint(s * (j + 1) - 2, (j + 1));
	//			T[k][2] = NumbPoint(s * (j + 1) - 1, (j + 1));
	//			k++;
	//			//T[k][0] = NumbPoint(i, j);
	//			//T[k][1] = NumbPoint((i + s), (j + 1));
	//			//T[k][2] = NumbPoint((i + 1), (j + 1));
	//			//k++;
	//		}
	//	}
	//}
	/*for (int s = 1; s < n; s++)
	{
		for (int j = 1; j < m; j++)
		{
			for (int i = (s * j) - 1; i <= (s * j) - 1; i++)
			{
				T[k][0] = NumbPoint(s * j - 1, j);
				T[k][1] = NumbPoint(s * (j + 1) - 2, (j + 1));
				T[k][2] = NumbPoint(s * (j + 1) - 1, (j + 1));
				k++;
				T[k][0] = NumbPoint(s * j - 1, j);
				T[k][1] = NumbPoint(s * (j + 1) - 1, (j + 1));
				T[k][2] = NumbPoint(0, j);
				k++;
			}
		}
	}
	for (int s = 1; s < n; s++)
	{
		for (int j = 1; j < m; j++)
		{
			for (int i = (s * j); i <= (s * j); i++)
			{
				T[k][0] = NumbPoint(i, j);
				T[k][1] = NumbPoint((i + s - 1), (j + 1));
				T[k][2] = NumbPoint((i + s), (j + 1));
				k++;
				T[k][0] = NumbPoint(i, j);
				T[k][1] = NumbPoint((i + s), (j + 1));
				T[k][2] = NumbPoint((i + 1), (j + 1));
				k++;
			}
		}
	}*/
	/*for (int j = 1; j < m; j++)
	{
		for (int i = (n - 1) * j; i <= (n * j) - 2; i++)
		{
			T[k][0] = NumbPoint(i, j);
			T[k][1] = NumbPoint((i + n - 1), (j + 1));
			T[k][2] = NumbPoint((i + n), (j + 1));
			k++;
			T[k][0] = NumbPoint(i, j);
			T[k][1] = NumbPoint((i + n), (j + 1));
			T[k][2] = NumbPoint((i + 1), j);
			k++;
		}
	}

	for (int j = 1; j < m; j++)
	{
		for (int i = (n * j) - 1; i <= (n * j) - 1; i++)
		{
			T[k][0] = NumbPoint(n * j - 1, j);
			T[k][1] = NumbPoint(n * (j + 1) - 2, (j + 1));
			T[k][2] = NumbPoint(n * (j + 1) - 1, (j + 1));
			k++;
			T[k][0] = NumbPoint(n * j - 1, j);
			T[k][1] = NumbPoint(n * (j + 1) - 1, (j + 1));
			T[k][2] = NumbPoint(0, j);
			k++;
		}
	}


	for (int j = 1; j < m; j++)
	{
		for (int i = (n * j); i <= (n * j); i++)
		{
			T[k][0] = NumbPoint(0, j);
			T[k][1] = NumbPoint(n * (j + 1) - 1, (j + 1));
			T[k][2] = NumbPoint(0, (j + 1));
			k++;
		}
	}*/

	//cout << k << endl;
	/*int i = n - 1;
	for (int j = 0; j < m; j++)
	{
		T[k][0] = (m + 1) * i + j;
		T[k][1] = (m + 1) * 0 + j;
		T[k][2] = (m + 1) * 0 + (j + 1);
		k++;
		T[k][0] = (m + 1) * i + j;
		T[k][1] = (m + 1) * 0 + (j + 1);
		T[k][2] = (m + 1) * i + (j + 1);
		k++;
	}*/
	savetostl();
}