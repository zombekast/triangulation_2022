#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;
int const n = 25, m = 10, p = n * (m + 1), N = 2 * n * m;
double X[p], Y[p], Z[p];
int T[N][3];
fstream TRstl;
#define PI 3.14159265358979323846
double pio(double);
double pio(double x)
{
	return cos(2*x)+2;
}

	int surfacesavetostl()
	{
		int k = 0, ia, ib, ic;
		double v, x1, y1, z1, x2, y2, z2, nx, ny, nz;
		fstream TRstl;
		remove("triangulation.stl");
		TRstl.open("triangulation.stl", ios::out | ios::app);
		TRstl << "solid <Triangulation>\n";
		for (k = 0; k < N; k++)
		{
			ia = T[k][0];
			ib = T[k][1];
			ic = T[k][2];
			x1 = X[ib] - X[ia];
			y1 = Y[ib] - Y[ia];
			z1 = Z[ib] - Z[ia];
			x2 = X[ic] - X[ia];
			y2 = Y[ic] - Y[ia];
			z2 = Z[ic] - Z[ia];
			nx = (y1 * z2 - y2 * z1);
			ny = (z1 * x2 - x1 * z2);
			nz = (x1 * y2 - x2 * y1);
			v = sqrt(nx * nx + ny * ny + nz * nz);
			nx = nx / v;
			ny = ny / v;
			nz = nz / v;
			TRstl << "facet normal " << nx << " " << ny << " " << nz << "\n";
			TRstl << "outer loop\n";
			TRstl << "vertex ";
			TRstl << X[ia] << " " << Y[ia] << " " << Z[ia] << "\n";
			TRstl << "vertex ";
			TRstl << X[ib] << " " << Y[ib] << " " << Z[ib] << "\n";
			TRstl << "vertex ";
			TRstl << X[ic] << " " << Y[ic] << " " << Z[ic] << "\n";
			TRstl << "endloop\n";
			TRstl << "endfacet\n";
		}
		TRstl << "endsolid";
		TRstl.close();
		return 0;
	}

int main()
{
	double alpha = 0, beta = 2 * PI, s, tau, rho, theta, H, h; // asin(1)
	int k = 0;
	H = 10;
	for (int i = 0; i <= n - 1; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			theta = i * (2 * PI) / n;
			X[k] = (pio(theta) * cos(theta));
			Y[k] = (pio(theta) * sin(theta));
			h = j * H / m;
			Z[k] = h;
			k++;
		}
	}

	k = 0;
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = 0; j < m; j++)
		{
			T[k][0] = (m + 1) * i + j;
			T[k][1] = (m + 1) * (i + 1) + j;
			T[k][2] = (m + 1) * i + (j + 1);
			k++;
			T[k][0] = (m + 1) * i + (j + 1);;
			T[k][1] = (m + 1) * (i + 1) + j;
			T[k][2] = (m + 1) * (i + 1) + (j + 1);;
			k++;
		}

	}
	int i = n - 1;
	for (int j = 0; j < m; j++)
	{
		T[k][0] = (m + 1) * i + j;
		T[k][1] = (m + 1) * 0 + j;
		T[k][2] = (m + 1) * i + (j + 1);
		k++;
		T[k][0] = (m + 1) * i + (j + 1);
		T[k][1] = (m + 1) * 0 + j;
		T[k][2] = (m + 1) * 0 + (j + 1);
		k++;
	}
	surfacesavetostl();
}