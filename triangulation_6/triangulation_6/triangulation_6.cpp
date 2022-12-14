#include <iostream>
#include <fstream>
#include <vector>
#include "math.h"
using namespace std;
int const n = 5, m = 5, p = 2 * (1 + n * m * (m + 1) / 2) - n * m, N = 2 * n * m * m;
double X[p], Y[p], Z[p];
int T[N][3];

int NumbPoint(int i, int j, int n)
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
	return Nij;
}

int surfacesavetostl(double X[p], double Y[p], double Z[p])
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
	int i, j, l, k = 0, kt, NSp;
	double R = 2.0, psi_j, rj, phi_ij;
	double PI = 2 * asin(1);
	double U[p], V[p], W[p];
	X[0] = 0.0;
	Y[0] = 0.0;
	Z[0] = R;
	U[0] = 2 * X[0];
	V[0] = Y[0];
	W[0] = Z[0] / 2;
	for (j = 1; j <= m; j++)
	{
		for (i = 0; i <= j * n - 1; i++)
		{
			k++;
			phi_ij = 2 * PI * i / (j * n);
			psi_j = PI * j / (2 * m);
			X[k] = R * sin(psi_j) * cos(phi_ij);
			Y[k] = R * sin(psi_j) * sin(phi_ij);
			Z[k] = R * cos(psi_j);
			U[k] = 2 * X[k];
			V[k] = Y[k];
			W[k] = Z[k] / 2;
			
		}
	}
	k++;
	NSp = k;
	X[k] = 0.0;
	Y[k] = 0.0;
	Z[k] = -R;
	U[k] = 2 * X[k];
	V[k] = Y[k];
	W[k] = Z[k] / 2;
	for (j = 1; j <= m - 1; j++)
	{
		for (i = 0; i <= j * n - 1; i++)
		{
			k++;
			phi_ij = 2 * PI * i / (j * n);
			psi_j = PI * j / (2 * m);
			X[k] = R * sin(psi_j) * cos(phi_ij);
			Y[k] = R * sin(psi_j) * sin(phi_ij);
			Z[k] = -R * cos(psi_j);
			U[k] = 2 * X[k];
			V[k] = Y[k];
			W[k] = Z[k] / 2;
		}
	}
	kt = 0;
	for (i = 1; i <= n - 1; i++)
	{
		T[kt][0] = 0;
		T[kt][1] = i;
		T[kt][2] = i + 1;
		kt++;
	}
	T[kt][0] = 0;
	T[kt][1] = n;
	T[kt][2] = 1;
	kt++;
	for (j = 1; j < m; j++)
	{
		for (l = 1; l <= n - 1; l++)
		{
			for (i = (l - 1) * j; i <= l * j; i++)
			{
				T[kt][0] = NumbPoint(i, j, n);
				T[kt][1] = NumbPoint(i + l - 1, j + 1, n);
				T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n);
				kt++;
				if (i < l * j)
				{
					T[kt][0] = NumbPoint(i, j, n);
					T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n);
					T[kt][2] = NumbPoint(i + 1, j, n);
					kt++;
				}
			}
		}
		//n sector
		l = n;
		for (i = (l - 1) * j; i <= l * j - 2; i++)
		{
			T[kt][0] = NumbPoint(i, j, n);
			T[kt][1] = NumbPoint(i + l - 1, j + 1, n);
			T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n);
			kt++;
			T[kt][0] = NumbPoint(i, j, n);
			T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n);
			T[kt][2] = NumbPoint(i + 1, j, n);
			kt++;
		}
		i = l * j - 1;
		T[kt][0] = NumbPoint(i, j, n);
		T[kt][1] = NumbPoint(i + l - 1, j + 1, n);
		T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n);
		kt++;
		T[kt][0] = NumbPoint(i, j, n);
		T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n);
		T[kt][2] = NumbPoint(0, j, n);
		kt++;
		i = l * j;
		T[kt][0] = NumbPoint(0, j, n);
		T[kt][1] = NumbPoint(i + l - 1, j + 1, n);
		T[kt][2] = NumbPoint(0, j + 1, n);
		kt++;
	}
	//2-nd semisphere
	//NSp=kt;
	for (i = 1; i <= n - 1; i++)
	{
		T[kt][0] = NSp;
		T[kt][2] = NSp + i;
		T[kt][1] = NSp + i + 1;
		kt++;
	}
	T[kt][0] = NSp;
	T[kt][2] = NSp + n;
	T[kt][1] = NSp + 1;
	kt++;
	for (j = 1; j < m - 1; j++)
	{
		for (l = 1; l <= n - 1; l++)
		{
			for (i = (l - 1) * j; i <= l * j; i++)
			{
				T[kt][0] = NumbPoint(i, j, n) + NSp;
				T[kt][2] = NumbPoint(i + l - 1, j + 1, n) + NSp;
				T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n) + NSp;
				kt++;
				if (i < l * j)
				{
					T[kt][0] = NumbPoint(i, j, n) + NSp;
					T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n) + NSp;
					T[kt][1] = NumbPoint(i + 1, j, n) + NSp;
					kt++;
				}
			}
		}
		//n sector
		l = n;
		for (i = (l - 1) * j; i <= l * j - 2; i++)
		{
			T[kt][0] = NumbPoint(i, j, n) + NSp;
			T[kt][2] = NumbPoint(i + l - 1, j + 1, n) + NSp;
			T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n) + NSp;
			kt++;
			T[kt][0] = NumbPoint(i, j, n) + NSp;
			T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n) + NSp;
			T[kt][1] = NumbPoint(i + 1, j, n) + NSp;
			kt++;
		}
		i = l * j - 1;
		T[kt][0] = NumbPoint(i, j, n) + NSp;
		T[kt][2] = NumbPoint(i + l - 1, j + 1, n) + NSp;
		T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n) + NSp;
		kt++;
		T[kt][0] = NumbPoint(i, j, n) + NSp;
		T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n) + NSp;
		T[kt][1] = NumbPoint(0, j, n) + NSp;
		kt++;
		i = l * j;
		T[kt][0] = NumbPoint(0, j, n) + NSp;
		T[kt][2] = NumbPoint(i + l - 1, j + 1, n) + NSp;
		T[kt][1] = NumbPoint(0, j + 1, n) + NSp;
		kt++;
	}
	j = m - 1;
	for (l = 1; l <= n - 1; l++)
	{
		for (i = (l - 1) * j; i <= l * j; i++)
		{
			T[kt][0] = NumbPoint(i, j, n) + NSp;
			T[kt][2] = NumbPoint(i + l - 1, j + 1, n);
			T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n);
			kt++;
			if (i < l * j)
			{
				T[kt][0] = NumbPoint(i, j, n) + NSp;
				T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n);
				T[kt][1] = NumbPoint(i + 1, j, n) + NSp;
				kt++;
			}
		}
	}
	//n sector
	l = n;
	for (i = (l - 1) * j; i <= l * j - 2; i++)
	{
		T[kt][0] = NumbPoint(i, j, n) + NSp;
		T[kt][2] = NumbPoint(i + l - 1, j + 1, n);
		T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n);
		kt++;
		T[kt][0] = NumbPoint(i, j, n) + NSp;
		T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n);
		T[kt][1] = NumbPoint(i + 1, j, n) + NSp;
		kt++;
	}
	i = l * j - 1;
	T[kt][0] = NumbPoint(i, j, n) + NSp;
	T[kt][2] = NumbPoint(i + l - 1, j + 1, n);
	T[kt][1] = NumbPoint(i + 1 + l - 1, j + 1, n);
	kt++;
	T[kt][0] = NumbPoint(i, j, n) + NSp;
	T[kt][2] = NumbPoint(i + 1 + l - 1, j + 1, n);
	T[kt][1] = NumbPoint(0, j, n) + NSp;
	kt++;
	i = l * j;
	T[kt][0] = NumbPoint(0, j, n) + NSp;
	T[kt][2] = NumbPoint(i + l - 1, j + 1, n);
	T[kt][1] = NumbPoint(0, j + 1, n);
	kt++;
	surfacesavetostl(U, V, W);
}
