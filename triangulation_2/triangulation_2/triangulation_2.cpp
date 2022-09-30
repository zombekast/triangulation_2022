#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;
/*vector <double> X, Y;
vector <int> TTA, TTB, TTC;*/
int const n = 10, m = 20, p = (n + 1) * (m + 1), N = 2 * n * m;
double X[p], Y[p];
int T[N][3];
fstream TRstl;
#define PI 3.14159265358979323846
double psi(double);
double psi(double x)
{
	return 2.0;
}
double phi(double);
double phi(double x)
{
	return 1.0;
}
/*int savetostl()
{
	int k = 0, ia, ib, ic;
	fstream TRstl;
	TRstl.open("triangulation.stl", ios::out | ios::in);
	TRstl.close();
	TRstl.open("triangulation.stl", ios::out | ios::app);
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
}*/

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

int main()
{
	double alpha = 0, beta = 2*PI, h, tau, rho, theta; // asin(1)
	int k = 0;
	for (int i = 0; i <= n -1; i++)
	{
		for (double j = 0; j <= m; j++)
		{
			tau = j / m;
			theta = alpha + i * (beta - alpha) / n;
			rho = tau * psi(theta) + (1 - tau) * phi(theta);
			//rho = phi(alpha + i * h) + j * tau * (psi(alpha + i * h) - phi(alpha + i * h));
			X[k] = (rho * cos(theta));
			Y[k] = (rho * sin(theta));
			
			//cout << X[k] << "    " << Y[k] << endl;
			k++;
		}
		//cout << X[i] << "    " << Y[i] << endl;
	}
	k = 0;
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = 0; j < m; j++)
		{
			T[k][0] = (m + 1) * i + j;
			T[k][1] = (m + 1) * (i + 1) + j;
			T[k][2] = (m + 1) * (i + 1) + (j + 1);
			k++;
			T[k][0] = (m + 1) * i + j;
			T[k][1] = (m + 1) * (i + 1) + (j + 1);
			T[k][2] = (m + 1) * i + (j + 1);
			k++;
		}
	
	}
	int i = n - 1;
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
	}
	savetostl();
}