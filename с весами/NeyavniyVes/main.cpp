#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include <sstream>

using namespace std;

double X = 0.5;
double T = 1000;
double p = 10950;
double c = 236;
double T0 = 323;
double Th = 373;
double Tc = 363;
double maxL = 6.3;

int N;
double h;

double* uOld;
double* uNew;
double* uLeft;
double* uRight;

double** U;

double t0;

int tmp = 0;

void PrintInFile(double time)
{
	std::ostringstream strs;
	strs << tmp++;
	std::string str = strs.str();
	std::string nameExample = "res_" + str + ".txt";

	char * fileName = new char[nameExample.length() + 1];

	strcpy(fileName, nameExample.c_str());

	FILE* file = fopen(fileName, "w");

	for (int i = 0; i < N; i++)
	{
		fprintf(file, "%0.3f %.3f\n", i * h, uNew[i]);
	}

	fclose(file);
}

double a(double T)
{
	return 5500 / (560 + T) + 0.942 * pow(10, -10) * pow(T, 3);
}

double func(int i, int j)
{
	return 0.0;
}

int main()
{

	double sigma = 0.5;
	cout << "N = "; cin >> N;

	h = X / (N - 1);
	double tau;

	tau = h * h * c / (2 * maxL);

	if (tau > ((h * h) / (4 * maxL * (0.5 - sigma))))
	{
		printf("incorrect tau or h!\n");
	}

	int M = (int)(T / tau) + 1;

	int n = N - 2;

	uOld = new double[N];
	uNew = new double[N];

	uLeft = new double[M];
	uRight = new double[M];

	U = new double*[n];

	for (int i = 0; i < M; i++)
	{
		uLeft[i] = Th;
	}

	for (int i = 0; i < M; i++)
	{
		uRight[i] = Tc;
	}

	for (int i = 0; i < N; i++)
	{
		uOld[i] = T0;
		uNew[i] = T0;
	}

	for (int i = 0; i < n; i++)
	{
		U[i] = new double[3];
	}

	double* B = new double[n];
	double* P = new double[n - 1];
	double* Q = new double[n];

	uOld[0] = Th;
	uOld[N - 1] = Tc;

	uNew[0] = Th;
	uNew[N - 1] = Tc;

	int t = 20;

	for (int j = 0; j < M - 1; j++)
	{
		U[0][0] = -(p * c / tau + (sigma) * (a(uOld[2]) + a(uOld[1])) / (h * h));
		U[0][1] = sigma * a(uOld[2]) / (h * h);
		U[0][2] = 0;

		for (int i = 1; i < n - 1; i++)
		{
			U[i][0] = (sigma) * a(uOld[i + 1]) / (h * h);
			U[i][1] = -(p * c / tau + (sigma) * (a(uOld[i + 2]) + a(uOld[i + 1])) / (h * h));
			U[i][2] = (sigma) * a(uOld[i + 2]) / (h * h);
		}

		U[n - 1][0] = (sigma)* a(uOld[n]) / (h * h);
		U[n - 1][1] = -(p * c / tau + (sigma) * (a(uOld[n + 1]) + a(uOld[n])) / (h * h));

		if (j == 0)
		{
			PrintInFile(0);
		}

		B[0] = -uOld[0] * (sigma) * a(uOld[1]) / (h * h) - p * c * uOld[1] / tau - (1 - sigma) * (a(uOld[2]) * (uOld[2] - uOld[1]) / (h * h) - a(uOld[1]) * (uOld[1] - uOld[0]) / (h * h));

		for (int i = 1; i < n - 1; i++)
		{
			B[i] = -p * c * uOld[i + 1] / tau - (1 - sigma) * (a(uOld[i + 2]) * (uOld[i + 2] - uOld[i+1]) / (h * h) - a(uOld[i+1]) * (uOld[i+1] - uOld[i]) / (h * h));
		}

		B[n - 1] = -uOld[n + 1] * (sigma) * a(uOld[n + 1]) / (h * h) - p * c * uOld[n] / tau - (1 - sigma) * (a(uOld[n + 1]) * (uOld[n + 1] - uOld[n]) / (h * h) - a(uOld[n]) * (uOld[n] - uOld[n - 1]) / (h * h));

		P[0] = -U[0][1] / U[0][0];
		Q[0] = B[0] / U[0][0];

		for (int i = 1; i < n - 1; i++)
		{
			P[i] = -U[i][2] / (U[i][1] + U[i][0] * P[i - 1]);
			Q[i] = (B[i] - U[i][0] * Q[i - 1]) / (U[i][1] + U[i][0] * P[i - 1]);
		}

		Q[n - 1] = (B[n - 1] - U[n - 1][0] * Q[n - 2]) / (U[n - 1][1] + U[n - 1][0] * P[n - 2]);

		uNew[n] = Q[n - 1];

		for (int i = n - 1; i >= 1; i--)
		{
			uNew[i] = P[i - 1] * uNew[i + 1] + Q[i - 1];
		}

		for (int i = 0; i < N; i++)
		{
			uOld[i] = uNew[i];
		}

		if (j == int(t / tau))
		{
			PrintInFile(int(t));
			t += 20;
		}
	}

	return 0;
}
