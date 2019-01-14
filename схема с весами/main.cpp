#include <iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>

using namespace std;

double lambda(double T) {
	return 5500/(560+T) + 0.942*pow(10.0,-10.0)*T*T*T;
}


int main()
{
    setlocale(LC_ALL, "rus");
	double u0 = 323, mu1 = 373, mu2 = 363, L = 0.5, T = 601, ro = 10950, c = 236;
	double sigma = 0.6;
	double h = 0.01;
    double maxLambda = lambda(323);
	double tau = 0.01;
	int n = (L / h) + 1;
	int stepT=T/(1.0*tau) ;

    char fileName[50];

    if ((tau> h*h/(2*maxLambda*(1-2*sigma))) && sigma<0.5 ) {

        cout << "Ќе выполн€етс€ критерий устойчивости" << endl;
    } else {
        double *uold = new double[n];
        double *u = new double[n];

        for (int i = 0; i < n; i++)
        {
            uold[i] = u0;
            u[i] = u0;
        }

        ofstream fOut("res00000.txt");
        for (int i = 0; i < n; i++)
        {
            fOut << h*i << " " << uold[i]-273 << "\n";
        }
        fOut.close();

            uold[0] = mu1;
            u[0] = mu1;

            double *P = new double[n];
            double *Q = new double[n];

            for (int j=1; j<=stepT; j++)
            {
                P[1] = 0.0;
                Q[1] = mu1;

                for (int i = 2; i<n; i++)
                {
                    uold[i-1]=u[i-1];

                    double gamma = 0.5*(lambda(u[i-2]) + lambda(u[i-1]))/h/h;
                    double alpha = 0.5*(lambda(u[i-2]) + lambda(u[i-1]))/h/h;
                    double beta = alpha+gamma+ro*c/tau;

                    P[i] = gamma / (beta - alpha * P[i-1]);
                    Q[i] = (alpha * Q[i-1] + ro*c/tau*u[i-1]) / (beta - alpha * P[i-1]);
                }

                uold[n-1] = mu2;
                u[n-1] = mu2;

                for (int i = n-2; i > 0; i--)
                {
                    u[i] = P[i+1] * u[i+1] + Q[i+1];
                    u[i]=sigma*u[i]+(1-sigma)*uold[i];
                }
               // cout<<j*tau<<endl;

                if ((abs(j*tau-10.0)<0.5*tau) || (abs(j*tau-60.0)<0.5*tau)  ||  (abs(j*tau-600.0)<0.5*tau)|| (abs(j*tau-300.0)<0.5*tau))
                {
                    sprintf(fileName, "res%05d.txt", (int)(j*tau));

                    ofstream fOut(fileName);
                    for (int i = 0; i < n; i++)
                    {
                        fOut << h*i << " " << u[i]-273 << "\n";
                    }
                    fOut.close();
                }
            }

            delete[] P;
            delete[] Q;
            delete[] uold;
            delete[] u;
    }

    return 0;
}
