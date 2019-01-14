#include <iostream>
#include <cmath>
#include <cstdio>

using namespace std;

double u0(double x)
{
    return cos(x * M_PI);
}
double v0(double t, double c)
{
   return cos( -c * t * M_PI);
}
double vL(double t, double c, double x)
{
    return M_PI * cos(M_PI * (x - c * t));
}
double ExactSolution(double t, double x, double c, double L)
{
    if(x > -1 && x < L+1)
        return u0(x - c*t);

}
double getNormL2(char *s, double *u, double h, int N, double t)
{
    double S = 0;

    for(int i = 0; i < N - 1; i++)
    {
        S+=fabs(u[i]-ExactSolution(t, i*h, M_PI, 20));
    }
    printf("%s : %f\n",s,S/N);
}


void fout( double* u,int N, double h, double t)
{
  char str[20];
 FILE* f;
    sprintf(str, "NC%d_%d.txt",N,(int)t);
   f = fopen(str,"w");

    for(int i = 0; i < N; i++)
        {fprintf(f,"%f %f\n", i*h, u[i]);}

    fclose(f);
}

double NeyavniyCenter(double T, double c,double h, double tau, int N)
{
    double* u = new double[N];
    for(int i = 0; i < N; i++)
        u[i] = u0(i*h);

    fout(u,N,h,0);
getNormL2("for t=0:",u,h,N,0);
    int n = N - 2;
     double** A = new double*[n];
    for(int i = 0; i < n; i++)
    {
        A[i] = new double[3];
    }

    A[0][0] = 1.0;
    A[0][1] = c*tau / (2*h);
    A[0][2] = 0;

    for(int i = 1; i < n - 1; i++)
     {
        A[i][0] = -c*tau / (2*h);
        A[i][1] = 1.0;
        A[i][2] = c*tau / (2*h);
     }

    A[n - 1][0] = -c*tau / (2*h);
    A[n - 1][1] = 1.0;

    double* B = new double[n];

    double* P = new double[n - 1];
    double* Q = new double[n];

    double t = 0;
    int k = 0;
    while(t < T) {
        t += tau;
        k++;

        for(int i = 1; i < n - 1; i++)
        {
            B[i] = u[i + 1];
        }

        B[0] = u[0]*c*tau/(2*h) + u[1];
        B[n - 1] = -u[N - 1]*c*tau/(2*h) + u[N - 2];

        u[0] = v0(t-tau,c);
        //u[N - 1] = vL(t-tau,c, k*h);
        u[N - 1] = u[N - 1] - c * tau * (u[N - 1] - u[N - 2]) / h;

        P[0] = - A[0][1]/A[0][0];
        Q[0] = B[0]/A[0][0];

        for(int i = 1; i < n - 1; i++)
        {
            P[i] = - A[i][2] / (A[i][1] + A[i][0] * P[i - 1]);
            Q[i] = (B[i] - A[i][0] * Q[i - 1]) / (A[i][1] + A[i][0] * P[i - 1]);
        }

        Q[n - 1] = (B[n - 1] - A[n - 1][0] * Q[n - 2]) / (A[n - 1][1] + A[n - 1][0] * P[n - 2]);

        u[n] = Q[n - 1];

        for(int i = n - 1; i >= 1; i--)
        {
            u[i] = P[i - 1] * u[i + 1] + Q[i - 1];
        }


        if(int(5.0/tau) == k)
        {
           fout(u,N,h,5);
 getNormL2("for t=5:",u,h,N,5);
        }

        if(int(10.0/tau) == k)
        {
           fout(u,N,h,10);
getNormL2("for t=10:",u,h,N,10);
        }

        if(int(15.0/tau) == k)
        {
           fout(u,N,h,15);
           getNormL2("for t=15:",u,h,N,15);

        }

        if(int(20.0/tau) == k)
        {
            fout(u,N,h,20);
 getNormL2("for t=20:",u,h,N,20);
        }
    }
   delete [] u;
}


int main()
{
    double L = 20;
    double T = 20;
    double c = M_PI;

    int N = 10000;
    double h = L / (N - 1);
    double tau = h /(c);
cout<<(pow(tau,2)+pow(h,2))<<endl;

//if(tau <= h/c)
    { NeyavniyCenter(T,c,h,tau,N);}
 //else {cout <<"Viberete drugoe tau";}


 char str[20];

    sprintf(str, "L_ExactSol_%d_%d.txt",N,1 );

    FILE* file = fopen(str,"w");

    for(int i = 0; i < N; i++)
        fprintf(file,"%f %f\n", i*h, ExactSolution(0, i*h, c,L));

    fclose(file);

    sprintf(str, "L_ExactSol_%d_%d.txt",N,2);
    file = fopen(str,"w");

    for(int i = 0; i < N; i++)
        fprintf(file,"%f %f\n", i*h, ExactSolution(5, i*h, c,L));

    fclose(file);

    sprintf(str, "L_ExactSol_%d_%d.txt",N,3);
    file = fopen(str,"w");

    for(int i = 0; i < N; i++)
        fprintf(file,"%f %f\n", i*h, ExactSolution(10, i*h, c,L));

    fclose(file);

    sprintf(str, "L_ExactSol_%d_%d.txt",N,4);
    file = fopen(str,"w");

    for(int i = 0; i < N; i++)
        fprintf(file,"%f %f\n", i*h, ExactSolution(15, i*h, c,L));

    fclose(file);

    sprintf(str, "L_ExactSol_%d_%d.txt",N,5);
    file = fopen(str,"w");

    for(int i = 0; i < N; i++)
        fprintf(file,"%f %f\n", i*h, ExactSolution(20, i*h, c,L));

    fclose(file);


    return 0;
}
