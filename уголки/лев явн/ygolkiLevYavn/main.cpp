#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cstring>

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
    //if(x > -1 && x < L+1)
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
    sprintf(str, "LY%d_%d.txt",N,(int)t);
   f = fopen(str,"w");

    for(int i = 0; i < N; i++)
        {fprintf(f,"%f %f\n", i*h, u[i]);}

    fclose(f);
}

double YavniyLeft(double T, double c,double h, double tau, int N)
{
    double* u = new double[N];
    double* x = new double[N];
    for(int i = 0; i < N; i++)
    {
        u[i] = u0(i*h);
        x[i] = i*h;
    }

    double* u_new = new double[N];

   fout(u,N,h,0);
getNormL2("for t=0:",u,h,N,0);
    double t = 0;
    int k = 0;
    while(t < T) {
        t += tau;
        k++;

        for(int i = 1; i < N ; i++)
           {
               u_new[i] = u[i] - c * tau * (u[i] - u[i - 1]) / h;
           }


        u_new[0] = v0(t,c);
        //u_new[N - 1] = vL(t,c, k*h);


        for(int i = 0; i < N; i++)
           {
               u[i] = u_new[i];
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
            cout<<endl;
            cout<<u[N-2]<<endl;
                        cout<<u[N-1]<<endl;
            fout(u,N,h,20);
            getNormL2("for t=20:",u,h,N,20);
        }

    }

    delete [] u;
    delete [] u_new;
}
int main()
{
    double L = 20;
    double T = 20;
    double c = M_PI;

    int N = 1000;
    double h = L / (N - 1);
    double tau = h /c ;
    cout<<h<<endl;
    cout<<tau<<endl;
cout<<(tau+h)<<endl;
    int k;

//if(tau <= h/c)
    { YavniyLeft(T,c,h,tau,N); }
//else {cout << "Viberete drugoe tau"; system("pause");}
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
