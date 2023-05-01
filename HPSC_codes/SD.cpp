#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
// #include </Users/utsavagarwal/Downloads/HPSC.ipynb>
#include <vector>
#include <fstream>
// #include <omp.h>

using namespace std;

int R;
const long double epsilon = 1e-6;
int N;

vector<long double> matrix_mult(vector<vector<long double> > A, vector<long double> x)
{
    int n = A.size();
    vector<long double> y(n, 0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            y[i] += A[i][j] * x[j];
        }
    }
    return y;
}

long double norm(vector<long double> x)
{
    long double sum = 0.0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}

int steepest_descent(vector<vector<long double> > A, vector<long double> b, vector<long double> &x)
{
    int n = A.size();
    vector<long double> r(n, 0);
    vector<long double> d(n, 0);
    vector<long double> Ad(n, 0);

    // Compute r = b - A * x
    for (int i = 0; i < n; i++)
    {
        r[i] = b[i];
        for (int j = 0; j < n; j++)
        {
            r[i] -= A[i][j] * x[j];
        }
    }

    // Initialize d and alpha
    d = r;
    long double alpha = 0.0;

    int count = 1;
    while (count)
    {
        // Compute alpha
        Ad = matrix_mult(A, d);
        long double Ad_dot_d = 0.0;
        for (int i = 0; i < n; i++)
        {
            Ad_dot_d += Ad[i] * d[i];
        }
        long double demo=norm(r);
        alpha = (demo*demo) / Ad_dot_d;

        // Update x and r
        for (int i = 0; i < n; i++)
        {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ad[i];
        }

        long double demo2=norm(r);
        // Compute beta and update d
        long double beta = (demo2*demo2) / (norm(d) * norm(d));
        for (int i = 0; i < n; i++)
        {
            d[i] = r[i] + beta * d[i];
        }

        // Check for convergence
        long double error = norm(r);
        if (error < epsilon)
        {
            return count;
        }
        count++;
    }
    return count;
}





int main(void)
{
    // Coefficient matrix
    // std::ios::sync_with_stdio(false);
    // std::cin.tie(NULL);
    srand(time(0));
    cout << "Enter the size of the interior mesh:\n";
    cin >> R;
    N = R * R;
    vector<long double> p(N * N);
    vector<vector<long double> > A(N, vector<long double>(N, 0));
    vector<long double> b(N);

    FILE *o1, *o2;
    o1 = fopen("Kmat.txt", "r"); // enter the name of the file in with matrix is stored a matrix

    for (int j = 0; j < (N * N); j++)
    {
        fscanf(o1, "%Lf\n", &p[j]);
        // cout << j << "-" << p[j] << endl;
    }
    cout << endl;
    fclose(o1);

    o2 = fopen("Fvec.txt", "r"); // enter the name of the file in with matrix is stored b matrix

    for (int j = 0; j < N; j++)
    {
        fscanf(o2, "%Lf\n", &b[j]);
        // cout<<j<<"-"<<b[j]<<endl;
    }
    fclose(o2);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = p[j * N + i];
        }
    }

    // the left hand side matrix a
    // long double a[N][N] = {{4, -1, -1}, {-2, 6, 1}, {-1, 1, 7}};

    // Right-hand side vector
    // long double b[N] = {676,2,3445,6,24,24534,55,423456,53432,435};

    // fdm_(A,b);
    // cout << endl
    //      << "the matrix a=" << endl
    //      << "{" << endl;
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j; j < N; j++)
    //     {
    //         cout << A[i][j] << " ";
    //     }
    //     // cout << endl;
    // }
    // cout << "}" << endl;
    // cout << "the b matrix is" << endl
    //      << "{" << endl;
    // // for (int i; i < N; i++)
    // // {
    // //     cout << b[i] << " ";
    // // }
    // cout << "}" << endl;
    vector<long double> x(N);
    int sd_count;
    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }
    cout << endl;
    sd_count = steepest_descent(A, b, x);

    cout << "+---------+------------------+" << endl
         << "| methode | no. of iteration |" << endl
         << "+---------+------------------+" << endl
         << "| steepest decent  | " << sd_count << " |" << endl;

    return 0;
}