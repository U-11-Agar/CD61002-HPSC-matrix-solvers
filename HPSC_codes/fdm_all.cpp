#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
// #include </Users/utsavagarwal/Downloads/HPSC.ipynb>
#include <vector>
#include <fstream>
#include <omp.h>

using namespace std;

int R;
const long double epsilon = 1e-6;
int N;

// void fdm_(long double &A[N][N], long double &b[N])
// {
//     // long double A[N][N]; long double b[N];
//     int MESH[R][R];
//     long double newnum;
//     newnum = (long double)R;
//     long double delx;
//     delx = 1.0 / (newnum + 2.0);
//     cout<<delx<<endl;
//     int ii, jj;
//     for (jj = 0; jj < R; jj++)
//     {
//         for (ii = 0; ii < R; ii++)
//         {
//             MESH[ii][jj] = jj * R + ii;
//         }
//     }
//     int n;
//     n = R * R;
//     // N = n;
//     //   long double A[n][n];long double b[n]; long double p[n*n];
//     int i, j, k;
//     for (i = 0; i < n; i++)
//     {
//         for (j = 0; j < n; j++)
//         {
//             A[i][j] = 0;
//         }
//         b[i] = delx * delx;
//     }
//     for (jj = 1; jj < (R - 1); jj++)
//     {
//         for (ii = 1; ii < (R - 1); ii++)
//         {
//             A[MESH[ii][jj]][MESH[ii + 1][jj]] = 1;
//             A[MESH[ii][jj]][MESH[ii - 1][jj]] = 1;
//             A[MESH[ii][jj]][MESH[ii][jj + 1]] = 1;
//             A[MESH[ii][jj]][MESH[ii][jj - 1]] = 1;
//             A[MESH[ii][jj]][MESH[ii][jj]] = -4;
//         }
//     }
//     long double Tb, Tt, Tl, Tr;
//     cout<<"Enter the value of tempertaure at bottom boundary:\n";
//     cin>>Tb;
//     cout<<"Enter the value of tempertaure at top boundary:\n";
//     cin>>Tt;
//     cout<<"Enter the value of tempertaure at left boundary:\n";
//     cin>>Tl;
//     cout<<"Enter the value of tempertaure at right boundary:\n";
//     cin>>Tr;
//     for (jj = 1; jj < R - 1; jj++)
//     {
//         A[MESH[0][jj]][MESH[1][jj]] = 1;
//         A[MESH[0][jj]][MESH[0][jj + 1]] = 1;
//         A[MESH[0][jj]][MESH[0][jj - 1]] = 1;
//         A[MESH[0][jj]][MESH[0][jj]] = -4;
//         b[MESH[0][jj]] = (delx * delx) - Tl;
//     }
//     for (jj = 1; jj < R - 1; jj++)
//     {
//         A[MESH[R - 1][jj]][MESH[R - 2][jj]] = 1;
//         A[MESH[R - 1][jj]][MESH[R - 1][jj + 1]] = 1;
//         A[MESH[R - 1][jj]][MESH[R - 1][jj - 1]] = 1;
//         A[MESH[R - 1][jj]][MESH[R - 1][jj]] = -4;
//         b[MESH[R - 1][jj]] = (delx * delx) - Tr;
//     }
//     for (ii = 1; ii < R - 1; ii++)
//     {
//         A[MESH[ii][0]][MESH[ii + 1][0]] = 1;
//         A[MESH[ii][0]][MESH[ii - 1][0]] = 1;
//         A[MESH[ii][0]][MESH[ii][1]] = 1;
//         A[MESH[ii][0]][MESH[ii][0]] = -4;
//         b[MESH[ii][0]] = (delx * delx) - Tb;
//     }
//     for (ii = 1; ii < R - 1; ii++)
//     {
//         A[MESH[ii][R - 1]][MESH[ii + 1][R - 1]] = 1;
//         A[MESH[ii][R - 1]][MESH[ii - 1][R - 1]] = 1;
//         A[MESH[ii][R - 1]][MESH[ii][R - 2]] = 1;
//         A[MESH[ii][R - 1]][MESH[ii][R - 1]] = -4;
//         b[MESH[ii][R - 1]] = (delx * delx) - Tt;
//     }
//     A[MESH[0][0]][MESH[0][1]] = 1;
//     A[MESH[0][0]][MESH[1][0]] = 1;
//     A[MESH[0][0]][MESH[0][0]] = -4;
//     b[MESH[0][0]] = (delx * delx) - Tl - Tb;
//     A[MESH[R - 1][0]][MESH[R - 2][0]] = 1;
//     A[MESH[R - 1][0]][MESH[R - 1][1]] = 1;
//     A[MESH[R - 1][0]][MESH[R - 1][0]] = -4;
//     b[MESH[R - 1][0]] = (delx * delx) - Tr - Tb;
//     A[MESH[R - 1][R - 1]][MESH[R - 1][R - 2]] = 1;
//     A[MESH[R - 1][R - 1]][MESH[R - 2][R - 1]] = 1;
//     A[MESH[R - 1][R - 1]][MESH[R - 1][R - 1]] = -4;
//     b[MESH[R - 1][R - 1]] = (delx * delx) - Tr - Tt;
//     A[MESH[0][R - 1]][MESH[0][R - 2]] = 1;
//     A[MESH[0][R - 1]][MESH[1][R - 1]] = 1;
//     A[MESH[0][R - 1]][MESH[0][R - 1]] = -4;
//     b[MESH[0][R - 1]] = (delx * delx) - Tl - Tt;
//     // for (i = 0; i < n; i++)
//     // {
//     //     cout<<"\n";
//     //     for (j = 0; j < n; j++)
//     //     {
//     //         cout<<"%Lf  ", A[i][j];
//     //     }
//     //     // cout<<"%Lf ", b[i]);
//     // }
//     // return
// }


int jacobi_solver(vector<vector<long double> > a, vector<long double>  b, vector<long double>  x)
{
    long double demo;
    long double demo_2;
    long double x_n[N];
    long long int count = 1;
    while (count)
    {
        demo = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                {
                    demo += a[i][j] * x[j];
                    // cout<<a[i][n+1];
                }
            }
            x_n[i] = (long double)((b[i] - demo) / a[i][i]);
            // cout << x_n[i] << "+______________________+" << x[i] << endl;
            demo = 0.0;
        }
        // ax=b
        // copy(x_n, x_n + N, x);
        // x=x_n;
        for (int i = 0; i < N; i++)
        {
            demo = x_n[i];
            x[i] = demo;
        }
        demo = 0.0;
        demo_2 = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                demo += x[j] * a[i][j];
            }
            demo = b[i] - demo;
            demo = fabs(demo);
            if (demo > demo_2)
            {
                demo_2 = demo;
            }
            demo = 0.0;
        }
        if (demo_2 < epsilon)
        {
            return count;
        }
        else
        {
            count += 1;
        }
    }
}

int gs_solver(vector<vector<long double> > a, vector<long double>  b, vector<long double>  x)
{
    long double demo;
    long double demo_2;
    long long int count = 1;
    while (count)
    {
        demo = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                {
                    demo += a[i][j] * x[j];
                    // cout<<a[i][n+1];
                }
            }
            x[i] = (long double)((b[i] - demo) / a[i][i]);
            demo = 0.0;
        }
        // ax=b
        // demo = 0.0;
        demo_2 = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                demo += x[j] * a[i][j];
            }
            demo = b[i] - demo;
            demo = fabs(demo);
            if (demo > demo_2)
            {
                demo_2 = demo;
            }
            demo = 0.0;
        }
        if (demo_2 < epsilon)
        {
            return count;
        }
        else
        {
            count += 1;
        }
    }
}

int sor_solver(vector<vector<long double> > a, vector<long double>  b, vector<long double>  x, float omega)
{

    long double demo;
    long double demo_2;
    long long int count = 1;

    while (count)
    {
        demo = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                {
                    demo += a[i][j] * x[j];
                    // cout<<a[i][n+1];
                }
            }
            x[i] = (long double)((omega * ((b[i] - demo) / a[i][i])) + (1 - omega) * x[i]);
            demo = 0.0;
        }
        // ax=b
        // demo = 0.0;
        demo_2 = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                demo += x[j] * a[i][j];
            }
            demo = b[i] - demo;
            demo = fabs(demo);
            if (demo > demo_2)
            {
                demo_2 = demo;
            }
            demo = 0.0;
        }
        if (demo_2 < epsilon)
        {
            return count;
        }
        else
        {
            count += 1;
        }
    }
}

// #include <stdio.h>

int jacobi_solver_omp(long double a[N][N], long double b[N], long double (&x)[N])
{
    long double demo;
    long double demo_2;
    long double x_n[N];
    long long int count = 1;
    // omp_set_num_threads(4);
    const int num_threads = 4;
    while (count)
    {
        demo = 0.0;
#pragma omp parallel for num_threads(num_threads) reduction(+ : demo)
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                {
                    demo += a[i][j] * x[j];
                }
            }
            x_n[i] = (long double)((b[i] - demo) / a[i][i]);
            demo = 0.0;
        }
#pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < N; i++)
        {
            x[i] = x_n[i];
        }
        demo = 0.0;
        demo_2 = 0.0;
#pragma omp parallel for num_threads(num_threads) reduction(max : demo_2)
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                demo += x[j] * a[i][j];
            }
            demo = b[i] - demo;
            demo = fabs(demo);
            if (demo > demo_2)
            {
                demo_2 = demo;
            }
            demo = 0.0;
        }
        if (demo_2 < epsilon)
        {
            return count;
        }
        else
        {
            count += 1;
        }
    }
}

int conjGrad(vector<vector<long double> > a, vector<long double>  b, vector<long double>  x)
{
    long double r[N], p[N], Ap[N];
    long double alpha, beta, rr0, rr1;
    // Initialize x and r
    for (int i = 0; i < N; i++)
    {
        x[i] = 0.0;
        r[i] = b[i];
        p[i] = r[i];
    }
    long long int count = 1;
    // Main iteration loop
    while (count)
    {
        rr0 = rr1 = 0.0;

        // Compute Ap
        for (int i = 0; i < N; i++)
        {
            Ap[i] = 0.0;
            for (int j = 0; j < N; j++)
            {
                Ap[i] += a[i][j] * p[j];
            }
        }

        // Compute alpha
        for (int i = 0; i < N; i++)
        {
            rr0 += r[i] * r[i];
            alpha += p[i] * Ap[i];
        }
        alpha = rr0 / alpha;

        // Update x and r
        for (int i = 0; i < N; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            rr1 += r[i] * r[i];
        }

        // Check for convergence
        if (sqrt(rr1) < epsilon)
        {
            return count;
        }
        else
        {
            count += 1;
        }

        // Compute beta
        beta = rr1 / rr0;

        // Update p
        for (int i = 0; i < N; i++)
        {
            p[i] = r[i] + beta * p[i];
        }
    }
}

// #define BLOCK_SIZE 32

// int gs_solver_omp(long double a[N][N], long double b[N], long double x[N])
// {
//     long double demo;
//     long double demo_2;
//     long long int count = 1;
//     int num_threads = omp_get_max_threads(); // get the maximum number of threads
//     while (count)
//     {
//         demo = 0.0;
//         #pragma omp parallel for num_threads(num_threads) reduction(+:demo)
//         for (int ii = 0; ii < N; ii += BLOCK_SIZE)
//         {
//             for (int jj = 0; jj < N; jj += BLOCK_SIZE)
//             {
//                 for (int i = ii; i < ii + BLOCK_SIZE; i++)
//                 {
//                     for (int j = jj; j < jj + BLOCK_SIZE; j += 4)
//                     {
//                         demo += a[i][j] * x[j];
//                         demo += a[i][j+1] * x[j+1];
//                         demo += a[i][j+2] * x[j+2];
//                         demo += a[i][j+3] * x[j+3];
//                     }
//                     x[i] = (long double)((b[i] - demo) / a[i][i]);
//                     demo = 0.0;
//                 }
//             }
//         }
//         // ax=b
//         // demo = 0.0;
//         demo_2 = 0.0;
//         #pragma omp parallel for num_threads(num_threads) reduction(max:demo_2)
//         for (int ii = 0; ii < N; ii += BLOCK_SIZE)
//         {
//             for (int jj = 0; jj < N; jj += BLOCK_SIZE)
//             {
//                 for (int i = ii; i < ii + BLOCK_SIZE; i++)
//                 {
//                     long double demo = 0.0;
//                     for (int j = jj; j < jj + BLOCK_SIZE; j++)
//                     {
//                         demo += x[j] * a[i][j];
//                     }
//                     demo = b[i] - demo;
//                     demo = fabs(demo);
//                     if (demo > demo_2)
//                     {
//                         demo_2 = demo;
//                     }
//                 }
//             }
//         }
//         if (demo_2 < epsilon)
//         {
//             return count;
//         }
//         else
//         {
//             count += 1;
//         }
//     }
// }

int main(void)
{// Coefficient matrix
    // std::ios::sync_with_stdio(false);
    // std::cin.tie(NULL);
    srand(time(0));
    cout << "Enter the size of the interior mesh:\n";
    cin >> R;
    N = R * R;
    vector<long double>  p(N*N);
    vector<vector<long double> > A(N,vector<long double>(N,0));
    vector<long double>  b(N);
    
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
    vector<long double>  x(N);
    cout << endl
         << "the matrix a=" << endl
         << "{" << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j; j < N; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "the b matrix is" << endl
         << "{" << endl;
    for (int i; i < N; i++)
    {
        cout << b[i] << " ";
    }
    cout << "}" << endl;
    long double x[N];
    int jacobi_count;
    int gs_count;
   
    cout << "initial x vector for jacobi" << endl;
    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }
    cout << endl;
    jacobi_count = jacobi_solver(A, b, x);
    cout << endl
         << "initial x vector for GS (Gauss Seidal)" << endl;
    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }
    cout << endl;
    gs_count = gs_solver(A, b, x);

    cout << "+---------+-----------------+---------------------+------------------+" << endl
         << "| methode | spartial radius | rate of convergence | no. of iteration |" << endl
         << "+---------+-----------------+---------------------+------------------+" << endl
         << "jacobi-null-null-" << jacobi_count << endl
         << "GS-null-null-" << gs_count << endl;
    //  << "sor-null-null-" << endl;
    vector<int> sor;
    int n = 1;
    vector<float> omega;
    omega.resize(n);
    omega[0] = 1;
    int demo = 1;
    while (demo)
    {
        cout << "initial x vector for SOR " << endl;
        for (int i = 0; i < N; i++)
        {
            // cin>>x[i];
            x[i] = rand() % 10;
            // cin>>demo_3;
            // cout << x[i] << " ";
        }
        cout << endl;
        omega.resize(n);
        cout << "enter the value of omega" << endl;
        cin >> omega[n - 1];
        sor.resize(n);
        sor[n - 1] = sor_solver(A, b, x, omega[n - 1]);
        cout << "do you want to try for a different value of omega if yes enter one else enter 0" << endl;
        cin >> demo;
        n += 1;
    }
    cout << "sor-null-null-" << endl;
    for (int i = 0; i < n - 1; i++)
    {
        cout << "omega-" << omega[i] << "-" << sor[i] << endl;
    }

    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }

    cout << endl;
    int jacobi_count_omp;
    jacobi_count_omp = jacobi_solver_omp(A, b, x);
    cout << "\n the iteration for openmp jacobi is " << jacobi_count_omp << endl;
    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }
    vector<vector<long double> > at(N,vector<long double>(N,0));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            at[j][i] = A[i][j];
        }
    }

    vector<vector<long double> >  new_a(N,vector<long double>(N,0));

    // Adding matrix A and its transpose AT
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            new_a[i][j] = A[i][j] + at[i][j];
        }
    }
    cout << endl;
    int conjGrad_count;
    conjGrad_count = conjGrad(new_a, b, x);
    cout << "\n the iteration for CG is " << conjGrad_count << endl;
    // int gs_omp=gs_solver_omp(a, b, x);
    // cout<<"\n the iteration for openmp gs is "<<gs_omp<<endl;
    return 0;
}