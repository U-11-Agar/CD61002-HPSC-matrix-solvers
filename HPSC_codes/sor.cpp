#include <iostream>
#include <cmath>
// #include <stdio.h>
#include <stdlib.h>
// #include </Users/utsavagarwal/Downloads/HPSC.ipynb>
#include <vector>
#include <fstream>
// #include <omp.h>
using namespace std;

int R;
const long double epsilon = 1e-6;

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

int N;

int sor_solver(vector<vector<long double> > a, vector<long double> b, vector<long double> x, float omega)
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
    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }
    cout << endl;

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

    cout << "+---------+------------------+" << endl
         << "| methode | no. of iteration |" << endl
         << "+---------+------------------+" << endl
         << "| SOR | ->" << endl;
    for (int i = 0; i < n - 1; i++)
    {
        cout << "omega=" << omega[i] << "-" << sor[i] << endl;
    }

    return 0;
}