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
#include <iostream>
#include <vector>

using namespace std;

// // Function to compute matrix-vector product A*x
// vector<long double> matrix_vector_product(const vector<vector<long double>>& A, const vector<long double>& x)
// {
//     int n = A.size();
//     vector<long double> y(n, 0.0);
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             y[i] += A[i][j] * x[j];
//         }
//     }
//     return y;
// }

// // Function to compute dot product of two vectors
// long double dot_product(const vector<long double>& x, const vector<long double>& y)
// {
//     int n = x.size();
//     double dot = 0.0;
//     for (int i = 0; i < n; i++) {
//         dot += x[i] * y[i];
//     }
//     return dot;
// }

// Conjugate gradient algorithm for solving Ax=b
int conjugate_gradient(vector<vector<long double> >A,  vector<long double> b, vector<long double>& x)
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
                Ap[i] += A[i][j] * p[j];
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
    
    return count;
}

int conjugate_gradient_omp(vector<vector<long double>> A, vector<long double> b, vector<long double>& x)
{
    long double r[N], p[N], Ap[N];
    long double alpha, beta, rr0, rr1;

    // Initialize x and r
    #pragma omp parallel for num_threads(num_threads)
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
        #pragma omp parallel for num_threads(num_threads) for private(j)
        for (int i = 0; i < N; i++)
        {
            Ap[i] = 0.0;
            for (int j = 0; j < N; j++)
            {
                Ap[i] += A[i][j] * p[j];
            }
        }

        // Compute alpha
        alpha = 0.0;
        #pragma omp parallel for num_threads(num_threads) reduction(+:alpha, rr0)
        for (int i = 0; i < N; i++)
        {
            rr0 += r[i] * r[i];
            alpha += p[i] * Ap[i];
        }
        alpha = rr0 / alpha;

        // Update x and r
        rr1 = 0.0;
        #pragma omp parallel for num_threads(num_threads) reduction(+:rr1)
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
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < N; i++)
        {
            p[i] = r[i] + beta * p[i];
        }
    }

    return count;
}

// // Example usage
// int main()
// {
//     // Define a 3x3 matrix A and vectors b and x0
//     vector<vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 3.0, -1.0}, {1.0, -1.0, 2.0}};
//     vector<double> b = {1.0, 2.0, 3.0};
//     vector<double> x0 = {0.0, 0.0, 0.0};

//     // Solve the system Ax=b using the conjugate gradient algorithm
//     double tol = 1e-6;
//     int max_iter = 1000;
//     vector<double> x = conjugate_gradient(A, b, x0, tol, max_iter);

//     // Print the solution
//     cout << "Solution: ";
//     for (int i = 0; i < x.size(); i++) {
//         cout << x[i] << " ";
//     }
//     cout << endl;

//     return 0


int main(void)
{
    // Coefficient matrix
    // std::ios::sync_with_stdio(false);
    // std::cin.tie(NULL);
    srand(time(0));
    cout << "Enter the size of the interior mesh:\n";
    cin >> R;
    N = R * R;
    vector<long double>  p(N*N);
    vector<vector<long double> > A(N,vector<long double>(N,0));
    // vector<vector<long double> > M(N,vector<long double>(N,0));
    vector<long double>  b(N);
    // for(int t=0;t<N;t++){
    //     M[t][t]= rand() % 10;
    // }
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
    int cg_count;
    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }
    cout << endl;
    cg_count = conjugate_gradient(A, b, x);
    
    int cg_count_omp;
    for (int i = 0; i < N; i++)
    {
        // cin>>x[i];
        x[i] = rand() % 10;
        // cin>>demo_3;
        // cout << x[i] << " ";
    }
    cout << endl;
    cg_count_omp = conjugate_gradient_omp(A, b, x);
    

    cout << "+---------+------------------+" << endl
         << "| methode | no. of iteration |" << endl
         << "+---------+------------------+" << endl
         << "| conjugate gradient with openmp | " << cg_count_omp <<"  |"<< endl;
    
    return 0;
}