#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>

using namespace std;

int main()
{
    const int n = 1000;        // Matrix size
    const int num_threads = 4; // Number of threads
    vector<vector<double> > A;  // Matrix
    vector<vector<double> > B;  // Matrix
    vector<vector<double> > C;  // Matrix
       // Result vector
    double sum;
    int i, j, k;
    A.resize(n, vector<double>(n, 1));
    B.resize(n, vector<double>(n, 1));
    C.resize(n, vector<double>(n, 0));
    // Initialize matrix and vector
    // for (i = 0; i < n; i++) {
    //     for (j = 0; j < n; j++) {
    //         A[i * n + j] = ;
    //     }
    //     x[i] = i + 1;
    // }
    // for (i = 0; i < n; i++)
    // {
    //     for (
    //         j = 0;
    //         j < n;
    //         j++)
    //     {
    //         cout << A[i][j] << endl;
    //     }
    // }
    // Compute matrix-vector product in parallel
    double start_time = omp_get_wtime();
#pragma omp parallel for private(j, sum) num_threads(num_threads)
    for (int i = 0; i < n; ++i)
    {
        for (int k = 0; k < n; ++k)
        {
            for (int j = 0; j < n; ++j)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
double end_time = omp_get_wtime();
double time_taken = end_time - start_time;

// // Print the result and time taken
// cout << "Result vector:" << endl;
// for (i = 0; i < n; i++) {
//     cout << y[i] << " ";
// }
cout << endl << "Time taken: " << time_taken << " seconds" << endl;

return 0;
}
