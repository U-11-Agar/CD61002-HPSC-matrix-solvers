/ Function to compute matrix-vector product A*x
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
