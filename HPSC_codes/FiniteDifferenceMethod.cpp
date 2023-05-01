#include<iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include<vector>
using namespace std;
int main()
{
  int R;
  printf("Enter the size of the interior mesh:\n");
  scanf("%d", &R);
  vector<vector<long double> >  MESH(R,vector<long double>(R,0));
  double newnum;
  newnum = (double)R;
  long double delx;
  delx = 1.0 / (newnum + 2.0);
  printf("%Lf", delx);
  int ii, jj;
  for (jj = 0; jj < R; jj++)
  {
    for (ii = 0; ii < R; ii++)
    {
      MESH[ii][jj] = jj * R + ii;
    }
  }
  int n;
  n = R * R;
  vector<vector<long double> > A(n,vector<long double>(n,0));
  vector<long double> b(n);
  vector<long double>p(n*n);
  int i, j, k;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      A[i][j] = 0;
    }
    b[i] = delx * delx;
  }
  for (jj = 1; jj < (R - 1); jj++)
  {
    for (ii = 1; ii < (R - 1); ii++)
    {
      A[MESH[ii][jj]][MESH[ii + 1][jj]] = 1;
      A[MESH[ii][jj]][MESH[ii - 1][jj]] = 1;
      A[MESH[ii][jj]][MESH[ii][jj + 1]] = 1;
      A[MESH[ii][jj]][MESH[ii][jj - 1]] = 1;
      A[MESH[ii][jj]][MESH[ii][jj]] = -4;
    }
  }
  double Tb, Tt, Tl, Tr;
  printf("Enter the value of tempertaure at bottom boundary:\n");
  scanf("%lf", &Tb);
  printf("Enter the value of tempertaure at top boundary:\n");
  scanf("%lf", &Tt);
  printf("Enter the value of tempertaure at left boundary:\n");
  scanf("%lf", &Tl);
  printf("Enter the value of tempertaure at right boundary:\n");
  scanf("%lf", &Tr);

  for (jj = 1; jj < R - 1; jj++)
  {
    A[MESH[0][jj]][MESH[1][jj]] = 1;
    A[MESH[0][jj]][MESH[0][jj + 1]] = 1;
    A[MESH[0][jj]][MESH[0][jj - 1]] = 1;
    A[MESH[0][jj]][MESH[0][jj]] = -4;
    b[MESH[0][jj]] = (delx * delx) - Tl;
  }
  for (jj = 1; jj < R - 1; jj++)
  {
    A[MESH[R - 1][jj]][MESH[R - 2][jj]] = 1;
    A[MESH[R - 1][jj]][MESH[R - 1][jj + 1]] = 1;
    A[MESH[R - 1][jj]][MESH[R - 1][jj - 1]] = 1;
    A[MESH[R - 1][jj]][MESH[R - 1][jj]] = -4;
    b[MESH[R - 1][jj]] = (delx * delx) - Tr;
  }
  for (ii = 1; ii < R - 1; ii++)
  {
    A[MESH[ii][0]][MESH[ii + 1][0]] = 1;
    A[MESH[ii][0]][MESH[ii - 1][0]] = 1;
    A[MESH[ii][0]][MESH[ii][1]] = 1;
    A[MESH[ii][0]][MESH[ii][0]] = -4;
    b[MESH[ii][0]] = (delx * delx) - Tb;
  }
  for (ii = 1; ii < R - 1; ii++)
  {
    A[MESH[ii][R - 1]][MESH[ii + 1][R - 1]] = 1;
    A[MESH[ii][R - 1]][MESH[ii - 1][R - 1]] = 1;
    A[MESH[ii][R - 1]][MESH[ii][R - 2]] = 1;
    A[MESH[ii][R - 1]][MESH[ii][R - 1]] = -4;
    b[MESH[ii][R - 1]] = (delx * delx) - Tt;
  }
  A[MESH[0][0]][MESH[0][1]] = 1;
  A[MESH[0][0]][MESH[1][0]] = 1;
  A[MESH[0][0]][MESH[0][0]] = -4;
  b[MESH[0][0]] = (delx * delx) - Tl - Tb;
  A[MESH[R - 1][0]][MESH[R - 2][0]] = 1;
  A[MESH[R - 1][0]][MESH[R - 1][1]] = 1;
  A[MESH[R - 1][0]][MESH[R - 1][0]] = -4;
  b[MESH[R - 1][0]] = (delx * delx) - Tr - Tb;
  A[MESH[R - 1][R - 1]][MESH[R - 1][R - 2]] = 1;
  A[MESH[R - 1][R - 1]][MESH[R - 2][R - 1]] = 1;
  A[MESH[R - 1][R - 1]][MESH[R - 1][R - 1]] = -4;
  b[MESH[R - 1][R - 1]] = (delx * delx) - Tr - Tt;
  A[MESH[0][R - 1]][MESH[0][R - 2]] = 1;
  A[MESH[0][R - 1]][MESH[1][R - 1]] = 1;
  A[MESH[0][R - 1]][MESH[0][R - 1]] = -4;
  b[MESH[0][R - 1]] = (delx * delx) - Tl - Tt;

  // for(i=0;i<n;i++){
  //    printf("\n");
  //    for(j=0;j<n;j++){
  //      printf("%lf  ",A[i][j]);
  //    }
  //    printf("%lf ",b[i]);
  //  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      p[j * n + i] = A[i][j];
    }
  }

  FILE *o1, *o2;
  o1 = fopen("Kmat.txt", "w");

  for (j = 0; j < n * n; j++)
  {
    fprintf(o1, "%Lf\n", p[j]);
  }

  fclose(o1);

  o2 = fopen("Fvec.txt", "w");

  for (j = 0; j < n; j++)
  {
    fprintf(o2, "%Lf\n", b[j]);
  }
  fclose(o2);
}
