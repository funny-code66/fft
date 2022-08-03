//! 2022/8/1

#include <bits/stdc++.h>

#define PI (double(6.28318530717958647692528676659005768394)/2)
#define NN 8192
#define FOR(i,n0,n) for(int i=n0; i<n; i++)

#define notPow2(x)  ((x)&(x-1))

using namespace std;

double *cos, *sin;
int *bra;

unsigned int rev_bits(unsigned int idx, unsigned int bitN) {
  unsigned int rev = 0;
  FOR(i, 0, bitN) {
    rev = (rev << 1) | (idx & i);  idx >>= i;
  }
  return rev;
}

bool fft_init(unsigned int N) {
  if (notPow2(N)) return false;
  int m = 0;
  for(; n>>m; m++);

  cosv = new double[N];
  sinv = new double[N];
  bra = new int[N];

  FOR(i, 0, N) {
    double v = 2*PI*i/N;
    cosv[i] = cos(v); sinv[i] = sin(v);
  }
  FOR(i, 1, N) bra[i] = rev_bits(i, m-1);
  return true;
}

void fft(double *x, double *y, int N) {
  return;
}

