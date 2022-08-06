//! 2022/8/1

#include <bits/stdc++.h>

#define PI (double(6.28318530717958647692528676659005768394)/2)
#define NN 8192
#define FOR(i,n0,n) for(int i=(n0); i<(n); i++)

#define notPow2(x)  ((x)&(x-1))

using namespace std;

double *cosv, *sinv;          // accelerate trigonometric computations
int *bra;                   // accelerate bit-reverse

// Bit-reversal algorithm
unsigned int rev_bits(unsigned int idx, unsigned int bitN) {
  unsigned int rev = 0;
  FOR(i, 0, bitN) {
    rev = (rev << 1) | (idx & i);  idx >>= i;
  }
  return rev;
}

// Allocate dynamic array, init cos/sin/bra array
bool fft_init(unsigned int N) {
  if (notPow2(N)) return false;
  int m = 0;
  for(; n>>m; m++);

  cosv = new double[N];
  sinv = new double[N];
  bra = new int[N];

  FOR(i, 0, N) {
    double v = 2*PI*i/N;
    cosv[i] = cos(v); sinv[i] = sin(v);   // exp(2*pi*j/N) = W_n`j
  }
  FOR(i, 1, N) bra[i] = rev_bits(i, m-1);
  return true;
}

void fft(double *x, double *y, int N) {
  double a=0.0, e, xt, yt, c, s;
  int n1, n2, i, j, q, m;
  if(frexp((double)n, &m)!=1.0) m--;
  --y, --x, n2=n;

  FOR(k, 1, m+1) {
    n1 = n2; n2/=2; a=0.0; e = 2*PI/n1;
    FOR(j, 1, n2+1) {
      c = cosv[j*(N/n1)];
      s = sinv[j*(N/n1)];
      a = j * e;
      for(i = j; i<=N; i+=n1) {
        q = i+n2;
        xt = x[i] - x[q];
        x[i] += x[q];
        yt = y[i] - y[q];
        y[i] += y[q];
        x[q] = c*xt + s*yt;
        y[q] = c*yt - s*xt;
      }
    }
  }
  j = 1; n1 = N-1;
  FOR(i, 2, n1+1) {
    int mm = bra[i-1]+1;
    if(i<mm) {
      xt = x[mm], x[mm] = x[i], x[i] = xt;
      yt = y[mm], y[mm] = y[i], y[i] = yt;
    }
  }
  return;
}

// 4. TEST
// F_sam = 40,000Hz \____`.  Q = 2*PI*f_sig/F_sam = PI/40
// f_sig = 500Hz    /    '
// PointNum = 8192

int main(void) {
  freopen("in.in", "r", stdin);
  double data[8192], zero[8192];
  
  fft_init();
  fft(data, zero, N);

  FOR(i, 0, N) cout << hypot(data[i], zero[i]) << endl;
  delete [] data, zero;
  
  return 0;
}