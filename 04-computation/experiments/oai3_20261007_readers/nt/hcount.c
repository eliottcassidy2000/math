/* class numbers h(-N) (primitive reduced forms) for N <= X, then min of L(1,chi) log N over
   fundamental N in [Nlo, X]; L(1,chi) = pi h / sqrt N (w=2, N>4). Lane nt sanity check. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static int gcd(int a,int b){while(b){int t=a%b;a=b;b=t;}return a;}
static int sqfree(long n){for(long p=2;p*p<=n;p++) if(n%(p*p)==0) return 0; return 1;}
static int fund(long N){ if(N%4==3) return sqfree(N); if(N%4==0){long n=N/4; return (n%4==1||n%4==2)&&sqfree(n);} return 0;}
int main(int argc,char**argv){ long X=atol(argv[1]), Nlo=atol(argv[2]);
  int *h=calloc(X+1,sizeof(int));
  for(long a=1;3*a*a<=X;a++) for(long b=-a+1;b<=a;b++) for(long c=a;;c++){ long N=4*a*c-b*b; if(N>X) break; if(N<=0) continue;
      if(b<0 && a==c) continue; if(gcd(gcd(a,labs(b)),c)!=1) continue; h[N]++; }
  double best=1e9; long bN=0; double best2=1e9; long bN2=0;
  for(long N=Nlo;N<=X;N++) if(fund(N)){ double L=M_PI*h[N]/sqrt((double)N); double v=L*log((double)N); if(v<best){best=v;bN=N;} if(L<best2){best2=L;bN2=N;} }
  printf("fundamental N in [%ld,%ld]: min L(1,chi)*log N = %.4f at N=%ld (h=%d); min L(1,chi) = %.4f at N=%ld (h=%d)\n",Nlo,X,best,bN,h[bN],best2,bN2,h[bN2]);
  return 0;}
