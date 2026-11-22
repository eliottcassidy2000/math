#include<bits/stdc++.h>
using namespace std;
using I=__int128_t;
long long P[65536],B[65536];
I p(I z){
 I V=0;
 while(z>=65536){
  int r=z&65535;z>>=16;
  V+=(65535-B[r])*z+P[r]+(r<32768&&(z&1));
 }
 return V+P[int(z)];
}
I mod(I n,I q){I r=n%q;return r<0?r+q:r;}
I floorq(I a,I b){return (a-mod(a,b))/b;}
I ceilq(I a,I b){return -floorq(-a,b);}
string out(I n){if(n<0)return "-"+out(-n);string r;do{r+='0'+n%10;n/=10;}while(n);reverse(r.begin(),r.end());return r;}
int main(int argc,char**argv){
 if(argc!=2 || atoi(argv[1])<1 || atoi(argv[1])>60) return 2;
 for(int i=0;i<65536;i++){
  P[i]=P[i>>1]+(i&1?0:(i>>1)+((i>>1)&1));
  for(int j=0;j<16;j++)B[i]|=((i>>j)&1)<<(15-j);
 }
 I E=(I(1)<<60)*((I(1)<<61)-1), N=E-1;
 // c, h, K_0, K_1, exact address width.
 const array<I,5> charts[]={{81,25,43,68,115},{411,75,367,392,113},{1021,125,915,940,112}};
 for(auto [c,h,k0,k1,k]:charts){
  for(int b=0;b<2;b++){
   int ell=atoi(argv[1]);I A=I(1)<<ell,P=E/A;
   I K=b?k1:k0,C=b?c:c+h,J=b?k1+1-c:k0+1;
   I tail=b*(A-1);
   I lo=ceilq((I(1)<<(k-1))-tail,A);
   I hi=floorq(min((I(1)<<k)-1,(N-K)/c)-tail,A);
   if(hi<lo){cout<<out(c)<<" "<<b<<" width "<<ell<<" 0 empty\n";continue;}
   I f0=b?J:J-h*hi,f1=b?J+h*hi:J;
   lo=max(lo,ceilq(P-floorq(f1,A),C)-b);
   hi=min(hi,floorq(P-ceilq(f0,A),C)-b);
   I residue=-1,step=2*h;
   for(I i=0;i<step;i++){
    if(mod(c*(A*i+tail)+K-N,step)==0){residue=i;break;}
   }
   if(residue<0){cout<<out(c)<<" "<<b<<" width "<<ell<<" 0 congruence\n";continue;}
   I g=gcd((long long)mod(c*A,step),(long long)step);step/=g;
   I tests=0,solutions=0;
   for(I z=lo+mod(residue-lo,step);z<=hi;z+=step){
    tests++;
    I v=p(z);
    I F=b?J+h*v+A*C*(z+1):J+h*(v-z+(z&1))+A*C*z;
    if(F==E){solutions++;cout<<"HIT "<<out(A*z+tail)<<"\n";}
   }
   cout<<out(c)<<" "<<b<<" width "<<ell<<" "<<out(tests)<<" tested "<<out(solutions)<<" hits"<<endl;
  }
 }
}
