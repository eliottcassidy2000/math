#include<bits/stdc++.h>
using namespace std;
using I=__int128_t;
I P(I a){I r=0;while(a){int b=a&1;a>>=1;if(!b)r+=a+(a&1);}return r;}
string str(I a){if(!a)return "0";if(a<0)return "-"+str(-a);string st;while(a){st+=char('0'+a%10);a/=10;}reverse(st.begin(),st.end());return st;}
I gcdI(I a,I b){return b?gcdI(b,a%b):a;}
I inv(I a,I m){if(m==1)return 0;I o=m,x=1,y=0;while(m){I d=a/m;I na=a-d*m;a=m;m=na;I nx=x-d*y;x=y;y=nx;}assert(a==1);x%=o;return x<0?x+o:x;}
I mul(I a,I b,I m){if(a<(I(1)<<62)&&b<(I(1)<<62))return a*b%m;I s=0;while(b){if(b&1)s=(s+a)%m;a=2*a%m;b>>=1;}return s;}
array<I,3> profile(uint64_t w){
 // Four DP states: has A (after the apex) started? and is the last B
 // (before the apex) contact legal? Start in state (A=no, B=valid).
 I S[2][2]={{1,0},{0,0}};I me=1;
 for(int i=63-__builtin_clzll(w|1);w&&i>=0;i--){
  int e=w>>i&1;
  if(!e)me+=I(1)<<(63-__builtin_clzll(w)-i);
  I T[2][2]{};
  for(int a=0;a<2;a++)for(int b=0;b<2;b++){
   T[a][b]+=S[a][b]; // other slot
   T[a][e]+=S[a][b]; // in B; legal iff contact bit 0
   if(a||e)T[1][b]+=S[a][b]; // A, starts at a contact 1
  }
  memcpy(S,T,sizeof S);
 }
 return {2*(S[0][0]+S[1][0]),me,P(w)};
}
int main(int argc,char**argv){
 int p=argc>1?atoi(argv[1]):61, maxw=argc>2?atoi(argv[2]):19, cap=argc>3?atoi(argv[3]):64;
 I E=(I(1)<<(p-1))*((I(1)<<p)-1),N=E-1;
 long long wcount=0,empty=0,scanned=0,skipped=0,cand=0,exact=0,mc=0;
 int minskip=130,maxskip=-1;
 I best=E,bw=-1,by=-1;int bt=-1;
 auto process=[&](uint64_t w){
  wcount++;
  auto [Q,me,pc]=profile(w);
  I A=25*(I(w)+1)+28*Q,F=25*(1+pc),D=10+28*me;
  auto check=[&](int t,int parity,I r,I mod,I B,I V){
   I fun_min=(V+B+F-1)/(B+F),fun_max=(V-F)/B;
   if(fun_min<1)fun_min=1;
   if(fun_max<fun_min){empty++;return;}
   if(parity!=-1){ // parity of y; u=y+1 has complementary parity.
    if((r&1)==parity)r+=mod;
    mod*=2;
   }
   I first=fun_min+(r-fun_min%mod+mod)%mod;
   I n=first>fun_max?0:1+(fun_max-first)/mod;
   if(n==0){empty++;return;}
   if(n>cap){skipped++;minskip=min(minskip,t);maxskip=max(maxskip,t);return;}
   scanned++;
   for(I u=first;u<=fun_max;u+=mod){
    cand++;
    I val=B*u+F*(1+P(u-1))-V;
    if(val==0){exact++;cout<<"EXACT w "<<w<<" t "<<t<<" y "<<str(u-1)<<endl;}
    if((val<0?-val:val)<best){best=val<0?-val:val;bw=w;by=u-1;bt=t;mc=val<0?-1:1;}
   }
  };
  auto base=[&](I V){I g=gcdI(A,F),m=F/g;return pair(g,m);};
  {
   auto [g,m]=base(N+D+25);
   for(int y=0;y<2;y++){
    I V=N+D+(!y?25:0);
    if(V%g){empty++;continue;}
    I R=mul(V/g%m,inv(A/g%m,m),m);
    check(0,y,R,m,A,V);
   }
  }
  I V=N+D;auto [g,m]=base(V);
  if(V%g){empty++;return;}
  I R=mul(V/g%m,inv(A/g%m,m),m);
  for(int t=1;t<125&&A<=V>>t;t++){
    R=(R+(R&1)*m)/2;
    check(t,-1,R,m,A<<t,V);
  }
 };
 process(0);
 for(uint64_t w=2;w<(1ULL<<maxw);w+=2)process(w);
 // Outside the complete short bank: all even words of width <= 60 and
 // Hamming weight <= 3.
 for(int a=maxw;a<60;a++){
  process(1ULL<<a);
  for(int b=1;b<a;b++){
   process((1ULL<<a)|(1ULL<<b));
   for(int c=1;c<b;c++)process((1ULL<<a)|(1ULL<<b)|(1ULL<<c));
  }
 }
 cout<<"width "<<maxw<<" w "<<wcount<<" empty "<<empty<<" scans "<<scanned<<" skipped "<<skipped
     <<" skips t-range "<<minskip<<".."<<maxskip<<" cand "<<cand<<" exact "<<exact
     <<" closest-grid-delta "<<mc<<"*"<<str(best)<<" at "<<str(bw)<<" "<<bt<<" "<<str(by)<<endl;
}
