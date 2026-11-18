#include<bits/stdc++.h>
using namespace std;using I=long long;
I F(I u){I z=1;while(u){I p=u/2;if(!(u&1))z+=p+(p&1);u=p;}return z;}
I cl(I x,I y){return x>=0?(x+y-1)/y:x/y;}
I count(I b,I e){return b*e+F(b-1)*F(e-1)-(b&1)*(e&1);}
int main(int argc,char**argv){for(int i=1;i<argc;i++){I N=atoll(argv[i]);bool hit=false;
for(I b=2;b*b<=N&&!hit;b++){
 I a=F(b-1),L=2;int k=1;
 while(L*L*b<N)L*=2,k++;
 for(I r=0;r<L;r++){
  I A=0;for(int j=0;j<k;j++)A=2*A+1-(r>>j&1);
  I B=b*L+a*A,M=N-b*(r+1)-a*(F(r)-1)+(b&1)*(r%2==0);
  for(I q=max(0LL,cl(M-2*a,B+a));q<=(M-a)/B;q++){
   I e=L*q+r+1;
   if(e>=b&&count(b,e)==N){cout<<N<<" "<<b<<" "<<e<<"\n";hit=true;break;}
  }
  if(hit)break;
 }
}
if(!hit)cout<<N<<" MISS\n";
}}
