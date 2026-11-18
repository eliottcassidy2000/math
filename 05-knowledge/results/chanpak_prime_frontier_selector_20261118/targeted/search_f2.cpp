#include<bits/stdc++.h>
using namespace std;using I=long long;
I F(I u){I z=1;while(u){I p=u/2;if(!(u&1))z+=p+(p&1);u=p;}return z;}
I ceildiv(I a,I b){return a>=0?(a+b-1)/b:a/b;}
int main(int argc,char**argv){
 ifstream f(argv[1]);vector<array<I,3>> bank;I a,b,c;
 while(f>>a>>b>>c)bank.push_back({a,b,c});
 for(int i=2;i<argc;i++){
  I N=atoll(argv[i]);bool found=0;
  for(auto [h,e,R]:bank){
   I L=2;int k=1;while((L*R)*(L*R)<N)L*=2,k++;
   for(I r=0;r<L;r++){
    I Z=0;for(int j=0;j<k;j++)Z=2*Z+1-(r>>j&1);
    I B=L*R+h*Z,C=N-R*(r+1)+e*(r%2==0)-h*(F(r)-1);
    for(I q=max(0LL,ceildiv(C-2*h,B+h));q<=(C-h)/B;q++){
     I v=L*q+r;
     if(h*F(v)-e*(v%2==0)+(v+1)*R==N){
      cout<<N<<" v,h,e,R "<<v<<" "<<h<<" "<<e<<" "<<R<<" f "<<F(v)<<endl;found=1;break;
     }
    }
    if(found)break;
   }
   if(found)break;
  }
  if(!found)cout<<N<<" MISS"<<endl;
 }
}
