#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
using I=int64_t;
I F(I u) { I z=1; while (u) { I p=u/2; if (!(u&1)) z+=p+(p&1); u=p; } return z; }
I fl(I x,I y){return x>=0?x/y:-( (-x+y-1)/y);}
I cl(I x,I y){return -fl(-x,y);}
struct W {I a=0,v=0;};
int main(int argc,char **argv){
 if(argc!=5 && argc!=6){cerr<<"usage: band first_odd last_odd F1_bmax output [F2_profiles]\n";return 2;}
 I lo=stoll(argv[1]),hi=stoll(argv[2]),cap=stoll(argv[3]);
 if(lo<1||hi<lo||hi>(1LL<<42)||!(lo&1)||!(hi&1)||cap<0||hi-lo>20000000){cerr<<"bad endpoints or scope\n";return 2;}
 vector<W> row((hi-lo)/2+1); I missed=row.size(),checked=0, max_b=0;
 auto seek=[&](I h,I end,I R,I vmin,I tag){
  I L=2; while((__int128)L*L*R<hi)L*=2;
  // Reverse complement via an incrementally bit-reversed counter.
  I Z=L-1;
  for(I r=0;r<L&&missed;r++){
   I B=L*R+h*Z;
   I K=R*(r+1)+h*(F(r)-1)-end*!(r&1);
   I q0=max({I(0),cl(lo-K-2*h,B+h),cl(vmin-r,L)});
   I q1=fl(hi-K-h,B);
   for(I q=q0;q<=q1;q++){
    ++checked;
    I v=L*q+r,N=K+B*q+h*(F(q)+(r<L/2)*(q&1));
    if(N>=lo&&N<=hi&&(N&1)&&!row[(N-lo)/2].a){
     if(h*F(v)-end*!(v&1)+R*(v+1)!=N)abort();
     row[(N-lo)/2]={tag,v}; --missed;
    }
   }
   I bit=L/2;
   while(bit&&(Z&bit)==0){Z^=bit;bit/=2;}
   Z^=bit;
  }
 };
 for(I b=1;b<=cap&&b<=hi/b&&missed;b++){seek(F(b-1),b&1,b,b-1,b);max_b=b;}
 I miss1=missed, j=0;
 if(argc==6){ifstream f(argv[5]);if(!f){cerr<<"missing profile file\n";return 2;} I id,h,e,R;while(missed&&f>>id>>h>>e>>R){++j;if(id<1||h<1||e<0||R<1||max({h,e,R})>1000000){cerr<<"bad profile\n";return 2;}seek(h,e,R,0,-id);}}
 ofstream out(argv[4]);if(!out)return 2;
 for(size_t i=0;i<row.size();i++)out<<lo+2*i<<" "<<row[i].a<<" "<<row[i].v<<"\n";
 cout<<"band "<<lo<<" "<<hi<<" size "<<row.size()<<" miss after F1 "<<miss1<<" final miss "<<missed<<" last b "<<max_b<<" F2 profiles "<<j<<" candidate q reads "<<checked<<"\n";
}
