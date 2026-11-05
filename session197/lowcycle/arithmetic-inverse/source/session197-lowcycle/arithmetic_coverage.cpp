// Exact positive image of a fixed finite prefix bank and one arbitrary final gadget.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <vector>
using namespace std;
struct S{int h,e,r;bool operator==(const S&o)const{return h==o.h&&e==o.e&&r==o.r;}};
struct Hash{size_t operator()(const S&s)const{return ((uint64_t)s.h*1000003+s.e)*1000003+s.r;}};
struct G{int h,q,b;};
int main(int argc,char**argv){
 int X=argc>1?atoi(argv[1]):10000001,B=argc>2?atoi(argv[2]):1000;
 if(X<3||B<1||B>X){cerr<<"bad bounds\n";return 2;}
 vector<G>g(X+1);g[1]={1,1,1};
 for(int L=2;L<=X;L++){int ell=(L+1)/2;auto t=g[ell];g[L]=L%2?G{t.h+ell-ell%2,t.q+t.b,2*t.b}:G{t.h,t.q,2*t.b};}
 unordered_set<S,Hash>seen;vector<S>q{{1,1,1}};seen.insert(q[0]);
 for(size_t at=0;at<q.size();at++){auto x=q[at];for(int L=1;L<=B/x.r;L++){auto t=g[L];int R=t.q*x.h+t.b*x.r;if(R>B)continue;S y{t.h*x.h-(L%2)*x.e+L*x.r,x.h-(L==1)*x.e,R};if(seen.insert(y).second)q.push_back(y);}}
 vector<unsigned char>bits((X+8)/8);auto mark=[&](int h){bits[h/8]|=1<<(h%8);};mark(1);long long trials=0;
 for(auto x:q)for(int L=1;L<=X/x.r;L++){trials++;long long h=1LL*g[L].h*x.h-(L%2)*x.e+1LL*L*x.r;if(h<=X)mark((int)h);}
 vector<int>missing;for(int h=3;h<=X;h+=2)if(!(bits[h/8]>>(h%8)&1))missing.push_back(h);
 cout<<"{\"bound\":"<<X<<",\"prefix_R_bound\":"<<B<<",\"prefix_states\":"<<q.size()<<",\"terminal_trials\":"<<trials<<",\"missing_odd\":[";
 for(size_t i=0;i<missing.size();i++){if(i)cout<<",";cout<<missing[i];}cout<<"]}\n";
 if(argc>3){ofstream f(argv[3],ios::binary);f.write((const char*)bits.data(),bits.size());if(!f)return 3;}
}
