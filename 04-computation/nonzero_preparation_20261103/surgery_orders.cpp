#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
uint64_t gates=0;
void need(bool yes,const char* text){++gates;if(!yes)throw runtime_error(text);}
vector<uint64_t> dp(const vector<unsigned>& rows) {
 unsigned n=rows.size(),size=1u<<n; vector<uint64_t>d(size*n);
 for(unsigned i=0;i<n;i++)d[(1u<<i)*n+i]=1;
 for(unsigned mask=1;mask<size;mask++)for(unsigned vs=mask;vs;vs&=vs-1){unsigned i=__builtin_ctz(vs);auto z=d[mask*n+i];if(!z)continue;
 for(unsigned q=rows[i]&(size-1-mask);q;q&=q-1){unsigned j=__builtin_ctz(q);d[(mask|(1u<<j))*n+j]+=z;}}
 return d;
}
int main(){try {unsigned n,p;uint64_t rows=0;while(cin>>n>>p){
 need(n>=3&&n<=14,"bounded exact core order");
 vector<unsigned>r(n),rev(n);for(auto&v:r)need(bool(cin>>v),"complete input row");
 unsigned size=1u<<n,full=size-1;
 bool typed=true;for(unsigned i=0;i<n;i++){typed=typed&&!(r[i]&~full)&&!(r[i]>>i&1);for(unsigned j=i+1;j<n;j++)typed=typed&&(((r[i]>>j&1)+(r[j]>>i&1))==1);}
 need(typed,"literal tournament adjacency");need(p&&p<full&&__builtin_popcount(p)<=2,"proper nonempty one/two entry set");
 vector<unsigned>reach=r;for(unsigned i=0;i<n;i++)reach[i]|=1u<<i;for(unsigned k=0;k<n;k++)for(unsigned i=0;i<n;i++)if(reach[i]>>k&1)reach[i]|=reach[k];bool connected=true;for(auto v:reach)connected=connected&&(v==full);need(connected,"actual strong core");
 for(unsigned i=0;i<n;i++)rev[i]=full^(1u<<i)^r[i];auto E=dp(r),S=dp(rev);vector<uint64_t>H(size),e(size),s(size);H[0]=1;
 for(unsigned m=1;m<size;m++)for(unsigned i=0;i<n;i++){H[m]+=E[m*n+i];if(p>>i&1){e[m]+=E[m*n+i];s[m]+=S[m*n+i];}}
 uint64_t p2=0,F=0,G=0,C=0;for(unsigned m=0;m<size;m++){p2+=H[m]*H[full^m];F+=e[m]*H[full^m];G+=s[m]*H[full^m];C+=e[m]*s[full^m];}
 uint64_t fact=1;for(unsigned i=2;i<=n;i++)fact*=i;
 need(H[full]<=fact&&p2<=(n+1)*fact&&F<=p2&&G<=p2&&C<=p2,"permutation/cut exact integer bounds");
 need(s[full]==2,"actual two-entry start budget");
 cout<<H[full]<<' '<<p2<<' '<<e[full]<<' '<<F<<' '<<G<<' '<<s[full]<<' '<<C<<'\n';++rows;}
 cerr<<"{\"rows\":"<<rows<<",\"gates\":"<<gates<<"}\n";
 }catch(const exception& e){cerr<<"REJECT "<<e.what()<<'\n';return 2;}
}
