#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;using I=long long;using W=__int128_t;
string dec(W v){if(v==0)return "0";bool neg=v<0;if(neg)v=-v;string s;while(v){s+=char('0'+v%10);v/=10;}if(neg)s+='-';reverse(s.begin(),s.end());return s;}
vector<I> endpoints(const vector<int>&rows){int n=rows.size(),N=1<<n;vector<I>d((size_t)N*n*n);for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(m>>v&1){int rest=m^(1<<v);I*out=&d[((size_t)m*n+v)*n];if(!rest){out[0]=1;continue;}int degree=__builtin_popcount((unsigned)rest)-1;for(int u=0;u<n;u++)if(rest>>u&1){const I*prev=&d[((size_t)rest*n+u)*n];bool edge=rows[u]>>v&1;for(int q=0;q<=degree;q++){I z=prev[q];if(edge)out[q]+=z;out[q+1]+=z;}}}return d;}
int main(){int count;cin>>count;cout<<"[";for(int caseid=0;caseid<count;caseid++){int n,z,G;cin>>n>>z>>G;if(n>17)throw runtime_error("current memory guard n<=17");vector<int>r(n),group(n),op(n);for(int&x:r)cin>>x;for(int&x:group)cin>>x;int N=1<<n,F=N-1;for(int i=0;i<n;i++)op[i]=F^r[i]^(1<<i);auto E=endpoints(r),S=endpoints(op);vector<vector<W>>ans(G,vector<W>(n-1));for(int m=1;m<F;m++)if(m>>z&1){int a=__builtin_popcount((unsigned)m)-1,b=n-a-2,c=F^m;const I*ez=&E[((size_t)m*n+z)*n],*sz=&S[((size_t)m*n+z)*n];for(int i=0;i<n;i++)if(c>>i&1){int g=group[i];const I*ei=&E[((size_t)c*n+i)*n],*si=&S[((size_t)c*n+i)*n];for(int p=0;p<=a;p++)for(int q=0;q<=b;q++)ans[g][p+q]+=W(ez[p])*si[q]-W(sz[p])*ei[q];}}
if(caseid)cout<<",";cout<<"{\"id\":"<<caseid<<",\"n\":"<<n<<",\"z\":"<<z<<",\"skews\":[";for(int g=0;g<G;g++){if(g)cout<<",";cout<<"[";for(int p=0;p<n-1;p++){if(p)cout<<",";cout<<dec(ans[g][p]);}cout<<"]";}cout<<"]}";cerr<<"case "<<caseid<<" done\n";
}cout<<"]\n";}
