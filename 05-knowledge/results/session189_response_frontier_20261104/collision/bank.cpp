#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
using I=long long;
I H(const vector<int>& old,int a,int c){
 int n=old.size(),N=n+2,full=(1<<N)-1; vector<int> rows=old; rows.resize(N); rows[n]=1<<(n+1);
 for(int i=0;i<n;++i){if(a>>i&1)rows[i]|=1<<n;else rows[n]|=1<<i;if(c>>i&1)rows[i]|=1<<(n+1);else rows[n+1]|=1<<i;}
 vector<I>d((1<<N)*N,0);for(int i=0;i<N;++i)d[(1<<i)*N+i]=1;
 for(int s=1;s<=full;++s){int v=s;while(v){int i=__builtin_ctz(v);v&=v-1;I z=d[s*N+i];if(!z)continue;int w=rows[i]&~s;while(w){int j=__builtin_ctz(w);w&=w-1;d[(s|(1<<j))*N+j]+=z;}}}
 I ans=0;for(int i=0;i<N;++i)ans+=d[full*N+i];return ans;
}
template<class T>void array_out(ostream&o,const vector<T>&a){o<<"[";for(size_t i=0;i<a.size();++i){if(i)o<<",";o<<a[i];}o<<"]";}
int main(int argc,char**argv){
 string out=argc>1?argv[1]:"bank.json";ofstream o(out);o<<"{\"format\":\"row-major masks a,c; complete Boolean bidegree at most (2,2)\",\"cases\":[";int caseid=0;
 for(int m:{5,6})for(int e: {0,2}){if(m==6&&e==0)continue;int qcount=e?8:1;if(m==6)qcount=2;
 for(int qi=0;qi<qcount;++qi){int qm=m==6?(qi?5:2):qi;vector<vector<int>> perms;vector<int> p(m);iota(p.begin(),p.end(),0);perms.push_back(p);reverse(p.begin(),p.end());perms.push_back(p);iota(p.begin(),p.end(),0);swap(p[1],p[2]);perms.push_back(p);iota(p.begin(),p.end(),0);rotate(p.begin(),p.begin()+1,p.end());perms.push_back(p);
 for(int pi=0;pi<(int)perms.size();++pi){p=perms[pi];int n=m+e;vector<int>rows(n);for(int i=0;i<m;++i)for(int j=i+1;j<m;++j)rows[p[i]]|=1<<p[j];if(e){for(int j=0;j<2;++j)for(int i=0;i<m;++i){if(qm>>j&1)rows[i]|=1<<(m+j);else rows[m+j]|=1<<i;}if(qm>>2&1)rows[m]|=1<<(m+1);else rows[m+1]|=1<<m;}
 vector<int> masks{0};for(int s=1;s<(1<<n);++s)if(__builtin_popcount((unsigned)s)==1)masks.push_back(s);for(int s=1;s<(1<<n);++s)if(__builtin_popcount((unsigned)s)==2)masks.push_back(s);int L=masks.size();vector<int>index(1<<n,-1);for(int i=0;i<L;++i)index[masks[i]]=i;vector<I>vals(L*L),coefs(L*L);for(int i=0;i<L;++i)for(int j=0;j<L;++j)vals[i*L+j]=H(rows,masks[i],masks[j]);
 for(int i=0;i<L;++i)for(int j=0;j<L;++j){int a=masks[i],c=masks[j];I z=0;for(int s=a;;s=(s-1)&a){for(int t=c;;t=(t-1)&c){int sign=(__builtin_popcount((unsigned)(a^s))+__builtin_popcount((unsigned)(c^t)))%2?-1:1;z+=sign*vals[index[s]*L+index[t]];if(!t)break;}if(!s)break;}coefs[i*L+j]=z;}
 if(caseid)o<<",";o<<"{\"id\":"<<caseid++<<",\"m\":"<<m<<",\"exterior_size\":"<<e<<",\"quotient_mask\":"<<qm<<",\"permutation_id\":"<<pi<<",\"order\":";array_out(o,p);o<<",\"rows\":";array_out(o,rows);o<<",\"masks\":";array_out(o,masks);o<<",\"values\":";array_out(o,vals);o<<",\"coefficients\":";array_out(o,coefs);o<<"}";
 cerr<<"case "<<caseid<<" n="<<n<<" done\n";
 }}}
 o<<"]}\n";
}
