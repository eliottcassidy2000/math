#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
void need(bool v,const char*s){if(!v)throw runtime_error(s);}
vector<unsigned> graph(unsigned n,unsigned code){vector<unsigned>a(n);unsigned k=0;for(unsigned i=0;i<n;i++)for(unsigned j=i+1;j<n;j++,k++){if(code>>k&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
vector<unsigned> op(vector<unsigned>a){unsigned all=(1u<<a.size())-1;for(unsigned i=0;i<a.size();i++)a[i]^=all^(1u<<i);return a;}
vector<unsigned> pairs(const vector<unsigned>&Q){unsigned n=Q.size();vector<unsigned>a(2*n);for(unsigned i=0;i<n;i++){a[2*i]|=1u<<(2*i+1);for(unsigned j=0;j<n;j++)if(Q[i]>>j&1){a[2*i]|=3u<<(2*j);a[2*i+1]|=3u<<(2*j);}}return a;}
vector<unsigned> attach(vector<unsigned>a,unsigned cm){unsigned n=a.size();a.push_back((1u<<n)-1+(1u<<(n+1)));a.push_back(0);for(unsigned i=0;i<n;i++){if(cm>>i&1)a[i]|=1u<<(n+1);else a[n+1]|=1u<<i;}return a;}
unsigned hp(const vector<unsigned>&a){unsigned n=a.size(),N=1u<<n;vector<unsigned char>d(size_t(N)*n);for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;for(unsigned s=1;s<N;s++)for(unsigned i=0;i<n;i++){unsigned z=d[size_t(s)*n+i];if(!z)continue;unsigned next=a[i]&~s;while(next){unsigned b=next&-next;next-=b;auto&v=d[size_t(s|b)*n+__builtin_ctz(b)];v=(v+z)&7;}}unsigned out=0;for(unsigned i=0;i<n;i++)out+=d[size_t(N-1)*n+i];return out&7;}
unsigned ends(const vector<unsigned>&a){unsigned n=a.size(),N=1u<<n;vector<unsigned>d(N);for(unsigned i=0;i<n;i++)d[1u<<i]=1u<<i;for(unsigned mask=1;mask<N;mask++){unsigned last=d[mask];while(last){unsigned b=last&-last;last-=b;unsigned next=a[__builtin_ctz(b)]&~mask;while(next){unsigned t=next&-next;next-=t;d[mask|t]^=t;}}}return d.back();}
int main(){unsigned checks=0;
 for(unsigned n=1;n<=4;n++)for(unsigned code=0;code<(1u<<(n*(n-1)/2));code++){
  auto q=graph(n,code),u=op(q);unsigned b=ends(q)^ends(u);
  for(unsigned i=0;i<n;i++){unsigned f=hp(attach(pairs(q),3u<<(2*i))),g=hp(attach(pairs(u),3u<<(2*i)));need(((f+8-g)&3)==2*((b>>i)&1),"endpoint inflation identity");checks+=2;}
 }
 auto q=graph(8,600653),u=op(q);unsigned f=hp(attach(pairs(q),3)),g=hp(attach(pairs(u),3));need(f==5&&g==7,"literal eighteen-name separator");
 cout<<"{\"status\":\"PASS\",\"small_literal_queries\":"<<checks<<",\"physical_old_names\":16,\"augmented_names\":18,\"quotient_code\":600653,\"a_mask\":0,\"c_mask\":3,\"source_mod8\":"<<f<<",\"target_mod8\":"<<g<<",\"quotient_start_parity\":"<<ends(u)<<",\"quotient_end_parity\":"<<ends(q)<<"}\n";
}
