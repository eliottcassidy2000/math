#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;void need(bool b,const char*s){if(!b)throw runtime_error(s);}
vector<int> opposite(vector<int>a){int all=(1<<a.size())-1;for(int i=0;i<(int)a.size();i++)a[i]^=all^(1<<i);return a;}
vector<int> ext(vector<int>a,int x,int c){int q=a.size();a.push_back(1<<(q+1));a.push_back(0);for(int i=0;i<q;i++){if(x>>i&1)a[i]|=1<<q;else a[q]|=1<<i;if(c>>i&1)a[i]|=1<<(q+1);else a[q+1]|=1<<i;}return a;}
int weighted(const vector<int>&Q,int a,int c){int q=Q.size(),n=q+2,N=1<<q;auto rows=ext(Q,a,c);vector<int>rs,idx(N,-1);for(int m=0;m<N;m++)if(__builtin_popcount((unsigned)m)<=2){idx[m]=rs.size();rs.push_back(m);}vector<unsigned char>d(rs.size()*N*4*n);auto at=[&](int r,int s,int p,int last)->unsigned char&{return d[((r*N+s)*4+p)*n+last];};for(int i=0;i<q;i++)at(0,1<<i,0,i)=1;at(0,0,1,q)=1;at(0,0,2,q+1)=1;
for(int ri=0;ri<(int)rs.size();ri++){int repeated=rs[ri];for(int visited=0;visited<N;visited++){if(repeated&~visited)continue;for(int ports=0;ports<4;ports++)for(int last=0;last<n;last++){int z=at(ri,visited,ports,last);if(!z)continue;int next=rows[last];while(next){int b=next&-next;next-=b;int v=__builtin_ctz((unsigned)b),rr=ri,ss=visited,pp=ports,weight=1;if(v<q){if(!(visited&b))ss|=b;else{if((repeated&b)||__builtin_popcount((unsigned)repeated)==2)continue;rr=idx[repeated|b];weight=2;}}else{int bit=1<<(v-q);if(ports&bit)continue;pp|=bit;}auto &y=at(rr,ss,pp,v);y=(y+weight*z)&7;}}}}
int total=0;for(int ri=0;ri<(int)rs.size();ri++)for(int last=0;last<n;last++)total+=at(ri,N-1,3,last);return total&7;}
vector<int> inflate(const vector<int>&Q){int q=Q.size();vector<int>T(2*q);for(int i=0;i<q;i++){T[2*i]|=1<<(2*i+1);for(int j=0;j<q;j++)if(Q[i]>>j&1){T[2*i]|=3<<(2*j);T[2*i+1]|=3<<(2*j);}}return T;}
int expanded(int m,int q){int x=0;for(int i=0;i<q;i++)if(m>>i&1)x|=3<<(2*i);return x;}
int literal(const vector<int>&rows){int n=rows.size(),N=1<<n;vector<unsigned char>d(N*n);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;for(int s=1;s<N;s++)for(int last=0;last<n;last++){int z=d[s*n+last];if(!z)continue;int next=rows[last]&~s;while(next){int b=next&-next;next-=b;auto &x=d[(s|b)*n+__builtin_ctz((unsigned)b)];x=(x+z)&7;}}int h=0;for(int i=0;i<n;i++)h+=d[(N-1)*n+i];return h&7;}
vector<int> graph(int n,int code){vector<int>Q(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(code>>k&1)Q[i]|=1<<j;else Q[j]|=1<<i;return Q;}
int main(){
long checks=0;vector<vector<int>>small;
for(int q=1;q<=4;q++){
 int count=1<<(q*(q-1)/2),N=1<<q;vector<vector<int>>bank(count,vector<int>(N*N));
 for(int code=0;code<count;code++){
  auto Q=graph(q,code),T=inflate(Q);
  for(int a=0;a<N;a++){for(int c=0;c<N;c++){
   int value=weighted(Q,a,c);need(value==literal(ext(T,expanded(a,q),expanded(c,q))),"literal run-DP control");bank[code][a+(c<<q)]=value;checks++;
  }}
 }
 if(q>=3){for(int code=0;code<count;code++){
  auto Q=graph(q,code);vector<int>record={q,code};
  auto coefficient=[&](int mask){int value=0;for(int s=mask;;s=(s-1)&mask){int delta=bank[code][s]-bank[count-1-code][s];value+=((__builtin_popcount((unsigned)(mask^s))&1)?-1:1)*delta;if(!s)break;}return (value%8+8)%8;};
  for(int j=0;j<q;j++){int L=0;for(int i=0;i<q;i++){if(i!=j){L+=coefficient((1<<i)|(1<<(q+j)));if(Q[j]>>i&1)L+=coefficient((1<<i)|(1<<j)|(1<<(q+j)));}}record.push_back(L%8);}
  for(int j=0;j<q;j++){for(int i=0;i<q;i++){record.push_back(i==j?0:coefficient((1<<i)|(1<<j)|(1<<(q+j))));}}
  small.push_back(record);
 }}
}
vector<int>Q={598,508,344,865,40,69,656,317,113,438};auto U=opposite(Q);int coef=0;
cout<<"{\"status\":\"PASS\",\"literal_small_query_checks\":"<<checks<<",\"small_trace_controls\":[";
for(int j=0;j<(int)small.size();j++){if(j)cout<<",";cout<<"[";for(int k=0;k<(int)small[j].size();k++){if(k)cout<<",";cout<<small[j][k];}cout<<"]";}
cout<<"],\"rows\":[";for(int i=0;i<10;i++){if(i)cout<<",";cout<<Q[i];}cout<<"],\"coarse_shared_cubic\":[3,0],\"actual_pair_queries\":[";
bool first=true;for(int bits=0;bits<8;bits++){
 int a=((bits&1)?1:0)|((bits&2)?8:0),c=(bits&4)?1:0;int f=weighted(Q,a,c),g=weighted(U,a,c),sgn=(__builtin_popcount((unsigned)bits)&1)?1:-1;coef+=sgn*(f-g);
 if(!first){cout<<",";} first=false;cout<<"{\"a\":"<<expanded(a,10)<<",\"c\":"<<expanded(c,10)<<",\"source_mod8\":"<<f<<",\"target_mod8\":"<<g<<",\"finite_difference_sign\":"<<sgn<<"}";
}
coef=(coef%8+8)%8;need(coef==4,"hostile cubic separator");cout<<"],\"coefficient_mod8\":"<<coef<<"}\n";
}
