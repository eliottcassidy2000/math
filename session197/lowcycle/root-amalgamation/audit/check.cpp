#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;using G=vector<unsigned>;using I=long long;
void need(bool x,const char*s){if(!x)throw runtime_error(s);}
struct P{I h,e,s,d,R,L,b,C;};
P profile(const G&A,int r=0){
 int n=A.size(),N=1<<n,full=N-1;vector<I> dp(N*n),root(N*n),H(N),S(N),E(N);
 for(int v=0;v<n;v++)dp[(1<<v)*n+v]=1;root[(1<<r)*n+r]=1;
 for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(m>>v&1){
  I x=dp[m*n+v],y=root[m*n+v];H[m]+=x;S[m]+=y;if(v==r)E[m]=x;
  unsigned out=A[v]&~m;while(out){int w=countr_zero(out);out&=out-1;dp[(m|(1<<w))*n+w]+=x;root[(m|(1<<w))*n+w]+=y;}
 }
 H[0]=1;int dmask=full^(1<<r);I R=0,L=0,b=0,C=0;
 for(int m=0;m<N;m++){R+=S[m]*H[full^m];L+=E[m]*H[full^m];if(!(m>>r&1))b+=H[m]*H[dmask^m];}
 for(int v=0;v<n;v++)if(A[v]>>r&1)C+=root[full*n+v];
 return {H[full],E[full],S[full],H[dmask],R,L,b,C};
}
G induce(const G&A,const vector<int>&v){G B(v.size());for(int i=0;i<(int)v.size();i++)for(int j=0;j<(int)v.size();j++)if(A[v[i]]>>v[j]&1)B[i]|=1<<j;return B;}
G normalize(const G&A,int r){vector<int>v{r};for(int i=0;i<(int)A.size();i++)if(i!=r)v.push_back(i);return induce(A,v);}
G opposite(G A){unsigned all=(1<<A.size())-1;for(int i=0;i<(int)A.size();i++)A[i]^=all^(1<<i);return A;}
G amalg(const G&A,const G&B){int n=A.size(),m=B.size();G T(n+m-1);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(A[i]>>j&1)T[i]|=1<<j;
 for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(B[i]>>j&1)T[i? n+i-1:0]|=1<<(j?n+j-1:0);
 for(int i=1;i<n;i++)for(int j=n;j<n+m-1;j++)T[i]|=1<<j;return T;}
G family(const G&V,const G&W,int j){G K(4);K[1]=1<<2;K[2]=1<<3;K[3]=1<<1;for(int i=0;i<3;i++){if(i<j)K[0]|=1<<(i+1);else K[i+1]|=1;}
 return amalg(amalg(V,K),opposite(W));}
G gadget(int L){int n=L-1,s=n?32-countl_zero((unsigned)n):0;G T(s+1);for(int i=0;i<s;i++){for(int j=i+1;j<s;j++)T[i]|=1<<j;if(n>>i&1)T[i]|=1<<s;else T[s]|=1<<i;}return T;}
pair<G,int> seam(const G&A,int a,const G&B,int b){int n=A.size(),m=B.size();G T=A;for(unsigned row:B)T.push_back(row<<n);
 for(int i=0;i<n;i++)for(int j=0;j<m;j++){if(i==a&&j==b)T[n+j]|=1<<i;else T[i]|=1<<(n+j);}return {T,n+b};}
const I M[4][4][4]={{{0,0,0,3},{0,0,0,-3},{0,0,12,-3},{3,0,-12,3}},{{0,0,0,3},{0,0,0,-2},{0,0,12,-6},{3,-1,-9,5}},{{0,0,0,3},{0,0,0,-1},{0,0,12,-9},{3,-2,-6,5}},{{0,0,0,3},{0,0,0,0},{0,0,12,-12},{3,-3,-3,3}}};
int main(){vector<pair<G,P>>bank;for(int n=1;n<=4;n++){int pairs=n*(n-1)/2;for(int mask=0;mask<(1<<pairs);mask++){G A(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(mask>>k&1)A[i]|=1<<j;else A[j]|=1<<i;
 for(int r=0;r<n;r++){G B=normalize(A,r);bank.push_back({B,profile(B)});}}}
 long general=0,families=0,boundaries=0,triangle_tails=0,closing_only=0,direct=0;array<long,4>counts{};
 for(auto&[A,p]:bank)for(auto&[B,q]:bank){P x=profile(amalg(A,B));need(x.h==(p.R-p.d)*(q.L-q.d)+(p.h-p.e)*q.d+(q.h-q.s)*p.d+p.e*q.s,"general H");need(x.R==(p.R-p.d)*q.b+p.d*q.R&&x.L==(q.L-q.d)*p.b+q.d*p.L,"general covers");need(x.d==p.d*q.d&&x.e==(B.size()==1?p.e:p.d*q.e)&&x.s==(A.size()==1?q.s:p.s*q.d),"general endpoints including singleton");general++;}
 for(auto&[V,p]:bank)if(p.s==1)for(auto&[W,q]:bank)if(q.s==1)for(int j=0;j<4;j++){
  G T=family(V,W,j);P x=profile(T);int expected=(V.size()>1&&W.size()>1)?3:(V.size()>1?3-j:(W.size()>1?j:(j==1||j==2)));
  need(x.C==expected,"all singleton cycle boundaries");boundaries++;
  bool allowed=(V.size()>1&&W.size()>1)||(V.size()>1&&W.size()==1&&j==0)||(V.size()==1&&W.size()>1&&j==3);
  if(!allowed)continue;array<I,4>a{p.h,p.e,p.R,p.d},b{q.h,q.e,q.R,q.d};I H=0;for(int s=0;s<4;s++)for(int t=0;t<4;t++)H+=a[s]*M[j][s][t]*b[t];need(x.h==H&&x.s==3*q.d&&x.e==3*p.d&&x.d==3*p.d*q.d,"physical bilinear family");families++;counts[j]++;
 }
 // Canonical Type-II initial graph, followed by actual one-entry suffixes and closure.
 vector<G> Ubank{{0},{2,4,1}};vector<vector<int>> suffixes{{},{1},{2},{3},{1,1}};
 for(const G&U:Ubank)for(int plen=0;plen<=1;plen++)for(int qlen=0;qlen<=1;qlen++)for(int pm=0;pm<(1<<plen);pm++)for(int j=0;j<4;j++){
  if(!qlen&&j==3)continue;int un=U.size(),tri=1+un+plen,qstart=tri+3,n=qstart+qlen;G X(n);
  for(int a=0;a<un;a++)for(int b=0;b<un;b++)if(U[a]>>b&1)X[1+a]|=1<<(1+b);
  for(int a=1;a<n;a++)for(int b=a+1;b<n;b++){if(a<1+un&&b<1+un)continue;if(a>=tri&&a<tri+3&&b>=tri&&b<tri+3)continue;X[a]|=1<<b;}
  X[tri]|=1<<(tri+1);X[tri+1]|=1<<(tri+2);X[tri+2]|=1<<tri;
  for(int a=1;a<n;a++){bool out=a==1;if(a>=1+un&&a<tri)out=pm>>(a-1-un)&1;if(a>=tri&&a<tri+3)out=a-tri<j;if(out)X[0]|=1<<a;else X[a]|=1;}
  need(profile(X).s==3&&profile(X).C>0,"canonical strong Type-II source");
  for(auto recipe:suffixes){G T=X;int r=0,first=-1;for(int L:recipe){G B=gadget(L);auto made=seam(T,r,B,B.size()-1);if(first<0)first=made.second;T=made.first;r=made.second;}
   auto closed=seam(T,r,G{0},0);if(first<0)first=closed.second;T=closed.first;
   vector<int>v{0},w{0},bundle;for(int a=1;a<tri;a++)v.push_back(a);for(int a=qstart;a<(int)T.size();a++)w.push_back(a);for(int a=n;a<(int)T.size();a++)bundle.push_back(a);
   G V=induce(T,v),Wop=induce(T,w),BG=induce(T,bundle);P vp=profile(V),wp=profile(Wop),bp=profile(BG,first-n);
   need(vp.s==1&&wp.e==1&&bp.e==1,"literal canonical bundling profiles");need(T==family(V,opposite(Wop),j),"exact canonical triangle-tail reconstruction");need(profile(T).C==3,"canonical final cycle count");triangle_tails++;if(recipe.empty()){need(BG.size()==1,"zero-suffix singleton bundle");closing_only++;}
  }
 }
 // Direct degree-one closure with its triangle in the deletion tail, including W singleton.
 for(const G&U:Ubank)for(int plen=0;plen<=2;plen++)for(int qlen=0;qlen<=2;qlen++){
  G V(U.size()+plen+1);for(int i=0;i<(int)U.size();i++)for(int j=0;j<(int)U.size();j++)if(U[i]>>j&1)V[1+i]|=1<<(1+j);
  V[0]|=2;for(int i=2;i<(int)V.size();i++)V[i]|=1;for(int i=1;i<(int)V.size();i++)for(int j=max(i+1,1+(int)U.size());j<(int)V.size();j++)V[i]|=1<<j;
  G W(qlen+1);for(int i=0;i<=qlen;i++)for(int j=i+1;j<=qlen;j++)W[i]|=1<<j;
  need(profile(V).s==1&&profile(W).s==1,"direct closure actual sources");G T=family(V,W,0);need(T[0]==2&&profile(T).C==3,"direct degree-one triangle closure");direct++;
 }
 cout<<"{\n  \"status\":\"PASS\",\n  \"rooted_labelled_inputs_through_four\":"<<bank.size()<<",\n  \"general_amalgamation_pairs\":"<<general<<",\n  \"physical_C3_pairs\":"<<families<<",\n  \"all_singleton_boundary_pairs\":"<<boundaries<<",\n  \"canonical_triangle_tail_reconstructions\":"<<triangle_tails<<",\n  \"closing_only_singleton_bundles\":"<<closing_only<<",\n  \"direct_degree_one_triangle_tail_closures\":"<<direct<<",\n  \"j_counts\":["<<counts[0]<<","<<counts[1]<<","<<counts[2]<<","<<counts[3]<<"]\n}\n";
}
