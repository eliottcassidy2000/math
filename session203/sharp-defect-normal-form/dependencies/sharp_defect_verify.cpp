#include <bits/stdc++.h>
using namespace std;
using I=long long; using G=vector<int>;
void need(bool q,const string&s){if(!q)throw runtime_error(s);}
struct Counts {
 G a; int n,N,F; vector<I> en,st,H;
 Counts(G g):a(g),n(g.size()),N(1<<n),F(N-1),en(N*n),st(N*n),H(N){
  for(int i=0;i<n;i++){need(a[i]>=0&&a[i]<N&&!(a[i]>>i&1),"row");for(int j=0;j<i;j++)need(((a[i]>>j&1)+(a[j]>>i&1))==1,"orientation");}
  H[0]=1;for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;
  for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if(A>>w&1){if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}
 }
 I h()const{return H[F];}
 I C(int S)const{I z=0;int A=S;while(true){z+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}return z;}
 I L(int v,int S)const{I z=0;for(int A=S;;A=(A-1)&S){z+=en[A*n+v]*H[S^A];if(!A)break;}return z;}
 I R(int v,int S)const{I z=0;for(int A=S;;A=(A-1)&S){z+=st[A*n+v]*H[S^A];if(!A)break;}return z;}
 array<I,9> profile(int v)const{int A=F^(1<<v);I h=H[F],s=st[F*n+v],e=en[F*n+v],l=L(v,F),r=R(v,F),d=H[A],b=C(A),B=h+l-s,D=r+b-e-B;return {h,s,e,l,r,d,b,B,D};}
};
G graph(int n,int mask){G a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){if(mask>>k&1)a[i]|=1<<j;else a[j]|=1<<i;}return a;}
G cyclic(const G&t){int m=t.size();G a(m+2);a[0]=2;for(int i=0;i<m;i++){a[1]|=1<<(i+2);a[i+2]=(t[i]<<2)|1;}return a;}
G zmove(G a,int v){int n=a.size();for(int i=0;i<n;i++)if(i!=v)a[i]|=1<<n;a.push_back(1<<v);return a;}
bool strong(G a){int n=a.size(),F=(1<<n)-1;for(int i=0;i<n;i++)a[i]|=1<<i;for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(a[i]>>k&1)a[i]|=a[k];return all_of(a.begin(),a.end(),[&](int x){return x==F;});}
int main(){I graphs=0,roots=0,positive=0,equality=0,physicalZ=0,cyclicGraphs=0;map<int,I> eqByOrder;set<I> eqStarts;for(int n=1;n<=6;n++)for(int code=0;code<(1<<(n*(n-1)/2));code++){
 G a=graph(n,code);Counts c(a);graphs++;
 for(int v=0;v<n;v++){
  roots++;auto[h,s,e,L,R,d,b,B,D]=c.profile(v);int A=c.F^(1<<v);
  need(s<=d,"start deletion");need(L>=h+d-s,"dual union inequality");if(L==h)need(s==d,"dual equality consequence");
  if(n==1){need(s==1&&D==0,"singleton hostile");continue;}if(!s)continue;positive++;
  I q=0;for(int i=0;i<n;i++)if(a[v]>>i&1)q+=c.L(i,A);
  I f=R-h+e-d,p=d-e;need(q==b-L,"exact outgoing L identity");need(f>=0&&p>=0&&q>=d,"nonnegative decomposition");need(D-s==f+2*p+q-d,"defect decomposition");need(D>=s,"defect floor");
  bool test=false;int w=-1;if(__builtin_popcount((unsigned)a[v])==1){w=__builtin_ctz((unsigned)a[v]);test=(c.L(w,A)==d&&c.en[A*n+w]==0);need(s==c.st[A*n+w]&&p==c.en[A*n+w]&&f==0,"unique outneighbor identities");}
  need((D==s)==test,"exact graph equality test");
  if(D==s){equality++;eqByOrder[n]++;eqStarts.insert(s);need(n>=3&&w>=0,"nonempty residual");need(c.H[A^(1<<w)]==s,"equality residual count");}
  if(n<=5){Counts z(zmove(a,v));auto zp=z.profile(v);need(zp[1]==s&&zp[7]==2*B+D,"physical Z readout");physicalZ++;}
 }
 G u=cyclic(a);need(strong(u),"cyclic source strong");Counts cu(u);auto[h,s,e,L,R,d,b,B,D]=cu.profile(0);I t=c.h(),C=c.C(c.F);
 need(h==t+C&&s==t&&e==t&&d==t&&L==2*C-t&&R==t+C&&b==2*C&&B==3*C-t&&D==t,"all cyclic coordinates");cyclicGraphs++;
 }
 G triangle={2,4,1};Counts ct(triangle);need(ct.h()==3&&ct.C(ct.F)==12,"triangle values");
 G source={14,4,8,2};Counts cs(source);need(cs.h()==3&&cs.L(0,cs.F)==3&&cs.st[cs.F*4]==3,"L=h does not imply start one");
 cout<<"{\n  \"status\": \"PASS\",\n  \"labeled_tournaments_through_6\": "<<graphs<<",\n  \"rooted_inputs\": "<<roots<<",\n  \"nontrivial_positive_start_inputs\": "<<positive<<",\n  \"equality_inputs\": "<<equality<<",\n  \"physical_Z_children_through_5\": "<<physicalZ<<",\n  \"cyclic_sources\": "<<cyclicGraphs<<",\n  \"equality_by_order\": {";bool first=true;for(auto[n,c]:eqByOrder){if(!first)cout<<", ";first=false;cout<<'"'<<n<<"\": "<<c;}cout<<"},\n  \"equality_start_values\": [";first=true;for(I s:eqStarts){if(!first)cout<<", ";first=false;cout<<s;}cout<<"]\n}\n";
}
