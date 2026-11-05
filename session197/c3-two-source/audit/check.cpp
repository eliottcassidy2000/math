#include <algorithm>
#include <array>
#include <bit>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>
using namespace std;using I=long long;using G=vector<unsigned>;using V=array<I,4>;
void need(bool b,const char*m){if(!b)throw runtime_error(m);}
G opposite(G A){unsigned full=(1<<A.size())-1;for(int i=0;i<(int)A.size();i++)A[i]^=full^(1<<i);return A;}
G induced(const G&A,const vector<int>&v){G B(v.size());for(int i=0;i<(int)v.size();i++)for(int j=0;j<(int)v.size();j++)if(A[v[i]]>>v[j]&1)B[i]|=1<<j;return B;}
struct Inv{int n,N,F;vector<I>end,start,H,R,L,e,s;I C;Inv(const G&A){n=A.size();N=1<<n;F=N-1;end.assign(N*n,0);start.assign(N*n,0);H.assign(N,0);R.assign(n,0);L.assign(n,0);e.assign(n,0);s.assign(n,0);G B=opposite(A);
 for(int i=0;i<n;i++)end[(1<<i)*n+i]=start[(1<<i)*n+i]=1;
 for(int m=1;m<N;m++)for(int i=0;i<n;i++)if(m>>i&1){H[m]+=end[m*n+i];for(int mode=0;mode<2;mode++){auto&d=mode?start:end;unsigned out=(mode?B[i]:A[i])&~m;while(out){int j=countr_zero(out);out&=out-1;d[(m|(1<<j))*n+j]+=d[m*n+i];}}}
 H[0]=1;for(int i=0;i<n;i++){e[i]=end[F*n+i];s[i]=start[F*n+i];for(int m=0;m<N;m++){R[i]+=start[m*n+i]*H[F^m];L[i]+=end[m*n+i]*H[F^m];}}
 vector<I>d(N*n);d[n]=1;for(int m=1;m<N;m++)for(int i=0;i<n;i++)if(d[m*n+i]){unsigned out=A[i]&~m;while(out){int j=countr_zero(out);out&=out-1;d[(m|(1<<j))*n+j]+=d[m*n+i];}}C=0;for(int i=0;i<n;i++)if(A[i]&1)C+=d[F*n+i];}
 V root(int r)const{return{H[F],e[r],R[r],H[F^(1<<r)]};}
 V response(unsigned S)const{I E=0,U=0,K=0;for(int a=0;a<n;a++)if(S>>a&1){E+=e[a];U+=R[a];for(int b=0;b<n;b++)if(b!=a&&(S>>b&1))for(int m=0;m<N;m++)K+=end[m*n+a]*start[(F^m)*n+b];}return{H[F],E,U,K};}
};
pair<G,int> birth(G A,unsigned selected){int n=A.size();for(int i=0;i<n;i++)if(!(selected>>i&1))A[i]|=1<<n;A.push_back(selected);return{A,n};}
pair<G,int> graft(G A,unsigned selected,const G&B,int r){int n=A.size(),m=B.size();for(unsigned row:B)A.push_back(row<<n);for(int i=0;i<n;i++)for(int j=0;j<m;j++)if(j==r&&(selected>>i&1))A[n+j]|=1<<i;else A[i]|=1<<(n+j);return{A,n+r};}
pair<G,unsigned> left(int kind,G A,int c){int x=c;if(kind==2){auto z=birth(A,1<<c);A=z.first;c=z.second;}auto z=birth(A,1<<c);A=z.first;int a=z.second;
 if(kind==3){z=birth(A,1<<a);return{z.first,(1u<<z.second)|(1u<<a)|(1u<<c)};}
 z=birth(A,(1<<a)|(1<<c));unsigned S=(1u<<z.second)|(1u<<(kind==0?a:kind==1?c:x));return{z.first,S};}
V response(int k,V v){I D=v[0],E=v[1],U=v[2],d=v[3];if(k==0)return{D-E+3*U-2*d,D-E+2*U-2*d,3*D-2*E+7*U-2*d,2*D-2*E+2*U};if(k==1)return{D-E+3*U-2*d,U,2*D-E+8*U-4*d,U+E};if(k==2)return{D+4*U,U+d,4*D-E+14*U-d,U+E+2*d};return{D+2*U,D+U+d,3*D-E+8*U-d,2*D+4*U};}
const I M[4][4][4]={{{0,0,0,3},{0,0,0,-3},{0,0,12,-3},{3,0,-12,3}},{{0,0,0,3},{0,0,0,-2},{0,0,12,-6},{3,-1,-9,5}},{{0,0,0,3},{0,0,0,-1},{0,0,12,-9},{3,-2,-6,5}},{{0,0,0,3},{0,0,0,0},{0,0,12,-12},{3,-3,-3,3}}};
G amalg(const G&A,int a,const G&W,int r,int j){vector<int>av,wv;for(int i=0;i<(int)A.size();i++)if(i!=a)av.push_back(i);for(int i=0;i<(int)W.size();i++)if(i!=r)wv.push_back(i);G B=opposite(W);int tri=1+av.size(),right=tri+3,n=right+wv.size();G T(n);vector<int>am(A.size()),bm(W.size());am[a]=bm[r]=0;for(int i=0;i<(int)av.size();i++)am[av[i]]=i+1;for(int i=0;i<(int)wv.size();i++)bm[wv[i]]=right+i;
 for(int i=0;i<(int)A.size();i++)for(int k=0;k<(int)A.size();k++)if(A[i]>>k&1)T[am[i]]|=1<<am[k];for(int i=0;i<(int)B.size();i++)for(int k=0;k<(int)B.size();k++)if(B[i]>>k&1)T[bm[i]]|=1<<bm[k];
 for(int i=1;i<tri;i++)for(int k=tri;k<n;k++)T[i]|=1<<k;for(int i=tri;i<right;i++)for(int k=right;k<n;k++)T[i]|=1<<k;
 for(int i=0;i<3;i++){T[tri+i]|=1<<(tri+(i+1)%3);if(i<j)T[0]|=1<<(tri+i);else T[tri+i]|=1;}return T;}
V one{1,1,1,1};bool allowed(int k,V v,V w){if(k<4)return w!=one&&(k>=2||v!=one);int j=k-4;return(v!=one&&w!=one)||(j==0&&v!=one&&w==one)||(j==3&&v==one&&w!=one);}
I value(int k,V v,V w){if(k<4){V a=response(k,v);return a[0]*w[0]-a[1]*w[1]+a[2]*w[2]-a[3]*w[3];}I z=0;for(int i=0;i<4;i++)for(int j=0;j<4;j++)z+=v[i]*M[k-4][i][j]*w[j];return z;}
G build(int k,const G&V,const G&W){if(k<4){auto [T,S]=left(k,V,0);return graft(T,S,opposite(W),0).first;}return amalg(V,0,W,0,k-4);}
G canonical(const G&A){int n=A.size(),F=(1<<n)-1;G best;vector<int>path{0};function<void(int,unsigned)>go=[&](int v,unsigned used){if(used==(unsigned)F){if(!(A[v]&1))return;for(int shift=0;shift<n;shift++){vector<int>p;for(int j=0;j<n;j++)p.push_back(path[(shift+j)%n]);G B=induced(A,p);if(best.empty()||B<best)best=B;}return;}unsigned out=A[v]&~used;while(out){int w=countr_zero(out);out&=out-1;path.push_back(w);go(w,used|(1<<w));path.pop_back();}};go(0,1);need(!best.empty(),"strong canonical graph");return best;}
vector<int> csv(string s){replace(s.begin(),s.end(),',',' ');istringstream q(s);vector<int>a;int x;while(q>>x)a.push_back(x);return a;}
void enumerate(int n,const function<void(const G&)>&f){int e=n*(n-1)/2;for(int mask=0;mask<(1<<e);mask++){G A(n);int p=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,p++)if(mask>>p&1)A[i]|=1<<j;else A[j]|=1<<i;f(A);}}
array<I,3> constants(int L){int n=L-1,g=1;for(int x=n;x;x/=2)if(!(x&1))g+=x/2+(x/2)%2;int bits=n?32-countl_zero((unsigned)n):0,Q=1;for(int i=0;i<bits;i++)if(!(n>>i&1))Q+=1<<(bits-1-i);return{g,Q,1<<bits};}
G gadget(int L){int n=L-1,s=n?32-countl_zero((unsigned)n):0;G G0(s+1);for(int i=0;i<s;i++){for(int j=i+1;j<s;j++)G0[i]|=1<<j;if(n>>i&1)G0[i]|=1<<s;else G0[s]|=1<<i;}return G0;}
int main(){long table=0,grafts=0,cone=0;for(int n=1;n<=5;n++)enumerate(n,[&](const G&A){Inv p(A);for(int r=0;r<n;r++){V v=p.root(r);if(n>1&&p.s[r]==1){need(v[0]-v[1]>=v[3]&&v[2]>=v[0]+v[3]-v[1]&&v[2]>=2*v[3],"physical cone");cone++;}for(int k=0;k<4;k++){auto [T,S]=left(k,A,r);need(Inv(T).response(S)==response(k,v),"universal response row");table++;}}});
 vector<pair<G,Inv>>small;for(int n=1;n<=3;n++)enumerate(n,[&](const G&A){small.emplace_back(A,Inv(A));});
 for(int n=1;n<=4;n++)enumerate(n,[&](const G&A){Inv p(A);for(unsigned S=1;S<(1u<<n);S++){V v=p.response(S);I sum=0;for(int r=0;r<n;r++)if(S>>r&1)sum+=p.s[r];for(auto&[B,q]:small)for(int r=0;r<(int)B.size();r++){G T=graft(A,S,B,r).first;Inv got(T);I h=v[0]*q.H[q.F]-v[1]*q.s[r]+v[2]*q.L[r]-v[3]*q.H[q.F^(1<<r)];need(got.H[got.F]==h,"arbitrary entry graft");if(B.size()>1)need(got.C==sum*q.e[r],"arbitrary entry cycles");grafts++;}}});
 map<G,V> sources;ifstream input("unlabelled-through7.csv");string line;while(getline(input,line)){auto a=csv(line);G A(a.begin()+1,a.end());need((int)A.size()==a[0],"bank row order");Inv p(A);for(int r=0;r<(int)A.size();r++)if(p.s[r]==1){vector<int>path{r};unsigned used=1<<r;while(used!=(unsigned)p.F){unsigned out=A[path.back()]&~used;int choice=-1;while(out){int w=countr_zero(out);out&=out-1;if(p.start[(p.F^used)*p.n+w]){need(choice<0,"unique starting path");choice=w;}}need(choice>=0,"source path exists");path.push_back(choice);used|=1<<choice;}sources[induced(A,path)]=p.root(r);}}
 array<int,8>sc{};for(auto&[A,v]:sources)sc[A.size()]++;need(sc==array<int,8>{0,1,1,2,5,13,34,89},"complete source graph inventory through7");
 map<int,set<G>> made,expected;long physical=0;for(int k=0;k<8;k++)for(auto&[V,v]:sources)for(auto&[W,w]:sources){if(!allowed(k,v,w))continue;int n=V.size()+W.size()+(k==2?3:2);if(n>10)continue;G T=build(k,V,W);Inv p(T);need(p.C==3&&p.H[p.F]==value(k,v,w),"physical template scalar and cycles");made[n].insert(canonical(T));made[n].insert(canonical(opposite(T)));physical++;}
 ifstream cb("c3-through10.csv");while(getline(cb,line)){auto a=csv(line);G T(a.begin()+2,a.end());need((int)T.size()==a[0]&&a[1]==3,"C3 census row");expected[a[0]].insert(canonical(T));}need(made==expected,"full exact graph-class agreement through10");
 long labelled_count=0;for(int n=5;n<=6;n++)enumerate(n,[&](const G&A){Inv p(A);if(p.C==3){need(made[n].count(canonical(A)),"complete labelled C3 universe");labelled_count++;}});
 // Unpruned arithmetic closure, with an independently reconstructed actual source per state.
 const int cap=167;map<V,pair<G,int>> bank;vector<V>queue{one};bank[one]={{0},0};
 for(size_t at=0;at<queue.size();at++){V v=queue[at];auto [A,r]=bank[v];auto add=[&](V w,const G&B,int root){if(w[2]<=cap&&!bank.count(w)){need(w[2]>v[2],"strict source increase");bank[w]={B,root};queue.push_back(w);}};
  G sink=A;for(unsigned&row:sink)row|=1<<A.size();sink.push_back(0);add(V{v[0],0,2*v[2],v[3]},sink,r);
  for(int L=1;L<=cap/v[2];L++){auto c=constants(L);V w{c[0]*v[0]-(L%2)*v[1]+L*v[2],v[0]-(L==1)*v[1],c[1]*v[0]+c[2]*v[2],v[0]};if(w[2]>cap)continue;G B=gadget(L);auto result=graft(A,1<<r,B,B.size()-1);add(w,result.first,result.second);}
 }
 need(bank.size()==811,"independent exact source bank");set<I> support;map<I,tuple<int,V,V>> witness;long tests=0;
 for(int k=0;k<8;k++)for(auto&[v,vr]:bank)for(auto&[w,wr]:bank)if(allowed(k,v,w)){I H=value(k,v,w);need(v[2]<=H/3&&w[2]<=H/3,"uniform target cap");tests++;if(H<=501){support.insert(H);witness.try_emplace(H,k,v,w);}}
 ifstream es("expected-support.txt");set<I>wanted;I h;while(es>>h)wanted.insert(h);need(support==wanted&&support.size()==217,"all unpruned scalar targets");int maxorder=0;
 for(auto&[H,rec]:witness){auto[k,v,w]=rec;auto[V,r]=bank[v];auto[W,a]=bank[w];vector<int>pv{r},pw{a};for(int i=0;i<(int)V.size();i++)if(i!=r)pv.push_back(i);for(int i=0;i<(int)W.size();i++)if(i!=a)pw.push_back(i);G T=build(k,induced(V,pv),induced(W,pw));Inv p(T);need(p.H[p.F]==H&&p.C==3,"independent bounded witness");maxorder=max(maxorder,(int)T.size());}
 cout<<"{\n\"status\":\"PASS\",\n\"universal_response_controls\":"<<table<<",\n\"arbitrary_grafts\":"<<grafts<<",\n\"physical_source_cone_controls\":"<<cone<<",\n\"source_classes_through7\":[1,1,2,5,13,34,89],\n\"literal_templates_through10\":"<<physical<<",\n\"C3_classes_through10\":[";for(int n=5;n<=10;n++){if(n>5)cout<<",";cout<<made[n].size();}cout<<"],\n\"labelled_C3_graphs_through6\":"<<labelled_count<<",\n\"independent_source_states\":"<<bank.size()<<",\n\"unpruned_pair_tests\":"<<tests<<",\n\"bounded_witnesses\":"<<support.size()<<",\n\"max_independent_witness_order\":"<<maxorder<<"\n}\n";
}
