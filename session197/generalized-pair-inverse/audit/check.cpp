#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using I=long long;using G=std::vector<unsigned>;using S=std::array<I,4>;
void need(bool x){if(!x)throw std::runtime_error("independent generalized inverse check failed");}
std::array<I,3> coef(int L){
 unsigned a=L-1;int m=0;for(unsigned x=a;x;x>>=1)++m;
 I d[3][3]={{1,0,0},{0,0,0},{0,0,0}};
 for(int i=0;i<m;++i){I z[3][3]={};bool toRoot=(a>>i)&1;
  for(int l=0;l<3;++l)for(int r=0;r<3;++r){z[toRoot?1:2][r]+=d[l][r];z[l][r?r:(toRoot?2:1)]+=d[l][r];}
  for(int l=0;l<3;++l)for(int r=0;r<3;++r)d[l][r]=z[l][r];
 }
 I g=d[0][0]+d[0][1]+d[1][0]+d[1][1],Q=1;
 // Reserved-start covers: choose the first base vertex after the root;
 // later selected vertices form its unique increasing continuation.
 for(int i=0;i<m;++i)if(!(a>>i&1))Q+=1LL<<(m-1-i);
 return {g,Q,1LL<<m};
}
G gadget(int L){unsigned a=L-1;int m=0;for(unsigned x=a;x;x>>=1)++m;G T(m+1);
 for(int i=0;i<m;++i){int v=i+1;if(a>>i&1)T[v]|=1;else T[0]|=1U<<v;
  for(int j=i+1;j<m;++j)T[v]|=1U<<(j+1);
 }return T;
}
int attach(G&T,int root,const G&D){int n=T.size(),m=D.size();
 for(auto&r:T)r|=((1U<<m)-1)<<n;for(auto r:D)T.push_back(r<<n);
 T[root]^=1U<<n;T[n]|=1U<<root;return n;
}
void sink(G&T){int n=T.size();for(auto&r:T)r|=1U<<n;T.push_back(0);}
struct Recipe{std::vector<int>Ls;int tail=0;};
std::map<S,Recipe> bank(int B){std::map<S,Recipe>A;if(!B)return A;S start={1,1,1,1};A[start]={};std::vector<S>q{start};
 for(size_t i=0;i<q.size();++i){auto s=q[i];for(int L=1;L<=B/s[2];++L){auto c=coef(L);
  S z={c[0]*s[0]-(L&1)*s[1]+L*s[2],s[0]-(L==1)*s[1],c[1]*s[0]+c[2]*s[2],s[0]};
  if(z[0]<=B&&!A.count(z)){auto r=A[s];r.Ls.push_back(L);A[z]=r;q.push_back(z);}
 }}return A;
}
std::pair<G,int> source(const Recipe&r){G T{0};int root=0;for(int L:r.Ls)root=attach(T,root,gadget(L));for(int k=0;k<r.tail;++k)sink(T);return{T,root};}
std::pair<G,int> seed(int kind,G V,int b){int n=V.size();
 if(kind==2){attach(V,b,G{0,5,1});int r=attach(V,n+1,G{0});return{V,r};}
 int a=attach(V,b,G{0}),y=attach(V,b,G{0});
 if(kind==1){attach(V,a,G{0});return{V,y};}return{V,a};
}
std::pair<G,int> build(int kind,G T,int root,int v){auto p=seed(kind,T,root);T=p.first;root=p.second;
 if(kind==3)root=attach(T,root,gadget(v));
 else for(unsigned a=v-1;a;a>>=1){if(a&1)attach(T,root,G{0});else sink(T);}
 return{T,root};
}
I formula(int k,S s,int v){I D=s[0],E=s[1],U=s[2],d=s[3],g=coef(v)[0],t=v&1;
 if(k==0)return(D-E+U)*g+(D+3*U-d)*v-U-E-t*(D-E+d);
 if(k==1)return(D+2*U)*g+(2*D-E+6*U-d)*v-E-U-t*(U+d);
 if(k==2)return(D+4*U)*g+(4*D-E+6*U)*v-t*U;
 return(D-E+3*U-2*d)*g+(2*D-E+4*U-d)*v-t*(U-d);
}
struct Record{int kind,v;S s;Recipe recipe;};
std::map<I,Record> image(int target,std::map<S,Recipe>&inputs){
 auto A=bank(target/7);for(auto [s,r]:A){if(s==S{1,0,2,1})continue;inputs.emplace(s,r);
  for(int k=1;s[0]+6*(s[2]<<k)<=target;++k){auto q=s;q[1]=0;q[2]<<=k;auto w=r;w.tail=k;inputs.emplace(q,w);}}
 std::map<I,Record>out;std::vector<I>gs(target+1);for(int v=1;v<=target;++v)gs[v]=coef(v)[0];
 // Deliberately use the larger v<=target rectangle, not the four source
 // parameter bounds. On physical sources every family has F(v)>=v.
 for(auto [s,r]:inputs)for(int k=0;k<4;++k){if(s==S{1,1,1,1}&&k!=2)continue;
  for(int v=(k==0?2:1);v<=target;++v){I D=s[0],E=s[1],U=s[2],d=s[3],g=gs[v],t=v&1,H;
   if(k==0)H=(D-E+U)*g+(D+3*U-d)*v-U-E-t*(D-E+d);
   else if(k==1)H=(D+2*U)*g+(2*D-E+6*U-d)*v-E-U-t*(U+d);
   else if(k==2)H=(D+4*U)*g+(4*D-E+6*U)*v-t*U;
   else H=(D-E+3*U-2*d)*g+(2*D-E+4*U-d)*v-t*(U-d);
   need(H>=v);if(H<=target)out.emplace(H,Record{k,v,s,r});
  }}return out;
}
std::uint64_t calls=0;
void words(const G&T,unsigned used,int last,int first,I&H,I&C){++calls;
 if(used+1==(1U<<T.size())){++H;if(first==0&&(T[last]&1))++C;return;}
 unsigned next=T[last]&~used;while(next){int v=__builtin_ctz(next);next&=next-1;words(T,used|(1U<<v),v,first,H,C);}
}
std::pair<I,I> literal(const G&T){I H=0,C=0;for(int v=0;v<(int)T.size();++v)words(T,1U<<v,v,v,H,C);return{H,C};}
struct Profile{S s;I start,C;};
Profile full(const G&T,int root){int n=T.size();std::vector<int>p(n);std::iota(p.begin(),p.end(),0);I H=0,e=0,R=0,start=0,C=0;
 do{int bad=0,rootpos=-1;for(int i=0;i<n;++i){if(p[i]==root)rootpos=i;if(i&&!(T[p[i-1]]>>p[i]&1))++bad;}
  if(!bad){++H;e+=p.back()==root;start+=p[0]==root;C+=p[0]==0&&(T[p.back()]&1);}
  int gapbad=rootpos>0&&!(T[p[rootpos-1]]>>root&1);R+=bad==gapbad;
 }while(std::next_permutation(p.begin(),p.end()));
 G D;for(int v=0;v<n;++v)if(v!=root){unsigned row=0;int x=0;for(int w=0;w<n;++w)if(w!=root){row|=((T[v]>>w)&1)<<x;++x;}D.push_back(row);}
 I d=D.empty()?1:literal(D).first;return{{H,e,R,d},start,C};
}
int main(int argc,char**argv){
 I seeds=0,outputs=0;for(int n=1;n<=4;++n){int pairs=n*(n-1)/2;
  for(int bits=0;bits<(1<<pairs);++bits){G V(n);int z=0;for(int a=0;a<n;++a)for(int b=a+1;b<n;++b,++z)if(bits>>z&1)V[a]|=1U<<b;else V[b]|=1U<<a;
   for(int r=0;r<n;++r){auto s=full(V,r);if(s.start!=1)continue;I D=s.s[0],E=s.s[1],U=s.s[2],d=s.s[3];
    std::array<S,3> want={S{D-E+3*U-2*d,U-d,2*D-E+4*U-d,D-E+U},S{3*D-2*E+6*U-2*d,D+U-d,3*D-E+8*U-d,D+2*U},S{5*D-E+9*U,D+3*U,5*D-E+10*U,D+4*U}};
    for(int k=0;k<3;++k){if(n==1&&k!=2)continue;auto A=seed(k,V,r);auto p=full(A.first,A.second);need(p.s==want[k]&&p.start==2&&p.C==(k==0?1:2));++seeds;}
    if(n<=3)for(int k=0;k<4;++k){if(n==1&&k!=2)continue;for(int v=(k==0?2:1);v<=5;++v){auto A=build(k,V,r,v);auto p=literal(A.first);need(p.first==formula(k,s.s,v)&&p.second==2);++outputs;}}
   }
  }
 }
 auto B=bank(1861);std::set<I>C1;for(auto [s,r]:B)C1.insert(s[0]);
 std::map<S,Recipe>inputs;auto Spt=image(1861,inputs);std::vector<I>missing;
 for(int n=1;n<=1861;n+=2)if(!C1.count(n)&&!Spt.count(n))missing.push_back(n);
 need(missing==std::vector<I>{7,21,67});need(B.size()==23432&&bank(265).size()==1316&&inputs.size()==2000&&Spt.size()==883);
 for(int d:{1,3}){S q={9,0,18,d};need(inputs.count(q));auto V=source(inputs[q]);need(full(V.first,V.second).s==q);
  auto N=seed(0,V.first,V.second);auto f=full(N.first,N.second);
  need(f.s[0]==(d==1?61:57)&&f.s[1]==(d==1?17:15)&&f.s[2]==(d==1?89:87));
  auto C=build(0,V.first,V.second,2);need(literal(C.first)==std::make_pair<I,I>(d==1?133:129,2));
 }
 if(argc==2){std::ofstream o(argv[1]);need(bool(o));for(auto[s,r]:B){o<<"C1";for(I x:s)o<<" "<<x;o<<"\n";}for(auto[s,r]:inputs){o<<"SRC";for(I x:s)o<<" "<<x;o<<"\n";}}
 int maxn=0;for(auto[H,w]:Spt){auto V=source(w.recipe);auto A=build(w.kind,V.first,V.second,w.v);auto p=literal(A.first);need(p.first==H&&p.second==2);maxn=std::max(maxn,(int)A.first.size());}
 std::cout<<"{\"status\":\"PASS\",\"literal_full_seed_profiles\":"<<seeds<<",\"literal_small_outputs\":"<<outputs<<",\"core_quadruples_265\":1316,\"candidate_sources\":2000,\"c1_quadruples_1861\":23432,\"complete_support_witnesses\":"<<Spt.size()<<",\"max_witness_order\":"<<maxn<<",\"literal_word_recursion_calls\":"<<calls<<",\"missing\":[7,21,67],\"support\":[";
 bool first=true;for(auto[H,w]:Spt){if(!first)std::cout<<",";first=false;std::cout<<H;}std::cout<<"]}\n";
}
