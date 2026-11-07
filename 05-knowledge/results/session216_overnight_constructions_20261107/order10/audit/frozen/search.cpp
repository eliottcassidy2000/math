#include <bits/stdc++.h>
using namespace std;using I=long long;using P=array<I,9>;using Key=array<I,9>;struct Target{int ix,mask,n;P p;};
Key key(P p){auto[h,s,e,L,R,d,C,b,T]=p;return{h,e,L,R,d,C,b,T,s+2*(s%3)};}
int main(int argc,char**argv){ifstream targets(argv[2]);int nt;targets>>nt;map<Key,vector<Target>>want;for(int k=0;k<nt;k++){Target t;targets>>t.ix>>t.mask>>t.n;for(auto&x:t.p)targets>>x;want[key(t.p)].push_back(t);}ifstream in(argv[1]);int total;in>>total;int limit=argc>3?atoi(argv[3]):total;
I graphs=0,early=0,linear=0,quad=0,candidates=0,hits=0;auto started=chrono::steady_clock::now();ofstream cand("CANDIDATES10.txt"),match("MATCHES.txt");
set<int>allowed{2,4,8,10,14,16,20,22};
for(int ix=0;ix<min(total,limit);ix++){int n;in>>n;vector<int>a(n);for(auto&v:a)in>>v;int N=1<<n,F=N-1;vector<I>en(N*n),st(N*n),H(N),C(N);H[0]=1;for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;
for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if(A>>w&1){if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}
graphs++;I d=H[F];if(d%3)continue;I minstart=LLONG_MAX;for(int i=0;i<n;i++)if(st[F*n+i])minstart=min(minstart,st[F*n+i]);if(minstart>22)continue;
for(int S=0;S<N;S++){int A=S;while(true){C[S]+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}}I b=C[F];if(b%12)continue;early++;
vector<I>EL(n),SR(n),ST(n);for(int A=1;A<N;A++)for(int u=0;u<n;u++){EL[u]+=en[A*n+u]*H[F^A];SR[u]+=st[A*n+u]*H[F^A];ST[u]+=st[A*n+u]*C[F^A];}
for(int mask=0;mask<N;mask++){I s=0,e=0,L=d,R=d,T=b;for(int u=0;u<n;u++)if(mask>>u&1){s+=st[F*n+u];R+=SR[u];T+=ST[u];}else{e+=en[F*n+u];L+=EL[u];}
linear++;if(!allowed.count(s)||e%3||L%12!=6||R%3||T%6)continue;quad++;
vector<I>E(N),S(N);for(int A=1;A<N;A++)for(int u=0;u<n;u++)if(mask>>u&1)S[A]+=st[A*n+u];else E[A]+=en[A*n+u];
I h=s+e;for(int A=1;A<N;A++)h+=E[A]*S[F^A];if(h%3)continue;
I ct=L+R-d;for(int A=1;A<N;A++)if(E[A])for(int B=F^A;B;B=(B-1)&(F^A))ct+=E[A]*S[B]*H[F^A^B];ct*=2;if(ct%12)continue;
P p{h,s,e,L,R,d,ct,b,T};candidates++;cand<<ix<<" "<<mask;for(auto v:p)cand<<" "<<v;cand<<"\n";
auto it=want.find(key(p));if(it!=want.end())for(auto t:it->second)if(t.p[1]%3!=s%3){hits++;match<<"MATCH "<<ix<<" "<<mask<<" target "<<t.ix<<" "<<t.mask<<"\n";for(auto v:p)match<<v<<" ";match<<"\n";for(auto v:t.p)match<<v<<" ";match<<"\n"<<flush;cerr<<"HIT "<<ix<<" "<<mask<<"\n";}
}
if(ix%10000==0)cerr<<"progress "<<ix<<" candidates "<<candidates<<"\n";
}
double elapsed=chrono::duration<double>(chrono::steady_clock::now()-started).count();cout<<"{\n\"graphs\":"<<graphs<<",\n\"cores_after_filters\":"<<early<<",\n\"linear_cut_checks\":"<<linear<<",\n\"quadratic_cut_checks\":"<<quad<<",\n\"candidates\":"<<candidates<<",\n\"hits\":"<<hits<<",\n\"seconds\":"<<elapsed<<"\n}\n";
}
