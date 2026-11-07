#include <bits/stdc++.h>
using namespace std;using I=long long;using P=array<I,9>;using Key=array<I,4>;
P mv(P p,char c){auto[h,s,e,L,R,d,C,b,T]=p;if(c=='O')return{h,s,0,L,2*R,d,2*C,2*b,3*T};if(c=='Z')return{h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b};return{h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C};}
struct Rec{int ix,root;string word;P p;};
int main(int argc,char**argv){ifstream in(argv[1]);int total;in>>total;vector<pair<Rec,int>> sources,rights;int roots=0;for(int ix=0;ix<total;ix++){int n;in>>n;vector<int>a(n);for(auto &v:a)in>>v;int N=1<<n,F=N-1;vector<I>en(N*n),st(N*n),H(N),C(N);H[0]=1;for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if(A>>w&1){if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}for(int S=0;S<N;S++){int A=S;while(true){C[S]+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}}
for(int r=0;r<n;r++){I s=st[F*n+r];if(n<2)continue;roots++;I L=0,R=0,T=0;for(int A=0;A<N;A++){L+=en[A*n+r]*H[F^A];R+=st[A*n+r]*H[F^A];T+=st[A*n+r]*C[F^A];}P p{H[F],s,en[F*n+r],L,R,H[F^(1<<r)],C[F],C[F^(1<<r)],T};if(n<=6)rights.push_back({{ix,r,"",p},n});if(s>0&&s<=26&&s%3)sources.push_back({{ix,r,"",p},n});}}
map<Key,array<pair<int,int>,2>> bank;
long long transitions=0;int matches=0;
for(int ai=0;ai<(int)sources.size();ai++){auto a=sources[ai].first;I s=a.p[1],r=s%3,j=(s-r)/3;
Key key{a.p[0]+a.p[3]-s+2*j,a.p[2],a.p[4]+a.p[7],a.p[5]};
auto it=bank.find(key);if(it==bank.end())it=bank.emplace(key,array<pair<int,int>,2>{make_pair(-1,-1),make_pair(-1,-1)}).first;
it->second[r-1]={ai,-1};}

for(int ai=0;ai<(int)sources.size();ai++){auto a=sources[ai].first;auto [h,s,e,L,R,d,C,b,T]=a.p;I B=h+L-s,V=R+b;
for(int wi=0;wi<(int)rights.size();wi++){auto w=rights[wi].first;auto[hw,sw,ew,Lw,Rw,q,Cw,bw,Tw]=w.p;
if(q==1||s*q>26||s*q%3==0)continue;I sp=s*q,rp=sp%3,j=(sp-rp)/3;
I BP=d*hw+(V-d)*Lw+(e-d)*sw+(B-e-V+d)*q;
Key key{BP+2*j,d*ew,bw*V+d*(Rw-bw),d*q};transitions++;
auto it=bank.find(key);if(it==bank.end())it=bank.emplace(key,array<pair<int,int>,2>{make_pair(-1,-1),make_pair(-1,-1)}).first;
it->second[rp-1]={ai,wi};}}
for(auto &[k,reps]:bank)if(reps[0].first>=0&&reps[1].first>=0&&(reps[0].second>=0||reps[1].second>=0)){matches++;cout<<"MATCH";for(auto x:k)cout<<" "<<x;cout<<"\n";for(auto [ai,wi]:reps){auto a=sources[ai].first;cout<<a.ix<<" "<<a.root<<" ";if(wi<0)cout<<"-1 -1\n";else{auto w=rights[wi].first;cout<<w.ix<<" "<<w.root<<"\n";}}}
cerr<<"sources "<<sources.size()<<" rights "<<rights.size()<<" transitions "<<transitions<<" keys "<<bank.size()<<" matches "<<matches<<"\n";
}
