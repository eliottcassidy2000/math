#include <bits/stdc++.h>
using namespace std;using I=long long;using P=array<I,9>;using Key=array<I,4>;
P mv(P p,char c){auto[h,s,e,L,R,d,C,b,T]=p;if(c=='O')return{h,s,0,L,2*R,d,2*C,2*b,3*T};if(c=='Z')return{h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b};return{h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C};}
struct Rec{int ix,root;string word;P p;};
int main(int argc,char**argv){ifstream in(argv[1]);int total;in>>total;map<Key,array<vector<Rec>,2>> bank;int roots=0;for(int ix=0;ix<total;ix++){int n;in>>n;vector<int>a(n);for(auto &v:a)in>>v;int N=1<<n,F=N-1;vector<I>en(N*n),st(N*n),H(N),C(N);H[0]=1;for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if(A>>w&1){if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}for(int S=0;S<N;S++){int A=S;while(true){C[S]+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}}
for(int r=0;r<n;r++){I s=st[F*n+r];if(s<=0||s>26||s%3==0||n<2)continue;roots++;I L=0,R=0,T=0;for(int A=0;A<N;A++){L+=en[A*n+r]*H[F^A];R+=st[A*n+r]*H[F^A];T+=st[A*n+r]*C[F^A];}P p{H[F],s,en[F*n+r],L,R,H[F^(1<<r)],C[F],C[F^(1<<r)],T};int residue=s%3,j=(s-residue)/3;
function<void(P,string)>walk=[&](P p,string w){I X=p[0]+p[3]-s+2*j;Key k{X,p[2],p[4]+p[7],p[5]};auto &slot=bank[k][residue-1];bool has=w.find('J')!=string::npos; if(slot.empty() || (has&&slot[0].word.find('J')==string::npos))slot={{ix,r,w,p}};if(w.size()==3)return;for(char c:string("OZJ"))walk(mv(p,c),w+c);};walk(p,"");}}
int count=0;for(auto &[k,z]:bank)if(!z[0].empty()&&!z[1].empty()){auto a=z[0][0],b=z[1][0];if(a.word.find('J')==string::npos&&b.word.find('J')==string::npos)continue;count++;cout<<"MATCH";for(auto x:k)cout<<" "<<x;cout<<"\n";for(auto v:{a,b}){cout<<v.ix<<" "<<v.root<<" "<<(v.word.empty()?"-":v.word);for(auto x:v.p)cout<<" "<<x;cout<<"\n";}}
cerr<<"roots "<<roots<<" states "<<bank.size()<<" matches "<<count<<"\n";
}
