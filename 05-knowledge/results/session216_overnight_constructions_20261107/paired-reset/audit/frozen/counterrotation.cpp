#include<bits/stdc++.h>
using namespace std;using I=long long;using P=array<I,7>;
P mv(P p,char c){auto[h,s,e,L,R,d,b]=p;if(c=='O')return{h,s,0,L,2*R,d,2*b};if(c=='Z')return{h+R-e,s,d,L+b,2*R-d,d,2*b};if(c=='P')return{h,0,e,2*L,R,d,2*b};return{h+L-s,d,e,2*L-d,R+b,d,2*b};}
I reversebits(I a,int m){I r=0;for(int i=0;i<m;i++){r=(r<<1)|(a&1);a>>=1;}return r;}
P transitive_root(I incoming,int m){P p{1,1,1,1,1,1,1};for(int i=0;i<m;i++){bool z=incoming>>i&1;p=mv(p,z?'Z':'O');if(i==0&&z)p[1]=0;}return p;}
int main(int argc,char**argv){int cap=argc>1?atoi(argv[1]):18;ifstream in("H5_CORE_STATES.txt");vector<pair<int,P>>cores;int mask;P p;while(in>>mask){for(auto&v:p)in>>v;cores.push_back({mask,p});}long long states=0,bitmatches=0,hits=0;
for(int N=0;N<=cap;N++){int m=N+5;for(int l=0;l<=N;l++){int r=N-l;for(auto [mask,p]:cores){function<void(P,int,string)>right=[&](P p,int left,string word){if(left){right(mv(p,'O'),left-1,word+'O');right(mv(p,'Z'),left-1,word+'Z');return;}states++;auto[h,s,e,L,R,d,b]=p;if(s!=5||e!=5)return;I a=L-d,rr=R-d;if(a<0||a>=b)return;if((b-1)-reversebits(a,m)!=rr)return;bitmatches++;P A=transitive_root(a,m);if(A[1]!=1||A[2]!=1||A[3]-1!=a||A[4]-1!=rr)throw runtime_error("bit compiler mismatch");if(A[0]!=h+2)return;hits++;cout<<"MATCH "<<N<<" "<<l<<" "<<r<<" "<<mask<<" "<<word<<" incoming "<<a<<"\n";for(auto q:{A,p}){for(auto z:q)cout<<z<<" ";cout<<"\n";}cout<<flush;};
function<void(P,int,string)>left=[&](P p,int rem,string word){if(rem){left(mv(p,'P'),rem-1,word+'P');left(mv(p,'Q'),rem-1,word+'Q');return;}right(p,r,word);};left(p,l,"");}}
cerr<<"padding "<<N<<" states "<<states<<" bitmatches "<<bitmatches<<" hits "<<hits<<"\n";if(hits)return 0;}
}
