// Exact augmented Hamilton counts for Q[blocks..., named child] with port bits uniform on blocks.
// Atoms: B block atoms (size s_i, ordered-cover profile q_i[1..s_i]) + S single atoms (named vertices + ports x,y).
// Count = sum over atom words (blocks visited h_i in [1,s_i] times, never consecutively; singles once) of prod q_i[h_i].
// Input: K jobs; each: B S; then B lines "s q1..qs"; then (B+S) rows of B+S chars arc matrix (atom a -> atom b).
// Output: one integer per job.
#include <bits/stdc++.h>
using namespace std; typedef unsigned long long u64;
int main(){int K;cin>>K;while(K--){int B,S;cin>>B>>S;vector<int>s(B);vector<vector<u64>>q(B);
 for(int i=0;i<B;i++){cin>>s[i];q[i].assign(s[i]+1,0);for(int h=1;h<=s[i];h++)cin>>q[i][h];}
 int A=B+S;vector<string>arc(A);for(auto&r:arc)cin>>r;
 vector<int>base(B+1,1);for(int i=0;i<B;i++)base[i+1]=base[i]*(s[i]+1);int H=base[B];
 int SM=1<<S; // state index: ((h*SM)+mask)*A+last
 vector<u64>dp((size_t)H*SM*A,0);
 auto idx=[&](int h,int m,int l){return ((size_t)h*SM+m)*A+l;};
 for(int a=0;a<A;a++){ if(a<B) dp[idx(base[a],0,a)]=1; else dp[idx(0,1<<(a-B),a)]=1; }
 // iterate states in order of total visits (h sum + popcount) to respect dependencies: process by layer
 int maxlen=0;for(int i=0;i<B;i++)maxlen+=s[i];maxlen+=S;
 vector<vector<pair<int,int>>> layer(maxlen+1);
 for(int h=0;h<H;h++){int tot=0,hh=h;for(int i=0;i<B;i++){tot+=hh%(s[i]+1);hh/=(s[i]+1);} for(int m=0;m<SM;m++){int L=tot+__builtin_popcount(m); if(L>=1&&L<=maxlen) layer[L].push_back({h,m});}}
 for(int L=1;L<maxlen;L++) for(auto [h,m]:layer[L]) for(int l=0;l<A;l++){u64 v=dp[idx(h,m,l)]; if(!v) continue;
   for(int t=0;t<A;t++){ if(arc[l][t]!='1') continue;
     if(t<B){ if(t==l) continue; int cur=(h/base[t])%(s[t]+1); if(cur>=s[t]) continue; dp[idx(h+base[t],m,t)]+=v; }
     else { int bit=1<<(t-B); if(m&bit) continue; dp[idx(h,m|bit,t)]+=v; } } }
 // final: all singles, all blocks h_i>=1
 unsigned __int128 tot=0;
 for(int h=0;h<H;h++){ u64 w=1; int hh=h; bool ok=true; for(int i=0;i<B;i++){int c=hh%(s[i]+1);hh/=(s[i]+1); if(c==0){ok=false;break;} w*=q[i][c];} if(!ok)continue;
   for(int l=0;l<A;l++) tot+=(unsigned __int128)dp[idx(h,SM-1,l)]*w; }
 u64 lo=(u64)tot; printf("%llu\n",(unsigned long long)lo);
}}
