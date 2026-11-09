#include <bits/stdc++.h>
using namespace std;using U=uint64_t;using W=uint32_t;using Sig=array<U,7>;
int order2(U h){U z=1;for(int l=1;l<12;l++){z=z*2%h;if(z==1)return l;}return 12;}
int main(int argc,char**argv){string base=argv[1];map<pair<U,U>,set<pair<U,U>>>seeds;U all=0,allcontent=0;
for(int n=1;n<=9;n++){ifstream f(base+"/classes9_"+to_string(n)+".txt");if(!f)throw runtime_error("missing class bank");string line;set<Sig>sigs;U classes=0,h1classes=0,satclasses=0;while(getline(f,line)){istringstream ss(line);U code,savedH;if(!(ss>>code>>savedH))continue;classes++;int N=1<<n,F=N-1,z=0;vector<int>a(n),rev(n);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,z++){if(code>>z&1){a[i]|=1<<j;rev[j]|=1<<i;}else{a[j]|=1<<i;rev[i]|=1<<j;}}
vector<U>en(N*n),st(N*n),h(N);h[0]=1;for(int v=0;v<n;v++)en[(1<<v)*n+v]=st[(1<<v)*n+v]=1;
for(int s=1;s<N;s++)for(int v=0;v<n;v++)if(s>>v&1){h[s]+=en[s*n+v];int av=a[v]&~s;while(av){int bit=av&-av;av-=bit;en[(s|bit)*n+__builtin_ctz((unsigned)bit)]+=en[s*n+v];}av=rev[v]&~s;while(av){int bit=av&-av;av-=bit;st[(s|bit)*n+__builtin_ctz((unsigned)bit)]+=st[s*n+v];}}
U H=h[F],Q=0;if(H!=savedH)throw runtime_error("H mismatch");for(int s=0;s<N;s++)Q+=h[s]*h[F^s];if(H==1){h1classes++;if(Q!=U(1)<<n)throw runtime_error("nontransitive H1");continue;}if(Q%H)continue;satclasses++;
vector<U>aw(n),gw(n),K(n*n);for(int s=0;s<N;s++){int r=F^s;for(int i=0;i<n;i++){aw[i]+=st[s*n+i]*h[r];gw[i]+=en[s*n+i]*h[r];for(int j=0;j<n;j++)K[i*n+j]+=en[s*n+i]*st[r*n+j];}}
vector<U>A(N),G(N),T(N),S(N),E(N);for(int m=0;m<N;m++){if(m){int bit=m&-m,v=__builtin_ctz((unsigned)bit),r=m^bit;A[m]=A[r]+aw[v];G[m]=G[r]+gw[v];S[m]=S[r]+st[F*n+v];E[m]=E[r]+en[F*n+v];T[m]=T[r];for(int u=0;u<n;u++)if(r>>u&1)T[m]+=K[u*n+v]+K[v*n+u];}sigs.insert({H,Q,A[m],G[m],T[m],S[m],E[m]});}}
U content=0;for(auto s:sigs){U H=s[0],Q=s[1],A=s[2],G=s[3];if(gcd(H,gcd(A,G))!=1)continue;content++;if(order2(H)<12)seeds[{H,Q}].insert({A,G});}all+=sigs.size();allcontent+=content;cerr<<"n "<<n<<" classes "<<classes<<" h1 "<<h1classes<<" saturated_classes "<<satclasses<<" signatures "<<sigs.size()<<" content_one "<<content<<endl;
}for(auto [hq,ps]:seeds){for(auto [a,g]:ps)cout<<hq.first<<" "<<hq.second<<" "<<a<<" "<<g<<"\n";}cerr<<"TOTAL saturated_signatures "<<all<<" content_one "<<allcontent<<"\n";
}
