#include <bits/stdc++.h>
using namespace std;
constexpr int n=10,P=101; int D, rr;
struct Forest {int sign; uint64_t mask; vector<pair<int,int>> e;};
using Key=pair<vector<int>,vector<array<int,3>>>;
vector<pair<int,int>> pairs;
map<Key,Forest> generate(){
 map<Key,Forest> reps; int ix[n][n]={};for(int i=0;i<(int)pairs.size();i++)ix[pairs[i].first][pairs[i].second]=i;
 vector<pair<int,int>> matching; vector<int> U;for(int i=2*rr;i<n;i++)U.push_back(i); vector<int> lens(rr);
 function<void(unsigned)> match=[&](unsigned left){if(left){int x=__builtin_ctz(left);left^=1u<<x;for(int y=x+1;y<2*rr;y++)if(left>>y&1){matching.push_back({x,y});match(left^(1u<<y));matching.pop_back();}return;}
  auto perm=U;do {function<void(int,int)> alloc=[&](int i,int rem){if(i==rr-1){lens[i]=rem;vector<int>S;vector<array<int,3>>T;Forest f{1,0,{}};int pos=0;for(int k=0;k<rr;k++){auto [a,b]=matching[k];int last=a;for(int z=0;z<lens[k];z++){int v=perm[pos++],lo=min(last,v),hi=max(last,v);f.mask|=1ull<<ix[lo][hi];f.sign*=last<v?1:-1;last=v;}int lo=min(last,b),hi=max(last,b);f.mask|=1ull<<ix[lo][hi];f.sign*=last<b?1:-1;T.push_back({a,b,lens[k]%2});}sort(T.begin(),T.end());auto key=make_pair(S,T);auto it=reps.find(key);if(it==reps.end()||f.mask<it->second.mask){for(int j=0;j<(int)pairs.size();j++)if(f.mask>>j&1)f.e.push_back(pairs[j]);reps[key]=f;}return;}for(int x=0;x<=rem;x++){lens[i]=x;alloc(i+1,rem-x);}};alloc(0,U.size());}while(next_permutation(perm.begin(),perm.end()));
 };match((1u<<(2*rr))-1);return reps;
}

int pw(int a,int b){int r=1;for(;b;b>>=1,a=a*a%P)if(b&1)r=r*a%P;return r;}
int coeff(int c,int r){int f=1;for(int i=1;i<=c;i++)f=f*i%P;return f*pw(2,r)%P*pw(pw(2,c-1),P-2)%P;}
int vals[100][100];
int coef(const Forest &f,const vector<pair<int,int>>&g,int N){
 int adj[100][2],deg[100]={},seen[100]={};
 auto add=[&](int a,int b){if(deg[a]>=2||deg[b]>=2)return false;adj[a][deg[a]++]=b;adj[b][deg[b]++]=a;return true;};
 for(auto[a,b]:f.e)if(!add(a,b))return 0;for(auto[a,b]:g)if(!add(a,b))return 0;
 int c=0,r=0,s=1;
 for(int i=0;i<N;i++)if(!seen[i]&&deg[i]<=1){c++;seen[i]=1;if(!deg[i])continue;r++;int at=i,prev=-1,len=0;while(true){int nx=-1;for(int k=0;k<deg[at];k++)if(adj[at][k]!=prev)nx=adj[at][k];if(nx<0)break;if(seen[nx])return 0;seen[nx]=1;if(at>nx)s=-s;len++;prev=at;at=nx;}if(len%2)return 0;}
 for(int i=0;i<N;i++)if(!seen[i])return 0;
 int v=vals[c][r];return s*f.sign<0?(P-v)%P:v;
}
vector<vector<int>> components(const Forest&f){vector<vector<int>> adj(n);for(auto[a,b]:f.e){adj[a].push_back(b);adj[b].push_back(a);}vector<vector<int>>cs;vector<int>seen(n);for(int i=0;i<n;i++)if(adj[i].size()<=1&&!seen[i]){vector<int>v;int at=i,prev=-1;while(at>=0){v.push_back(at);seen[at]=1;int nx=-1;for(int b:adj[at])if(b!=prev)nx=b;prev=at;at=nx;}cs.push_back(v);}return cs;}
int main(int argc,char**argv){
 if(argc!=2)return 5; rr=atoi(argv[1]);if(rr<1||rr>5)return 6; D=(int[]){0,1,6,60,735,945}[rr]; string base="block"+to_string(rr);
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)pairs.emplace_back(i,j);auto reps=generate();cerr<<"representatives "<<reps.size()<<"\n";if(reps.size()!=D)return 2;
 vector<Forest> fs;for(auto &z:reps)fs.push_back(z.second);
 ofstream rep(base+".reps");for(auto &f:fs)rep<<f.mask<<' '<<f.sign<<'\n';rep.close();
 for(int c=1;c<100;c++)for(int r=0;r<=c;r++)vals[c][r]=coeff(c,r);
 vector<vector<uint8_t>> basis(D);uint8_t mult[P][P];for(int i=0;i<P;i++)for(int j=0;j<P;j++)mult[i][j]=i*j%P;
 mt19937 rng(146160); int rank=0;ofstream cert(base+".contexts");
 for(int trial=0;trial<100000&&rank<D;trial++){
 auto cs=components(fs[trial%D]);int N=n;vector<pair<int,int>>g;
 shuffle(cs.begin(),cs.end(),rng);
 while(cs.size()>1&&rng()%4!=0){auto a=cs.back();cs.pop_back();auto b=cs.back();cs.pop_back();if(rng()%2)reverse(a.begin(),a.end());if(rng()%2)reverse(b.begin(),b.end());int last=a.back();int count=1+rng()%3;for(int t=0;t<count;t++){g.push_back({last,N});a.push_back(N);last=N++;}g.push_back({last,b.front()});a.insert(a.end(),b.begin(),b.end());cs.push_back(a);shuffle(cs.begin(),cs.end(),rng);}
 for(auto &v:cs){int left=1+rng()%3,right=1+rng()%3;if((v.size()-1+left+right)%2)right++;int at=v.front();while(left--){g.push_back({at,N});at=N++;}at=v.back();while(right--){g.push_back({at,N});at=N++;}}
 int paths=rng()%3;while(paths--){int a=N++,b=N++,c=N++;g.push_back({a,b});g.push_back({b,c});}N+=rng()%3;
 vector<uint8_t>v(D);for(int j=0;j<D;j++)v[j]=coef(fs[j],g,N);
 for(int j=0;j<D;j++)if(v[j]){if(basis[j].empty()){int inv=pw(v[j],P-2);for(int k=j;k<D;k++)v[k]=mult[inv][v[k]];basis[j]=v;rank++;cert<<N<<' '<<g.size();for(auto[a,b]:g)cert<<' '<<a<<' '<<b;cert<<'\n';if(rank%100==0||rank==D){cert.flush();cerr<<"trial "<<trial<<" rank "<<rank<<"\n";}break;}else{auto m=mult[v[j]];auto &b=basis[j];for(int k=j;k<D;k++){int a=v[k]-m[b[k]];v[k]=a+(a<0?P:0);}}}
 if(trial%5000==0)cerr<<"progress "<<trial<<' '<<rank<<'\n';
 }
 cerr<<"FINAL "<<rank<<'\n';return rank==D?0:3;
}
