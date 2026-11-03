#include <bits/stdc++.h>
using namespace std;
constexpr int n=7, D=2668,P=101;
struct Forest {int sign; uint32_t mask; vector<pair<int,int>> e;};
using Key=pair<vector<int>,vector<array<int,3>>>;
vector<pair<int,int>> pairs;
map<Key,Forest> generate(){
 map<Key,Forest> reps;
 for(uint32_t mask=0; mask<(1u<<21); mask++){
  int deg[n]={},adj[n][2]={};bool ok=true; Forest f{1,mask,{}};
  for(int i=0;i<21;i++)if(mask>>i&1){auto [a,b]=pairs[i];if(deg[a]==2||deg[b]==2){ok=false;break;}adj[a][deg[a]++]=b;adj[b][deg[b]++]=a;f.e.emplace_back(a,b);} if(!ok)continue;
  bool seen[n]={}; vector<int>S;vector<array<int,3>> T;
  for(int i=0;i<n;i++)if(!seen[i]&&deg[i]<=1){if(!deg[i]){seen[i]=1;S.push_back(i);continue;}int at=i,prev=-1,len=0,sg=1;seen[at]=1;while(true){int nx=-1;for(int k=0;k<deg[at];k++)if(adj[at][k]!=prev)nx=adj[at][k];if(nx<0)break;if(seen[nx]){ok=false;break;}sg*=at<nx?1:-1;len++;seen[nx]=1;prev=at;at=nx;}if(!ok)break; T.push_back({i,at,(len-1)%2});f.sign*=sg;}
  if(!ok||count(seen,seen+n,true)!=n)continue;sort(T.begin(),T.end());reps.emplace(make_pair(S,T),f);
 }
 return reps;
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
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)pairs.emplace_back(i,j);auto reps=generate();cerr<<"representatives "<<reps.size()<<"\n";if(reps.size()!=D)return 2;
 vector<Forest> fs;for(auto &z:reps)fs.push_back(z.second);
 ofstream rep("reps7.txt");for(auto &f:fs)rep<<f.mask<<' '<<f.sign<<'\n';rep.close();
 for(int c=1;c<100;c++)for(int r=0;r<=c;r++)vals[c][r]=coeff(c,r);
 vector<vector<uint8_t>> basis(D);uint8_t mult[P][P];for(int i=0;i<P;i++)for(int j=0;j<P;j++)mult[i][j]=i*j%P;
 mt19937 rng(146146); int rank=0;ofstream cert("contexts7.generated.txt");
 for(int trial=0;trial<100000&&rank<D;trial++){
 auto cs=components(fs[trial%D]);int N=n;vector<pair<int,int>>g;
 shuffle(cs.begin(),cs.end(),rng);
 while(cs.size()>1&&rng()%4!=0){auto a=cs.back();cs.pop_back();auto b=cs.back();cs.pop_back();if(rng()%2)reverse(a.begin(),a.end());if(rng()%2)reverse(b.begin(),b.end());int last=a.back();int count=1+rng()%3;for(int t=0;t<count;t++){g.push_back({last,N});a.push_back(N);last=N++;}g.push_back({last,b.front()});a.insert(a.end(),b.begin(),b.end());cs.push_back(a);shuffle(cs.begin(),cs.end(),rng);}
 for(auto &v:cs){if(rng()%2)reverse(v.begin(),v.end());int count=(v.size()-1)%2;if(rng()%3==0)count+=2;while(count--){g.push_back({v.back(),N});v.push_back(N++);}}N+=rng()%3;
 vector<uint8_t>v(D);for(int j=0;j<D;j++)v[j]=coef(fs[j],g,N);
 for(int j=0;j<D;j++)if(v[j]){if(basis[j].empty()){int inv=pw(v[j],P-2);for(int k=j;k<D;k++)v[k]=mult[inv][v[k]];basis[j]=v;rank++;cert<<N<<' '<<g.size();for(auto[a,b]:g)cert<<' '<<a<<' '<<b;cert<<'\n';if(rank%100==0||rank==D){cert.flush();cerr<<"trial "<<trial<<" rank "<<rank<<"\n";}break;}else{auto m=mult[v[j]];auto &b=basis[j];for(int k=j;k<D;k++){int a=v[k]-m[b[k]];v[k]=a+(a<0?P:0);}}}
 if(trial%5000==0)cerr<<"progress "<<trial<<' '<<rank<<'\n';
 }
 cerr<<"FINAL "<<rank<<'\n';return rank==D?0:3;
}
