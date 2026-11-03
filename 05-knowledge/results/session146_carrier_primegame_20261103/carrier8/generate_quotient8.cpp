#include <bits/stdc++.h>
using namespace std;
constexpr int n=8, D=2233,P=101;
struct Forest {int sign; uint32_t mask; vector<pair<int,int>> e;};
using Key=pair<vector<int>,vector<array<int,3>>>;
vector<pair<int,int>> pairs;
map<Key,Forest> generate(){
 map<Key,Forest> reps; int ix[n][n]={};for(int i=0;i<(int)pairs.size();i++)ix[pairs[i].first][pairs[i].second]=i;
 vector<int>w(n);iota(w.begin(),w.end(),0);
 do {for(int cut=0;cut<(1<<(n-1));cut++) {vector<array<int,3>> T;vector<int>S; Forest f{1,0,{}};bool valid=true;
 for(int a=0;a<n;){int b=a;while(b<n-1&&!(cut>>b&1))b++;if(a==b){valid=false;break;} int l=w[a],r=w[b],sign=1;T.push_back({min(l,r),max(l,r),(b-a-1)%2});for(int j=a;j<b;j++){int x=w[j],y=w[j+1];int lo=min(x,y),hi=max(x,y);f.mask |= (1u<<ix[lo][hi]);sign*=x<y?1:-1;}if(l>r&&(b-a)%2)sign=-sign;f.sign*=sign;a=b+1;}
 if(!valid)continue;sort(T.begin(),T.end());auto key=make_pair(S,T);auto it=reps.find(key);if(it==reps.end()||f.mask<it->second.mask){f.e.clear();for(int j=0;j<(int)pairs.size();j++)if(f.mask>>j&1)f.e.push_back(pairs[j]);reps[key]=f;}
 }}while(next_permutation(w.begin(),w.end()));return reps;
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
 ofstream rep("reps8.generated.txt");for(auto &f:fs)rep<<f.mask<<' '<<f.sign<<'\n';rep.close();
 for(int c=1;c<100;c++)for(int r=0;r<=c;r++)vals[c][r]=coeff(c,r);
 vector<vector<uint8_t>> basis(D);uint8_t mult[P][P];for(int i=0;i<P;i++)for(int j=0;j<P;j++)mult[i][j]=i*j%P;
 mt19937 rng(146150); int rank=0;ofstream cert("contexts8.generated.txt");
 ifstream start("start.txt"); string line;
 while(getline(start,line)) {istringstream in(line);int N,k;in>>N>>k;vector<pair<int,int>> g(k);for(auto &e:g)in>>e.first>>e.second;vector<uint8_t>v(D);for(int j=0;j<D;j++)v[j]=coef(fs[j],g,N);bool added=false;for(int j=0;j<D;j++)if(v[j]){if(basis[j].empty()){int inv=pw(v[j],P-2);for(int k=j;k<D;k++)v[k]=mult[inv][v[k]];basis[j]=v;rank++;cert<<line<<'\n';added=true;break;}else{auto m=mult[v[j]];auto &b=basis[j];for(int k=j;k<D;k++){int a=v[k]-m[b[k]];v[k]=a+(a<0?P:0);}}}if(!added){cerr<<"bad start"<<endl;return 4;}}
 cerr<<"START rank "<<rank<<endl; cert.flush();
 for(int trial=0;trial<100000&&rank<D;trial++){
 vector<int> missing;for(int j=0;j<D;j++)if(basis[j].empty())missing.push_back(j);auto cs=components(fs[missing[rng()%missing.size()]]);int N=n;vector<pair<int,int>>g;
 shuffle(cs.begin(),cs.end(),rng);
 while(cs.size()>1&&rng()%4!=0){auto a=cs.back();cs.pop_back();auto b=cs.back();cs.pop_back();if(rng()%2)reverse(a.begin(),a.end());if(rng()%2)reverse(b.begin(),b.end());int last=a.back();int count=1+rng()%3;for(int t=0;t<count;t++){g.push_back({last,N});a.push_back(N);last=N++;}g.push_back({last,b.front()});a.insert(a.end(),b.begin(),b.end());cs.push_back(a);shuffle(cs.begin(),cs.end(),rng);}
 for(auto &v:cs){int left=rng()%4,right=rng()%4;if((v.size()-1+left+right)%2)right++;int at=v.front();while(left--){g.push_back({at,N});at=N++;}at=v.back();while(right--){g.push_back({at,N});at=N++;}}
 int paths=rng()%3;while(paths--){int a=N++,b=N++,c=N++;g.push_back({a,b});g.push_back({b,c});}N+=rng()%3;
 vector<uint8_t>v(D);for(int j=0;j<D;j++)v[j]=coef(fs[j],g,N);
 for(int j=0;j<D;j++)if(v[j]){if(basis[j].empty()){int inv=pw(v[j],P-2);for(int k=j;k<D;k++)v[k]=mult[inv][v[k]];basis[j]=v;rank++;cert<<N<<' '<<g.size();for(auto[a,b]:g)cert<<' '<<a<<' '<<b;cert<<'\n';if(rank%100==0||rank==D){cert.flush();cerr<<"trial "<<trial<<" rank "<<rank<<"\n";}break;}else{auto m=mult[v[j]];auto &b=basis[j];for(int k=j;k<D;k++){int a=v[k]-m[b[k]];v[k]=a+(a<0?P:0);}}}
 if(trial%5000==0)cerr<<"progress "<<trial<<' '<<rank<<'\n';
 }
 cerr<<"FINAL "<<rank<<'\n';return rank==D?0:3;
}
