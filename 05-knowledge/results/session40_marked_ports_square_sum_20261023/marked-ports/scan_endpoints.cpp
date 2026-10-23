#include <bits/stdc++.h>
using namespace std;
struct Case{int c,s,t;};
string signature(const array<int,8>&A,int s,int t,vector<vector<int>>&cells){
 array<int,8> col;for(int v=0;v<8;v++)col[v]=(v==s?0:v==t?1:2);
 string sig;int prev=3;
 for(int r=0;r<8;r++){
  vector<vector<int>> keys(8);for(int v=0;v<8;v++){keys[v].push_back(col[v]);vector<int> z(prev);for(int w=0;w<8;w++)if(A[v]>>w&1)z[col[w]]++;keys[v].insert(keys[v].end(),z.begin(),z.end());}
  auto uni=keys;sort(uni.begin(),uni.end());uni.erase(unique(uni.begin(),uni.end()),uni.end());
  auto hist=keys;sort(hist.begin(),hist.end());sig.push_back(prev);for(auto &x:hist)for(int v:x)sig.push_back(v);
  for(int v=0;v<8;v++)col[v]=lower_bound(uni.begin(),uni.end(),keys[v])-uni.begin();
  int now=uni.size();if(now==prev){cells.assign(now,{});for(int v=0;v<8;v++)cells[col[v]].push_back(v);break;}prev=now;
 }
 return sig;
}
int main(int argc,char**argv){if(argc!=2){cerr<<"usage: scan_endpoints classes_input_n8.tsv\n";return 2;}ifstream f(argv[1]);string line;getline(f,line);vector<unsigned>codes;while(getline(f,line)){stringstream ss(line);unsigned code,id;ss>>id>>code;codes.push_back(code);} sort(codes.begin(),codes.end());if(codes.size()!=6880||adjacent_find(codes.begin(),codes.end())!=codes.end()){cerr<<"invalid ordinary class inventory\n";return 2;}
 unordered_map<string,vector<Case>> groups,augmented;map<vector<int>,int> shapes;long long pairs=0,discrete=0;vector<vector<int>>cells;
 for(int ci=0;ci<(int)codes.size();ci++){
  auto code=codes[ci];array<int,8>A{};int k=27;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++,k--){if(code>>k&1)A[i]|=1<<j;else A[j]|=1<<i;}
  for(int s=0;s<8;s++){
   int dp[256]{};dp[1<<s]=1<<s;for(int m=1;m<256;m++)for(int v=0;v<8;v++)if(dp[m]>>v&1){int q=A[v]&~m;while(q){int b=q&-q; q-=b;dp[m|b]|=b;}}
   for(int t=0;t<8;t++)if(dp[255]>>t&1){pairs++;auto sig=signature(A,s,t,cells);vector<int>shape;for(auto&c:cells)shape.push_back(c.size());sort(shape.begin(),shape.end());shapes[shape]++;if(shape.size()==8)discrete++;groups[sig].push_back({ci,s,t}); vector<vector<int>>tr;for(auto&b:cells)if(b.size()==3)tr.push_back(b); int chir=0;if(tr.size()==2){auto X=tr[0],Y=tr[1];int q=0;for(int v:Y)q+=(A[X[0]]>>v)&1;if(q==1||q==2){int u=X[0],w=-1,pu=-1,pw=-1;for(int v:X)if(A[u]>>v&1)w=v;for(int v:Y){if(((A[u]>>v)&1)==(q==1))pu=v;if(((A[w]>>v)&1)==(q==1))pw=v;}chir=(A[pu]>>pw&1)?1:2;}}sig.push_back(chir);augmented[sig].push_back({ci,s,t});}
  }
 }
 cout<<"classes "<<codes.size()<<" actual_endpoint_pairs "<<pairs<<" discrete "<<discrete<<" signatures "<<groups.size()<<"\n";
 for(auto &[sh,n]:shapes){cout<<"shape";for(int x:sh)cout<<" "<<x;cout<<" : "<<n<<"\n";}
 vector<vector<Case>> bad;for(auto &[key,g]:groups){set<int>cs;for(auto z:g)cs.insert(z.c);if(cs.size()>1)bad.push_back(g);} sort(bad.begin(),bad.end(),[](auto&a,auto&b){return tie(a[0].c,a[0].s,a[0].t)<tie(b[0].c,b[0].s,b[0].t);});
 int augbad=0;for(auto &[key,g]:augmented){set<int>cs;for(auto z:g)cs.insert(z.c);if(cs.size()>1)augbad++;}cout<<"augmented_signatures "<<augmented.size()<<" augmented_collision_groups "<<augbad<<"\n";
 cout<<"collision_groups "<<bad.size()<<"\n";for(auto&g:bad){cout<<"collision";for(auto z:g)cout<<" "<<codes[z.c]<<":"<<z.s<<":"<<z.t;cout<<"\n";}
return 0;
}
