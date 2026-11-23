#include <bits/stdc++.h>
using namespace std;using U=uint64_t;
int n,F;vector<pair<int,int>> edges;vector<vector<int>> idx;vector<array<unsigned,12>> d;unordered_map<U,int> memo;
int ht(U m){
 auto a=memo.find(m);if(a!=memo.end())return a->second;
 array<int,12> rows{};
 for(int k=0;k<edges.size();k++){auto [i,j]=edges[k];if(m>>k&1)rows[i]|=1<<j;else rows[j]|=1<<i;}
 for(int s=1;s<=F;s++)for(unsigned w=s;w;w&=w-1){
  int i=__builtin_ctz(w),t=s^(1<<i);unsigned z=0;
  if(!t)z=1;else for(unsigned b=t&rows[i];b;b&=b-1)z+=d[t][__builtin_ctz(b)];
  d[s][i]=z;
 }
 int h=0;for(int i=0;i<n;i++)h+=d[F][i];memo.emplace(m,h);return h;
}
U start;vector<int> script, heights;unsigned long long hits=0, labelcount=0;int E,extra; vector<int> cuts(5); unsigned long long progress[12]{};
bool arc(U m,int i,int j){return bool(m>>idx[i][j]&1)==(i<j);}
bool iso(U A,U B){
 int da[12]{},db[12]{},p[12]{};
 for(auto [i,j]:edges){if(arc(A,i,j))da[i]++;else da[j]++;if(arc(B,i,j))db[i]++;else db[j]++;}
 auto rec=[&](auto&&f,int i,int s)->bool{
  if(i==n)return true;
  for(int j=0;j<n;j++)if(!(s>>j&1)&&da[i]==db[j]){
   bool ok=true;for(int t=0;t<i;t++)if(arc(A,t,i)!=arc(B,p[t],j)){ok=false;break;}
   if(ok){p[i]=j;if(f(f,i+1,s|1<<j))return true;}
  }
  return false;
 };
 return rec(rec,0,0);
}
void walk(U m,int step){
 ++progress[step];
 if(step==heights.size()){
  if(m==start||iso(start,m)){
   if(hits++<2 || (m==start&&!labelcount)){
    cout<<"HIT order "<<n<<" cuts ";for(int t:cuts)cout<<t; cout<<" labelled "<<(m==start)<<" start "<<start<<" ops ";
    for(int k:script)cout<<edges[k].first<<","<<edges[k].second<<" ";cout<<" heights ";
    U z=start;cout<<9<<" ";for(int k:script){z^=1ull<<k;cout<<ht(z)<<" ";}cout<<endl;
   }
   if(m==start)++labelcount;
  }
  return;
 }
 int goal=heights[step];
 for(int k=0;k<E;k++){
  U z=m^(1ull<<k);
  if(goal&&ht(z)!=goal)continue;
  script.push_back(k); walk(z,step+1); script.pop_back();
 }
}
U player(vector<vector<int>> C){
 U ans=0;array<int,12> score{};int s=0;
 for(auto block:C){
  for(int i=0;i<block.size();i++)for(int j=i+1;j<block.size();j++)if(block[i]>>j&1)ans|=1ull<<idx[s+i][s+j];
  for(int i=s;i<s+block.size();i++)for(int j=s+block.size();j<n;j++)ans|=1ull<<idx[i][j];
  s+=block.size();
 }
 if(s!=n)abort();
 return ans;
}
int main(int argc,char**argv){
 n=atoi(argv[1]);F=(1<<n)-1;d.resize(F+1);idx.assign(n,vector<int>(n));
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){idx[i][j]=idx[j][i]=edges.size();edges.emplace_back(i,j);}E=edges.size();
 vector<U> S;
 for(int l=0;l<=n-5;l++)for(auto core:{vector<int>{6,4,8,19,7},vector<int>{2,4,9,19,7}}){
  vector<vector<int>> C(l,{0});C.push_back(core);C.insert(C.end(),n-5-l,{0});S.push_back(player(C));
 }
 for(int l=0;l<=n-6;l++)for(int r=0;r<=n-6-l;r++){
  vector<vector<int>> C(l,{0});C.push_back({2,4,1});C.insert(C.end(),r,{0});C.push_back({2,4,1});C.insert(C.end(),n-6-l-r,{0});S.push_back(player(C));
 }
 extra=argc>2?atoi(argv[2]):1;
 auto c=[&](auto&&f,int i,int a)->void{
  if(i==5){
   if(a)return;
   heights.clear();int h[]={81,63,27,45,9};
   for(int j=0;j<5;j++){heights.insert(heights.end(),cuts[j],0);heights.push_back(h[j]);}
   for(U v:S){start=v;if(ht(v)!=9)abort();walk(v,0);}
   cout<<"PART ";for(int t:cuts)cout<<t;cout<<" total hits "<<hits<<" progress";for(int j=0;j<=5+extra;j++)cout<<" "<<progress[j];cout<<endl;
   fill(begin(progress),end(progress),0);return;
  }
  for(int d=0;d<=a;d++){cuts[i]=d;f(f,i+1,a-d);}
 };c(c,0,extra);
 cout<<"SUMMARY "<<n<<" starters "<<S.size()<<" hits "<<hits<<" labelled "<<labelcount<<" cache "<<memo.size()<<endl;
}
