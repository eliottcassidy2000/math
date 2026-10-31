// Input: one line per gentourng representative, "n upper-triangle-bits".
// Preserve every small-H representative through order 9, then all one-vertex
// extensions at order 10. Moon's strong-ear theorem proves this cover complete.
#define main inherited_profile_main
#include "strong_profiles.cpp"
#undef main
#include <string>

I count_h(const vector<int>& rows,I cap){
 int n=rows.size(),full=(1<<n)-1;
 vector<array<I,11>> dp(1<<n);
 for(int v=0;v<n;v++)dp[1<<v][v]=1;
 for(int m=1;m<=full;m++){
  I H=0;for(int v=0;v<n;v++)if(m>>v&1){
   int a=m^(1<<v);
   if(a)for(int u=0;u<n;u++)if((a>>u&1) && (rows[u]>>v&1))dp[m][v]+=dp[a][u];
   H+=dp[m][v];
  }
  // Every induced Hamilton path has a deterministic insertion extension.
  if(H>cap)return cap+1;
 }
 I h=0;for(int v=0;v<n;v++)h+=dp[full][v];return h;
}

void retain_profiles(const vector<int>& rows){
 int n=rows.size(),full=(1<<n)-1;
 vector<vector<I>> end(1<<n,vector<I>(n)),start=end;
 vector<I> H(1<<n);H[0]=1;
 for(int v=0;v<n;v++)end[1<<v][v]=start[1<<v][v]=1;
 for(int m=1;m<=full;m++)for(int v=0;v<n;v++)if(m>>v&1){
  int a=m^(1<<v);
  if(a)for(int u=0;u<n;u++)if(a>>u&1){
   if(rows[u]>>v&1)end[m][v]+=end[a][u];
   if(rows[v]>>u&1)start[m][v]+=start[a][u];
  }
  H[m]+=end[m][v];
 }
 I E=0;vector<I> pv(n),qv(n);vector<vector<I>> W(n,vector<I>(n));
 for(int m=0;m<=full;m++){
  E+=H[m]*H[full^m];
  for(int v=0;v<n;v++)if(m>>v&1){pv[v]+=end[m][v]*H[full^m];qv[v]+=start[m][v]*H[full^m];}
  if(m && m!=full)for(int u=0;u<n;u++)if(m>>u&1)for(int v=0;v<n;v++)if(!(m>>v&1))W[u][v]+=end[m][u]*start[full^m][v];
 }
 for(int mask=0;mask<=full;mask++){
  I p=0,q=0,s=0,t=0,C=0;
  for(int v=0;v<n;v++){
   if(mask>>v&1){q+=qv[v];s+=start[full][v];}
   else {p+=pv[v];t+=end[full][v];}
   if(mask>>v&1)for(int u=0;u<n;u++)if(!(mask>>u&1) && (rows[v]>>u&1)){
    if((W[u][v]+W[v][u])%2)throw 1;
    C+=(W[u][v]+W[v][u])/2;
   }
  }
  Profile b={H[full],E,p,q,s,t,H[full]+2*C};
  if(!profiles.count(b))profiles[b]={rows,mask};
 }
}
int main(int argc,char**argv){
 int cap=argc>1?stoi(argv[1]):127,seedmax=argc>2?stoi(argv[2]):9;
 map<int,uint64_t> total,small;vector<vector<int>> last;
 int n;string bits;
 while(cin>>n>>bits){
  vector<int> rows(n);int at=0;
  if(n>1 && bits.size()!=n*(n-1)/2)throw 2;
  for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){
   if(bits[at++]=='1')rows[u]|=1<<v;else rows[v]|=1<<u;
  }
  total[n]++;if(!strong(rows))throw 3;
  if(count_h(rows,cap)>cap)continue;
  small[n]++;retain_profiles(rows);if(n==seedmax)last.push_back(rows);
 }
 cerr<<"order"<<seedmax<<" retained="<<last.size()<<" profiles="<<profiles.size()<<"\n";
 for(int order=seedmax+1;order<=10;order++){
  vector<vector<int>> next;
  for(const auto &row:last)for(int mask=0;mask<(1<<(order-1));mask++){
   vector<int> rows=row;rows.push_back(mask);
   for(int v=0;v<order-1;v++)if(!(mask>>v&1))rows[v]|=1<<(order-1);
   total[order]++;if(!strong(rows) || count_h(rows,cap)>cap)continue;
   small[order]++;retain_profiles(rows);next.push_back(rows);
  }
  last=next;cerr<<"order"<<order<<" retained="<<last.size()<<" profiles="<<profiles.size()<<"\n";
 }
 int Q=2*(cap+1);
 cout<<"{\"cap_H\":"<<cap<<",\"max_order\":10,\"gauge\":\"gentourng_then_strong_ears\",\"seed_max_order\":"<<seedmax<<",\"profile_filter\":{\"q_max\":"<<Q<<",\"F_max\":"<<2*Q<<"},\"raw_profiles\":"<<profiles.size()<<",\"levels\":[";bool first=true;
 for(auto [n,t]:total){if(!first)cout<<",";first=false;cout<<"{\"n\":"<<n<<",\"tested\":"<<t<<",\"retained\":"<<small[n]<<"}";}
 cout<<"],\"profiles\":[";first=true;
 for(auto [b,w]:profiles){if(b[3]>Q || b[6]>2*Q)continue;if(!first)cout<<",";first=false;
  cout<<"{\"profile\":[";for(int i=0;i<7;i++){if(i)cout<<",";cout<<b[i];}
  cout<<"],\"mask\":"<<w.mask<<",\"rows\":[";for(int i=0;i<w.rows.size();i++){if(i)cout<<",";cout<<w.rows[i];}cout<<"]}";
 }
 cout<<"]}\n";
 return 0;
}
