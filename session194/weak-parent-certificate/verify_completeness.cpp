#include <bits/stdc++.h>
using namespace std; using Code=uint64_t;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
Code canon(const vector<int>&r){
 int n=r.size();if(n==0)return 0;
 vector<pair<int,int>>deg;for(int v=0;v<n;v++)deg.push_back({__builtin_popcount((unsigned)r[v]),v});sort(deg.begin(),deg.end());
 vector<vector<int>>groups;for(auto [d,v]:deg){if(groups.empty()||__builtin_popcount((unsigned)r[groups.back()[0]])!=d)groups.push_back({});groups.back().push_back(v);}
 vector<int>p(n);Code best=~Code(0);
 function<void(int,int)>rec=[&](int g,int at){if(g==(int)groups.size()){
   Code z=0;for(int a:p){int row=0;for(int j=0;j<n;j++)if(r[a]>>p[j]&1)row|=1<<j;z=(z<<n)|row;}best=min(best,z);return;
  }auto q=groups[g];do{copy(q.begin(),q.end(),p.begin()+at);rec(g+1,at+q.size());}while(next_permutation(q.begin(),q.end()));};
 rec(0,0);return best;
}
int main(int argc,char**argv){
 try{
  need(argc==2,"pass bank path");ifstream f(argv[1]);need(f.good(),"missing bank");int lines;f>>lines;vector<vector<vector<int>>>bank(9);int n;
  for(int l=0;l<lines;l++){f>>n;need(n>=1&&n<=8,"bad order");vector<int>r(n);for(int&i:r)f>>i;for(int i=0;i<n;i++){need(r[i]>=0&&r[i]<(1<<n)&&!(r[i]>>i&1),"bad row");for(int j=0;j<i;j++)need(((r[i]>>j)&1)+((r[j]>>i)&1)==1,"bad pair");}bank[n].push_back(r);}
  int extra;need(!(f>>extra),"extra bank data");need(lines==7412,"wrong row count");
  vector<set<Code>>codes(9);codes[0].insert(0);
  for(int n=1;n<=8;n++){for(auto&r:bank[n])codes[n].insert(canon(r));need(codes[n].size()==bank[n].size(),"duplicate isomorphism class");}
  bank[0].push_back({});cout<<"{\"status\":\"PASS\",\"growth_steps\":[";long long total=0;
  for(int n=0;n<8;n++){
   set<Code>generated;long long presentations=0;
   for(auto&r:bank[n])for(int mask=0;mask<(1<<n);mask++){
    vector<int>rr=r;for(int i=0;i<n;i++)if(!(mask>>i&1))rr[i]|=1<<n;rr.push_back(mask);generated.insert(canon(rr));presentations++;
   }
   need(generated==codes[n+1],"growth closure mismatch at order "+to_string(n+1));total+=presentations;
   if(n)cout<<",";cout<<"{\"child_order\":"<<n+1<<",\"all_extensions\":"<<presentations<<",\"classes\":"<<generated.size()<<"}";cout.flush();
  }
  cout<<"],\"total_extensions\":"<<total<<",\"canonicalization\":\"all permutations inside outdegree classes, with full adjacency retained\"}"<<endl;
 }catch(exception&e){cerr<<e.what()<<endl;return 1;}
}
