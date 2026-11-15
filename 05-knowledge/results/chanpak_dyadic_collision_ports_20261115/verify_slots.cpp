#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;

uint32_t arcs[30];
int n, path[30], place[30][30], before[30][30], after[30][30], count_h;
long long words = 0, marked_slots = 0, profiles = 0;

void dfs(int depth, uint32_t remaining) {
  if (!remaining) {
    ++count_h;
    for (int j=0;j<n;++j) {
      int root=path[j]; uint32_t left=0,right=0;
      for(int i=0;i+1<n;++i) if(i!=j && i+1!=j) {
        int u=place[root][path[i]], v=place[root][path[i+1]];
        if(u>v) (i<j?left:right) |= ((1u<<u)-1)^((1u<<v)-1);
      }
      for(int i=0;i<n-2;++i) {before[root][i]+=(left>>i)&1; after[root][i]+=(right>>i)&1;}
    }
    return;
  }
  for(uint32_t s=depth?remaining&arcs[path[depth-1]]:remaining;s;s&=s-1) {
    path[depth]=__builtin_ctz(s);
    dfs(depth+1,remaining^(s&-s));
  }
}

void test(const vector<int>& links) {
  vector<int> junctions{0}, interior;
  int h=1;
  for (int s:links) {
    int a=junctions.back(),b=a+s-1;
    for(int j=a+1;j<b;++j) interior.push_back(j);
    junctions.push_back(b);
    h*=1+(1<<(s-2));
  }
  n=junctions.back()+1; if(n>14)return;
  for(int i=0;i<n;++i) {
    arcs[i]=0;
    for(int j=i+1;j<n;++j) arcs[i]|=1u<<j;
    for(int j=0;j<n;++j) before[i][j]=after[i][j]=0;
  }
  for(int i=1;i<(int)junctions.size();++i) {
    int a=junctions[i-1],b=junctions[i];
    arcs[a]^=1u<<b; arcs[b]^=1u<<a;
  }
  vector<int> c(junctions.rbegin(),junctions.rend());c.insert(c.end(),interior.begin(),interior.end());
  for(int i=0;i<n;++i) {
    assert((arcs[c[i]]>>c[(i+1)%n])&1);
    for (int j=1;j<n;++j) place[c[i]][c[(i+j)%n]]=j-1;
  }
  count_h=0;dfs(0,(1u<<n)-1);assert(count_h==h);
  for(int root:c) {
    int mark=-1;
    for(int i=0;i<(int)junctions.size();++i) if(junctions[i]==root) mark=i;
    for(int q=0;q<n-2;++q) {
      int charge=before[root][q]+after[root][q];
      if (charge) assert((mark>=0?45:225)*charge>=h);
      else {
        for(int a=0;a<n;++a) if(a!=root && place[root][a]<=q)
          for(int b=0;b<n;++b) if(b!=root && place[root][b]>q)
            assert((arcs[a]>>b)&1);
        if(mark<0) assert(links.size()==1 && n==3);
        else {
          if(mark>1)
            for(int a=0;a<=junctions[mark-1];++a) assert(place[root][a]<=q);
          if(mark+2<(int)junctions.size())
            for(int a=junctions[mark+1];a<n;++a) assert(place[root][a]>q);
        }
      }
      ++marked_slots;
    }
    ++profiles;
  }
  ++words;
}

void generate(vector<int>& links,int cost,int threes,int fives){
  if(!links.empty()) test(links);
  if(links.size()==4)return;
  for(int s:{3,4,5}){
    int a=threes+(s==3),b=fives+(s==5);
    if(a+2*b>2 || cost+s-1>13)continue;
    links.push_back(s);generate(links,cost+s-1,a,b);links.pop_back();
  }
}

int main(){
  vector<int> links;generate(links,0,0,0);
  cout<<"chain_words="<<words<<" rooted_profiles="<<profiles<<" marked_slots="<<marked_slots<<"\n";
}
