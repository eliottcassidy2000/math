#include <iostream>
#include <vector>
#include <set>
using namespace std;
int n, rows[8],full;
int count_dfs(int used,int last){
 if(used==full)return n>=3 && (rows[last]&1);
 int todo=rows[last]&~used,total=0;
 while(todo){int bit=todo&-todo;todo-=bit;total+=count_dfs(used|bit,__builtin_ctz((unsigned)bit));}
 return total;
}
int main(){
 for(n=1;n<=7;n++){
  full=(1<<n)-1;vector<pair<int,int>> arcs;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)arcs.emplace_back(i,j);
  unsigned total=1u<<arcs.size();set<int> values;
  for(unsigned code=0;code<total;code++){
   fill(rows,rows+n,0);
   for(unsigned e=0;e<arcs.size();e++){auto [i,j]=arcs[e];if(code>>e&1)rows[i]|=1<<j;else rows[j]|=1<<i;}
   values.insert(count_dfs(1,0));
  }
  cout<<n<<'\t'<<total<<'\t';bool first=true;for(int c:values){cout<<(first?"":",")<<c;first=false;}cout<<'\n';
 }
}
