#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;using I=long long;
int E=9,q=1953125;I visited=0;vector<int>w;
int height(I x){if(!x)return E;int h=0;while(x%5==0&&h<E){x/=5;h++;}return h;}
void dfs(I P,I B,I Q,int A,int S,bool nonone){
 for(int a=1;A+a<=21;a++){
  I pp=3*P,bb=3*B+Q,qq=Q*(1LL<<a),z=bb-pp+qq;
  int charge=min(height(pp-qq),height(bb));w.push_back(a);visited++;
  int previous_precision=max(1,height(B-P+Q));
  if((nonone||a>1)&&z%q==0&&height(pp-qq)>=previous_precision&&S+E==10&&A+a==21){
   cout<<"{\"word\":[";for(int i=0;i<(int)w.size();i++)cout<<(i?",":"")<<w[i];
   cout<<"],\"P\":"<<pp<<",\"B\":"<<bb<<",\"Q\":"<<qq<<",\"clock\":21,\"actual_charge\":"<<S+charge<<",\"fake_terminal_charge\":10,\"previous_precision\":"<<previous_precision<<",\"actual_terminal_content\":"<<charge<<"}\n";
  }
  dfs(pp,bb,qq,A+a,S+charge,nonone||a>1);w.pop_back();
 }
}
int main(){dfs(1,0,1,0,0,false);cerr<<"visited "<<visited<<"\n";}
