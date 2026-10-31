// Exact hostile probe: all compositions of25 into15 positive letters.
// Tests actual prefixes followed by one allowed artificial terminal lift.
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;using I=long long;
const int E=9,q=1953125;I leaves=0,hits=0;vector<int>w;
int height(I x){if(!x)return E;int h=0;while(x%5==0&&h<E){x/=5;h++;}return h;}
void dfs(I P,I B,I Q,int A,int S){
 int left=15-w.size();if(left==0)return;
 for(int a=1;a<=25-A-(left-1);a++){
  if(left==1&&A+a!=25)continue;
  I pp=3*P,bb=3*B+Q,qq=Q*(1LL<<a),z=bb-pp+qq;
  int charge=min(height(pp-qq),height(bb));w.push_back(a);
  if(left==1){
   ++leaves;int precision=max(2,height(B-P+Q));
   if(z%q==0&&height(pp-qq)>=precision&&S+E==12){
    ++hits;cout<<"{\"word\":[";for(int i=0;i<(int)w.size();i++)cout<<(i?",":"")<<w[i];
    cout<<"],\"P\":"<<pp<<",\"B\":"<<bb<<",\"Q\":"<<qq<<",\"clock\":25,\"actual_charge\":"<<S+charge<<",\"fake_total_charge\":12,\"previous_precision\":"<<precision<<",\"actual_terminal_content\":"<<charge<<"}\n";
   }
  }else dfs(pp,bb,qq,A+a,S+charge);
  w.pop_back();
 }
}
int main(){dfs(1,0,1,0,0);cerr<<"leaves "<<leaves<<" hits "<<hits<<"\n";}
