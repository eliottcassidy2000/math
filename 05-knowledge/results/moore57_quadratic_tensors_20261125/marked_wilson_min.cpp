// Exact minimum of the one-source route-collision action through a small
// number of reference edits. Reads "nested_core_check.py --board" output.
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
using namespace std;
struct Contact {int y;array<int,3> p;};
vector<Contact>B[53];
Contact chosen[53];
int old[53],use[3][56],best,score;
array<int,53> ans;
vector<int> holes;
void push(Contact a){for(int g=0;g<3;g++)score+=use[g][a.p[g]]++;}
void pop(Contact a){for(int g=0;g<3;g++)score-=--use[g][a.p[g]];}
int price(Contact a){return use[0][a.p[0]]+use[1][a.p[1]]+use[2][a.p[2]];}
void finish(int t){
 if(score>=best)return;
 if(t==2){best=score;for(int h=0;h<53;h++)ans[h]=chosen[h].y;return;}
 int h=holes[t];
 for(auto a:B[h])if(score+price(a)<best){chosen[h]=a;push(a);finish(t+1);pop(a);}
}
void change(vector<int>& H,int t){
 if(score>=best)return;
 if(t==(int)H.size()){finish(0);return;}
 int h=H[t];
 for(auto a:B[h])if(a.y!=old[h]&&score+price(a)<best){chosen[h]=a;push(a);change(H,t+1);pop(a);}
}
void pick(vector<int>&H,int t,int left){
 if(!left){change(H,0);return;}
 for(int h=t;h<53;h++)if(old[h]>=0){
   H.push_back(h);
   auto a=chosen[h];
   pop(a);pick(H,h+1,left-1);push(a);
   chosen[h]=a;H.pop_back();
 }
}
int main(int argc,char**argv){
 int n;cin>>n;if(n!=53)return 2;
 Contact marker;
 for(int h=0;h<n;h++){
  int k;cin>>old[h]>>k;
  while(k--){
   Contact a;cin>>a.y>>a.p[0]>>a.p[1]>>a.p[2];
   if(!h){marker=a;chosen[h]=a;continue;}
   if(equal(a.p.begin(),a.p.end(),marker.p.begin(),[](int x,int y){return x!=y;})){
    B[h].push_back(a);
    if(a.y==old[h]){chosen[h]=a;push(a);}
   }
  }
  if(old[h]<0)holes.push_back(h);
 }
 if(holes.size()!=2||score!=0||marker.y!=7)return 3;
 int lim=argc>1?stoi(argv[1]):2;
 best=100000;
 for(int k=0;k<=lim;k++){
   vector<int> H;pick(H,1,k);
   cout<<k<<" "<<best;for(int y:ans)cout<<" "<<y;
   cout<<endl;
 }
}
