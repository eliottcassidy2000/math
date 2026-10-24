// Independent explicit section-skeleton audit. These are actual representative states,
// not independently minimized pair costs. All intermediate representatives concatenate.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>
using namespace std;
int plusq(int a,int b){return (((a/16+b/16)%128)*16)|((a^b)&15);}
int minusq(int a){return ((128-a/16)%128)*16|(a&15);}
struct Core{uint32_t mask=0;int q=0,w=999;};
struct Helper{int a,b,q;};
int main(int argc,char**argv){
 if(argc!=2)throw runtime_error("input filename required");ifstream input(argv[1]);int nc,nh;input>>nc>>nh;
 if(nc!=28||nh!=51)throw runtime_error("declared inventory");
 vector<int> id(nc+nh),lab(nc+nh),cost(nc+nh);for(int i=0;i<nc+nh;++i)input>>id[i]>>lab[i]>>cost[i];if(!input)throw runtime_error("complete input");
 vector<Helper> helpers;int helperIndex[51][51];for(int a=0;a<51;++a)for(int b=0;b<51;++b)helperIndex[a][b]=-1;
 for(int a=0;a<nh;++a)for(int b=a+1;b<nh;++b){helperIndex[a][b]=helperIndex[b][a]=helpers.size();helpers.push_back({a,b,plusq(lab[nc+a],lab[nc+b])});}
 cout<<"[";bool first=true;
 for(int grade=2;grade<=4;++grade){
  vector<Core> section(2048);int states=0;
  function<void(int,int,uint32_t,int,int)> enumerate=[&](int next,int left,uint32_t mask,int q,int w){
   if(!left){++states;auto &r=section[q];if(w<r.w||(w==r.w&&mask<r.mask))r={mask,q,w};return;}
   for(int i=next;i<=nc-left;++i)enumerate(i+1,left-1,mask|(1u<<i),plusq(q,lab[i]),w+cost[i]);};
  enumerate(0,grade,0,0,0);int charges=0;for(auto c:section)charges+=c.w!=999;
  vector<uint8_t> support(2048*2048,255),work(2048*2048,255);
  for(int a=0;a<2048;++a)if(section[a].w!=999)for(int b=0;b<2048;++b)if(section[b].w!=999){
   uint32_t m=section[a].mask^section[b].mask;int n=__builtin_popcount(m),w=0;for(int i=0;i<nc;++i)if(m>>i&1)w+=cost[i];support[a*2048+b]=n;work[a*2048+b]=w;}
  for(int sign:{1,-1}){
   int minSheets=1276,maxSheets=0,maxUsedClock=0,maxUsedSupport=0;uint64_t checksum=1469598103934665603ULL;long long totalSheets=0,acceptedEdges=0;int maxAcceptedClock=0;
   for(int total=0;total<2048;++total){
    vector<int> residual(helpers.size(),-1);int root=-1,count=0;
    for(int j=0;j<(int)helpers.size();++j){int q=plusq(total,minusq(helpers[j].q));if(sign<0)q=minusq(q);if(section[q].w!=999){residual[j]=q;if(root<0)root=j;++count;}}
    if(root<0)throw runtime_error("nonempty joint charge missing");minSheets=min(minSheets,count);maxSheets=max(maxSheets,count);totalSheets+=count;
    vector<int> queue={root};vector<bool> visited(helpers.size());visited[root]=true;
    for(size_t head=0;head<queue.size();++head){int j=queue[head];auto h=helpers[j];
     for(int retain:{h.a,h.b}){int drop=retain==h.a?h.b:h.a;
      for(int insert=0;insert<nh;++insert){if(insert==retain||insert==drop)continue;int k=helperIndex[retain][insert];if(residual[k]<0)continue;
       int ix=residual[j]*2048+residual[k],n=2+support[ix],w=cost[nc+drop]+cost[nc+insert]+work[ix];
       if(n<=8&&w<=33){if(j<k){++acceptedEdges;maxAcceptedClock=max(maxAcceptedClock,w);}if(visited[k])continue;visited[k]=true;queue.push_back(k);maxUsedClock=max(maxUsedClock,w);maxUsedSupport=max(maxUsedSupport,n);
        for(uint64_t z:{(uint64_t)total,(uint64_t)j,(uint64_t)k,(uint64_t)w,(uint64_t)n}){checksum^=z;checksum*=1099511628211ULL;}}
      }
     }
    }
    if((int)queue.size()!=count){cerr<<"grade="<<grade<<" sign="<<sign<<" q="<<total<<" reached="<<queue.size()<<" sheets="<<count<<"\n";throw runtime_error("section skeleton disconnected");}
   }
   if(!first)cout<<",";first=false;cout<<"{\"virtual_grade\":"<<grade<<",\"sign\":"<<sign<<",\"enumerated_core_subsets\":"<<states<<",\"core_charge_count\":"<<charges<<",\"total_charges_verified\":2048,\"total_helper_sheets\":"<<totalSheets<<",\"min_helper_sheets\":"<<minSheets<<",\"max_helper_sheets\":"<<maxSheets<<",\"max_selected_edge_clock\":"<<maxUsedClock<<",\"max_selected_edge_support\":"<<maxUsedSupport<<",\"accepted_edges\":"<<acceptedEdges<<",\"max_accepted_edge_clock\":"<<maxAcceptedClock<<",\"forest_fnv64\":\""<<checksum<<"\",\"section_choice\":\"minimum selected weight then minimum binary mask\"}";
  }
 }
 cout<<"]\n";
}
