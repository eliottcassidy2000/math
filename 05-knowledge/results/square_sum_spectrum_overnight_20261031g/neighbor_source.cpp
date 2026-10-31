#define main exhaustive_source_main
#include "power_sources.cpp"
#undef main
int main(){
 vector<vector<int>> seeds={{284,25,154,496,128,23,311,99,182},{420,389,72,451,463,478,259,68,132},{228,193,394,451,143,158,52,64,243}};
 int tested=0,seed_index=-1,first_edge=-1,second_edge=-1;
 for(int si=0;si<seeds.size()&&!hits.count(237);si++){
  auto seed=seeds[si];vector<pair<int,int>> es;for(int u=0;u<9;u++)for(int v=u+1;v<9;v++)es.push_back({u,v});
  for(int a=0;a<es.size()&&!hits.count(237);a++)for(int b=a;b<es.size()&&!hits.count(237);b++){
   auto rows=seed;auto flip=[&](int i){auto [u,v]=es[i];rows[u]^=1<<v;rows[v]^=1<<u;};flip(a);if(b!=a)flip(b);test(rows);tested++;
   if(hits.count(237)){seed_index=si;first_edge=a;second_edge=b;}
  }
 }
 cout<<"{\"target\":237,\"tested_neighbors\":"<<tested<<",\"seed_index\":"<<seed_index<<",\"flipped_edge_indices\":["<<first_edge<<","<<second_edge<<"],\"hits\":[";
 if(hits.count(237)){auto h=hits[237];cout<<"{\"deficit\":237,\"mask\":"<<h.mask<<",\"state\":["<<h.F<<","<<h.q<<","<<h.g<<","<<h.h<<"],\"rows\":[";for(int i=0;i<h.rows.size();i++){if(i)cout<<",";cout<<h.rows[i];}cout<<"]}";}
 cout<<"]}\n";return 0;
}
