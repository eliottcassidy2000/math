#define main exhaustive_source_main
#include "power_sources.cpp"
#undef main
int main(){
 vector<vector<int>> seeds={{284,25,154,496,128,23,311,99,182},{420,389,72,451,463,478,259,68,132},{228,193,394,451,143,158,52,64,243}};
 int tested=0,si_found=-1,mask_found=-1;
 for(int si=0;si<seeds.size()&&!hits.count(237);si++)for(int mask=0;mask<512&&!hits.count(237);mask++){
  auto rows=seeds[si];rows.push_back(mask);for(int v=0;v<9;v++)if(!(mask>>v&1))rows[v]|=1<<9;
  test(rows);tested++;if(hits.count(237)){si_found=si;mask_found=mask;}
 }
 cout<<"{\"target\":237,\"tested_insertions\":"<<tested<<",\"seed_index\":"<<si_found<<",\"insertion_mask\":"<<mask_found<<",\"hits\":[";
 if(hits.count(237)){auto h=hits[237];cout<<"{\"deficit\":237,\"mask\":"<<h.mask<<",\"state\":["<<h.F<<","<<h.q<<","<<h.g<<","<<h.h<<"],\"rows\":[";for(int i=0;i<h.rows.size();i++){if(i)cout<<",";cout<<h.rows[i];}cout<<"]}";}
 cout<<"]}\n";return 0;
}
