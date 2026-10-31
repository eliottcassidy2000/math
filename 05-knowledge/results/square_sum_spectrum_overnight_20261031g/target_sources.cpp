#define main exhaustive_source_main
#include "power_sources.cpp"
#undef main
#include <set>
int main(int argc,char**argv){
 int cap=argc>1?stoi(argv[1]):50000;set<I>wanted={117,169,223,237};int used=0;
 for(int s=0;s<cap;s++){
  int n=9;vector<int>rows(n);
  for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)if(next()&1)rows[u]|=1<<v;else rows[v]|=1<<u;
  test(rows);used++;
  bool done=true;for(I c:wanted)if(!hits.count(c))done=false;
  if(done)break;
 }
 cout<<"{\"seed\":123456789,\"parent_order\":9,\"sampled_parents\":"<<used<<",\"sample_cap\":"<<cap<<",\"hits\":[";bool first=true;
 for(I c:wanted)if(hits.count(c)){auto h=hits[c];if(!first)cout<<",";first=false;cout<<"{\"deficit\":"<<c<<",\"mask\":"<<h.mask<<",\"state\":["<<h.F<<","<<h.q<<","<<h.g<<","<<h.h<<"],\"rows\":[";
  for(int i=0;i<h.rows.size();i++){if(i)cout<<",";cout<<h.rows[i];}cout<<"]}";
 }
 cout<<"],\"scope\":\"Deterministic source discovery; no absence claim outside the sample.\"}\n";
}
