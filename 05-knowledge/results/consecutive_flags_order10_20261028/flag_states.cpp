// Exact fixed-path states with one distinguished vertex, through order seven.
#define main generation_entry_unused
#include "generate.cpp"
#undef main
void write32(ofstream&f,uint32_t x){for(int i=0;i<4;++i)f.put(char(x>>(8*i)));}
int main(int argc,char**argv){try{
 need(argc==2,"usage: flag_states output_directory");filesystem::path dir=argv[1];filesystem::create_directories(dir);
 production_bound=-1;
 for(int n=2;n<=7;++n){
  int m=(n-1)*(n-2)/2;uint64_t total=uint64_t(1)<<m;
  map<pair<uint64_t,U>,uint32_t>ids;vector<array<uint64_t,3>>meta;
  ofstream f(dir/("flag_states"+to_string(n)+".u32"),ios::binary);
  for(uint64_t x=0;x<total;++x){
   Rows r{};int k=0;
   for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
    bool forward=j==i+1||((x>>k)&1);
    if(j!=i+1)++k;
    r[forward?i:j]|=U(1)<<(forward?j:i);
   }
   auto child=canonical(r,n).code;
   for(int v=0;v<n;++v){
    auto q=pointed(r,n,v);auto it=ids.find(q);uint32_t id;
    if(it==ids.end()){
     id=uint32_t(ids.size());ids.emplace(q,id);meta.push_back({q.first,q.second,uint64_t(child)});
    }else{id=it->second;need(meta[id][2]==uint64_t(child),"pointed key determines child");}
    write32(f,id);
   }
  }
  f.close();ofstream out(dir/("flag_meta"+to_string(n)+".json"));out<<"[";
  for(size_t i=0;i<meta.size();++i){if(i)out<<",";out<<"["<<meta[i][0]<<","<<meta[i][1]<<","<<meta[i][2]<<"]";}out<<"]\n";
  cout<<"order="<<n<<" path_root_states="<<total*n<<" pointed_types="<<ids.size()<<"\n";
 }
 return 0;
}catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
