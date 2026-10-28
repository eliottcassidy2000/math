// Alternate full-bank verifier: stable colors followed by ALL permutations
// inside color cells, without individualize-and-refine branching.
// Shared low-level row encoding is not an independent component.
#define main generation_entry_unused
#include "generate.cpp"
#undef main
#include <functional>
#include <numeric>
struct AuditCanon{uint64_t code=~uint64_t(0),aut=0;};
AuditCanon alternate(const Rows&r,int n){
 array<int,12>color{},next{},order{};for(int i=0;i<n;++i){color[i]=pc(r[i]);order[i]=i;}
 for(;;){
  array<array<int,13>,12>sig{};array<int,12>size{};
  for(int i=0;i<n;++i)++size[color[i]];
  int before=0;for(auto x:size)before+=(x!=0);
  for(int i=0;i<n;++i){
   sig[i][0]=color[i];
   for(int j=0;j<n;++j)if(r[i]>>j&1)++sig[i][color[j]+1];
  }
  sort(order.begin(),order.begin()+n,[&](int a,int b){return sig[a]<sig[b];});
  int current=0;next[order[0]]=0;
  for(int k=1;k<n;++k){current+=(sig[order[k]]!=sig[order[k-1]]);next[order[k]]=current;}
  color=next;if(current+1==before)break;
 }
 vector<vector<int>>cells(*max_element(color.begin(),color.begin()+n)+1);
 for(int i=0;i<n;++i)cells[color[i]].push_back(i);
 Perm p{};AuditCanon best;
 function<void(int,int)>walk=[&](int cell,int offset){
  if(cell==int(cells.size())){
   uint64_t code=uint64_t(encode(r,n,p));
   if(code<best.code){best.code=code;best.aut=0;}
   if(code==best.code)++best.aut;return;
  }
  auto block=cells[cell];
  do{
   for(int j=0;j<int(block.size());++j)p.at(offset+j)=block[j];
   walk(cell+1,offset+int(block.size()));
  }while(next_permutation(block.begin(),block.end()));
 };
 walk(0,0);return best;
}
vector<uint64_t>read64(const filesystem::path&path){
 ifstream f(path,ios::binary);need(bool(f),"binary dictionary open");
 vector<uint64_t>v;for(;;){uint64_t x=0;int c=f.get();if(c==EOF)break;x=uint8_t(c);
  for(int j=1;j<8;++j){c=f.get();need(c!=EOF,"complete dictionary word");x|=uint64_t(uint8_t(c))<<(8*j);}v.push_back(x);}
 return v;
}
int main(int argc,char**argv){try{
 need(argc==4,"usage: verify_bank n directory receipt");
 int n=stoi(argv[1]);need(n>=2&&n<=10,"audit order2..10");filesystem::path dir=argv[2];
 auto started=chrono::steady_clock::now();auto parents=read64(dir/("parents"+to_string(n-1)+".u64"));
 ifstream flags(dir/("flags"+to_string(n)+".varint"),ios::binary);need(bool(flags),"flags opened");
 vector<uint64_t>codes;map<uint64_t,uint64_t>hist;uint64_t address=0,previous=0,mass=0,fac=1;
 for(int i=2;i<=n;++i)fac*=i;
 for(;;){
  int c=flags.get();if(c==EOF)break;uint64_t delta=0;int shift=0;
  for(;;){need(shift<63,"varint length");delta|=uint64_t(c&127)<<shift;if(!(c&128))break;
   shift+=7;c=flags.get();need(c!=EOF,"complete varint");}
  address+=delta;need(codes.empty()||address>previous,"strict flag order");previous=address;
  uint64_t pid=address>>(n-1);U mask=address&((U(1)<<(n-1))-1);
  need(pid<parents.size(),"parent id bound");Rows r=extend(decode(parents[pid],n-1),n-1,mask);
  auto a=alternate(r,n);need(a.aut&&fac%a.aut==0,"alternate group order");codes.push_back(a.code);
  ++hist[a.aut];mass+=fac/a.aut;
  if(codes.size()%1000000==0)cerr<<"audit_classes="<<codes.size()<<"\n";
 }
 sort(codes.begin(),codes.end());need(adjacent_find(codes.begin(),codes.end())==codes.end(),"no isomorphic duplicate");
 need(mass==(uint64_t(1)<<(n*(n-1)/2)),"complete labelled orbit mass");
 ofstream f(argv[3]);f<<"{\"status\":\"PASS\",\"order\":"<<n<<",\"classes\":"<<codes.size()
  <<",\"labelled_mass\":"<<mass<<",\"full_class_uniqueness\":true,\"automorphism_order_histogram\":{";
 bool comma=false;for(auto[g,k]:hist){if(comma)f<<",";comma=true;f<<"\""<<g<<"\":"<<k;}
 f<<"},\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";
 cout<<"PASS full distinct class bank; count="<<codes.size()<<" mass="<<mass<<"\n";
 return 0;
}catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
