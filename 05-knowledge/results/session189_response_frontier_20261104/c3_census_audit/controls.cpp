// Auxiliary audit controls: does not include producer code.
#define main audit_program_main
#include "audit.cpp"
#undef main
#include <random>
unsigned subset_module(int n,const R&r){
 unsigned F=(1u<<n)-1;
 for(unsigned m=1;m<F;++m)if(pc(m)>=2){
  bool module=true;for(int x=0;x<n;++x)if(!(m>>x&1u)){
   bool first=true,orientation=false;
   for(int y=0;y<n;++y)if(m>>y&1u){bool edge=(r[x]>>y)&1u;if(first){orientation=edge;first=false;}else if(edge!=orientation)module=false;}
  }
  if(module)return m;
 }
 return 0;
}
int main(int argc,char**argv){try{
 bool compact=(argc==3&&std::string(argv[1])=="--compact");
 require(argc==2||compact,"usage: controls DATA_DIR OR controls --compact CENSUS_TSV");
 U all_labeled=0;
 for(int n=2;n<=6;++n)for(U code=0;code<(U(1)<<(n*(n-1)/2));++code){
  R r=decode(n,code);require(bool(proper_module(n,r))==bool(subset_module(n,r)),"closure/subset mismatch exhaustive order <=6");++all_labeled;
 }
 std::ifstream in(compact?argv[2]:std::string(argv[1])+"/tower_9.txt");require(bool(in),"input missing");std::string line;std::getline(in,line);
 U rows_checked=0,relabel_checks=0;std::mt19937 rng(189);
 while(std::getline(in,line)){
  U code=0;if(compact){std::istringstream row(line);unsigned node;row>>node>>code;require(bool(row)&&node==rows_checked,"invalid compact row");}else code=std::stoull(line.substr(line.rfind('|')+1));R r=decode(9,code);
  require(bool(proper_module(9,r))==bool(subset_module(9,r)),"closure/subset mismatch at order 9");
  if(rows_checked%97==0){
   Canon c=canonical(9,r);std::array<int,9> p;std::iota(p.begin(),p.end(),0);std::shuffle(p.begin(),p.end(),rng);R t{};
   for(int i=0;i<9;++i)for(int j=0;j<9;++j)if(r[p[i]]>>p[j]&1u)t[i]|=1u<<j;
   Canon d=canonical(9,t);require(c.key==d.key&&c.aut==d.aut,"relabel invariance failure");++relabel_checks;
  }
  ++rows_checked;
 }
 std::cout<<"PASS: module closure equals full subset predicate on "<<all_labeled<<" labeled tournaments of orders 2..6 and all "<<rows_checked<<" order-9 merged representatives; "<<relabel_checks<<" deterministic shuffled-label canonicalization controls.\n";
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
