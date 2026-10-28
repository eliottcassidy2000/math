// Consecutive aligned tournament pairs and normalized anti-involution doubles.
#define main generation_entry_unused
#include "generate.cpp"
#undef main
#include <numeric>
Rows double_pair(const Rows&T,const Rows&W,U mask,int n,bool odd){
 Rows D{};
 for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
  bool a=T[i]>>j&1,b=W[i]>>j&1,agree=a==b;
  D[a?i:j]|=U(1)<<(a?j:i);
  D[a?n+j:n+i]|=U(1)<<(a?n+i:n+j);
  D[agree?i:n+j]|=U(1)<<(agree?n+j:i);
  D[agree?j:n+i]|=U(1)<<(agree?n+i:j);
 }
 for(int i=0;i<n;++i){
  bool d=!odd||((mask>>i)&1);
  D[d?i:n+i]|=U(1)<<(d?n+i:i);
  if(odd){D[2*n]|=U(1)<<i;D[n+i]|=U(1)<<(2*n);}
 }
 return D;
}
int main(int argc,char**argv){try{
 need(argc==2,"usage: pair_bridge receipt.json");ofstream out(argv[1]);out<<"{\"status\":\"FINITE-EXACT\",\"rows\":[";
 bool comma=false;
 for(int n=1;n<=4;++n)for(int odd=0;odd<=1;++odd){
  int e=n*(n-1)/2,N=2*n+odd;map<uint64_t,uint64_t>forward,backward;
  uint64_t examined=0;vector<uint64_t>aut_gap;Perm p{};for(int i=0;i<n;++i)p[i]=i;
  uint64_t universe=uint64_t(1)<<(2*e+(odd?n:0));
  for(uint64_t x=0;x<universe;++x){
   Rows T=decode(x&((uint64_t(1)<<e)-1),n),W=decode((x>>e)&((uint64_t(1)<<e)-1),n);
   U mask=odd?(x>>(2*e)):0;Rows D=double_pair(T,W,mask,n,odd);
   Perm involution{};for(int i=0;i<2*n;++i)involution[i]=(i+n)%(2*n);if(odd)involution[2*n]=2*n;
   for(int i=0;i<N;++i)for(int j=0;j<N;++j)if(i!=j)
    need(((D[i]>>j)&1)==((D[involution[j]]>>involution[i])&1),"literal anti-involution");
   uint64_t source=~uint64_t(0),source_aut=0;
   for(int i=0;i<n;++i)p[i]=i;
   do{
    uint64_t y=uint64_t(encode(T,n,p))|(uint64_t(encode(W,n,p))<<e)|(uint64_t(move_mask(mask,p,n))<<(2*e));
    if(y<source){source=y;source_aut=0;}if(y==source)++source_aut;
   }while(next_permutation(p.begin(),p.begin()+n));
   Canon can=canonical(D,N);uint64_t target=uint64_t(can.code);
   auto a=forward.emplace(source,target);need(a.second||a.first->second==target,"well-defined orbit map");
   auto b=backward.emplace(target,source);need(b.second||b.first->second==source,"injective orbit map");
   // Compare full target automorphisms with the centralizer of the chosen involution.
   Perm invbase{};for(int i=0;i<N;++i)invbase[can.orders[0][i]]=i;
   uint64_t centralizer=0;
   for(auto order:can.orders){
    Perm g{};for(int i=0;i<N;++i)g[i]=order[invbase[i]];
    bool commutes=true;
    for(int i=0;i<N;++i)commutes&=g[involution[i]]==involution[g[i]];
    centralizer+=commutes;
   }
   need(centralizer==source_aut,"source automorphisms equal target centralizer");
   need(can.orders.size()%centralizer==0,"anti-involution conjugacy index");
   if(aut_gap.empty()&&centralizer<can.orders.size())aut_gap={x,source_aut,can.orders.size(),target};
   ++examined;
  }
  if(comma)out<<",";comma=true;
  out<<"{\"n\":"<<n<<",\"odd\":"<<(odd?"true":"false")<<",\"target_order\":"<<N
     <<",\"labelled_pair_inputs\":"<<examined<<",\"pair_orbits\":"<<forward.size()
     <<",\"distinct_self_converse_classes\":"<<backward.size()<<",\"group_loss_witness\":[";
  for(size_t i=0;i<aut_gap.size();++i){if(i)out<<",";out<<aut_gap[i];}out<<"]}";
  cout<<"n="<<n<<" odd="<<odd<<" inputs="<<examined<<" classes="<<forward.size()<<"\n";
 }
 out<<"],\"checks\":"<<checks<<"}\n";return 0;
}catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
