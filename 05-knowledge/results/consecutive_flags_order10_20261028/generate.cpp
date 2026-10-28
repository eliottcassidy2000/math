// Parent-only canonical augmentation with stable refinement.
// Inherits the exact mechanism of root_ports_closed_carriers_20261028.
// Generation contains no target count and never canonicalizes a target child.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using U=uint16_t;
using Rows=array<U,12>;
using Perm=array<uint8_t,12>;
using Code=__uint128_t;
using Cells=vector<U>;
uint64_t checks=0,leaves=0,canons=0;
int production_bound=-1;
void need(bool b,const char*m){++checks;if(!b)throw runtime_error(m);}
int pc(U x){return __builtin_popcount(unsigned(x));}
int first(U x){return __builtin_ctz(unsigned(x));}
Cells initial(const Rows&r,int n){
 array<U,12>a{};for(int v=0;v<n;++v)a[pc(r[v])]|=U(1)<<v;
 Cells c;for(auto x:a)if(x)c.push_back(x);return c;
}
Cells stable(const Rows&r,Cells cells){
 for(;;){
  Cells out;
  for(U cell:cells){
   if(pc(cell)==1){out.push_back(cell);continue;}
   array<pair<uint64_t,int>,12>sig{};int count=0;
   for(U rem=cell;rem;rem&=rem-1){
    int v=first(rem);uint64_t s=0;
    for(U mask:cells)s=(s<<4)|pc(r[v]&mask);
    sig[count++]={s,v};
   }
   sort(sig.begin(),sig.begin()+count);
   uint64_t previous=~uint64_t(0);U part=0;
   for(int k=0;k<count;++k){
    if(k&&sig[k].first!=previous){out.push_back(part);part=0;}
    previous=sig[k].first;part|=U(1)<<sig[k].second;
   }
   out.push_back(part);
  }
  if(out.size()==cells.size())return out;
  cells=move(out);
 }
}
Code encode(const Rows&r,int n,const Perm&p){
 Code code=0;int bit=0;
 for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit)
  code|=Code((r[p[i]]>>p[j])&1)<<bit;
 return code;
}
uint64_t encode(const Rows&r,int n){
 Perm p{};for(int i=0;i<n;++i)p[i]=i;
 need(n<=10,"raw export width");return uint64_t(encode(r,n,p));
}
Rows decode(uint64_t code,int n){
 Rows r{};int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit)
  r[(code>>bit)&1?i:j]|=U(1)<<((code>>bit)&1?j:i);
 return r;
}
struct Canon{Code code=~Code(0);vector<Perm>orders;};
void visit(const Rows&r,int n,Cells cells,Canon&answer){
 cells=stable(r,move(cells));int split=-1;
 for(int i=0;i<int(cells.size());++i)if(pc(cells[i])>1){split=i;break;}
 if(split<0){
  ++leaves;Perm p{};for(int i=0;i<n;++i)p[i]=first(cells[i]);
  Code code=encode(r,n,p);
  if(code<answer.code){answer.code=code;answer.orders.clear();}
  if(code==answer.code)answer.orders.push_back(p);
  return;
 }
 U part=cells[split];
 for(U rem=part;rem;rem&=rem-1){
  U bit=rem&-rem;Cells next=cells;
  next[split]=bit;next.insert(next.begin()+split+1,part^bit);
  visit(r,n,move(next),answer);
 }
}
Canon canonical(const Rows&r,int n){
 need(production_bound<0||n<=production_bound,"parent-only canonicalization");
 ++canons;Canon a;visit(r,n,initial(r,n),a);need(!a.orders.empty(),"nonempty transports");return a;
}
U move_mask(U mask,const Perm&p,int n){
 U out=0;for(int i=0;i<n;++i)out|=U((mask>>p[i])&1)<<i;return out;
}
Rows extend(const Rows&r,int n,U mask){
 Rows out=r;out[n]=mask;
 for(int v=0;v<n;++v)out[v]|=U(1-((mask>>v)&1))<<n;
 return out;
}
pair<Rows,U> deletion(const Rows&r,int N,int root){
 Rows out{};U mask=0;Perm p{};int k=0;
 for(int i=0;i<N;++i)if(i!=root)p[k++]=i;
 for(int i=0;i<N-1;++i)out[i]=move_mask(r[p[i]],p,N-1);
 mask=move_mask(r[root],p,N-1);return {out,mask};
}
pair<uint64_t,U> pointed(const Rows&r,int N,int root){
 auto d=deletion(r,N,root);Canon c=canonical(d.first,N-1);U mask=(U(1)<<(N-1))-1;
 for(auto&p:c.orders)mask=min(mask,move_mask(d.second,p,N-1));
 return {uint64_t(c.code),mask};
}
bool strong(const Rows&r,int n){
 for(int reversed=0;reversed<2;++reversed){
  U done=0,todo=1;
  while(todo){
   int v=first(todo);todo&=todo-1;done|=U(1)<<v;U next=0;
   if(!reversed)next=r[v];
   else for(int w=0;w<n;++w)if(r[w]>>v&1)next|=U(1)<<w;
   todo|=next&~done;
  }
  if(done!=(U(1)<<n)-1)return false;
 }
 return true;
}
void write64(ofstream&f,uint64_t x){for(int i=0;i<8;++i)f.put(char(x>>(8*i)));}
void varint(ofstream&f,uint64_t x){
 while(x>=128){f.put(char((x&127)|128));x>>=7;}f.put(char(x));
}
struct Parent{uint64_t code;Rows rows;vector<Perm>group;};
void controls(){
 for(int n:{3,5,7,12}){
  Rows r{};int expected=0;
  if(n<12){
   vector<int>steps=n==3?vector<int>{1}:n==5?vector<int>{1,2}:vector<int>{1,2,4};
   for(int v=0;v<n;++v)for(int d:steps)r[v]|=U(1)<<((v+d)%n);
   expected=n==3?3:n==5?5:21;
  }else{
   for(int block=0;block<3;++block)for(int v=0;v<4;++v)
    r[4*block+v]=((v<3?U((1<<((v+1)%3))|8):0)<<(4*block))|(U(15)<<(4*((block+1)%3)));
   expected=81;
  }
  need(int(canonical(r,n).orders.size())==expected,"strong atom group control");
  if(n==12){
   for(int i=0;i<4;++i)for(int j=0;j<4;++j)if(i!=j)r[i]^=U(1)<<j;
   need(canonical(r,n).orders.size()==27,"81/27 native group control");
  }
 }
}
int main(int argc,char**argv){try{
 need(argc==3,"usage: generate max_order output_directory");
 int maximum=stoi(argv[1]);need(maximum>=1&&maximum<=10,"explicit scope1..10");
 filesystem::path dir=argv[2];filesystem::create_directories(dir);controls();
 vector<uint64_t>parents{0};ofstream report(dir/"generation.json");report<<"{\n\"results\":[\n";
 uint64_t fact=1;
 for(int n=0;n<maximum;++n){
  auto start=chrono::steady_clock::now();int N=n+1;fact*=N;production_bound=n;
  uint64_t leaf_start=leaves,canon_start=canons;
  vector<Parent>normalized;normalized.reserve(parents.size());
  for(auto raw:parents){
   Rows r=decode(raw,n);Canon c=canonical(r,n);Rows nr=decode(uint64_t(c.code),n);
   Perm inverse{};for(int i=0;i<n;++i)inverse[c.orders[0][i]]=i;
   vector<Perm>group;
   for(auto p:c.orders){Perm g{};for(int i=0;i<n;++i)g[i]=inverse[p[i]];
    need(encode(nr,n,g)==c.code,"actual parent action");group.push_back(g);}
   normalized.push_back({uint64_t(c.code),nr,move(group)});
  }
  sort(normalized.begin(),normalized.end(),[](auto&a,auto&b){return a.code<b.code;});
  for(size_t i=1;i<normalized.size();++i)need(normalized[i-1].code!=normalized[i].code,"distinct canonical parents");
  ofstream dict(dir/("parents"+to_string(n)+".u64"),ios::binary);
  for(auto&p:normalized)write64(dict,p.code);dict.close();need(bool(dict),"parent dictionary output");
  uint64_t prep_leaves=leaves-leaf_start;
  uint64_t states=0,degree_reject=0,stable_reject=0,key_reject=0,queries=0,accepted=0,mass=0,ns=0,early_single=0;
  map<uint64_t,uint64_t>aut_hist;map<int,uint64_t>selected_sizes;
  vector<uint64_t>children;
  ofstream output(dir/("flags"+to_string(N)+".varint"),ios::binary);
  uint64_t previous=0;
  for(size_t id=0;id<normalized.size();++id){
   auto&p=normalized[id];array<bool,1024>seen{};
   for(U mask=0;mask<(U(1)<<n);++mask){
    if(seen[mask])continue;
    for(auto&g:p.group){U m=move_mask(mask,g,n);need(m>=mask,"minimum mask orbit");seen[m]=true;}
    ++states;int d=pc(mask);bool reject=false;
    for(int v=0;v<n;++v)if(pc(p.rows[v])+1-((mask>>v)&1)>d){reject=true;break;}
    if(reject){++degree_reject;continue;}
    Rows child=extend(p.rows,n,mask);Cells cells=initial(child,N);U selected=cells.back();
    if(pc(selected)==1)++early_single;
    else selected=stable(child,move(cells)).back();
    if(!(selected>>n&1)){++stable_reject;continue;}
    ++selected_sizes[pc(selected)];uint64_t orbit=1;
    auto own=make_pair(p.code,mask);
    for(U rem=selected&~(U(1)<<n);rem;rem&=rem-1){
     ++queries;auto key=pointed(child,N,first(rem));
     if(key<own){++key_reject;reject=true;break;}
     orbit+=(key==own);
    }
    if(reject)continue;
    uint64_t stabilizer=0;for(auto&g:p.group)stabilizer+=(move_mask(mask,g,n)==mask);
    uint64_t aut=orbit*stabilizer;need(aut&&fact%aut==0&&aut%2==1,"child orbit-stabilizer");
    ++aut_hist[aut];mass+=fact/aut;++accepted;
    if(N>=2)ns+=strong(child,N);else ++ns;
    uint64_t address=(uint64_t(id)<<n)|mask;
    need(accepted==1||address>previous,"ordered flag address");
    varint(output,address-previous);previous=address;
    if(N<maximum||N<=9)children.push_back(encode(child,N));
   }
   if(N==10&&(id+1)%20000==0)
    cerr<<"parents="<<id+1<<"/"<<normalized.size()<<" accepted="<<accepted<<" queries="<<queries<<"\n";
  }
  output.close();need(bool(output),"complete flag output");
  need(accepted+degree_reject+stable_reject+key_reject==states,"candidate partition");
  need(mass==(uint64_t(1)<<(N*(N-1)/2)),"exact labelled mass");
  if(N<=9){ofstream raw(dir/("raw"+to_string(N)+".u64"),ios::binary);for(auto c:children)write64(raw,c);raw.close();need(bool(raw),"raw output");}
  double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();
  if(n)report<<",\n";
  report<<"{\"order\":"<<N<<",\"classes\":"<<accepted<<",\"pointed_states\":"<<states
   <<",\"strong_classes\":"<<ns<<",\"degree_rejections\":"<<degree_reject<<",\"stable_rejections\":"<<stable_reject
   <<",\"key_rejections\":"<<key_reject<<",\"competing_root_queries\":"<<queries
   <<",\"early_degree_singletons\":"<<early_single<<",\"preparation_leaves\":"<<prep_leaves
   <<",\"competitor_leaves\":"<<leaves-leaf_start-prep_leaves<<",\"canonical_searches\":"<<canons-canon_start
   <<",\"maximum_canonicalized_order\":"<<n<<",\"labelled_mass\":"<<mass<<",\"seconds\":"<<seconds
   <<",\"automorphism_order_histogram\":{";bool sep=false;
  for(auto[a,c]:aut_hist){if(sep)report<<",";sep=true;report<<"\""<<a<<"\":"<<c;}
  report<<"},\"selected_stable_cell_histogram\":{";sep=false;
  for(auto[a,c]:selected_sizes){if(sep)report<<",";sep=true;report<<"\""<<a<<"\":"<<c;}report<<"}}";report.flush();
  cout<<"order="<<N<<" classes="<<accepted<<" states="<<states<<" queries="<<queries<<" seconds="<<seconds<<"\n"<<flush;
  parents=move(children);
 }
 report<<"\n],\"status\":\"FINITE-EXACT\",\"checks\":"<<checks<<",\"hostile_automorphism_orders\":[81,27]}\n";
 return 0;
}catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
