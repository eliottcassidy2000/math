// Fixed finite order12 families; no all-order12 completeness claim.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using U=std::uint64_t;using Rows=std::vector<unsigned>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
Rows decode(int n,U code){Rows a(n);int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit)if((code>>bit)&1)a[i]|=1u<<j;else a[j]|=1u<<i;return a;}
Rows opposite(const Rows&a){int n=a.size();Rows b(n);for(int i=0;i<n;++i)b[i]=((1u<<n)-1)^(1u<<i)^a[i];return b;}
std::array<unsigned char,4096> parity{};
unsigned parity_ends(const Rows&a){int n=a.size();unsigned full=(1u<<n)-1;std::array<unsigned,4096> dp{};Rows pred=opposite(a);
 for(unsigned mask=1;mask<=full;++mask){unsigned bits=mask,value=0;while(bits){int v=__builtin_ctz(bits);bits&=bits-1;unsigned rest=mask^(1u<<v);if(!rest||parity[dp[rest]&pred[v]])value|=1u<<v;}dp[mask]=value;}return dp[full];}
template<class V>std::vector<V> endpoints(const Rows&a,bool mod4){int n=a.size();unsigned full=(1u<<n)-1;std::vector<V> dp((full+1)*n);Rows pred=opposite(a);
 for(unsigned mask=1;mask<=full;++mask){unsigned bits=mask;while(bits){int v=__builtin_ctz(bits);bits&=bits-1;unsigned rest=mask^(1u<<v),before=rest&pred[v];U x=!rest;while(before){int u=__builtin_ctz(before);before&=before-1;x+=dp[rest*n+u];}dp[mask*n+v]=mod4?(x&3):x;}}return dp;}
struct Data{U h=0;std::vector<U>s,e,c;};
template<class V>Data data(const Rows&a,bool mod4){int n=a.size();unsigned full=(1u<<n)-1;auto e=endpoints<V>(a,mod4),s=endpoints<V>(opposite(a),mod4);Data d;d.s.resize(n);d.e.resize(n);d.c.resize(n*n);
 for(int i=0;i<n;++i){d.e[i]=e[full*n+i];d.s[i]=s[full*n+i];d.h+=d.e[i];}
 for(unsigned A=1;A<full;++A){unsigned us=A;while(us){int u=__builtin_ctz(us);us&=us-1;if(!e[A*n+u])continue;unsigned vs=full^A;while(vs){int v=__builtin_ctz(vs);vs&=vs-1;d.c[u*n+v]+=U(e[A*n+u])*s[(full^A)*n+v];}}}
 if(mod4){d.h&=3;for(auto&v:d.c)v&=3;}return d;}
bool symmetric4(const Data&d){int n=d.s.size();for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if((d.c[i*n+j]-d.c[j*n+i])&3)return false;return true;}
void compare(const Data&a,const Data&b){need((a.h&3)==b.h,"H mod4");for(unsigned i=0;i<a.s.size();++i)need((a.s[i]&3)==b.s[i]&&(a.e[i]&3)==b.e[i],"endpoints mod4");for(unsigned i=0;i<a.c.size();++i)need((a.c[i]&3)==b.c[i],"C mod4");}
Rows inflate(const Rows&q,const std::vector<Rows>&modules){int n=0;std::vector<int>offset;for(auto&m:modules){offset.push_back(n);n+=m.size();}Rows a(n);for(unsigned i=0;i<modules.size();++i){for(unsigned u=0;u<modules[i].size();++u)a[offset[i]+u]|=modules[i][u]<<offset[i];for(unsigned j=0;j<modules.size();++j)if((q[i]>>j)&1)for(unsigned u=0;u<modules[i].size();++u)for(unsigned v=0;v<modules[j].size();++v)a[offset[i]+u]|=1u<<(offset[j]+v);}return a;}
Rows birth(const Rows&old,unsigned outmask){int n=old.size();Rows a=old;a.push_back(outmask);for(int i=0;i<n;++i)if(!((outmask>>i)&1))a[i]|=1u<<n;return a;}
U next(U&state){state^=state<<13;state^=state>>7;state^=state<<17;return state;}
Rows random_rows(int n,U&state){Rows a(n);for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(next(state)&1)a[i]|=1u<<j;else a[j]|=1u<<i;return a;}
void flip(Rows&a,int i,int j){a[i]^=1u<<j;a[j]^=1u<<i;}
Rows self_converse12(U&state){Rows a(12);auto arc=[&](int u,int v){a[u]|=1u<<v;};for(int i=0;i<6;++i){if(next(state)&1)arc(i,i+6);else arc(i+6,i);for(int j=i+1;j<6;++j){if(next(state)&1){arc(i,j);arc(j+6,i+6);}else{arc(j,i);arc(i+6,j+6);}if(next(state)&1){arc(i,j+6);arc(j,i+6);}else{arc(j+6,i);arc(i+6,j);}}}return a;}
struct Stats{U tested=0,balanced=0,joint=0,integer_audits=0;};
struct Witness{std::string family;U index;Rows a;Data d;bool joint;};
std::map<std::string,Stats> stats;std::vector<Witness>witnesses,small_witnesses;U all=0,min_skew=67;Witness best;
void test(const std::string&family,const Rows&a,U index){need(a.size()==12,"order12 candidate");auto&z=stats[family];++z.tested;++all;unsigned e=parity_ends(a),s=parity_ends(opposite(a));bool balance=e==s;Data exact;
 if(balance||z.tested<=16||z.tested%4096==0){exact=data<U>(a,false);++z.integer_audits;for(int i=0;i<12;++i)need((exact.s[i]&1)==((s>>i)&1)&&(exact.e[i]&1)==((e>>i)&1),"fast predicate exact audit");}
 if(balance){++z.balanced;auto modular=data<unsigned char>(a,true);compare(exact,modular);bool joint=symmetric4(modular);z.joint+=joint;U skew=0;for(int i=0;i<12;++i)for(int j=i+1;j<12;++j)skew+=((modular.c[i*12+j]-modular.c[j*12+i])&3)!=0;if(skew<min_skew){min_skew=skew;best={family,index,a,exact,joint};}if(joint||z.balanced<=2)witnesses.push_back({family,index,a,exact,joint});}
 if(all%32768==0)std::cerr<<"tested="<<all<<" family="<<family<<" balanced="<<z.balanced<<" joint="<<z.joint<<"\n";
}
template<class V>void array(const std::vector<V>&a){std::cout<<"[";for(unsigned i=0;i<a.size();++i){if(i)std::cout<<",";std::cout<<U(a[i]);}std::cout<<"]";}
int main(int argc,char**argv){try{need(argc==2,"bank controls");for(int i=1;i<4096;++i)parity[i]=parity[i>>1]^(i&1);std::ifstream f(argv[1]);need(bool(f),"bank open");std::map<int,std::vector<Rows>>bank;for(int n=1;n<=8;++n){unsigned count;f>>count;for(unsigned i=0;i<count;++i){U code;f>>code;need(bool(f),"bank code");bank[n].push_back(decode(n,code));}}
 std::map<int,U>csym;std::map<int,U>jointsmall;U smallcomp=0;for(int n=1;n<=8;++n)for(const auto&a:bank[n]){auto d=data<unsigned char>(a,true);bool bal=true;for(int i=0;i<n;++i)bal&=((d.s[i]-d.e[i])&1)==0;csym[n]+=symmetric4(d);jointsmall[n]+=symmetric4(d)&&bal;if(n<=6||smallcomp%256==0)compare(data<U>(a,false),d);if(symmetric4(d))small_witnesses.push_back({"small_core_order"+std::to_string(n),0,a,data<U>(a,false),bal});++smallcomp;}
 // All labelled six-vertex quotients, each vertex replaced by a directed edge.
 for(U code=0;code<(U(1)<<15);++code)test("all_TT2_six_vertex_quotients",inflate(decode(6,code),std::vector<Rows>(6,Rows{2,0})),code);
 U ix=0;for(U q=0;q<64;++q)for(const auto&a:bank[4])for(const auto&b:bank[4])test("quotient4_modules_2_2_4_4",inflate(decode(4,q),{{2,0},{2,0},a,b}),ix++);
 ix=0;for(U q=0;q<8;++q)for(const auto&a:bank[4])for(const auto&b:bank[4])for(const auto&c:bank[4])test("quotient3_modules_4_4_4",inflate(decode(3,q),{a,b,c}),ix++);
 ix=0;for(const auto&a:bank[6])for(const auto&b:bank[6])test("ordered_pairs_of_six_vertex_classes",inflate({2,0},{a,b}),ix++);
 Rows transitive(11),cyclic(11),paley(11);for(int i=0;i<11;++i)for(int j=0;j<11;++j){if(i<j)transitive[i]|=1u<<j;int diff=(j-i+11)%11;if(diff>=1&&diff<=5)cyclic[i]|=1u<<j;if(diff==1||diff==3||diff==4||diff==5||diff==9)paley[i]|=1u<<j;}
 for(unsigned mask=0;mask<2048;++mask){test("TT11_all_births",birth(transitive,mask),mask);test("cyclic11_all_births",birth(cyclic,mask),mask);test("Paley11_all_births",birth(paley,mask),mask);}
 U paired_state=0x13579bdf2468ace1ULL;for(U k=0;k<2048;++k){Rows a=self_converse12(paired_state);test("self_converse12_fixed_involution",a,k);Rows b=a;int i=next(paired_state)%12,j=next(paired_state)%11;if(j>=i)++j;flip(b,i,j);test("self_converse12_one_flip",b,k);int u=next(paired_state)%12,v=next(paired_state)%11;if(v>=u)++v;flip(b,u,v);test("self_converse12_two_flip",b,k);}
 U state=0x123456789abcdefULL;for(U k=0;k<262144;++k)test("specified_random12",random_rows(12,state),k);
 if(min_skew<67)witnesses.push_back(best);witnesses.insert(witnesses.end(),small_witnesses.begin(),small_witnesses.end());
 std::cout<<"{\"small_core_C_symmetric_mod4\":{";bool first=true;for(auto[n,c]:csym){if(!first)std::cout<<",";first=false;std::cout<<"\""<<n<<"\":"<<c;}std::cout<<"},\"small_core_joint_survivors\":{";first=true;for(auto[n,c]:jointsmall){if(!first)std::cout<<",";first=false;std::cout<<"\""<<n<<"\":"<<c;}
 std::cout<<"},\"families\":{";first=true;for(auto&[name,z]:stats){if(!first)std::cout<<",";first=false;std::cout<<"\""<<name<<"\":{\"tested\":"<<z.tested<<",\"endpoint_parity_balanced\":"<<z.balanced<<",\"joint_survivors\":"<<z.joint<<",\"integer_audits\":"<<z.integer_audits<<"}";}
 std::cout<<"},\"witnesses\":[";first=true;for(auto&w:witnesses){if(!first)std::cout<<",";first=false;std::cout<<"{\"family\":\""<<w.family<<"\",\"index\":"<<w.index<<",\"joint_survivor\":"<<(w.joint?"true":"false")<<",\"rows\":";array(w.a);std::cout<<",\"H\":"<<w.d.h<<",\"S\":";array(w.d.s);std::cout<<",\"E\":";array(w.d.e);std::cout<<",\"C\":";array(w.d.c);std::cout<<"}";}
 std::cout<<"],\"total_candidates\":"<<all<<",\"minimum_nonzero_skew_pairs_order12\":"<<min_skew<<",\"random_final_state\":"<<state<<",\"paired_final_state\":"<<paired_state<<",\"scope\":\"Specified labelled families and deterministic streams; not all order12 classes\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
