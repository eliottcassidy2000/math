// Exhaust all fixed-Hamilton-path presentations, not a sample of classes.
// Canonical forms permute every equal-degree cell; ports, when present, are fixed.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using U=uint64_t; using T=std::array<unsigned,8>;
int n;
void need(bool v,const char* s){if(!v)throw std::runtime_error(s);}
T presentation(unsigned m){T t{};int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
 bool f=j==i+1 || ((m>>k)&1);if(j>i+1)k++;if(f)t[i]|=1u<<j;else t[j]|=1u<<i;}return t;}
U encode(const T&t,const std::vector<int>&p){U c=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)c|=U((t[p[i]]>>p[j])&1)<<(8*i+j);return c;}
U canonical(const T&t,bool roots=false){std::map<int,std::vector<int>> gr;std::vector<int> p;
 if(roots){p.push_back(0);p.push_back(n-1);}
 for(int i=0;i<n;i++)if(!roots||(i!=0&&i!=n-1))gr[__builtin_popcount(t[i])].push_back(i);
 std::vector<std::vector<int>> gs;for(auto&[d,v]:gr)gs.push_back(v);
 U best=~U(0);auto visit=[&](auto&&self,int k)->void{if(k==int(gs.size())){best=std::min(best,encode(t,p));return;}
 auto v=gs[k];do{auto old=p.size();p.insert(p.end(),v.begin(),v.end());self(self,k+1);p.resize(old);}while(std::next_permutation(v.begin(),v.end()));};visit(visit,0);return best;}
bool strong(const T&t){for(int rev=0;rev<2;rev++){unsigned seen=1;for(;;){unsigned next=seen;
 for(int i=0;i<n;i++)if(seen>>i&1)for(int j=0;j<n;j++)if((rev?t[j]>>i:t[i]>>j)&1)next|=1u<<j;
 if(next==seen)break;seen=next;}if(seen!=((1u<<n)-1))return false;}return true;}
bool blue(const T&t){for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(((t[i]>>j)&1)!=((t[n-1-j]>>(n-1-i))&1))return false;return true;}
struct Record{std::map<U,unsigned> joint,free,bl;unsigned mass=0;};
int main(int argc,char**argv){int top=argc>1?std::stoi(argv[1]):8;need(top>=3&&top<=8,"range3..8");std::cout<<"[\n";
 for(n=3;n<=top;n++){unsigned M=(n-1)*(n-2)/2,N=1u<<M;std::vector<U> fc(N),gc(N);std::set<U> all;
 for(unsigned m=0;m<N;m++){auto t=presentation(m);fc[m]=canonical(t);gc[m]=canonical(t,true);all.insert(fc[m]);}
 std::set<std::array<U,2>> fp,gp;std::set<std::array<U,3>> fpc,gpc;
 for(unsigned m=0;m<N;m++){U b=blue(presentation(m));unsigned k=m^(N-1);fp.insert({fc[m],fc[k]});gp.insert({gc[m],gc[k]});fpc.insert({fc[m],fc[k],b});gpc.insert({gc[m],gc[k],b});}
 std::map<U,Record> rec;std::set<U> target;unsigned total=0,maxf=0,blueatoms=0;
 for(unsigned m=0;m<N;m++){T t=presentation(m);if(!strong(t))continue;auto&r=rec[fc[m]];U h=gc[m^(N-1)];
 r.joint[h]++;r.free[fc[m^(N-1)]]++;r.mass++;if(blue(t))r.bl[h]++;target.insert(h);total++;}
 std::set<std::map<U,unsigned>> js,fs,bs;unsigned maxcollapse=0;U witness=0;
 for(auto&[a,r]:rec){js.insert(r.joint);fs.insert(r.free);if(!r.bl.empty())blueatoms++;
 unsigned mass=0;for(auto[h,c]:r.free)mass+=c;need(mass==r.mass,"free pushforward conserves mass");
 for(auto[h,c]:r.joint)maxf=std::max(maxf,c);
 // Compare actual pushforward with the old dictionary-overwrite bug.
 std::map<U,unsigned> wrong;for(auto[h,c]:r.joint){T t{};for(int i=0;i<n;i++)t[i]=(h>>(8*i))&255;wrong[canonical(t)]=c;}
 unsigned wm=0;for(auto[h,c]:wrong)wm+=c;if(r.mass-wm>maxcollapse){maxcollapse=r.mass-wm;witness=a;}}
 std::cout<<"{\"n\":"<<n<<",\"presentations\":"<<N<<",\"classes\":"<<all.size()<<",\"strong_classes\":"<<rec.size()<<",\"strong_presentations\":"<<total<<",\"joint_targets\":"<<target.size()<<",\"joint_signatures\":"<<js.size()<<",\"free_signatures\":"<<fs.size()<<",\"blue_atoms\":"<<blueatoms<<",\"max_joint_fibre\":"<<maxf<<",\"old_overwrite_max_mass_loss\":"<<maxcollapse<<",\"mass_loss_source_code\":"<<witness<<",\"free_pair_cells\":"<<fp.size()<<",\"rooted_pair_cells\":"<<gp.size()<<",\"free_pair_colour_cells\":"<<fpc.size()<<",\"rooted_pair_colour_cells\":"<<gpc.size()<<"}"<<(n<top?",\n":"\n");std::cout.flush();
 }std::cout<<"]\n";
}
