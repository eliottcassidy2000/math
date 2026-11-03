// Complete radius-two native return scan with a requested full mod-seven port state.
#define NATIVE_LIBRARY
#include "native_probe.cpp"
std::vector<std::array<I,3>> observers(const G&T){int n=T.size();unsigned full=(1u<<n)-1;auto e=ends(T),s=ends(converse(T));std::vector<I>hs(1u<<n);hs[0]=1;
 for(unsigned mask=1;mask<=full;++mask)for(int r=0;r<n;++r)hs[mask]+=e[mask*n+r];std::vector<std::array<I,3>>v;
 for(int r=0;r<n;++r){I C=0;for(unsigned mask=1;mask<=full;++mask)C+=s[mask*n+r]*hs[full^mask];need((C-e[full*n+r])%2==0,"integral observer");v.push_back({hs[full],(C-e[full*n+r])/2,(C+e[full*n+r])/2});}return v;}
int main(int argc,char**argv){need(argc>=4,"desired x y; source rows");I x=std::stoll(argv[1]),y=std::stoll(argv[2]);G T;for(int i=3;i<argc;++i)T.push_back(std::stoul(argv[i]));I h=H(T);auto ps=paths(T);std::set<G>seen;I checked=0,eligible=0,hits=0;
 for(auto&p:ps){G U=native(T,p);for(auto&q:paths(U)){G V=native(U,q);if(!seen.insert(V).second)continue;++checked;if(H(V)!=h)continue;auto b=branches(V,h);if(b[0].empty()||b[1].empty())continue;++eligible;auto vv=observers(V);
  for(int r=0;r<(int)V.size();++r)if(vv[r][1]%7==x%7&&vv[r][2]%7==y%7){++hits;std::cout<<"HIT root "<<r<<" state "<<vv[r][0]<<' '<<vv[r][1]<<' '<<vv[r][2]<<'\n';printlist("path",p);printlist("path",q);std::cout<<"target";for(auto z:V)std::cout<<' '<<z;std::cout<<'\n';printlist("plus2",b[0]);printlist("plus4",b[1]);}
 }}std::cout<<"TOTAL endpoints "<<checked<<" branching "<<eligible<<" phase_hits "<<hits<<'\n';
}
