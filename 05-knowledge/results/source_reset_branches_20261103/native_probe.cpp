// Find literal short native resets; stops at the first successful endpoint.
#define EAR_LIBRARY
#include "ear_probe.cpp"
#include <set>
using Path=std::vector<int>;
std::vector<Path> paths(const G&T){std::vector<Path>out;Path p;unsigned full=(1u<<T.size())-1;
#ifdef PERMUTATION_PATHS
 p.resize(T.size());std::iota(p.begin(),p.end(),0);
 do{bool ok=true;for(int i=1;i<(int)p.size();++i)if(!(T[p[i-1]]>>p[i]&1)){ok=false;break;}if(ok)out.push_back(p);}while(std::next_permutation(p.begin(),p.end()));
 return out;
#else
 auto dfs=[&](auto&&self,unsigned mask,int last)->void{if(mask==full){out.push_back(p);return;}unsigned next=(last<0?full:T[last])&~mask;while(next){unsigned b=next&-next;next-=b;int v=__builtin_ctz(b);p.push_back(v);self(self,mask|b,v);p.pop_back();}};
 dfs(dfs,0,-1);return out;
#endif
}
G native(const G&T,const Path&p){need(p.size()==T.size(),"spine length");G U=converse(T);unsigned mask=0;for(int v:p){need(!(mask>>v&1),"spine distinctness");mask|=1u<<v;}
 for(int j=1;j<(int)p.size();++j){int a=p[j-1],b=p[j];need(T[a]>>b&1,"Hamilton spine");U[a]|=1u<<b;U[b]&=~(1u<<a);}return U;}
std::array<std::vector<int>,2> branches(const G&T,I h){auto c=allcuts(T);std::array<std::vector<int>,2>b;for(int k=0;k<(int)c.size();++k)for(int j=0;j<2;++j)if(c[k]==3*h+2+2*j)b[j].push_back(k);return b;}
void printlist(const char*name,const std::vector<int>&v){std::cout<<name;for(auto x:v)std::cout<<' '<<x;std::cout<<'\n';}
bool strong(const G&T){for(auto&U:std::vector<G>{T,converse(T)}){unsigned mask=1;while(true){unsigned old=mask;for(int i=0;i<(int)U.size();++i)if(mask>>i&1)mask|=U[i];if(old==mask)break;}if(mask!=(1u<<T.size())-1)return false;}return true;}
#ifndef NATIVE_LIBRARY
int main(int argc,char**argv){need(argc>=3,"[--strong] maximum steps, rows");bool strong_only=std::string(argv[1])=="--strong";int pos=strong_only?2:1;int maxsteps=std::stoi(argv[pos]);G T;for(int i=pos+1;i<argc;++i)T.push_back(std::stoul(argv[i]));I h=H(T);auto ps=paths(T);need((I)ps.size()==h,"path count");long long count=0;std::set<G>seen;std::vector<G>first;std::vector<Path>fp;
 auto test=[&](const G&U,const std::vector<Path>&word)->bool{if(!seen.insert(U).second)return false;++count;if((strong_only&&!strong(U))||H(U)!=h)return false;auto b=branches(U,h);if(b[0].empty()||b[1].empty())return false;
  std::cout<<"SUCCESS H "<<h<<" probes "<<count<<" steps "<<word.size()<<'\n';for(auto&p:word)printlist("path",p);std::cout<<"target";for(auto r:U)std::cout<<' '<<r;std::cout<<'\n';printlist("plus2",b[0]);printlist("plus4",b[1]);return true;};
 if(test(T,{}))return 0;
 for(auto&p:ps){G U=native(T,p);first.push_back(U);fp.push_back(p);if(test(U,{p}))return 0;}
 if(maxsteps>=2)for(int i=0;i<(int)first.size();++i)for(auto&q:paths(first[i]))if(test(native(first[i],q),{fp[i],q}))return 0;
 std::cout<<"NO_SUCCESS H "<<h<<" probes "<<count<<" radius "<<maxsteps<<'\n';
}
#endif
