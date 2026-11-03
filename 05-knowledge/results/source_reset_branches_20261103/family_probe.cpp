// Search a fully specified two-native reset template for uniform clone lifts.
#define NATIVE_LIBRARY
#include "native_probe.cpp"
std::vector<I> clone_coeff(const G&T,int v){Path p;for(int i=0;i<(int)T.size();++i)if(i!=v)p.push_back(i);std::vector<I>c(T.size()+1);
 do{int t=((T[v]>>p.front())&1)+((T[p.back()]>>v)&1),q=0;bool good=true;
  for(int i=1;i<(int)p.size();++i){int x=p[i-1],y=p[i];bool slot=(T[x]>>v&1)&&(T[v]>>y&1);t+=slot;if(!(T[x]>>y&1)){++q;if(!slot){good=false;break;}}}
  if(good){I b=1;for(int j=0;j<=q;++j){c[t-j]+=(j%2?-1:1)*b;b=b*(q-j)/(j+1);}}
 }while(std::next_permutation(p.begin(),p.end()));c[0]=0;while(c.size()>1&&c.back()==0)c.pop_back();return c;
}
std::pair<int,int> neighbours(const Path&p,int v){int i=std::find(p.begin(),p.end(),v)-p.begin();return{i? p[i-1]:-1,i+1<(int)p.size()?p[i+1]:-1};}
bool affine(const std::vector<I>&a,const std::vector<I>&b,int d){for(int j=1;j<(int)std::max(a.size(),b.size());++j)if((j<(int)b.size()?b[j]:0)!=3*(j<(int)a.size()?a[j]:0)+(j==1?d:0))return false;return true;}
int main(int argc,char**argv){need(argc>=2,"rows");G T;for(int i=1;i<argc;++i)T.push_back(std::stoul(argv[i]));I h=H(T);auto ps=paths(T);std::vector<std::vector<I>>cs;for(int v=0;v<(int)T.size();++v)cs.push_back(clone_coeff(T,v));
 long long pairs=0,neutral=0,branching=0,inactive=0,lifts=0;std::set<std::pair<G,int>>seen;
 for(auto&p:ps){G U=native(T,p);for(auto&q:paths(U)){++pairs;G V=native(U,q);if(H(V)!=h)continue;++neutral;auto b=branches(V,h);if(b[0].empty()||b[1].empty())continue;++branching;
  for(int v=0;v<(int)T.size();++v)if(neighbours(p,v)==neighbours(q,v)&&seen.insert({V,v}).second){++inactive;if(cs[v]!=clone_coeff(V,v))continue;++lifts;std::array<std::vector<int>,2>hits;
   for(int j=0;j<2;++j)for(int cut:b[j])if(affine(cs[v],clone_coeff(ear(V,cut),v),2+2*j))hits[j].push_back(cut);
   std::cout<<"LIFT vertex "<<v<<" branches "<<hits[0].size()<<' '<<hits[1].size()<<" coeff";for(auto x:cs[v])std::cout<<' '<<x;std::cout<<'\n';printlist("path",p);printlist("path",q);std::cout<<"target";for(auto x:V)std::cout<<' '<<x;std::cout<<'\n';printlist("plus2",hits[0]);printlist("plus4",hits[1]);
  }
 }}
 std::cout<<"TOTAL pairs "<<pairs<<" neutral "<<neutral<<" branching "<<branching<<" inactive "<<inactive<<" clone_lifts "<<lifts<<'\n';
}
