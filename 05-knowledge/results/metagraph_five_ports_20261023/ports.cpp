// Focused COMPLETE target fibres: all open Hamiltonian presentations of every
// class in the unresolved order-nine one-round-v4 buckets. This is not a
// complete order-nine metagraph enumeration. Canonicalization is inherited;
// every output row carries literal target-path and source-isomorphism witnesses.
#define main inherited_blue_census_main
#include "../blue_quotient_20261020b/census.cpp"
#undef main
#include <unordered_set>
static uint64_t checks2=0;
static void require2(bool x,string why){checks2++;if(!x)throw runtime_error(why);}
static void write32(ofstream&f,uint32_t z){for(int i=0;i<4;i++)f.put((z>>(8*i))&255);}
static void write64(ofstream&f,uint64_t z){for(int i=0;i<8;i++)f.put((z>>(8*i))&255);}
static int selected_lut[1024][5];
static vector<array<int,5>> subsets5;
static vector<pair<int,int>> selected_profile(const Rows&a){
 array<pair<int,int>,9>count{};for(auto S:subsets5){int code=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)code=2*code+((a[S[i]]>>S[j])&1);for(int i=0;i<5;i++){int selected=selected_lut[code][i];count[S[i]].first+=(selected&1)!=0;count[S[i]].second+=(selected&2)!=0;}}
 vector<pair<int,int>>v(count.begin(),count.end());sort(v.begin(),v.end());return v;
}
int main(int argc,char**argv){try{
 require2(argc==4,"usage ports all_code_aut.tsv targets.tsv output_directory");string dir=argv[3];ifstream universe(argv[1]);Code code;int aut;map<Code,int> auts;uint64_t mass=0;
 while(universe>>code>>aut){require2(auts.emplace(code,aut).second,"duplicate class");require2(362880%aut==0,"aut divisor");mass+=362880/aut;}
 require2(auts.size()==191536&&mass==(1ULL<<36),"n9 universe");
 for(int code=0;code<1024;code++){Rows a=decode(code,5);for(int root=0;root<5;root++){vector<int>p{root};for(int v=0;v<5;v++)if(v!=root)p.push_back(v);uint64_t best=~0ULL;do{best=min(best,encode(a,p));}while(next_permutation(p.begin()+1,p.end()));selected_lut[code][root]=(best==66?1:0)|(best==204?2:0);}}
 for(int a=0;a<9;a++)for(int b=a+1;b<9;b++)for(int c=b+1;c<9;c++)for(int d=c+1;d<9;d++)for(int e=d+1;e<9;e++)subsets5.push_back({a,b,c,d,e});require2(subsets5.size()==126,"five subset universe");
 ifstream select(argv[2]);vector<Code>targets;while(select>>code){require2(auts.count(code),"target absent");targets.push_back(code);}require2(targets.size()==201,"target selection count");
 vector<pair<int,int>>cells;for(int i=0;i<9;i++)for(int j=i+2;j<9;j++)cells.push_back({i,j});array<int,28>rev{};for(int k=0;k<28;k++)rev[k]=find(cells.begin(),cells.end(),make_pair(8-cells[k].second,8-cells[k].first))-cells.begin();
 ofstream binary(dir+"/focused_ports.bin",ios::binary),summary(dir+"/focused_targets.tsv",ios::binary);summary<<"target\taut\tH\topen_paths\topen_presentations\tblue_presentations\n";
 uint64_t totalhp=0,totalopen=0,totalmasks=0,totalblue=0;unordered_set<uint32_t>seen;set<Code>sources;
 for(Code target:targets){Rows a=decode(target,9);auto expected_profile=selected_profile(a);uint64_t H=ham(a,9),paths=0,opened=0;map<uint32_t,pair<int,Code>> masks;vector<int>P;
  function<void(unsigned)>dfs=[&](unsigned used){if(P.size()==9){paths++;if(!(a[P[0]]>>P[8]&1))return;opened++;uint32_t targetmask=0;for(int k=0;k<28;k++)if(a[P[cells[k].first]]>>P[cells[k].second]&1)targetmask|=1u<<k;uint32_t mask=targetmask^((1u<<28)-1);auto&v=masks[mask];v.first++;if(v.first==1)v.second=packed(P);return;}unsigned avail=(P.empty()?511:a[P.back()])&~used;while(avail){int v=__builtin_ctz(avail);avail&=avail-1;P.push_back(v);dfs(used|1u<<v);P.pop_back();}};dfs(0);
  require2(paths==H,"independent DFS/Held Karp count");uint64_t blue=0;
  for(auto&[mask,multiplicity_path]:masks){require2(multiplicity_path.first==auts.at(target),"Aut orbit multiplicity");require2(seen.insert(mask).second,"target fibre collision");Rows s{};for(int i=0;i<8;i++)s[i]|=1u<<(i+1);for(int k=0;k<28;k++){auto[i,j]=cells[k];if(mask>>k&1)s[i]|=1u<<j;else s[j]|=1u<<i;}require2(!(s[0]>>8&1),"source not closed");auto c=canonical(s,9);require2(auts.count(c.code)&&c.aut==auts.at(c.code),"source owner outside universe");require2(encode(s,unpack(c.order,9))==c.code,"source witness");sources.insert(c.code);uint32_t reflected=0;for(int k=0;k<28;k++)if(mask>>k&1)reflected|=1u<<rev[k];blue+=(mask==reflected);
   Rows transported=converse(s,9);for(int i=0;i<8;i++){transported[i]^=1u<<(i+1);transported[i+1]^=1u<<i;}require2(selected_profile(transported)==expected_profile,"literal transported rooted-five pair");
   write32(binary,mask);write64(binary,c.code);write64(binary,target);write64(binary,c.order);write64(binary,multiplicity_path.second);
  }
  totalhp+=paths;totalopen+=opened;totalmasks+=masks.size();totalblue+=blue;summary<<target<<'\t'<<auts.at(target)<<'\t'<<H<<'\t'<<opened<<'\t'<<masks.size()<<'\t'<<blue<<'\n';
  if(totalmasks%25000<masks.size())cerr<<"focused target fibres "<<totalmasks<<" literal ports\n";
 }
 ofstream out(dir+"/port_producer.json",ios::binary);out<<"{\"target_classes\":"<<targets.size()<<",\"source_classes\":"<<sources.size()<<",\"paths\":"<<totalhp<<",\"open_paths\":"<<totalopen<<",\"closed_source_presentations\":"<<totalmasks<<",\"blue_presentations\":"<<totalblue<<",\"checks\":"<<checks2<<"}\n";cout<<"target classes "<<targets.size()<<" source classes "<<sources.size()<<" literal ports "<<totalmasks<<" checks "<<checks2<<'\n';
 }catch(exception&e){cerr<<e.what()<<'\n';return 1;}return 0;}
