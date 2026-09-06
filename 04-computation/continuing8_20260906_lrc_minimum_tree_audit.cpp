#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;
using I=long long;
using Word=array<int,7>;
const I INF=10000000000LL;
uint64_t gates=0, unpruned=0, wordchecks=0, evaluations=0, compatible=0;
void need(bool x,const string& s){++gates;if(!x){cerr<<"FAILED "<<s<<"\n";exit(1);}}
I mod(I a,I b){I r=a%b;return r<0?r+b:r;}
I floorq(I a,I b){return (a-mod(a,b))/b;}
I ceilq(I a,I b){return -floorq(-a,b);}
int gd[91][91];
struct Geo{int p,q;I L;vector<pair<I,I>> intervals;};
bool danger(I x,I den,int p){I z=mod(p*x,den);return 14*min(z,den-z)<den;}
Geo geometry(int p,int q){
    Geo G{p,q,14LL*p*q,{}};
    vector<I> walls={0,G.L};
    for(int k=0;k<=p;k++)for(int sign:{-1,1}){I x=(14LL*k+sign)*q;if(0<x&&x<G.L)walls.push_back(x);}
    for(int k=0;k<=q;k++)for(int sign:{-1,1}){I x=(14LL*k+sign)*p;if(0<x&&x<G.L)walls.push_back(x);}
    sort(walls.begin(),walls.end());walls.erase(unique(walls.begin(),walls.end()),walls.end());
    vector<pair<I,I>> cells;
    for(size_t i=1;i<walls.size();i++){
        I a=walls[i-1],b=walls[i];
        if(danger(a+b,2*G.L,p)&&danger(a+b,2*G.L,q)){
            if(!cells.empty()&&cells.back().second==a&&danger(a,G.L,p)&&danger(a,G.L,q))cells.back().second=b;
            else cells.push_back({a,b});
        }
    }
    need(cells.size()>=2&&cells.front().first==0&&cells.back().second==G.L,"raw wall origin components");
    G.intervals.push_back({cells.back().first-G.L,cells.front().second});
    G.intervals.insert(G.intervals.end(),cells.begin()+1,cells.end()-1);
    return G;
}
I count_at(const Geo&G,int n,I phase2){
    I ans=0;for(auto [a,b]:G.intervals)ans+=ceilq(2LL*n*b-phase2,2*G.L)-floorq(2LL*n*a-phase2,2*G.L)-1;return ans;
}
pair<I,I> capacity(const Geo&G,int n){
    map<I,pair<int,int>> events;
    for(auto [a,b]:G.intervals){events[mod(n*a,G.L)].first++;events[mod(n*b,G.L)].second++;}
    I best=count_at(G,n,0),phase=0,cur=count_at(G,n,1);
    if(cur<best){best=cur;phase=1;}
    for(auto [wall,se]:events){
        if(wall==0)continue;
        I at=cur-se.second;
        if(at<best){best=at;phase=2*wall;}
        cur=at+se.first;
        if(cur<best){best=cur;phase=2*wall+1;}
    }
    return {best,phase};
}
I literal(const Geo&G,int n,I phase2){
    I ans=0,den=2*G.L*n;
    for(int j=0;j<n;j++){I x=2*G.L*j+phase2;if(danger(x,den,G.p)&&danger(x,den,G.q))ans++;}
    return ans;
}
vector<Geo> atlas(){
    vector<int> primes;
    for(int p=2;p<=356;p++){bool prime=true;for(int d=2;d*d<=p;d++)if(p%d==0){prime=false;break;}if(prime&&p%3==2)primes.push_back(p);}
    set<int> sums={1};
    for(int p:primes){vector<int> old(sums.begin(),sums.end());for(int s:old){if(s*p<=356)sums.insert(s*p);if((I)s*p*p<=356)sums.insert(s*p*p);}}
    vector<Geo> out;
    for(int s:sums)if(s>=3)for(int p=1;2*p<s;p++)if(gcd(p,s-p)==1)out.push_back(geometry(p,s-p));
    need(out.size()==5855,"multiplicatively generated strict atlas");return out;
}
array<unordered_set<uint64_t>,7> profiles;
vector<int> masks;
uint64_t key(int c,const vector<int>& w){uint64_t k=c;for(int x:w)k=(k<<7)|x;return k;}
bool valid(const Word&w){
    int all=0;for(int x:w)all=gd[all][x];if(all!=1)return false;
    int cache[128]={0};
    auto get=[&](auto&&self,int m)->int{
        if(cache[m])return cache[m];int bit=__builtin_ctz((unsigned)m),rest=m&(m-1);
        return cache[m]=rest?gd[w[bit]][self(self,rest)]:w[bit];
    };
    for(int m:masks){
        int c=get(get,m),k=7-__builtin_popcount((unsigned)m);array<int,6> a{};int z=0;
        for(int j=0;j<7;j++)if(!(m&(1<<j)))a[z++]=gd[c][w[j]];
        sort(a.begin(),a.begin()+k);uint64_t v=c;for(int j=0;j<k;j++)v=(v<<7)|a[j];
        if(!profiles[k].count(v))return false;
    }
    return true;
}
void enumerate(const vector<int>&D,Word&w,int pos,int first,vector<Word>&out){
    if(pos==7){unpruned++;if(valid(w))out.push_back(w);return;}
    for(int j=first;j<(int)D.size();j++){w[pos]=D[j];enumerate(D,w,pos+1,j,out);}
}
void write_words(const vector<Word>&W,const string&path){
    ofstream f(path,ios::binary);f<<'[';bool first=true;
    for(auto&w:W){if(!first)f<<',';first=false;f<<'[';for(int j=0;j<7;j++){if(j)f<<',';f<<w[j];}f<<']';}f<<']';
}
I mst(const Word&w,const array<array<I,91>,91>& weights){
    array<array<I,3>,21> edges;int z=0;
    for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)edges[z++]={weights[w[i]][w[j]],i,j};
    sort(edges.begin(),edges.end());array<int,7> parent; iota(parent.begin(),parent.end(),0);
    auto find=[&](int x){while(parent[x]!=x)x=parent[x];return x;};
    I total=0;int used=0;
    for(auto ed:edges){int a=find(ed[1]),b=find(ed[2]);if(a!=b){parent[a]=b;total+=ed[0];if(++used==6)break;}}
    return total;
}
struct Expected{int t,domain;I words,minimum,ownerE,ownerTree,expected_events,expected_compatible;Word owner;vector<array<I,9>> survivors;map<pair<int,int>,vector<I>> weights;};
int main(int argc,char**argv){
    if(argc!=3){cerr<<"input and output directory required\n";return 2;}
    ifstream in(argv[1]);string dir=argv[2];need(bool(in),"native input open");
    for(int a=0;a<=90;a++)for(int b=0;b<=90;b++)gd[a][b]=gcd(a,b);
    int np;in>>np;for(int j=0;j<np;j++){int k,c;in>>k>>c;vector<int>w(k);for(int&x:w)in>>x;profiles[k].insert(key(c,w));}
    for(int size=6;size>=1;size--)for(int m=1;m<127;m++)if(__builtin_popcount((unsigned)m)==size)masks.push_back(m);
    need(masks.size()==126,"all126 proper nonempty positional subsets");
    int nd;in>>nd;vector<vector<int>> domains(nd);for(auto&D:domains){int k;in>>k;D.resize(k);for(int&x:D)in>>x;}
    int nt;in>>nt;vector<Expected> ex(nt);
    for(auto&E:ex){
        int ns,nw;in>>E.t>>E.domain>>E.words>>E.minimum>>E.ownerE>>E.ownerTree>>E.expected_events>>E.expected_compatible;for(int&x:E.owner)in>>x;
        in>>ns;E.survivors.resize(ns);for(auto&s:E.survivors)for(I&x:s)in>>x;
        in>>nw;for(int j=0;j<nw;j++){int a,b,k;in>>a>>b>>k;vector<I>v(k);for(I&x:v)in>>x;E.weights[{a,b}]=v;}
    }
    need(bool(in),"complete native input");
    auto AT=atlas();
    set<pair<int,int>> atlas_pairs;for(auto&G:AT)atlas_pairs.insert({G.p,G.q});
    for(auto pq:vector<pair<int,int>>{{1,1},{1,13},{1,112},{1,355},{11,263},{23,323},{33,278},{66,289}}){
        auto G=geometry(pq.first,pq.second);
        for(int n:{1,2,7,13,14,28,29,720,968,2485}){
            auto [best,phase]=capacity(G,n);set<I>walls={0};for(auto ab:G.intervals){walls.insert(mod(n*ab.first,G.L));walls.insert(mod(n*ab.second,G.L));}
            I brute=INF;for(I x:walls){brute=min(brute,count_at(G,n,2*x));brute=min(brute,count_at(G,n,2*x+1));}
            need(best==brute,"all-wall direct ceil-floor control");need(best==literal(G,n,phase),"literal strict original-grid at minimizing phase");
        }
    }
    need(capacity(geometry(1,355),2485).first==0,"exact near-resonance endpoint hostile");
    vector<vector<Word>> banks(nd);I distinct=0;
    for(int d=0;d<nd;d++){
        Word w{};enumerate(domains[d],w,0,0,banks[d]);distinct+=banks[d].size();
        write_words(banks[d],dir+"/words_"+to_string(d)+".json");
        cerr<<"DOMAIN "<<d+1<<'/'<<nd<<" valid="<<banks[d].size()<<"\n";
    }
    I minupper=INF,totalwords=0;int excluded=0,hostiles=0;I pairtables=0;
    for(auto&E:ex){
        auto&D=domains[E.domain];auto&bank=banks[E.domain];need((I)bank.size()==E.words,"unpruned full word count at "+to_string(E.t));
        array<array<I,91>,91> W;for(auto&row:W)row.fill(INF);set<int>es;set<int>states(D.begin(),D.end());
        for(int a:D)for(int b:D)es.insert(gcd(a,b));
        I localcompatible=0,localevents=0;
        for(int e:es){int n=E.t/e;for(auto&G:AT){
            int a=e*gcd(n,G.p),b=e*gcd(n,G.q);if(a>b)swap(a,b);if(!states.count(a)||!states.count(b))continue;
            localcompatible++;I sep=0;for(auto [lo,hi]:G.intervals)sep+=floorq((I)n*(hi-lo)-1,G.L);sep*=e;
            if(sep>=W[a][b])continue;
            I v=e*capacity(G,n).first;evaluations++;localevents++;need(sep<=v&&v<=E.t,"separate open-arc floor lower bound");
            if(v<W[a][b])W[a][b]=W[b][a]=v;
        }}
        compatible+=localcompatible;
        if(E.expected_events>=0)need(localevents==E.expected_events,"complete event-evaluation count");
        if(E.expected_compatible>=0)need(localcompatible==E.expected_compatible,"complete compatible-edge count");
        need(E.weights.size()==D.size()*(D.size()+1)/2,"complete forced-pair table");
        for(auto [ab,v]:E.weights){
            int a=ab.first,b=ab.second;pairtables++;need(W[a][b]==v[0],"independent raw-wall minimum pair credit at "+to_string(E.t));
            if(v.size()>1){int e=v[1],p=v[2],q=v[3];auto G=geometry(p,q);int n=E.t/e;
                need(e==gcd(a,b)&&atlas_pairs.count({p,q}),"attainer gcd and exact strict-atlas membership");
                pair<int,int> got={e*gcd(n,p),e*gcd(n,q)};if(got.first>got.second)swap(got.first,got.second);
                need(got==ab,"attainer actual state pair");need(e*capacity(G,n).first==v[0],"attainer count reconstruction");
                I sep=0;for(auto [lo,hi]:G.intervals)sep+=floorq((I)n*(hi-lo)-1,G.L);need(e*sep==v[4],"attainer separate-open-interval column");
                if(E.t==7200||E.t==12672){auto cp=capacity(G,n);need(cp.first==literal(G,n,cp.second),"complete hostile/minimum-clock attainers on literal native grids");}
            }else need(v[0]==INF,"impossible-pair infinity convention");
        }
        I worst=INF;vector<array<I,9>> bad;
        for(auto&w:bank){
            I Ecost=-E.t,res=0;for(int d:w){Ecost+=(I)d*ceilq(E.t,7*d);res+=(I)d*mod(-E.t/d,7);}
            need(7*Ecost==res&&Ecost>=0,"native marginal residue identity");
            I tree=mst(w,W),margin=tree-Ecost;worst=min(worst,margin);wordchecks++;totalwords++;
            if(margin<=0){array<I,9>row{};for(int j=0;j<7;j++)row[j]=w[j];row[7]=Ecost;row[8]=tree;bad.push_back(row);}
        }
        need(worst==E.minimum,"all-word Kruskal minimum at "+to_string(E.t));need(bad==E.survivors,"complete compatible cost survivors at "+to_string(E.t));
        I oe=-E.t;for(int d:E.owner)oe+=(I)d*ceilq(E.t,7*d);
        need(valid(E.owner)&&oe==E.ownerE&&mst(E.owner,W)==E.ownerTree&&E.ownerTree-oe==worst,"producer extremizer independently valid and attaining");
        if(E.t>=12000){need(bad.empty(),"upper clock closes");minupper=min(minupper,worst);excluded++;}else{need(E.t==7200&&bad.size()==15&&worst==-74,"7200 lower-envelope stopping control");hostiles++;}
        cerr<<"CLOCK "<<E.t<<" words="<<bank.size()<<" margin="<<worst<<" survivors="<<bad.size()<<"\n";
    }
    cout<<"ATLAS 5855; PROPER_SUBSETS 126; DOMAIN_COUNT "<<nd<<"\n";
    cout<<"UNPRUNED_WORDS "<<unpruned<<"; DISTINCT_VALID_WORDS "<<distinct<<"; WORD_CLOCK_EVALUATIONS "<<totalwords<<"\n";
    cout<<"PAIR_TABLE_ENTRIES "<<pairtables<<"; COMPATIBLE_ATLAS_EDGES "<<compatible<<"; EXACT_PROFILE_EVALUATIONS "<<evaluations<<"\n";
    cout<<"UPPER_CLOCKS_EXCLUDED "<<excluded<<"; MINIMUM_MARGIN "<<minupper<<"; HOSTILE_7200_COUNT 15; HOSTILE_MARGIN -74\n";
    need(excluded==549&&hostiles==1,"declared complete selected scope");
    cout<<"PASS "<<gates<<" always-active native exact gates\n";
}
