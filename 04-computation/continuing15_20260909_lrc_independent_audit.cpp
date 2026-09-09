// Self-contained adaptation of the independently audited continuing14 consumer.
// New bridge census uses literal edge deletion; large composite products use int128.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#include <tuple>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif
using namespace std;
using I=long long;
using Word=array<int,7>;
uint64_t gates=0,raw_words=0,capacity_evaluations=0,zero_ratios=0,open_zero_ratios=0,literal_points=0;
void need(bool ok,const string& why){++gates;if(!ok){cerr<<"FAILED "<<why<<'\n';exit(1);}}
I mod(I a,I b){I r=a%b;return r<0?r+b:r;}
I floorq(I a,I b){return (a-mod(a,b))/b;}
I ceilq(I a,I b){return -floorq(-a,b);}
bool danger(I x,I denominator,int speed){I r=(I)(((__int128)speed*x)%denominator);if(r<0)r+=denominator;return (__int128)14*min(r,denominator-r)<denominator;}
struct Geometry{int p,q;I denominator;vector<pair<I,I>> intervals;};
Geometry geometry(int p,int q){
    Geometry G{p,q,14LL*p*q,{}};
    vector<I> boundaries={0,G.denominator};
    for(int j=0;j<=p;j++)for(int sign:{-1,1}){I x=(14LL*j+sign)*q;if(0<x&&x<G.denominator)boundaries.push_back(x);}
    for(int j=0;j<=q;j++)for(int sign:{-1,1}){I x=(14LL*j+sign)*p;if(0<x&&x<G.denominator)boundaries.push_back(x);}
    sort(boundaries.begin(),boundaries.end());boundaries.erase(unique(boundaries.begin(),boundaries.end()),boundaries.end());
    vector<pair<I,I>> cells;
    for(size_t j=1;j<boundaries.size();j++){
        I a=boundaries[j-1],b=boundaries[j];
        if(danger(a+b,2*G.denominator,p)&&danger(a+b,2*G.denominator,q)){
            if(!cells.empty()&&cells.back().second==a&&danger(a,G.denominator,p)&&danger(a,G.denominator,q))cells.back().second=b;
            else cells.push_back({a,b});
        }
    }
    need(cells.size()>=2&&cells.front().first==0&&cells.back().second==G.denominator,"raw midpoint geometry retains circle wrap");
    G.intervals.push_back({cells.back().first-G.denominator,cells.front().second});
    G.intervals.insert(G.intervals.end(),cells.begin()+1,cells.end()-1);
    return G;
}
I interval_count(const Geometry&G,int n,I phase2){
    I value=0;
    for(auto[a,b]:G.intervals)value+=ceilq(2LL*n*b-phase2,2*G.denominator)-floorq(2LL*n*a-phase2,2*G.denominator)-1;
    return value;
}
struct Capacity{I ordinary,chamber,ordinary_phase,chamber_phase;};
Capacity capacity(const Geometry&G,int n){
    vector<pair<I,int>> events;
    for(auto[a,b]:G.intervals){events.push_back({mod(n*a,G.denominator),1});events.push_back({mod(n*b,G.denominator),-1});}
    sort(events.begin(),events.end());
    I initial=interval_count(G,n,0),cur=interval_count(G,n,1);
    Capacity out{min(initial,cur),cur,initial<=cur?0:1,1};
    size_t j=0;
    while(j<events.size()){
        I wall=events[j].first;int starts=0,ends=0;
        while(j<events.size()&&events[j].first==wall){if(events[j].second==1)starts++;else ends++;j++;}
        if(wall==0)continue;
        I at=cur-ends;
        if(at<out.ordinary){out.ordinary=at;out.ordinary_phase=2*wall;}
        cur=at+starts;
        if(cur<out.ordinary){out.ordinary=cur;out.ordinary_phase=2*wall+1;}
        if(cur<out.chamber){out.chamber=cur;out.chamber_phase=2*wall+1;}
    }
    return out;
}
array<I,3> literal(const Geometry&G,int n,I phase2){
    array<I,3> counts{};I denominator=2*G.denominator*n;
    for(int j=0;j<n;j++){
        I x=2*G.denominator*j+phase2;bool a=danger(x,denominator,G.p),b=danger(x,denominator,G.q);
        counts[0]+=a;counts[1]+=b;counts[2]+=a&&b;literal_points++;
    }
    return counts;
}
vector<Geometry> atlas(){
    vector<int> primes;
    for(int p=2;p<=356;p++){
        bool prime=true;for(int d=2;d*d<=p;d++)if(p%d==0){prime=false;break;}
        if(prime&&p%3==2)primes.push_back(p);
    }
    set<int> sums={1};
    for(int p:primes){vector<int> before(sums.begin(),sums.end());for(int v:before){if(v*p<=356)sums.insert(v*p);if((I)v*p*p<=356)sums.insert(v*p*p);}}
    vector<Geometry> out;
    for(int v:sums)if(v>=3)for(int p=1;2*p<v;p++)if(gcd(p,v-p)==1)out.push_back(geometry(p,v-p));
    need(out.size()==5855,"independently multiplicative full strict atlas");return out;
}
void controls(){
    for(auto[p,q]:vector<pair<int,int>>{{1,1},{1,13},{1,112},{1,355},{11,263},{23,323},{33,320},{48,307}}){
        auto G=geometry(p,q);
        for(int n:{7,14,29,651,1260,1890,1939,1953,2485}){
            set<I> walls={0};for(auto[a,b]:G.intervals){walls.insert(mod(n*a,G.denominator));walls.insert(mod(n*b,G.denominator));}
            I ordinary=1000000000,chamber=1000000000;
            for(I w:walls){ordinary=min(ordinary,interval_count(G,n,2*w));I value=interval_count(G,n,2*w+1);ordinary=min(ordinary,value);chamber=min(chamber,value);}
            auto found=capacity(G,n);
            need(found.ordinary==ordinary&&found.chamber==chamber,"every direct wall and chamber versus independent signed sweep");
            need(literal(G,n,found.ordinary_phase)[2]==ordinary,"literal control native intersection");
        }
    }
}
int common[91][91];
array<unordered_set<uint64_t>,7> profiles;
vector<int> inside_masks;
bool valid(const Word&w){
    int total=0;for(int x:w)total=common[total][x];if(total!=1)return false;
    int gcd_subset[128]={};
    for(int mask=1;mask<128;mask++){
        int index=__builtin_ctz((unsigned)mask);gcd_subset[mask]=common[gcd_subset[mask&(mask-1)]][w[index]];
    }
    for(int selected:inside_masks){
        int c=gcd_subset[127^selected],k=__builtin_popcount((unsigned)selected),j=0;
        array<int,6> restricted{};
        for(int i=0;i<7;i++)if(selected&(1<<i))restricted[j++]=common[c][w[i]];
        sort(restricted.begin(),restricted.begin()+k);uint64_t signature=c;
        for(int i=0;i<k;i++)signature=(signature<<7)|restricted[i];
        if(!profiles[k].count(signature))return false;
    }
    return true;
}
void enumerate(const vector<int>& alphabet,int position,int first,Word&w,vector<Word>&bank){
    if(position==7){raw_words++;if(valid(w))bank.push_back(w);return;}
    for(int i=first;i<(int)alphabet.size();i++){w[position]=alphabet[i];enumerate(alphabet,position+1,i,w,bank);}
}
using Table=array<array<int,91>,91>;
bool connected(const Word&w,const Table& counts){
    array<int,7> parent; iota(parent.begin(),parent.end(),0);
    auto find=[&](int a){while(parent[a]!=a)a=parent[a];return a;};
    for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)if(counts[w[i]][w[j]]>0)parent[find(i)]=find(j);
    for(int i=1;i<7;i++)if(find(i)!=find(0))return false;return true;
}
void writeword(ostream&out,const Word&w){for(int x:w)out<<' '<<x;}
bool connected_bits(const array<int,7>&adj,int remove_a=-1,int remove_b=-1){
    int reached=1,front=1;
    while(front){int next=0;for(int i=0;i<7;i++)if(front&(1<<i)){
        int row=adj[i];if(i==remove_a)row&=~(1<<remove_b);if(i==remove_b)row&=~(1<<remove_a);next|=row;
    }next&=~reached;reached|=next;front=next;}return reached==127;
}
int main(int argc,char**argv){
#ifdef _WIN32
    _setmode(_fileno(stdout),_O_BINARY);
#endif
    need(argc==4,"mode, independent input, output directory");
    string mode=argv[1],directory=argv[3];ifstream input(argv[2]);
    if(mode=="words"){
        for(int a=0;a<=90;a++)for(int b=0;b<=90;b++)common[a][b]=gcd(a,b);
        int np;input>>np;for(int j=0;j<np;j++){
            int k,c;input>>k>>c;uint64_t signature=c;for(int i=0;i<k;i++){int a;input>>a;signature=(signature<<7)|a;}profiles[k].insert(signature);
        }
        for(int mask=1;mask<127;mask++)inside_masks.push_back(mask);
        need(inside_masks.size()==126,"all proper subset profiles without pruning");
        int nd;input>>nd;map<int,vector<int>> alphabets;
        for(int i=0;i<nd;i++){int d,k;input>>d>>k;vector<int> row(k);for(int&a:row)input>>a;alphabets.emplace(d,row);}
        int nt;input>>nt;vector<pair<int,int>> clocks(nt);for(auto&[t,d]:clocks)input>>t>>d;
        int nq;input>>nq;map<tuple<int,int,int>,int> inherited;
        for(int i=0;i<nq;i++){int n,a,b,count;input>>n>>a>>b>>count;inherited[{n,a,b}]=count;}
        need(bool(input),"complete independent census input");
        map<int,vector<Word>> words;ofstream domain_file(directory+"/domains.txt",ios::binary);
        uint64_t accepted=0,residual=0,evaluated=0,single_eligible=0,cheap_eligible=0;
        for(auto&[d,alphabet]:alphabets){
            uint64_t before=raw_words;Word w{};enumerate(alphabet,0,0,w,words[d]);accepted+=words[d].size();
            domain_file<<d<<' '<<raw_words-before<<' '<<words[d].size()<<'\n';
            ofstream file(directory+"/words_"+to_string(d)+".json",ios::binary);file<<'[';bool first=true;
            for(const auto&w:words[d]){if(!first)file<<',';first=false;file<<'[';for(int i=0;i<7;i++){if(i)file<<',';file<<w[i];}file<<']';}file<<']';
        }
        ofstream residual_file(directory+"/residuals.txt",ios::binary),clock_file(directory+"/clocks.txt",ios::binary);
        for(auto[t,d]:clocks){
            Table counts{};for(int a:alphabets[d])for(int b:alphabets[d]){
                int e=gcd(a,b);counts[a][b]=inherited[{t/e,min(a,b)/e,max(a,b)/e}];
            }
            uint64_t rc=0,se=0,ce=0;
            for(const auto&w:words[d]){
                evaluated++;array<int,7> adj{};
                for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)if(counts[w[i]][w[j]]){adj[i]|=1<<j;adj[j]|=1<<i;}
                if(!connected_bits(adj))continue;
                map<int,int> mult;for(int a:w)mult[a]++;
                bool slots=true;for(auto[a,m]:mult)if(counts[a][a]==0){I bound=0;for(auto[b,n]:mult)if(a!=b)bound+=(I)n*counts[a][b];if(bound<m)slots=false;}
                if(!slots)continue;
                residual++;rc++;array<int,7> bridges{};
                for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)if((adj[i]&(1<<j))&&!connected_bits(adj,i,j)){bridges[i]|=1<<j;bridges[j]|=1<<i;}
                set<array<int,3>> single,cheap;
                for(int middle=0;middle<7;middle++)for(int i=0;i<7;i++)if(bridges[middle]&(1<<i))for(int j=i+1;j<7;j++)if(bridges[middle]&(1<<j)){
                    int a=w[i],b=w[middle],c=w[j];array<int,3> key{min(a,c),b,max(a,c)};
                    if(counts[a][b]==1||counts[c][b]==1)single.insert(key);
                    if((I)counts[a][b]*counts[c][b]<=256)cheap.insert(key);
                }
                se+=!single.empty();ce+=!cheap.empty();residual_file<<t;writeword(residual_file,w);
                residual_file<<' '<<single.size();for(auto k:single)for(int a:k)residual_file<<' '<<a;
                residual_file<<' '<<cheap.size();for(auto k:cheap)for(int a:k)residual_file<<' '<<a;residual_file<<'\n';
            }
            single_eligible+=se;cheap_eligible+=ce;clock_file<<t<<' '<<d<<' '<<words[d].size()<<' '<<rc<<' '<<se<<' '<<ce<<'\n';
        }
        cout<<"DOMAINS "<<nd<<" UNPRUNED "<<raw_words<<" ACCEPTED "<<accepted<<" WORD_CLOCKS "<<evaluated<<'\n';
        cout<<"CLOCKS "<<nt<<" RESIDUAL_WORDS "<<residual<<" SINGLE_ELIGIBLE "<<single_eligible<<" CHEAP_ELIGIBLE "<<cheap_eligible<<'\n';
    }else if(mode=="banks"){
        auto ratios=atlas();int nk;input>>nk;
        ofstream file(directory+"/banks.txt",ios::binary),meta(directory+"/bank_counts.txt",ios::binary);uint64_t typed=0,zero=0;
        for(int i=0;i<nk;i++){
            int t,a,b;input>>t>>a>>b;int e=gcd(a,b),n=t/e;uint64_t before_typed=typed,before_zero=zero;
            for(const auto&G:ratios){
                int da=e*gcd(n,G.p),db=e*gcd(n,G.q);
                if(!((da==a&&db==b)||(da==b&&db==a)))continue;
                auto cap=capacity(G,n);typed++;
                need(cap.ordinary>=0&&cap.chamber>=cap.ordinary,"native bank ordered capacities");
                if(cap.chamber)continue;
                auto actual=literal(G,n,cap.chamber_phase);
                need(actual[2]==0&&actual[0]==n/7&&actual[1]==n/7,"every independently regenerated chamber zero has literal attained marginals");
                if(da==a&&db==b){file<<t<<' '<<a<<' '<<b<<' '<<G.p<<' '<<G.q<<'\n';zero++;}
                if(db==a&&da==b){file<<t<<' '<<a<<' '<<b<<' '<<G.q<<' '<<G.p<<'\n';zero++;}
            }
            meta<<t<<' '<<a<<' '<<b<<' '<<typed-before_typed<<' '<<zero-before_zero<<'\n';
        }
        need(bool(input),"complete ratio-bank input");
        cout<<"BANKS "<<nk<<" TYPED_NO_SKIP "<<typed<<" ORIENTED_ZEROS "<<zero<<" LITERAL_POINTS "<<literal_points<<'\n';
    }else if(mode=="capacities"){
        controls();
        int nk;input>>nk;ofstream file(directory+"/capacities.txt",ios::binary);
        for(int i=0;i<nk;i++){
            int n,p,q;input>>n>>p>>q;need(p>0&&p<q&&gcd(p,q)==1,"arbitrary composite coprime ratio");
            auto G=geometry(p,q);auto cap=capacity(G,n);
            need(cap.ordinary>=0&&cap.chamber>=cap.ordinary,"endpoint capacities ordered");
            auto ordinary=literal(G,n,cap.ordinary_phase),chamber=literal(G,n,cap.chamber_phase);
            need(ordinary[2]==cap.ordinary,"literal all-phase minimizing endpoint grid");
            need(chamber[2]==cap.chamber,"literal chamber minimizing endpoint grid");
            file<<n<<' '<<p<<' '<<q<<' '<<cap.ordinary<<' '<<cap.chamber<<' '<<cap.ordinary_phase<<' '<<cap.chamber_phase<<'\n';
        }
        need(bool(input),"complete arbitrary product-capacity input");
        cout<<"COMPOSITE_CAPACITIES "<<nk<<" LITERAL_POINTS "<<literal_points<<'\n';
    }else need(false,"known independent consumer mode");
    cout<<"NATIVE_ALWAYS_ACTIVE_GATES "<<gates<<'\n';
}
