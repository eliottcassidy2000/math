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
using namespace std;
using I=long long;
using Word=array<int,7>;
uint64_t gates=0,raw_words=0,capacity_evaluations=0,zero_ratios=0,open_zero_ratios=0,literal_points=0;
void need(bool ok,const string& why){++gates;if(!ok){cerr<<"FAILED "<<why<<'\n';exit(1);}}
I mod(I a,I b){I r=a%b;return r<0?r+b:r;}
I floorq(I a,I b){return (a-mod(a,b))/b;}
I ceilq(I a,I b){return -floorq(-a,b);}
bool danger(I x,I denominator,int speed){I r=mod(speed*x,denominator);return 14*min(r,denominator-r)<denominator;}
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
int main(int argc,char**argv){
    need(argc==3,"independent input and output paths");ifstream input(argv[1]);string directory=argv[2];
    for(int a=0;a<=90;a++)for(int b=0;b<=90;b++)common[a][b]=gcd(a,b);
    int count;input>>count;
    for(int j=0;j<count;j++){int k,c;input>>k>>c;uint64_t signature=c;for(int i=0;i<k;i++){int x;input>>x;signature=(signature<<7)|x;}profiles[k].insert(signature);}
    for(int k=1;k<=6;k++)for(int mask=1;mask<127;mask++)if(__builtin_popcount((unsigned)mask)==k)inside_masks.push_back(mask);
    need(inside_masks.size()==126,"all proper selected-position profile tests");
    int nd;input>>nd;vector<vector<int>> alphabets(nd);
    for(auto&row:alphabets){int size;input>>size;row.resize(size);for(int&x:row)input>>x;}
    int nt;input>>nt;vector<pair<int,int>> clocks(nt);set<int> quotients,declared;
    for(auto&[t,d]:clocks){input>>t>>d;declared.insert(t);for(int a:alphabets[d])for(int b:alphabets[d]){need(t%(7*a)==0,"whole original margin alphabet zero budget");quotients.insert(t/gcd(a,b));}}
    need(bool(input)&&declared==quotients,"all actual quotient clocks covered by exact declared set");
    controls();auto ratios=atlas();
    map<int,pair<Table,Table>> banks;
    ofstream table_file(directory+"/zero_pairs.txt",ios::binary),critical_file(directory+"/critical_banks.txt",ios::binary),all_zero_file(directory+"/all_zero_ratios.txt",ios::binary);
    set<int> critical={1939,1953};for(int a:{4,8,9,16,18,24})for(int b:{4,8,9,16,18,24})critical.insert(7056/gcd(a,b));
    uint64_t examined=0;
    for(int n:quotients){
        Table ordinary{},chamber{};
        for(const auto&G:ratios){
            int a=gcd(n,G.p),b=gcd(n,G.q);examined++;
            if(a>90||b>90||a%7==0||b%7==0)continue;
            // Deliberately no lower-bound skip and no first-witness/ratio exit.
            auto value=capacity(G,n);capacity_evaluations++;
            need(value.ordinary>=0&&value.chamber>=value.ordinary,"unpruned native capacities ordered");
            int copies=a==b?2:1,lo=min(a,b),hi=max(a,b);
            if(value.ordinary==0){
                ordinary[lo][hi]+=copies;zero_ratios++;
                need(literal(G,n,value.ordinary_phase)[2]==0,"every ordinary zero has literal full-grid witness");
            }
            if(value.chamber==0){
                chamber[lo][hi]+=copies;open_zero_ratios++;
                auto physical=literal(G,n,value.chamber_phase);
                need(physical[2]==0&&physical[0]==n/7&&physical[1]==n/7,"every chamber zero attains both actual marginal bounds");
            }
            int bits=(value.ordinary==0?1:0)|(value.chamber==0?2:0);
            if(bits)all_zero_file<<n<<' '<<G.p<<' '<<G.q<<' '<<a<<' '<<b<<' '<<bits<<'\n';
            if(bits&&critical.count(n))critical_file<<n<<' '<<G.p<<' '<<G.q<<' '<<a<<' '<<b<<' '<<bits<<'\n';
        }
        for(int a=1;a<=90;a++)for(int b=a;b<=90;b++){
            if(ordinary[a][b])table_file<<n<<' '<<a<<' '<<b<<' '<<(1|(chamber[a][b]?2:0))<<' '<<ordinary[a][b]<<' '<<chamber[a][b]<<'\n';
            ordinary[b][a]=ordinary[a][b];chamber[b][a]=chamber[a][b];
        }
        banks.emplace(n,make_pair(ordinary,chamber));
    }
    cerr<<"UNPRUNED_NATIVE ratios="<<examined<<" typed="<<capacity_evaluations<<" zeros="<<zero_ratios<<'\n';
    vector<vector<Word>> words(nd);ofstream domain_file(directory+"/domains.txt",ios::binary);
    uint64_t valid_words=0;
    for(int d=0;d<nd;d++){
        uint64_t before=raw_words;Word w{};enumerate(alphabets[d],0,0,w,words[d]);valid_words+=words[d].size();
        domain_file<<d<<' '<<words[d].size()<<' '<<raw_words-before<<'\n';
        ofstream out(directory+"/words_"+to_string(d)+".json",ios::binary);out<<'[';bool first=true;
        for(const auto&row:words[d]){if(!first)out<<',';first=false;out<<'[';for(int j=0;j<7;j++){if(j)out<<',';out<<row[j];}out<<']';}out<<']';
    }
    cerr<<"UNPRUNED_DOMAINS raw="<<raw_words<<" valid="<<valid_words<<'\n';
    ofstream survey_file(directory+"/survey.txt",ios::binary),stop_file(directory+"/stop.txt",ios::binary);
    I evaluated_words=0,closed0=0,closed1=0,closed2=0;
    for(auto[t,d]:clocks){
        Table ordinary{},chamber{};
        for(int a:alphabets[d])for(int b:alphabets[d]){int e=gcd(a,b);ordinary[a][b]=banks.at(t/e).first[a/e][b/e];chamber[a][b]=banks.at(t/e).second[a/e][b/e];}
        I n0=0,n1=0,n2=0;Word first0{},first1{},first2{};
        for(const auto&w:words[d]){
            evaluated_words++;
            if(connected(w,ordinary)){if(!n0)first0=w;n0++;}
            if(!connected(w,chamber))continue;
            if(!n1)first1=w;n1++;
            map<int,int> multiplicity;for(int a:w)multiplicity[a]++;
            bool possible=true;
            for(auto[a,needed]:multiplicity)if(chamber[a][a]==0){
                I slots=0;for(auto[b,number]:multiplicity)if(a!=b)slots+=(I)number*chamber[a][b];
                if(slots<needed)possible=false;
            }
            if(possible){if(!n2)first2=w;n2++;if(t==7056){for(int a:w)stop_file<<a<<' ';stop_file<<'\n';}}
        }
        need(n2<=n1&&n1<=n0,"all independent clock-stage residuals nested");
        closed0+=n0==0;closed1+=n1==0;closed2+=n2==0;
        survey_file<<t<<' '<<d<<' '<<words[d].size()<<' '<<n0<<' '<<n1<<' '<<n2;
        writeword(survey_file,first0);writeword(survey_file,first1);writeword(survey_file,first2);survey_file<<'\n';
    }
    cout<<"ATLAS 5855; QUOTIENTS "<<quotients.size()<<"; DOMAINS "<<nd<<'\n';
    cout<<"NATIVE_RATIOS_EXAMINED "<<examined<<"; TYPED_CAPACITIES_NO_SKIPS "<<capacity_evaluations<<'\n';
    cout<<"ORDINARY_ZERO_RATIOS "<<zero_ratios<<"; CHAMBER_ZERO_RATIOS "<<open_zero_ratios<<"; LITERAL_GRID_POINTS "<<literal_points<<'\n';
    cout<<"UNPRUNED_MULTISETS "<<raw_words<<"; DISTINCT_VALID_WORDS "<<valid_words<<"; WORD_CLOCKS "<<evaluated_words<<'\n';
    cout<<"CUMULATIVE_CLOSURES "<<closed0<<' '<<closed1<<' '<<closed2<<'\n';
    cout<<"NATIVE_ALWAYS_ACTIVE_GATES "<<gates<<'\n';
}
