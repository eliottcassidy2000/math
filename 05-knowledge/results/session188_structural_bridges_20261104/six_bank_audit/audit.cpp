#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
static void need(bool p,const char *s){if(!p)throw runtime_error(s);}
static int edge_index(int i,int j){if(i>j)swap(i,j);return i*(11-i)/2+j-i-1;}
static int hp(const array<unsigned,8>&r){
    int d[256][8]={}; for(int i=0;i<8;i++)d[1<<i][i]=1;
    // Backward-ending recurrence, distinct from forward path extension.
    for(unsigned mask=1;mask<256;mask++)for(int v=0;v<8;v++)if(mask&(1u<<v)){
        unsigned rest=mask^(1u<<v); if(!rest)continue;
        int z=0;for(int u=0;u<8;u++)if((rest&(1u<<u))&&(r[u]&(1u<<v)))z+=d[rest][u];
        d[mask][v]=z;
    }
    return accumulate(d[255],d[255]+8,0);
}
static array<unsigned,8> graph(int code,int a,int c){
    array<unsigned,8> r{};for(int i=0;i<6;i++)for(int j=i+1;j<6;j++){
        if(code&(1<<edge_index(i,j)))r[i]|=1u<<j;else r[j]|=1u<<i;
    }
    r[6]|=1u<<7;for(int u=0;u<6;u++)for(int p=0;p<2;p++){
        if((p?c:a)&(1<<u))r[u]|=1u<<(p+6);else r[p+6]|=1u<<u;
    }return r;
}
int main(int argc,char**argv){
    need(argc==2,"supply binary bank path");const int N=32768,K=484;
    vector<unsigned char> saved(N*K);ifstream f(argv[1],ios::binary);
    f.read(reinterpret_cast<char*>(saved.data()),saved.size());need(bool(f)&&f.peek()==EOF,"exact raw byte length");
    for(auto x:saved)need(x<4,"canonical residue byte");
    vector<int> basis,index(4096,-1);
    // Enumerate the independently stated 22-by-22 named basis directly.
    vector<int> masks={0};for(int u=0;u<6;u++)masks.push_back(1<<u);
    for(int u=0;u<6;u++)for(int v=u+1;v<6;v++)masks.push_back((1<<u)|(1<<v));
    sort(masks.begin(),masks.end());need(masks.size()==22,"single-port basis count");
    for(int c:masks)for(int a:masks){index[a|(c<<6)]=basis.size();basis.push_back(a|(c<<6));}
    need(basis.size()==K&&is_sorted(basis.begin(),basis.end()),"named coordinate order");
    vector<unsigned char> rebuilt(N*K);array<int,8>w; iota(w.begin(),w.end(),0);
    uint64_t words=0,compatible_words=0,updates=0;
    do{
        unsigned old_one=0,old_zero=0,port_one=0,port_zero=0;bool reject=false;
        for(int t=1;t<8;t++){
            int u=w[t-1],v=w[t];
            if(u<6&&v<6){unsigned e=1u<<edge_index(u,v);if(u<v)old_one|=e;else old_zero|=e;}
            else if(u>=6&&v>=6){if(u==7)reject=true;}
            else if(u<6)port_one|=1u<<(6*(v-6)+u);
            else port_zero|=1u<<(6*(u-6)+v);
        }
        if(reject)continue;
        words++;
        need(!(old_one&old_zero)&&!(port_one&port_zero),"literal word conflicting incidence");
        unsigned free=32767u^(old_one|old_zero);
        compatible_words+=uint64_t(1)<<__builtin_popcount(free);
        // Direct physical old-edge completions: no old-edge monomial expansion
        // and no zeta transform. Port coefficients are their literal finite differences.
        for(unsigned s=port_zero;;s=(s-1)&port_zero){
            int k=index[port_one|s];need(k>=0,"port bidegree bound");
            unsigned char add=__builtin_popcount(s)%2?3:1;
            auto *row=rebuilt.data()+size_t(k)*N;
            for(unsigned z=free;;z=(z-1)&free){auto &x=row[old_one|z];x=(x+add)&3;updates++;if(!z)break;}
            if(!s)break;
        }
    }while(next_permutation(w.begin(),w.end()));
    need(rebuilt==saved,"complete direct-completion bank differs");
    set<vector<unsigned char>> vectors;
    for(int code=0;code<N;code++){vector<unsigned char> v(K);for(int k=0;k<K;k++)v[k]=saved[size_t(k)*N+code];need(vectors.insert(v).second,"duplicate complete named vector");}
    uint64_t interpolation_queries=0,interpolation_coefficients=0,full_queries=0;
    // Independent controls differ from the author's chosen source list.
    set<int> controls={0,32767,408,15462,17,12345};for(int j=0;j<24;j++)controls.insert((j*4051+1013)&32767);
    for(int code:controls){vector<int> values(4096);
        for(int m:basis){values[m]=hp(graph(code,m&63,m>>6));interpolation_queries++;}
        for(int k=0;k<K;k++){unsigned m=basis[k];int z=0;
            for(unsigned s=m;;s=(s-1)&m){z+=(__builtin_popcount(m^s)%2?-1:1)*values[s];if(!s)break;}
            need(((z%4)+4)%4==saved[size_t(k)*N+code],"independent literal interpolation");interpolation_coefficients++;
        }
    }
    for(int code:{408,15462})for(int m=0;m<4096;m++){
        int want=hp(graph(code,m&63,m>>6)),got=0;
        for(int k=0;k<K;k++)if((basis[k]&m)==basis[k])got+=saved[size_t(k)*N+code];
        need((want&3)==(got&3),"complete literal truth table");full_queries++;
    }
    need(hp(graph(408,6,1))==455&&hp(graph(15462,6,1))==477,"literal named higher-degree witness");
    // Independent full ordered-cover counts on every old tournament. In each
    // permutation every backwards adjacency forces a cut, all other cuts optional.
    uint64_t cover_profiles=0,cover_coordinates=0;array<uint64_t,7> cover_sums{};
    int choose[6][6]={};for(int i=0;i<=5;i++){choose[i][0]=choose[i][i]=1;for(int j=1;j<i;j++)choose[i][j]=choose[i-1][j-1]+choose[i-1][j];}
    for(int code=0;code<N;code++){
        array<int,6> p; iota(p.begin(),p.end(),0);array<int,7> q{};
        do{int backwards=0;for(int i=1;i<6;i++){int u=p[i-1],v=p[i];bool forward=bool(code&(1<<edge_index(u,v)))==(u<v);backwards+=!forward;}
            for(int k=backwards+1;k<=6;k++)q[k]+=choose[5-backwards][k-1-backwards];
        }while(next_permutation(p.begin(),p.end()));
        need(q[1]>0&&(q[1]&1),"old Hamilton count odd positive");
        need((q[1]&3)==saved[code],"independent q1 equals all-zero port query");
        need(q[6]==720&&q[5]==1800,"independent universal high-cover counts");
        for(int k=1;k<=6;k++){cover_sums[k]+=q[k];cover_coordinates++;}cover_profiles++;
    }
    cout<<"{\"status\":\"PASS\",\"binary_bytes\":"<<saved.size()<<",\"all_compared_residues\":"<<saved.size()<<",\"literal_words\":"<<words<<",\"compatible_word_tournament_pairs\":"<<compatible_words<<",\"direct_residue_updates\":"<<updates<<",\"exact_distinct_full_vectors\":"<<vectors.size()<<",\"independent_interpolation_sources\":"<<controls.size()<<",\"independent_interpolation_queries\":"<<interpolation_queries<<",\"independent_interpolated_coefficients\":"<<interpolation_coefficients<<",\"complete_truth_table_queries\":"<<full_queries<<",\"ordered_cover_profiles\":"<<cover_profiles<<",\"ordered_cover_coordinates\":"<<cover_coordinates<<",\"cover_coordinate_sums\":[";
    for(int k=1;k<=6;k++){if(k>1)cout<<",";cout<<cover_sums[k];}cout<<"]}\n";
}
