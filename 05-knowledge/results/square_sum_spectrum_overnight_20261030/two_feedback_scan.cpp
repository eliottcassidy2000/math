// Exact ordinary Hamilton-path spectrum of the two-feedback word grammar.
// Gap states are empty / valid-last / invalid-last; the final gap has two states.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
using Row = std::array<uint64_t,36>;
std::vector<std::pair<int,int>> arcs[4];
std::vector<int> accept;
std::vector<uint64_t> witness;
std::vector<int> length;
uint64_t nodes=0;
int depth_limit, limit;
int index(int b,int a,int c,int d){return b*18+a*6+c*2+d;}
void visit(const Row &row,int depth,uint64_t word){
    ++nodes;uint64_t h=0;for(int i:accept)h+=row[i];
    if(h<=uint64_t(limit)&&length[h]<0){length[h]=depth;witness[h]=word;}
    if(depth==depth_limit||h>uint64_t(limit))return;
    for(int letter=0;letter<4;++letter){
        // A leading 00 is a source. Omitting it preserves H and all later orientations.
        if(depth==0&&letter==0)continue;
        Row next{};for(auto [i,j]:arcs[letter])next[j]+=row[i];
        visit(next,depth+1,(word<<2)|letter);
    }
}
int main(int argc,char **argv){
    depth_limit=argc>1?std::atoi(argv[1]):10;
    limit=argc>2?std::atoi(argv[2]):624;
    if(depth_limit<0||depth_limit>20||limit<1)return 2;
    witness.resize(limit+1);length.assign(limit+1,-1);
    for(int b=0;b<2;++b){
        int p[2]={b,1-b};
        for(int s0=0;s0<3;++s0)for(int s1=0;s1<3;++s1)for(int s2=0;s2<2;++s2){
            int s[3]={s0,s1,s2},i=index(b,s0,s1,s2);
            if(s0!=2&&s1!=2&&(s1!=0||p[0]==0))accept.push_back(i);
            for(int a=0;a<4;++a)for(int g=0;g<3;++g){
                if(s[g]==0&&g>0&&!(a>>p[g-1]&1))continue;
                int t[3]={s0,s1,s2};
                t[g]=(g==2||!(a>>p[g]&1))?1:2;
                arcs[a].emplace_back(i,index(b,t[0],t[1],t[2]));
            }
        }
    }
    Row row{};row[0]=row[18]=1;visit(row,0,0);
    std::cout<<"{\"word_depth\":"<<depth_limit<<",\"limit\":"<<limit<<",\"visited_prefixes\":"<<nodes<<",\"missing_odds\":[";
    bool comma=false;for(int h=1;h<=limit;h+=2)if(length[h]<0){if(comma)std::cout<<',';comma=true;std::cout<<h;}
    std::cout<<"],\"witnesses\":[";comma=false;
    for(int h=1;h<=limit;h+=2)if(length[h]>=0){if(comma)std::cout<<',';comma=true;std::cout<<'['<<h<<','<<length[h]<<','<<witness[h]<<']';}
    std::cout<<"]}\n";
}
