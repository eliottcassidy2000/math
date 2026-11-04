#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using Rows = std::array<unsigned,12>;
using DP = std::array<unsigned,4096>;

DP endpoints(const Rows& a) {
    DP d{};
    for(unsigned m=1;m<4096;++m) {
        unsigned todo=m;
        while(todo) {
            int v=__builtin_ctz(todo); todo&=todo-1;
            unsigned rest=m^(1u<<v), pred=4095u^(1u<<v)^a[v];
            if(!rest || __builtin_parity(d[rest]&pred)) d[m]|=1u<<v;
        }
    }
    return d;
}

Rows reverse(const Rows& a) {
    Rows b{}; for(int u=0;u<12;++u)b[u]=4095u^(1u<<u)^a[u]; return b;
}

bool cover_gate(const Rows& a) {
    Rows odd{};
    for(int start=0;start<12;++start) {
        DP d{}; d[1u<<start]=1u<<start;
        for(unsigned m=1;m<4096;++m) {
            unsigned vs=d[m];
            if(__builtin_parity(m)) odd[start]^=vs;
            while(vs) {
                int u=__builtin_ctz(vs);vs&=vs-1;
                unsigned todo=a[u]&~m;
                while(todo) {unsigned bit=todo&-todo;todo-=bit;d[m|bit]^=bit;}
            }
        }
    }
    for(int u=0;u<12;++u)for(int v=u+1;v<12;++v)
        if(((odd[u]>>v)^(odd[v]>>u))&1)return false;
    return true;
}

struct Stats {
    unsigned total=0, endpoint=0, cubic=0, combined_old=0, combined_all=0;
    std::vector<Rows> old_survivors, cubic_survivors;
};

void inspect(const Rows& a, Stats& z) {
    for(int u=0;u<12;++u) {
        if(a[u]>>u&1)throw std::runtime_error("loop");
        for(int v=u+1;v<12;++v)
            if(((a[u]>>v)&1)+((a[v]>>u)&1)!=1)throw std::runtime_error("arc");
    }
    DP e=endpoints(a),s=endpoints(reverse(a));
    bool ep=e[4095]==s[4095];
    Rows columns{};
    for(unsigned m=1;m<4096;++m) {
        int size=__builtin_popcount(m);
        if((size&1)||size>8)continue;
        unsigned outside=4095^m;
        while(outside) {int v=__builtin_ctz(outside);outside&=outside-1;columns[v]^=s[m];}
    }
    bool cubic=true;for(unsigned col:columns)cubic&=col==0;
    ++z.total; z.endpoint+=ep; z.cubic+=cubic;
    if(cubic)z.cubic_survivors.push_back(a);
    if(ep && cover_gate(a)) {
        ++z.combined_old;z.combined_all+=cubic;z.old_survivors.push_back(a);
    }
}

void print_rows(const std::vector<Rows>& rows) {
    std::cout<<"[";for(size_t i=0;i<rows.size();++i){if(i)std::cout<<",";
    std::cout<<"[";for(int j=0;j<12;++j){if(j)std::cout<<",";std::cout<<rows[i][j];}std::cout<<"]";}std::cout<<"]";
}

int main(int argc,char**argv) {try {
    if(argc!=2)throw std::runtime_error("expected quotient bank input");
    std::cout<<"{\"families\":[";
    for(int family=0;family<4;++family) {
        Stats z;
        if(family<3) {
            Rows core{};
            for(int u=0;u<11;++u)for(int v=u+1;v<11;++v) {
                int d=(v-u+11)%11;
                bool edge=family==0 || (family==1 ? d<=5 : (d==1||d==3||d==4||d==5||d==9));
                if(edge)core[u]|=1u<<v;else core[v]|=1u<<u;
            }
            for(unsigned mask=0;mask<2048;++mask) {
                Rows a=core;
                for(int u=0;u<11;++u)if(mask>>u&1)a[u]|=1u<<11;else a[11]|=1u<<u;
                inspect(a,z);
            }
        } else {
            std::ifstream f(argv[1]);unsigned count;f>>count;
            if(count!=56)throw std::runtime_error("complete order6 quotient count");
            for(unsigned i=0;i<count;++i) {
                unsigned code;f>>code;if(!f)throw std::runtime_error("truncated quotient bank");
                Rows a{};for(int u=0;u<6;++u)a[2*u]|=1u<<(2*u+1);
                int bit=0;for(int u=0;u<6;++u)for(int v=u+1;v<6;++v,++bit) {
                    int from=(code>>bit&1)?u:v,to=(code>>bit&1)?v:u;
                    a[2*from]|=3u<<(2*to);a[2*from+1]|=3u<<(2*to);
                }
                inspect(a,z);
            }
            unsigned extra;if(f>>extra)throw std::runtime_error("quotient bank extra data");
        }
        if(family)std::cout<<",";
        std::cout<<"{\"family\":"<<family<<",\"named_presentations\":"<<z.total
                 <<",\"endpoint_parity_survivors\":"<<z.endpoint
                 <<",\"cubic_survivors\":"<<z.cubic
                 <<",\"old_combined_survivors\":"<<z.combined_old
                 <<",\"all_three_survivors\":"<<z.combined_all
                 <<",\"old_combined_rows\":";
        print_rows(z.old_survivors);std::cout<<",\"cubic_rows\":";print_rows(z.cubic_survivors);std::cout<<"}";
    }
    std::cout<<"]}\n";
} catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}}
