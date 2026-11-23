// Read canonical digraph6 tournaments from standard input, one per line.
// "extend" emits every strong one-vertex extension of path count at most 699.
// "read" emits the endpoint and one-flip data for divisors of the target heights.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Count = std::int64_t;
using Rows = std::vector<int>;
constexpr int max_order = 13;
using Endpoints = std::array<Count, max_order>;

int first(int bits) { return __builtin_ctz(static_cast<unsigned>(bits)); }

Rows decode(const std::string& line) {
    int n = line.at(1) - 63;
    Rows rows(n);
    for (int i=0, k=0; i<n; ++i) {
        for (int j=0; j<n; ++j, ++k) {
            if (((line.at(2+k/6)-63) >> (5-k%6)) & 1) rows[i] |= 1<<j;
        }
    }
    return rows;
}

std::string encode(const Rows& rows) {
    int n=rows.size(), bits=0, len=0;
    std::string out="&";
    out.push_back(63+n);
    for (int row : rows) for (int j=0; j<n; ++j) {
        bits=2*bits+((row>>j)&1);
        if (++len==6) {
            out.push_back(63+bits);
            bits=len=0;
        }
    }
    if (len) out.push_back(63+(bits<<(6-len)));
    return out;
}

void process(const Rows& rows, bool extend) {
    int n=rows.size(), size=1<<n, full=size-1;
    std::vector<Endpoints> start(size), end(size);
    std::vector<Count> h(size);
    h[0]=1;
    for (int m=1; m<size; ++m) for (int bits=m; bits; bits&=bits-1) {
        int i=first(bits), rest=m^(1<<i);
        if (!rest) start[m][i]=end[m][i]=1;
        for (int b=rest&rows[i]; b; b&=b-1) start[m][i]+=start[rest][first(b)];
        for (int b=rest&~rows[i]; b; b&=b-1) end[m][i]+=end[rest][first(b)];
        h[m]+=start[m][i];
    }
    if (!extend) {
        bool needed=false;
        for (int target : {99,101,199,201,299,301,399,401,499,501,599,601,699})
            if (target%h[full]==0) needed=true;
        if (!needed) return;
    }
    Count q=0;
    Endpoints starts{}, ends{};
    std::array<Endpoints,max_order> joined{};
    for (int m=0; m<size; ++m) {
        int rest=full^m;
        q+=h[m]*h[rest];
        for (int bits=m; bits; bits&=bits-1) {
            int i=first(bits);
            starts[i]+=start[m][i]*h[rest];
            ends[i]+=end[m][i]*h[rest];
            if (end[m][i]) for (int b=rest; b; b&=b-1) {
                int j=first(b);
                joined[i][j]+=end[m][i]*start[rest][j];
            }
        }
    }
    if (!extend) {
        std::cout<<encode(rows)<<" "<<h[full]<<" "<<q;
        for (const Endpoints* data : {&start[full], &end[full], &starts, &ends})
            for (int i=0; i<n; ++i) std::cout<<" "<<(*data)[i];
        for (int i=0; i<n; ++i) for (int j=0; j<i; ++j) {
            int direction=(rows[j]>>i&1) ? 1 : -1;
            std::cout<<" "<<h[full]+direction*(joined[i][j]-joined[j][i]);
        }
        std::cout<<"\n";
        return;
    }
    // Splitting an apex path before and after the apex gives a quadratic
    // polynomial in its outgoing contact bits. Update it one added bit at a time.
    Endpoints linear{};
    for (int i=0; i<n; ++i) {
        linear[i]=start[full][i]-end[full][i];
        for (int j=0; j<n; ++j) linear[i]+=joined[j][i];
    }
    std::vector<Count> apex(size);
    apex[0]=h[full];
    // The old tournament is strong. A mixed contact set is necessary and
    // sufficient for its extension to stay strong.
    for (int m=1; m<full; ++m) {
        int i=first(m), rest=m&(m-1);
        apex[m]=apex[rest]+linear[i];
        for (int b=rest; b; b&=b-1) {
            int j=first(b);
            apex[m]-=joined[i][j]+joined[j][i];
        }
        if (apex[m]>699) continue;
        Rows grown=rows;
        for (int j=0; j<n; ++j) if (!(m&(1<<j))) grown[j]|=1<<n;
        grown.push_back(m);
        std::cout<<encode(grown)<<"\n";
    }
}

int main(int argc,char** argv) {
    if (argc!=2 || (std::string(argv[1])!="extend" && std::string(argv[1])!="read")) return 2;
    std::string line;
    while (std::getline(std::cin,line)) {
        if (line.empty() || line[0]!='&') return 2;
        process(decode(line), std::string(argv[1])=="extend");
    }
}
