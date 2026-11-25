// All multilinear chain-inflation profiles certified by ordinary triangle covers.
// Input: nauty gentourng -cq upper triangles, one per strong isomorphism class.
#include <array>
#include <iostream>
#include <set>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>

using Rows = std::array<int,10>;
using Coeff = std::vector<long long>;

Rows decode(const std::string& line, int n) {
    if (line.size() != n*(n-1)/2) throw std::runtime_error("invalid record length");
    Rows rows{};
    int at=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) {
        if(line[at]=='1') rows[i]|=1<<j;
        else if(line[at]=='0') rows[j]|=1<<i;
        else throw std::runtime_error("invalid triangle");
        ++at;
    }
    return rows;
}

Rows guards(const Rows& rows, int n) {
    Rows result{};
    for(int i=0;i<n;i++) {
        result[i]=(1<<n)-1-(1<<i);
        for(int out=rows[i];out;out&=out-1) {
            int j=__builtin_ctz(out);
            for(int back=rows[j]&~rows[i]&~(1<<i);back;back&=back-1)
                result[i]&=(1<<j)|(back&-back);
        }
    }
    return result;
}

Coeff profile(const Rows& rows, int n, const std::vector<int>& ports) {
    Rows step{}, bit{};
    int length=1, size=1;
    for (int p:ports) { step[p]=length; length*=3; }
    for (int i=0; i<n; ++i) if (!step[i]) { bit[i]=size; size*=2; }
    std::vector<std::array<long long,10>> paths(length*size);
    for (int i=0; i<n; i++) paths[length*bit[i]+step[i]][i]=1;
    int full=(1<<ports.size())-1;
    Coeff f(full+1);
    for (int mask=0;mask<size;mask++) {
        int free=0;
        for (int j=0; j<n; ++j) if (bit[j] && !(bit[j]&mask)) free|=1<<j;
        for (int state=0; state<length; ++state) {
            int choices=free, repeat=0, read=0;
            for (int j=0; j<int(ports.size()); ++j) {
                int count=state/step[ports[j]]%3;
                if (count<2) choices|=1<<ports[j];
                if (count==2) repeat|=1<<j;
                if (count) read|=1<<j;
            }
            int at=mask*length+state;
            for (int last=0;last<n;last++) {
                auto num=paths[at][last];
                if (!num) continue;
                if (mask==size-1 && read==full) f[repeat]+=num;
                for (int s=rows[last]&choices;s;s&=s-1) {
                    int next=__builtin_ctz(s);
                    paths[at+length*bit[next]+step[next]][next]+=num;
                }
            }
        }
    }
    for (int j=0;j<int(ports.size());++j) for (int i=0;i<=full;++i)
        if (!(i>>j&1)) f[i]-=2*f[i|(1<<j)];
    return f;
}

int main(int argc, char **argv) {
    if (argc<2 || argc>5) return 2;
    int n=std::stoi(argv[1]), k=(argc>=3 ? std::stoi(argv[2]):3);
    if(n<3 || n>10 || k<1 || k>=n) return 2;
    std::vector<int> candidates;
    for(int mask=1;mask<(1<<n);mask++)
        if(__builtin_popcount(mask)==k) candidates.push_back(mask);
    std::set<Coeff> seen, wanted;
    bool reduced=false, select=false;
    std::string wanted_path;
    for(int i=3;i<argc;++i) {
        if(std::string(argv[i])=="--reduced") reduced=true;
        else { wanted_path=argv[i]; select=true; }
    }
    if(select) {
        std::ifstream input(wanted_path);
        if(!input) return 2;
        long long value;
        while(input>>value) {
            Coeff f(1<<k);
            f[0]=value;
            for(int i=1;i<(1<<k);i++) if(!(input>>f[i])) return 2;
            wanted.insert(f);
        }
    }
    std::string line;
    long long count=0,sets=0;
    while (std::cin>>line) {
        ++count;
        auto rows=decode(line,n), guard=guards(rows,n);
        int twin_sites=0;
        if(reduced) for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
            if(!((rows[i]^rows[j]) & ~((1<<i)|(1<<j)))) twin_sites|=(1<<i)|(1<<j);
        for (int mask:candidates) {
            if(mask&twin_sites) continue;
            std::vector<int> ports;
            for(int v=mask;v;v&=v-1) {
                int i=__builtin_ctz(v);
                if(!(guard[i]&~mask)) break;
                ports.push_back(i);
            }
            if(ports.size()!=k) continue;
            ++sets;
            auto f=profile(rows,n,ports);
            if (select ? !wanted.count(f) : !seen.insert(f).second) continue;
            std::cout<<line;
            for(auto i:ports) std::cout<<' '<<i;
            for(auto x:f) std::cout<<' '<<x;
            std::cout<<'\n';
        }
    }
    std::cerr<<n<<" graphs "<<count<<" eligible sets "<<sets<<" polynomials "<<seen.size()<<'\n';
}
