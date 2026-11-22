// All positive independent transitive widths on a marked 27->29 seed. C++17.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Row = std::vector<int>;
using Small = std::uint8_t;

struct Spectrum {
    int prime, d, size;
    Row powers;
    std::vector<Small> lower,upper;
    Spectrum(int p, Row rows, int a, int b):prime(p),d(rows.size()-2),size(1) {
        for(int i=0;i<d;++i) { powers.push_back(size);size*=p; }
        lower=values(rows,a,b);
        rows[a]^=1<<b; rows[b]^=1<<a;
        upper=values(rows,a,b);
    }

    // Each state is a port multiplicity vector and a two-bit mask of ordinary vertices.
    std::vector<Small> values(const Row& rows,int a,int b) const {
        int n=rows.size(),p=prime;
        Row port;
        for(int i=0;i<n;++i) if(i!=a&&i!=b) port.push_back(i);
        Row pred(n);
        for(int v=0;v<n;++v) for(int j=0;j<n;++j) if(rows[j]>>v&1) pred[v]|=1<<j;
        std::vector<Small> dp(size*4*n),spectrum(size);
        Row multiplicities(d);
        for(int code=0;code<size;++code) {
            int digits=code,used=0;
            for(int j=0;j<d;++j) {
                multiplicities[j]=digits%p;digits/=p;used+=multiplicities[j];
            }
            for(int mask=0;mask<4;++mask) {
                int length=used+(mask&1)+(mask>>1);
                auto write=[&](int v,int previous) {
                    int sum=length==1;
                    for(int bits=pred[v];bits;bits&=bits-1) sum+=dp[previous*n+__builtin_ctz(bits)];
                    dp[(code*4+mask)*n+v]=sum%p;
                };
                for(int j=0;j<d;++j) if(multiplicities[j]) write(port[j],(code-powers[j])*4+mask);
                if(mask&1) write(a,code*4+mask-1);
                if(mask&2) write(b,code*4+mask-2);
            }
            if(std::find(multiplicities.begin(),multiplicities.end(),0)==multiplicities.end()) {
                int at=(code*4+3)*n;
                spectrum[code]=std::accumulate(dp.begin()+at,dp.begin()+at+n,0)%p;
            }
        }
        std::vector<Row> sur(p,Row(p));
        sur[0][0]=1;
        for(int w=1;w<p;++w)
            for(int r=1;r<p;++r) sur[w][r]=r*(sur[w-1][r]+sur[w-1][r-1])%p;
        // Apply a one-dimensional change of basis along each tensor coordinate.
        for(int stride:powers) {
            std::vector<Small> next(size);
            for(int block=0;block<size;block+=stride*p) for(int off=0;off<stride;++off) {
                int base=block+off;
                for(int w=1;w<p;++w) {
                    int result=0;
                    for(int r=1;r<=w;++r) result+=sur[w][r]*spectrum[base+r*stride];
                    next[base+(w%(p-1))*stride]=result%p;
                }
            }
            spectrum.swap(next);
        }
        return spectrum;
    }

    std::pair<int,int> lookup(const Row& widths) const {
        int index=0;
        for(int j=0;j<d;++j) index+=(widths[j]%(prime-1))*powers[j];
        return {lower[index],upper[index]};
    }
};

int literal(const Row& rows,const Row& widths,int a,int b) {
    Row names;
    int next=0;
    for(int v=0;v<int(rows.size());++v) {
        int w=(v==a||v==b)?1:widths[next++];
        names.insert(names.end(),w,v);
    }
    int n=names.size(),full=(1<<n)-1;
    Row in(n);
    for(int i=0;i<n;++i) for(int j=0;j<n;++j)
        if((names[i]==names[j] && j<i) || (rows[names[j]]>>names[i]&1)) in[i]|=1<<j;
    Row dp((full+1)*n);
    for(int mask=1;mask<=full;++mask) for(int x=mask;x;x&=x-1) {
        int v=__builtin_ctz(x),prev=mask^(1<<v),num=!prev;
        for(int y=prev&in[v];y;y&=y-1) num+=dp[prev*n+__builtin_ctz(y)];
        dp[mask*n+v]=num;
    }
    return std::accumulate(dp.begin()+full*n,dp.end(),0);
}

int e61_flags(const Spectrum& table,const Row& widths) {
    int p=table.prime, u=1;
    for(int i=0;i<60;++i) u=2*u%p;
    int c=u*(2*u-1)%p;
    const auto [l,h]=table.lookup(widths);
    return int(l==(c+p-1)%p && h==(c+1)%p) + 2*int(h==(c+p-1)%p && l==(c+1)%p);
}

void lifts(Row& w, int i, int current, std::vector<Row>& out) {
    if(i==int(w.size())) { out.push_back(w);return; }
    int base=w[i];
    for(;current<=155;w[i]+=12,current+=12) lifts(w,i+1,current,out);
    w[i]=base;
}

int bank_search() {
    int n=0,bank=0,retained=0;
    std::int64_t boxes=0,total_lifts=0,left=0;
    while(std::cin>>n) {
        if(n<3||n>7) return 2;
        Row rows(n);
        for(int& r:rows) if(!(std::cin>>r)) return 2;
        ++bank;
        std::vector<Spectrum> tables;
        for(int p:{3,5,7,13}) tables.emplace_back(p,rows,0,1);
        int count=1,d=n-2;
        for(int i=0;i<d;++i) count*=12;
        std::vector<Row> possible;
        for(int code=0;code<count;++code) {
            int tmp=code, flags=3;Row w(d);
            for(int& x:w) {x=tmp%12+1;tmp/=12;}
            for(auto& tab:tables) { flags&=e61_flags(tab,w);if(!flags) break; }
            if(flags) {
                ++boxes;
                lifts(w,0,std::accumulate(w.begin(),w.end(),0),possible);
            }
        }
        if(possible.empty()) continue;
        ++retained;total_lifts+=possible.size();
        for(int q:{11,17,19,23}) {
            tables.emplace_back(q,rows,0,1);
            possible.erase(std::remove_if(possible.begin(),possible.end(),[&](const Row& w) {
                int f=3;
                for(const auto& tab:tables) f&=e61_flags(tab,w);
                return !f;
            }),possible.end());
            if(possible.empty()) break;
        }
        left+=possible.size();
        for(auto& w:possible) {
            std::cout << "unresolved rows";
            for(int x:rows) std::cout<<' '<<x;
            std::cout<<" widths"; for(int x:w) std::cout<<' '<<x;
            std::cout<<'\n';
        }
    }
    std::cout<<"quotients "<<bank<<", with passing phase "<<retained<<", boxes "<<boxes
             <<", width vectors within bound "<<total_lifts<<", unresolved "<<left<<'\n';
    return 0;
}

int main(int argc,char** argv) {
    if(argc==2 && std::string(argv[1])=="--bank") return bank_search();
    if(argc<5) return 2;
    int a=std::stoi(argv[1]),b=std::stoi(argv[2]);
    Row rows;
    for(int i=3;i<argc;++i) rows.push_back(std::stoi(argv[i]));
    Row primes={3,5,7,13};
    std::vector<Spectrum> tables;
    for(int p:primes) tables.emplace_back(p,rows,a,b);
    int d=rows.size()-2, boxes=1;
    for(int i=0;i<d;++i) boxes*=12;
    std::array<std::array<int,16>,12> counts{};
    Row other=rows; other[a]^=1<<b; other[b]^=1<<a;
    for(int code=0;code<boxes;++code) {
        Row w(d); int digits=code;
        for(int &v:w) {v=1+digits%12;digits/=12;}
        bool control=std::accumulate(w.begin(),w.end(),0)<=d+3;
        int l=control?literal(rows,w,a,b):0, h=control?literal(other,w,a,b):0;
        std::array<int,12> forward{},backward{};
        for(int i=0;i<4;++i) {
            const auto [x,y]=tables[i].lookup(w);
            int p=primes[i];
            if(control && (l%p!=x || h%p!=y)) throw std::logic_error("literal control failed");
            int power=1;
            for(int r=0;r<12;++r) {
                int center=(power*(2*power-1))%p;
                power=2*power%p;
                if(x==(center+p-1)%p && y==(center+1)%p) forward[r]|=1<<i;
                if(y==(center+p-1)%p && x==(center+1)%p) backward[r]|=1<<i;
            }
        }
        for(int r=0;r<12;++r)
            for(int mask=1;mask<16;++mask)
                counts[r][mask]+=((forward[r]&mask)==mask || (backward[r]&mask)==mask);
    }
    std::cout << "seed on " << rows.size() << " vertices, arc " << a << ',' << b << ", " << boxes << " width boxes\n";
    for(int r=0;r<12;++r) {
        std::cout << "exponent " << r+1 << " mod 12: " << counts[r][15] << " boxes; minimal rejecting masks";
        for(int mask=1;mask<16;++mask) {
            if(counts[r][mask]) continue;
            bool minimal=true;
            for(int i=0;i<4;++i) if((mask>>i&1) && (mask^(1<<i)) && !counts[r][mask^(1<<i)]) minimal=false;
            if(minimal) std::cout << ' ' << mask;
        }
        std::cout << '\n';
    }
}
