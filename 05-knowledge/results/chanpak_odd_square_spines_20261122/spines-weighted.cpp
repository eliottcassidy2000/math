// Exact search on |a*x-b*y| <= width in the two canonical E61 charts.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Int = __int128_t;
using Coord = int64_t;
using Pair = std::pair<Coord, Coord>;
constexpr Int target61 = (Int(1) << 121) - (Int(1) << 60) - 1;
struct Chart { Coord h, q, c, t, e, r, k; };
constexpr Chart four{2465,2293760,1344270,712944,1740,986,423385};
constexpr Chart eight{37265,144834560,84839910,57602544,24360,14906,33838015};

Coord mod(Int value, Coord period) {
    Coord result = value % period;
    return result < 0 ? result + period : result;
}

Coord inverse(Coord value, Coord modulus) {
    Coord r=value, s=modulus, x=1, y=0;
    while (s) {
        Coord q = r/s;
        Coord nr=r-q*s, nx=x-q*y;
        r=s; s=nr; x=y; y=nx;
    }
    assert(r==1);
    return mod(x,modulus);
}

Int read_word(Coord x) {
    Int total=0;
    while (x) {
        Coord half=x/2;
        if ((x&1)==0) total+=half+(half&1);
        x=half;
    }
    return total;
}

Int height(const Chart& ch, Coord x, Coord y) {
    return ch.k+Int(ch.q)*x*y+Int(ch.t)*x+Int(ch.c)*y+
           ch.e*(y&1)+ch.r*(x&1)+ch.h*(read_word(x)+read_word(y));
}

// Largest n >= 0 with q*n^2+v*n <= bound, or -1 if there is none.
// |v| <= q makes the polynomial nondecreasing on integer n >= 0.
Coord root(Int bound, Coord q=1, Coord v=0) {
    if (bound<0) return -1;
    long double discriminant=(long double)v*v+4.0L*q*(long double)bound;
    Coord n=static_cast<Coord>((std::sqrt(discriminant)-v)/(2.0L*q));
    auto value=[&](Coord j) {return (Int(q)*j+v)*j;};
    while (value(n)>bound) --n;
    while (value(n+1)<=bound) ++n;
    return n;
}

// floor(4ab * value / q), avoiding the larger unscaled product.
Int scaled(Int value, Coord a, Coord b, Coord q) {
    assert(value>=0);
    Coord scale=4*a*b;
    return (value/q)*scale+((value%q)*scale)/q;
}

struct Stats {
    uint64_t sums=0, occupied=0, reads=0;
    std::set<Pair> hits;
};

void branch(const Chart& ch, Int target, Coord a, Coord b, Coord width,
            Coord i, Coord j, Stats& stats) {
    const Coord h=ch.h, Q=4*ch.q;
    const Coord T=2*(ch.t+ch.q*j), C=2*(ch.c+ch.q*i);
    const Coord K=ch.k+ch.q*i*j+ch.t*i+ch.c*j+ch.r*i+ch.e*j;
    if (target<K) return;
    // Physical x=2X+i and y=2Y+j. S=aX+bY.
    const Coord delta=a*i-b*j, bound=(width+std::abs(delta))/2;
    const Coord upper=root(scaled(target-K,a,b,Q)+Int(bound)*bound);
    const Coord linear=std::max((T+2*h+a-1)/a,(C+2*h+b-1)/b);
    Int minimum=target-K-h*(i+j)-Int(linear)*upper;
    Coord lower=(minimum<=0) ? 0 : root(scaled(minimum,a,b,Q));
    // A downward rounding at the lower endpoint can only add a slice.
    Coord inv_a=inverse(a%b,b), inv_h=inverse(h%b,b);

    // For every X mod h, solve a single LINEAR congruence for S mod h.
    // Multiplication by b is invertible mod h for the selected directions.
    std::vector<std::vector<Coord>> allowed(h);
    for (Coord x=0;x<h;++x) {
        Coord slope=mod(Int(Q)*x+C,h);
        Coord rhs=mod(Int(a)*Q*x*x-(b*T-a*C)*Int(x)-b*(K-target),h);
        Coord g=std::gcd(slope,h);
        if (rhs%g) continue;
        Coord period=h/g;
        Coord s0=period==1 ? 0 : mod(Int(rhs/g)*inverse(slope/g,period),period);
        for (Coord s=s0;s<h;s+=period) allowed[s].push_back(x);
    }

    const Coord curvature=Q*a*b;
    for (Coord s=lower;s<=upper;++s) {
        ++stats.sums;
        const auto& residues=allowed[s%h];
        if (residues.empty()) continue;
        ++stats.occupied;
        Coord residue=mod(inv_a*(s%b),b);
        Int vertex=Int(Q)*s+b*T-a*C;
        // First lattice point closest to the continuous vertex; either tie works.
        Int shift=vertex-Int(2*a)*Q*residue+curvature;
        Coord offset=shift>=0 ? shift/(2*curvature) : -((-shift+2*curvature-1)/(2*curvature));
        Coord mid_x=residue+b*offset, mid_y=(s-a*mid_x)/b;
        Coord beta=static_cast<Coord>(Int(b)*T-Int(a)*C+
                                      Int(Q)*(Int(b)*mid_y-Int(a)*mid_x));
        assert(std::abs(beta)<=curvature);
        Int excess=K+Int(Q)*mid_x*mid_y+Int(T)*mid_x+Int(C)*mid_y-target;
        Int allowance=Int(h)*(2*(s/std::min(a,b))+i+j);

        for (Coord sign : {-1,1}) {
            Coord v=-sign*beta;
            Coord low=std::max(Coord(sign==1),root(excess-1,curvature,v)+1);
            Coord high=root(excess+allowance,curvature,v);
            if (low>high) continue;
            Coord left=sign>0 ? mid_x+b*low : mid_x-b*high;
            Coord right=sign>0 ? mid_x+b*high : mid_x-b*low;
            left=std::max(Coord(0),left);
            right=std::min(s/a,right);
            for (Coord xr : residues) {
                Coord initial=xr+h*mod(Int(residue-xr)*inv_h,b);
                Coord first=left+mod(initial-left,h*b);
                for (Coord X=first;X<=right;X+=h*b) {
                    Coord Y=(s-a*X)/b;
                    Coord x=2*X+i, y=2*Y+j;
                    if (std::abs(a*x-b*y)>width) continue;
                    ++stats.reads;
                    if (height(ch,x,y)==target) stats.hits.emplace(x,y);
                }
            }
        }
    }
}

Stats scan(const Chart& ch, Int target, Coord a, Coord b, Coord width) {
    Stats stats;
    for (Coord i : {0,1})
        for (Coord j : {0,1})
            branch(ch,target,a,b,width,i,j,stats);
    return stats;
}

void check() {
    for (const Chart& ch : {four,eight}) for (Pair dir : {Pair{2,3},{11,2},{3,11}}) {
        for (Pair planted : {Pair{0,0},{0,17},{29,0},{16,27},{91,62},{97,100}}) {
            auto [a,b]=dir;
            auto [u,v]=planted;
            auto target=height(ch,u,v);
            Coord gap=std::abs(a*u-b*v);
            for (Coord width : {std::max(Coord(0),gap-1),gap,gap+1}) {
                std::set<Pair> expected;
                // Independent exhaustive lattice control from the positive linear terms.
                for (Coord x=0;ch.k+Int(ch.t)*x<=target;++x) {
                    Coord ylo=std::max(Coord(0),(a*x-width)/b-1);
                    for (Coord y=ylo;b*y<=a*x+width;++y)
                        if (std::abs(a*x-b*y)<=width && height(ch,x,y)==target)
                            expected.emplace(x,y);
                }
                if (scan(ch,target,a,b,width).hits!=expected)
                    throw std::logic_error("planted brute-force comparison failed");
            }
        }
    }
    int long_tests=0;
    for (const Chart& ch : {four,eight}) for (Coord a : {2,3,11})
        for (Coord b : {2,3,11}) if (a!=b) for (Coord phase : {0,1}) {
            Coord k=(Coord(1)<<43)+phase;
            Coord x=b*k+12345678+phase, y=a*k+1-phase;
            Coord gap=std::abs(a*x-b*y);
            Int target=height(ch,x,y);
            if (!scan(ch,target,a,b,gap).hits.count({x,y}))
                throw std::logic_error("large planted tube boundary was lost");
            if (scan(ch,target,a,b,gap-1).hits.count({x,y}))
                throw std::logic_error("point outside a large tube was admitted");
            long_tests+=2;
        }
    std::cout<<"108 exact brute-force comparisons and "<<long_tests
             <<" large boundary membership tests passed\n";
}

int main(int argc,char** argv) {
    if (argc==2 && std::string(argv[1])=="--check") {check();return 0;}
    if (argc!=5) {std::cerr<<"usage: spine a b half-width\n";return 2;}
    int n=std::stoi(argv[1]);
    Coord a=std::stoll(argv[2]), b=std::stoll(argv[3]), width=std::stoll(argv[4]);
    const std::array<Coord,3> directions={2,3,11};
    auto member=[&](Coord k) {return std::find(directions.begin(),directions.end(),k)!=directions.end();};
    if ((n!=4 && n!=8) || !member(a) || !member(b) || a==b || width<0 || width>(Coord(1)<<43)) return 2;
    auto stats=scan(n==4 ? four : eight,target61,a,b,width);
    std::cout<<n<<" "<<a<<" "<<b<<" "<<width<<" "<<stats.sums<<" "<<stats.occupied
             <<" "<<stats.reads<<" "<<stats.hits.size()<<"\n";
    for (auto [x,y] : stats.hits) std::cout<<"HIT "<<x<<" "<<y<<"\n";
}
