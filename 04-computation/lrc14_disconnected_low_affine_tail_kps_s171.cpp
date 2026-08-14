// Exact affine-tail compiler for the disconnected-low LRC14 pair floor.
//
// The Python wrapper supplies the 79 small-ruler contexts not already closed
// by the monotone midpoint envelope.  This program enumerates the canonical
// 22,890 nonzero Dirichlet rays, computes their exact homogenized limits,
// applies the proved O(1/N) averaging bound and the many-turn exit, and
// evaluates every remaining gcd<=3 physical row by integer floor moments.

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

namespace {

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }
void require(bool condition, const std::string& message) { if (!condition) fail(message); }

std::string text128(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    u128 magnitude = negative ? u128(-(value + 1)) + 1 : u128(value);
    std::string answer;
    while (magnitude) { answer.push_back(char('0' + magnitude % 10)); magnitude /= 10; }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

i64 narrow(i128 value, const char* label) {
    if (value < std::numeric_limits<i64>::min() || value > std::numeric_limits<i64>::max())
        fail(std::string("int64 overflow in ") + label + ": " + text128(value));
    return i64(value);
}

i64 floor_div(i64 a, i64 b) {
    require(b > 0, "nonpositive floor divisor");
    i64 q = a / b, r = a % b;
    if (r < 0) --q;
    return q;
}
i64 ceil_div(i64 a, i64 b) { return -floor_div(-a, b); }
i64 positive_mod(i64 a, i64 m) { i64 r = a % m; return r < 0 ? r + m : r; }

struct Fraction {
    i64 n = 0, d = 1;
    Fraction() = default;
    Fraction(i64 numerator, i64 denominator) : n(numerator), d(denominator) {
        require(d != 0, "zero denominator");
        if (d < 0) { n = -n; d = -d; }
        i64 g = std::gcd(n < 0 ? -n : n, d); n /= g; d /= g;
    }
};

Fraction subtract(Fraction a, Fraction b) {
    return Fraction(narrow(i128(a.n) * b.d - i128(b.n) * a.d, "fraction subtraction"),
                    narrow(i128(a.d) * b.d, "fraction denominator"));
}
bool less(Fraction a, Fraction b) { return i128(a.n) * b.d < i128(b.n) * a.d; }
std::string qtext(Fraction x) { return x.d == 1 ? std::to_string(x.n) : std::to_string(x.n) + "/" + std::to_string(x.d); }

struct Moments { i64 s0, s1, s2; };

Moments floor_moments(i64 n, i64 m, i64 a, i64 b) {
    if (n == 0) return {0, 0, 0};
    require(n >= 0 && m > 0, "floor-moment domain");
    i128 S1 = i128(n) * (n - 1) / 2;
    i128 S2 = i128(n) * (n - 1) * (2 * i128(n) - 1) / 6;
    i64 qa = floor_div(a, m), ar = positive_mod(a, m);
    i64 qb = floor_div(b, m), br = positive_mod(b, m);
    i128 b0 = i128(qa) * S1 + i128(qb) * n;
    i128 b1 = i128(qa) * S2 + i128(qb) * S1;
    i128 b2 = i128(qa) * qa * S2 + 2 * i128(qa) * qb * S1 + i128(qb) * qb * n;
    if (ar == 0) return {narrow(b0,"fm0"),narrow(b1,"fm1"),narrow(b2,"fm2")};
    i64 height = narrow((i128(ar) * (n - 1) + br) / m, "floor height");
    if (height == 0) return {narrow(b0,"fmz0"),narrow(b1,"fmz1"),narrow(b2,"fmz2")};
    Moments u = floor_moments(height, ar, m, m - br + ar - 1);
    i128 r0 = i128(n) * height - u.s0;
    i128 r1 = i128(height) * S1 - (i128(u.s2) - u.s0) / 2;
    i128 r2 = i128(n) * height * height - 2 * i128(u.s1) - u.s0;
    return {narrow(b0+r0,"fmr0"),narrow(b1+r1,"fmr1"),
            narrow(b2+2*i128(qa)*r1+2*i128(qb)*r0+r2,"fmr2")};
}

struct Prefix { i64 count, sum; };
Prefix residue_prefix(i64 n,i64 m,i64 a,i64 b,i64 threshold,const Moments& base,i64 total) {
    if (threshold <= 0) return {0,0};
    if (threshold >= m) return {n,total};
    Moments shifted = floor_moments(n,m,a,b+m-threshold);
    i64 d0=shifted.s0-base.s0,d1=shifted.s1-base.s1;
    i64 y=(shifted.s2-base.s2-d0)/2;
    i64 high=narrow(i128(a)*d1+i128(b)*d0-i128(m)*y,"prefix sum");
    return {n-d0,total-high};
}

i64 triangle_sum(i64 n,i64 m,i64 a,i64 b,i64 peak,i64 L,const Moments& base,i64 total) {
    if (peak <= 0) return 0;
    i64 radius=(peak-1)/L,turns=radius/m,tail=radius%m;
    i128 answer=i128(n)*(2*i128(turns)*peak-i128(L)*m*turns*turns);
    Prefix low=residue_prefix(n,m,a,b,tail+1,base,total);
    answer+=i128(peak-L*turns*m)*low.count-i128(L)*low.sum;
    Prefix before=residue_prefix(n,m,a,b,m-tail,base,total);
    i64 hc=n-before.count,hs=total-before.sum;
    answer+=i128(peak-L*(turns+1)*m)*hc+i128(L)*hs;
    return narrow(answer,"triangle sum");
}

struct Context { int L,cell,e,f,mask; };
struct Ray { int d,a,c,p0,q0; };

Fraction physical_mass(const Context& x,i64 p,i64 q) {
    i64 z=i64(x.L)*p-x.e,w=i64(x.L)*q-x.f;
    Context y=x;
    if (z>w) { std::swap(z,w); std::swap(p,q); std::swap(y.e,y.f); }
    i64 r=positive_mod(i64(y.e)*y.cell,y.L),s=positive_mod(i64(y.f)*y.cell,y.L);
    i64 determinant=narrow(i128(r)*w-i128(s)*z,"phase determinant");
    require(determinant%y.L==0,"nonintegral phase");
    i64 b=positive_mod(determinant/y.L,z),a=w%z;
    Moments base=floor_moments(p,z,a,b);
    i64 total=narrow(i128(a)*p*(p-1)/2+i128(b)*p-i128(z)*base.s0,"residue total");
    i64 unit=y.L/14;
    i64 outer=narrow(i128(unit)*(z+w),"outer peak"),inner=narrow(i128(unit)*(w-z),"inner peak");
    i64 numerator=triangle_sum(p,z,a,b,outer,y.L,base,total)-triangle_sum(p,z,a,b,inner,y.L,base,total);
    i64 denominator=narrow(i128(z)*w,"mass denominator");
    require(numerator>=0&&denominator>0,"mass sign");
    return Fraction(numerator,denominator);
}

i128 triangle_primitive2(i64 radius,i64 u) {
    if (u<=-radius) return 0;
    if (u<=0) return i128(u+radius)*(u+radius);
    if (u<radius) return i128(radius)*radius+2*i128(radius)*u-i128(u)*u;
    return 2*i128(radius)*radius;
}

i128 integrated_tent2(i64 radius,i64 left,i64 right,i64 period) {
    require(left<=right&&period>0,"integrated tent domain");
    i64 lo=floor_div(-radius-right,period)-2;
    i64 hi=floor_div(radius-left,period)+2;
    i128 answer=0;
    for(i64 shift=lo;shift<=hi;++shift)
        answer+=triangle_primitive2(radius,right+shift*period)-triangle_primitive2(radius,left+shift*period);
    return answer;
}

i64 tent_value_numerator(i64 radius,i64 point,i64 period) {
    i64 lo=floor_div(-radius-point,period)-2;
    i64 hi=floor_div(radius-point,period)+2;
    i128 answer=0;
    for(i64 shift=lo;shift<=hi;++shift)
        answer+=std::max<i64>(0,radius-std::llabs(point+shift*period));
    return narrow(answer,"tent value");
}

Fraction ray_limit(const Context& x,int d,int a,int c) {
    int D=d+a,k=std::gcd(d,D),P=d/k,Q=D/k;
    i64 M=i64(k)*x.L;
    i64 r=positive_mod(i64(x.e)*x.cell,x.L),s=positive_mod(i64(x.f)*x.cell,x.L);
    i64 A=i64(x.L)*c+i64(D)*x.e-i64(d)*x.f;
    i64 origin=-i64(D)*r+i64(d)*s;
    i64 RA=i64(P+Q)*M/14,RB=i64(std::abs(P-Q))*M/14;
    require((i64(P+Q)*M)%14==0&&RA>=RB,"tent scaling");
    if(A==0) {
        i64 numerator=tent_value_numerator(RA,origin,M)-tent_value_numerator(RB,origin,M);
        return Fraction(numerator,narrow(i128(M)*P*Q,"constant limit denominator"));
    }
    i64 other=origin-A,left=std::min(origin,other),right=std::max(origin,other);
    i128 numerator=integrated_tent2(RA,left,right,M)-integrated_tent2(RB,left,right,M);
    i128 denominator=2*i128(M)*std::llabs(A)*P*Q;
    return Fraction(narrow(numerator,"limit numerator"),narrow(denominator,"limit denominator"));
}

std::vector<Ray> make_rays() {
    std::vector<Ray> rays;
    for(int d=1;d<=8;++d) for(int a=0;a<=7*d;++a) for(int c=-46;c<=46;++c) {
        if(c==0||(a==0&&c<0)||(a==7*d&&c>0)) continue;
        for(int p0=1;p0<=d;++p0) if((a*p0+c)%d==0)
            rays.push_back({d,a,c,p0,p0+(a*p0+c)/d});
    }
    require(rays.size()==22890,"ray count");
    return rays;
}

std::vector<Context> read_contexts(const std::string& path) {
    std::ifstream input(path); require(bool(input),"cannot open contexts");
    std::vector<Context> rows; Context x;
    while(input>>x.L>>x.cell>>x.e>>x.f>>x.mask) rows.push_back(x);
    require(input.eof(),"malformed contexts");
    require(rows.size()==79,"hostile context count");
    for(const auto& r:rows) {
        require(r.mask>=1&&r.mask<=7,"context mask");
        int unit=r.L/14;
        int re=int(positive_mod(i64(r.e)*r.cell,r.L));
        int rf=int(positive_mod(i64(r.f)*r.cell,r.L));
        require(unit<=re&&re+r.e<=r.L-unit&&unit<=rf&&rf+r.f<=r.L-unit,"unsafe context");
    }
    return rows;
}

u64 mix(u64 x) {
    x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);
}
void hash_value(u64& state,u64 value){state^=mix(value+state);state*=0x100000001b3ULL;}
std::string hex64(u64 x){std::ostringstream out;out<<std::hex<<std::setw(16)<<std::setfill('0')<<x;return out.str();}

struct Minimum {
    bool valid=false; Fraction value; Ray ray{}; Context context{}; i64 n=0,p=0,q=0;
};
bool better(const Minimum& x,const Minimum& y) {
    if(!y.valid) return true;
    if(less(x.value,y.value)) return true;
    if(less(y.value,x.value)) return false;
    return std::tie(x.ray.d,x.ray.a,x.ray.c,x.ray.p0,x.context.L,x.context.cell,x.context.e,x.context.f,x.n)
         < std::tie(y.ray.d,y.ray.a,y.ray.c,y.ray.p0,y.context.L,y.context.cell,y.context.e,y.context.f,y.n);
}

struct Result {
    u64 periodic_candidates=0,scope_skips=0,envelope_skips=0,midpoint_skips=0,literal_checks=0,failures=0;
    u64 hash=0xcbf29ce484222325ULL; Minimum minimum;
    i64 largest_head=0; u64 convergence_exits=0,turn_exits=0;
};

i64 cone_start(const Ray& ray) {
    int D=ray.d+ray.a;
    // A formal affine ray covers a physical Dirichlet incidence only once
    // its inherited approximation condition |c|<=p/9 holds.
    i64 n=std::max<i64>({0,ceil_div(264-ray.p0,ray.d),
                         ceil_div(9*i64(std::abs(ray.c))-ray.p0,ray.d),
                         ceil_div(1-ray.q0,D)});
    while(i64(ray.p0)+i64(ray.d)*n>=i64(ray.q0)+i64(D)*n) ++n;
    while(i64(ray.q0)+i64(D)*n>=8*(i64(ray.p0)+i64(ray.d)*n)) ++n;
    return n;
}

i64 convergence_start(const Context& x,const Ray& ray,Fraction limit) {
    Fraction gap=subtract(limit,Fraction(1,294)); require(gap.n>0,"nonpositive limit gap");
    int D=ray.d+ray.a;
    i64 A=i64(x.L)*ray.c+i64(D)*x.e-i64(ray.d)*x.f;
    i64 z0=i64(x.L)*ray.p0-x.e;
    // With alpha=p0-e/L, beta=q0-f/L, k=gcd(d,d+a), P=d/k,
    // put N=kn, epsilon=alpha/P and C=P beta-Q alpha=A/(kL).
    // Reparametrizing the first clock exactly and comparing the N full
    // cycles with the limiting phase integral gives
    // |I_n-J| <= [2eps/7+6/(49P)+|C|(3+eps)/P]/(N+eps)
    //          = [2z0/7+6L/49+3|A|/k+|A|z0/(dL)]/z.
    // The 1/7 and 6/49 terms use the exact mass and primitive oscillation
    // of the first circular interval, rather than bounding its tail by one.
    int k=std::gcd(ray.d,D);
    i64 absA=std::llabs(A);
    i128 loss_numerator=14*i128(z0)*k*ray.d*x.L
                         +6*i128(x.L)*k*ray.d*x.L
                         +147*i128(absA)*ray.d*x.L
                         +49*i128(absA)*z0*k;
    i128 loss_denominator=49*i128(k)*ray.d*x.L;
    // Require gap*z >= loss_numerator/loss_denominator.
    i128 required_numerator=loss_numerator*gap.d;
    i128 required_denominator=loss_denominator*gap.n;
    require(required_numerator>=0&&required_denominator>0,"convergence ratio");
    i64 required_z=narrow((required_numerator+required_denominator-1)/required_denominator,
                          "required convergence z");
    i64 n=std::max<i64>(0,ceil_div(required_z-z0,i64(x.L)*ray.d));
    i64 z=z0+i64(x.L)*ray.d*n;
    require(i128(gap.n)*z*loss_denominator>=i128(gap.d)*loss_numerator,
            "convergence threshold");
    return n;
}

struct Rational128 { i128 n=0,d=1; };
i128 gcd128(i128 a,i128 b){if(a<0)a=-a;if(b<0)b=-b;while(b){i128 r=a%b;a=b;b=r;}return a;}
Rational128 reduced128(i128 n,i128 d){require(d!=0,"rational128 zero denominator");if(d<0){n=-n;d=-d;}i128 g=gcd128(n,d);return{n/g,d/g};}
Rational128 add128(Rational128 x,Rational128 y){
    i128 g=gcd128(x.d,y.d),xd=x.d/g,yd=y.d/g;
    return reduced128(x.n*yd+y.n*xd,xd*y.d);
}
Rational128 gamma_error(int K,int g,int L,int endpoint){
    if(K%2==0)return reduced128(i128(endpoint)*K,2*i128(g*L)*K-2*endpoint);
    return reduced128(i128(endpoint)*(i128(K)*K+1),2*i128(K)*(i128(g*L)*K-endpoint));
}
bool midpoint_closes(const Context& x,i64 p,i64 q,i64 g){
    int P=int(p/g),Q=int(q/g);require(std::gcd(P,Q)==1,"midpoint primitive");
    Rational128 phase1{1,105};
    Rational128 phase2{i128(P)*Q-12,49*i128(P)*Q};
    Rational128 phase=(phase1.n*phase2.d>=phase2.n*phase1.d)?phase1:phase2;
    Rational128 error=add128(gamma_error(P,int(g),x.L,x.e),gamma_error(Q,int(g),x.L,x.f));
    i64 C=std::llabs(i64(Q)*x.e-i64(P)*x.f);
    error=add128(error,reduced128(i128(C)*(C/x.L+1),2*i128(g)*g*x.L*P*Q));
    // phase-error >= 1/294, compared without forming the subtraction.
    return (phase.n*error.d-error.n*phase.d)*294>=phase.d*error.d;
}

i64 turn_start(const Context& x,const Ray& ray) {
    int D=ray.d+ray.a;
    i64 A=std::llabs(i64(x.L)*ray.c+i64(D)*x.e-i64(ray.d)*x.f);
    i64 coefficient=A-5*i64(x.L)*ray.d;
    if(coefficient<=0) return std::numeric_limits<i64>::max()/4;
    i64 b0=ray.p0/ray.d,z0=i64(x.L)*ray.p0-x.e;
    // The generalized many-turn lemma beats Dmax/5 from p>=264, but its
    // stronger edgewise 1/294 conclusion starts at p>=273.  Keep the nine
    // intervening raw levels in the literal/convergence certificate so this
    // compiler really proves the advertised universal pair floor.
    i64 n=std::max<i64>({0,ceil_div(5*z0-b0*A,coefficient),
                         ceil_div(273-ray.p0,ray.d),
                         ceil_div(9*i64(std::abs(ray.c))-ray.p0,ray.d)});
    i64 p=ray.p0+i64(ray.d)*n,z=i64(x.L)*p-x.e;
    require(p>=273,"edgewise turn threshold scope");
    require(p>=9*i64(std::abs(ray.c)),"Dirichlet turn threshold scope");
    require(i128(p/ray.d)*A>=5*i128(z),"turn threshold");
    return n;
}

Result prove_ray(const Ray& ray,const std::vector<Context>& contexts) {
    Result out;
    i64 start=cone_start(ray);
    for(const Context& x:contexts) {
        Fraction limit=ray_limit(x,ray.d,ray.a,ray.c);
        i64 conv=convergence_start(x,ray,limit),turn=turn_start(x,ray);
        i64 stop=std::min(conv,turn);
        if(conv<=turn) ++out.convergence_exits; else ++out.turn_exits;
        if(stop<start) stop=start;
        out.largest_head=std::max(out.largest_head,stop-start);
        int D=ray.d+ray.a;
        i64 period=std::llabs(ray.c);
        for(i64 residue=0;residue<period;++residue) {
            i64 witness_n=start+positive_mod(residue-start,period);
            i64 witness_p=ray.p0+i64(ray.d)*witness_n,witness_q=ray.q0+i64(D)*witness_n;
            i64 residue_g=std::gcd(witness_p,witness_q);
            // gcd(p,q) divides c and therefore is constant modulo |c|.
            require(residue_g==std::gcd(std::gcd(witness_p, witness_q),period),"periodic gcd");
            i64 count=witness_n>=stop?0:1+(stop-1-witness_n)/period;
            if(residue_g>3||witness_p/residue_g+witness_q/residue_g<8
               ||(witness_p/residue_g==3&&witness_q/residue_g==5)) {
                out.scope_skips+=count;continue;
            }
            if(!(x.mask&(1<<(residue_g-1)))) {out.envelope_skips+=count;continue;}
            for(i64 n=witness_n;n<stop;n+=period) {
                ++out.periodic_candidates;
                i64 p=ray.p0+i64(ray.d)*n,q=ray.q0+i64(D)*n;
                require(std::gcd(p,q)==residue_g,"gcd period changed");
                if(midpoint_closes(x,p,q,residue_g)){++out.midpoint_skips;continue;}
                Fraction value=physical_mass(x,p,q); ++out.literal_checks;
                for(i64 field:std::array<i64,14>{ray.d,ray.a,ray.c,ray.p0,ray.q0,x.L,x.cell,x.e,x.f,n,p,q,value.n,value.d})
                    hash_value(out.hash,u64(field));
                Minimum candidate{true,value,ray,x,n,p,q};
                if(better(candidate,out.minimum)) out.minimum=candidate;
                if(i128(value.n)*294<value.d) ++out.failures;
            }
        }
    }
    return out;
}

std::string minimum_text(const Minimum& m) {
    if(!m.valid) return "none";
    std::ostringstream out;out<<qtext(m.value)<<"@(d,a,c,p0,q0;L,cell,e,f;n,p,q)=("
        <<m.ray.d<<','<<m.ray.a<<','<<m.ray.c<<','<<m.ray.p0<<','<<m.ray.q0<<';'
        <<m.context.L<<','<<m.context.cell<<','<<m.context.e<<','<<m.context.f<<';'
        <<m.n<<','<<m.p<<','<<m.q<<')';return out.str();
}

int probe(const std::string& path) {
    std::ifstream input(path);require(bool(input),"cannot open probes");
    Context x;Ray r;i64 n;
    while(input>>x.L>>x.cell>>x.e>>x.f>>r.d>>r.a>>r.c>>r.p0>>r.q0>>n) {
        int D=r.d+r.a;i64 p=r.p0+i64(r.d)*n,q=r.q0+i64(D)*n;
        Fraction mass=physical_mass(x,p,q),limit=ray_limit(x,r.d,r.a,r.c);
        std::cout<<x.L<<' '<<x.cell<<' '<<x.e<<' '<<x.f<<' '<<r.d<<' '<<r.a<<' '<<r.c<<' '
                 <<r.p0<<' '<<r.q0<<' '<<n<<' '<<mass.n<<' '<<mass.d<<' '<<limit.n<<' '<<limit.d<<'\n';
    }
    require(input.eof(),"malformed probes");return 0;
}

int scan(const std::string& context_path,const std::string& summary_path,int threads) {
    auto contexts=read_contexts(context_path);auto rays=make_rays();
    std::set<std::tuple<int,int,int>> triples;
    for(const auto&r:rays) triples.insert({r.d,r.a,r.c});
    require(triples.size()==17206,"unique triple count");

    Minimum limit_minimum;u64 limit_checks=0;
    for(auto [d,a,c]:triples) for(const auto&x:contexts) {
        Fraction value=ray_limit(x,d,a,c);++limit_checks;
        Ray r{d,a,c,0,0};Minimum candidate{true,value,r,x,0,0,0};
        if(better(candidate,limit_minimum))limit_minimum=candidate;
    }
    require(limit_checks==1359274,"limit check count");
    require(limit_minimum.value.n==709&&limit_minimum.value.d==48048,"limit minimum");

    std::vector<Result> results(rays.size());std::atomic<std::size_t> next{0};
    std::vector<std::thread> workers;
    for(int worker=0;worker<threads;++worker) workers.emplace_back([&](){
        for(;;){std::size_t i=next.fetch_add(1,std::memory_order_relaxed);if(i>=rays.size())return;
            results[i]=prove_ray(rays[i],contexts);}
    });
    for(auto&worker:workers)worker.join();

    std::ofstream summary(summary_path,std::ios::binary);require(bool(summary),"cannot create summary");
    u64 candidates=0,scope=0,envelope=0,midpoint=0,checks=0,failures=0,conv=0,turn=0; i64 maxhead=0;Minimum minimum;
    for(std::size_t i=0;i<rays.size();++i){const auto&r=rays[i];const auto&z=results[i];
        candidates+=z.periodic_candidates;scope+=z.scope_skips;envelope+=z.envelope_skips;midpoint+=z.midpoint_skips;checks+=z.literal_checks;
        failures+=z.failures;conv+=z.convergence_exits;turn+=z.turn_exits;maxhead=std::max(maxhead,z.largest_head);
        if(z.minimum.valid&&better(z.minimum,minimum))minimum=z.minimum;
        summary<<r.d<<' '<<r.a<<' '<<r.c<<' '<<r.p0<<' '<<r.q0<<' '<<z.periodic_candidates<<' '
               <<z.scope_skips<<' '<<z.envelope_skips<<' '<<z.midpoint_skips<<' '<<z.literal_checks<<' '<<z.failures<<' '
               <<z.convergence_exits<<' '<<z.turn_exits<<' '<<z.largest_head<<' '<<hex64(z.hash)<<' '
               <<minimum_text(z.minimum)<<'\n';
    }
    require(bool(summary),"summary write");require(failures==0,"literal failure");
    Fraction limit_gap=subtract(limit_minimum.value,Fraction(1,294));
    std::cout<<"LRC14 DISCONNECTED-LOW AFFINE TAIL CERTIFICATE\n";
    std::cout<<"status=PROVED analytic averaging/turn exits + FINITE-EXACT residual\n";
    std::cout<<"contexts="<<contexts.size()<<";rays="<<rays.size()<<";unique_limit_triples="<<triples.size()<<"\n";
    std::cout<<"limit_checks="<<limit_checks<<";limit_minimum="<<minimum_text(limit_minimum)
             <<";limit_gap="<<qtext(limit_gap)<<"\n";
    std::cout<<"ray_contexts="<<u64(rays.size())*contexts.size()<<";convergence_exits="<<conv<<";turn_exits="<<turn
             <<";maximum_finite_head="<<maxhead<<"\n";
    std::cout<<"periodic_candidates="<<candidates<<";scope_skips="<<scope<<";envelope_skips="<<envelope
             <<";midpoint_skips="<<midpoint<<";literal_checks="<<checks<<"\n";
    std::cout<<"literal_minimum="<<minimum_text(minimum)<<";failures_below_1/294="<<failures<<"\n";
    std::cout<<"summary_rows="<<rays.size()<<"\n";
    std::cout<<"conclusion=every_Dirichlet-admissible_affine-ray_incidence_in_the_remaining_gcd<=3_scope_has_mass>=1/294\n";
    return 0;
}

} // namespace

int main(int argc,char**argv){
    try{
        if(argc==3&&std::string(argv[1])=="--probe")return probe(argv[2]);
        if(argc<3||argc>4){std::cerr<<"usage: affine-tail CONTEXTS SUMMARY [THREADS] or --probe PROBES\n";return 2;}
        int threads=argc==4?std::stoi(argv[3]):int(std::thread::hardware_concurrency());if(threads<1)threads=1;
        return scan(argv[1],argv[2],threads);
    }catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
