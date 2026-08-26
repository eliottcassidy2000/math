// Independent endpoint-toggle/direct-body audit for THM-4231.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <omp.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr std::array<int, 30> P = {
    8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
    120,126,132,143,145,168,170,176,190,193,240,252,264,286,290
};

struct Cell { u64 width; u32 failed; };
struct Row { u32 deletion; u64 activation; u64 mass; int components; u64 surplus; };
struct BodyTables {
    std::vector<std::pair<u32,u64>> mass;
    std::vector<std::pair<u32,int>> component;
    std::size_t raw_mass=0,raw_component=0;
};

struct Fnv1a64 {
    u64 value=0xcbf29ce484222325ULL;
    void add(u64 word) {
        for(int shift=0;shift<64;shift+=8) {
            value^=(word>>shift)&0xffULL;
            value*=0x100000001b3ULL;
        }
    }
};

[[noreturn]] static void fail(const std::string& message) {
    std::cerr << "FAIL " << message << "\n";
    std::exit(2);
}
static void require(bool ok, const std::string& message) { if (!ok) fail(message); }

static u64 gcd64(u64 a, u64 b) { return std::gcd(a,b); }
static u64 lcm64(u64 a, u64 b) { return a / gcd64(a,b) * b; }
static u64 ceil128(u128 n, u128 d) { return static_cast<u64>((n+d-1)/d); }

static std::string dec128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x) { s.push_back(char('0' + x%10)); x/=10; }
    std::reverse(s.begin(),s.end()); return s;
}

static std::string labels(u32 mask) {
    std::string s="{"; bool first=true;
    for(int i=0;i<30;i++) if(mask&(u32(1)<<i)) {
        if(!first) s+=","; first=false; s+=std::to_string(P[i]);
    }
    return s+"}";
}

static u32 mask_of(std::initializer_list<int> xs) {
    u32 ans=0;
    for(int x:xs) {
        auto it=std::find(P.begin(),P.end(),x);
        require(it!=P.end(),"unknown label");
        ans |= u32(1) << (it-P.begin());
    }
    return ans;
}

/* Independent endpoint-state sweep.  The state on (0, first wall) is all
   failures.  + events enter a bad tooth; - events leave one. */
static std::pair<u64,std::vector<Cell>> build_cells() {
    u64 D=1;
    for(int s:P) D=lcm64(D,u64(14*s));
    struct Event { u32 enter=0, leave=0; };
    std::map<u64,Event> events;
    events[0]; events[D];
    for(int i=0;i<30;i++) {
        int s=P[i]; u64 unit=D/u64(14*s); u32 bit=u32(1)<<i;
        for(int k=0;k<s;k++) {
            events[u64(14*k+1)*unit].leave |= bit;
            events[u64(14*k+13)*unit].enter |= bit;
        }
    }
    std::vector<Cell> cells; cells.reserve(events.size()-1);
    u32 state=(u32(1)<<30)-1; u64 previous=0;
    for(auto const& [wall,event]:events) {
        if(wall>previous) cells.push_back({wall-previous,state});
        state &= ~event.leave;
        state |= event.enter;
        previous=wall;
    }
    require(D==18241159416480ULL,"common denominator");
    require(cells.size()==7133,"cell count");
    return {D,cells};
}

template<class T> static T subset_sum(u32 mask, const std::unordered_map<u32,T>& table) {
    T total{}; u32 sub=mask;
    for(;;) {
        auto it=table.find(sub); if(it!=table.end()) total += it->second;
        if(sub==0) return total;
        sub=(sub-1)&mask;
    }
}

static std::vector<Row> build_rows(u64 D, const std::vector<Cell>& cells) {
    std::unordered_map<u32,u64> widths;
    std::unordered_map<u32,int> cell_counts, joins;
    for(auto const& c:cells) if(std::popcount(c.failed)<=6) {
        widths[c.failed]+=c.width; cell_counts[c.failed]++;
    }
    for(size_t i=0;i<cells.size();i++) {
        u32 joined=cells[i].failed | cells[(i+1)%cells.size()].failed;
        if(std::popcount(joined)<=6) joins[joined]++;
    }
    std::vector<Row> rows; rows.reserve(150000);
    int nonpositive=0,equalities=0;
    u32 x=(u32(1)<<6)-1, limit=u32(1)<<30;
    while(x<limit) {
        u64 mass=subset_sum(x,widths);
        int components=subset_sum(x,cell_counts)-subset_sum(x,joins);
        require(mass>0 && components>0,"row geometry");
        __int128 signed_surplus=__int128(45)*mass-__int128(4)*D;
        if(signed_surplus==0) equalities++;
        if(signed_surplus<=0) nonpositive++;
        else {
            u64 surplus=u64(signed_surplus);
            u64 activation=ceil128(u128(108)*u64(components)*D,u128(7)*surplus);
            rows.push_back({x,activation,mass,components,surplus});
        }
        u32 low=x&-x, ripple=x+low;
        if(ripple>=limit || ripple==0) break;
        x=ripple | (((ripple^x)>>2)/low);
    }
    require(rows.size()==140082,"strict row count");
    require(nonpositive==453693 && equalities==0,"nonpositive row count");
    std::sort(rows.begin(),rows.end(),[](auto const&a,auto const&b){
        return std::tie(a.activation,a.deletion)<std::tie(b.activation,b.deletion);
    });
    return rows;
}

static std::pair<u64,int> body_geometry(u32 body, const std::vector<Cell>& cells) {
    u64 mass=0; int components=0;
    bool previous=(cells.back().failed&body)==0;
    for(auto const& c:cells) {
        bool safe=(c.failed&body)==0;
        if(safe) { mass+=c.width; if(!previous) components++; }
        previous=safe;
    }
    require(mass>0 && components>0,"body geometry");
    return {mass,components};
}

static BodyTables build_body_tables(const std::vector<Cell>& cells) {
    std::unordered_map<u32,u64> mass;
    std::unordered_map<u32,int> component;
    for(auto const& c:cells) { mass[c.failed]+=c.width; component[c.failed]++; }
    for(size_t i=0;i<cells.size();i++)
        component[cells[i].failed|cells[(i+1)%cells.size()].failed]--;
    BodyTables out;
    out.raw_mass=mass.size();
    for(auto const& [mask,value]:component) if(value) out.raw_component++;
    // A nine-body can be disjoint only from masks of size at most 21.
    for(auto const& [mask,value]:mass) if(std::popcount(mask)<=21) out.mass.push_back({mask,value});
    for(auto const& [mask,value]:component)
        if(value && std::popcount(mask)<=21) out.component.push_back({mask,value});
    std::sort(out.mass.begin(),out.mass.end());
    std::sort(out.component.begin(),out.component.end());
    require(out.raw_mass==2950 && out.raw_component==1459,"raw body coefficient census");
    require(out.mass.size()==2939 && out.component.size()==1457,"effective body coefficient census");
    return out;
}

static u64 body_activation_coeff(u32 body,u64 D,const BodyTables& tables,
                                 u64* mass_out=nullptr,int* comp_out=nullptr,
                                 u64* surplus_out=nullptr) {
    u64 mass=0; int components=0;
    for(auto const& [failed,width]:tables.mass) if((failed&body)==0) mass+=width;
    for(auto const& [failed,coefficient]:tables.component)
        if((failed&body)==0) components+=coefficient;
    require(mass>0 && components>0,"coefficient body geometry");
    __int128 surplus_signed=__int128(45)*mass-__int128(4)*D;
    if(mass_out) *mass_out=mass; if(comp_out) *comp_out=components;
    if(surplus_signed<=0) { if(surplus_out)*surplus_out=0; return std::numeric_limits<u64>::max(); }
    u64 surplus=u64(surplus_signed); if(surplus_out)*surplus_out=surplus;
    return ceil128(u128(108)*u64(components)*D,u128(7)*surplus);
}

static u64 body_activation(u32 body,u64 D,const std::vector<Cell>&cells,
                           u64* mass_out=nullptr,int* comp_out=nullptr,u64* surplus_out=nullptr) {
    auto [mass,components]=body_geometry(body,cells);
    __int128 surplus_signed=__int128(45)*mass-__int128(4)*D;
    if(mass_out) *mass_out=mass; if(comp_out) *comp_out=components;
    if(surplus_signed<=0) { if(surplus_out)*surplus_out=0; return std::numeric_limits<u64>::max(); }
    u64 surplus=u64(surplus_signed); if(surplus_out)*surplus_out=surplus;
    return ceil128(u128(108)*u64(components)*D,u128(7)*surplus);
}

static u64 mix64(u64 x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct Scan { u64 q; size_t edges=0; u64 bodies=0,checks=0,covers=0,bad=0;
              u64 positive=0,nonpositive=0,max_k=0,max_k_count=0,min_k=std::numeric_limits<u64>::max();
              u32 max_body=0,min_body=0,min_mass_body=0,max_comp_body=0;
              u64 max_mass=0,max_surplus=0,min_mass=std::numeric_limits<u64>::max();
              int max_comp=0,min_mass_comp=0,max_component_value=0;
              u64 min_mass_k=0,min_mass_surplus=0,right_equalities=0;
              u64 digest_xor=0,digest_sum=0,fnv7=0;
              u64 repair_max=0,repair_max_count=0,max_body_repair=0,theta_max=0,theta_max_count=0;
              u64 second_k=0,second_k_count=0,mass_max=0;
              u32 mass_max_body=0,min_comp_body=0;
              int min_component_value=std::numeric_limits<int>::max(); };

static Scan scan_cutoff(u64 q,u64 D,const BodyTables&tables,const std::vector<Row>&rows,
                        const std::vector<u32>& bodies,bool verbose=true) {
    std::vector<u32> edges;
    for(auto const&r:rows) if(r.activation<=q) edges.push_back(r.deletion); else break;
    std::sort(edges.begin(),edges.end(),[&](u32 a,u32 b){
        u64 ka=mix64(u64(a)^q*0xd6e8feb86659fd93ULL),kb=mix64(u64(b)^q*0xd6e8feb86659fd93ULL);
        return std::tie(ka,a)<std::tie(kb,b);
    });
    Scan out;out.q=q;out.edges=edges.size();out.bodies=bodies.size();
    const bool full_ledger=edges.empty();
    std::vector<u64> ledger_mass,ledger_rho;
    std::vector<u32> ledger_k,ledger_component;
    if(full_ledger) {
        ledger_mass.resize(bodies.size());ledger_rho.resize(bodies.size());
        ledger_k.resize(bodies.size());ledger_component.resize(bodies.size());
    }
    std::vector<Scan> locals(omp_get_max_threads());
#pragma omp parallel
    {
      Scan& local=locals[omp_get_thread_num()];
#pragma omp for schedule(static)
      for(std::int64_t bi=0;bi<std::int64_t(bodies.size());bi++) {
        u32 body=bodies[bi];
        bool cover=true;
        for(u32 e:edges) { local.checks++; if((body&e)==0){cover=false;break;} }
        if(cover) {
            local.covers++;
            u64 mass=0,surplus=0;int comp=0;
            u64 k=body_activation_coeff(body,D,tables,&mass,&comp,&surplus);
            if(full_ledger) {
                require(k<=std::numeric_limits<u32>::max(),"body kappa u32");
                u64 rho=std::numeric_limits<u64>::max();
                for(auto const& row:rows) if((body&row.deletion)==0) {rho=row.activation;break;}
                require(rho!=std::numeric_limits<u64>::max(),"missing disjoint repair");
                ledger_mass[bi]=mass;ledger_component[bi]=u32(comp);
                ledger_k[bi]=u32(k);ledger_rho[bi]=rho;
            }
            if(k>q) local.bad++;
            if(k==std::numeric_limits<u64>::max()) local.nonpositive++;
            else {
                local.positive++;
                if(k<local.min_k || (k==local.min_k && body<local.min_body)) {
                    local.min_k=k;local.min_body=body;
                }
                u128 numerator=u128(108)*u64(comp)*D;
                u128 denominator=u128(7)*surplus;
                require(denominator*u128(k-1)<numerator && numerator<=denominator*u128(k),
                        "body ceiling audit");
                if(numerator==denominator*u128(k)) local.right_equalities++;
            }
            if(k>local.max_k) {
                local.max_k=k;local.max_body=body;local.max_mass=mass;
                local.max_comp=comp;local.max_surplus=surplus;
                local.max_k_count=1;
            } else if(k==local.max_k) {
                local.max_k_count++;
                if(body<local.max_body) {
                    local.max_body=body;local.max_mass=mass;
                    local.max_comp=comp;local.max_surplus=surplus;
                }
            }
            if(mass<local.min_mass || (mass==local.min_mass && body<local.min_mass_body)) {
                local.min_mass=mass;local.min_mass_body=body;local.min_mass_comp=comp;
                local.min_mass_k=k;local.min_mass_surplus=surplus;
            }
            if(comp>local.max_component_value ||
               (comp==local.max_component_value && body<local.max_comp_body)) {
                local.max_component_value=comp;local.max_comp_body=body;
            }
            u64 h=mix64(body);
            h=mix64(h^k);h=mix64(h^mass);h=mix64(h^u64(comp));h=mix64(h^surplus);
            local.digest_xor^=h;local.digest_sum+=h;
        }
      }
    }
    for(auto const& local:locals) {
        out.checks+=local.checks;out.covers+=local.covers;out.bad+=local.bad;
        out.positive+=local.positive;out.nonpositive+=local.nonpositive;
        out.right_equalities+=local.right_equalities;
        out.digest_xor^=local.digest_xor;out.digest_sum+=local.digest_sum;
        if(local.min_k<out.min_k || (local.min_k==out.min_k && local.min_body<out.min_body)) {
            out.min_k=local.min_k;out.min_body=local.min_body;
        }
        if(local.max_k>out.max_k) {
            out.max_k=local.max_k;out.max_body=local.max_body;out.max_mass=local.max_mass;
            out.max_comp=local.max_comp;out.max_surplus=local.max_surplus;
            out.max_k_count=local.max_k_count;
        } else if(local.max_k==out.max_k) {
            out.max_k_count+=local.max_k_count;
            if(local.max_body<out.max_body) {
                out.max_body=local.max_body;out.max_mass=local.max_mass;
                out.max_comp=local.max_comp;out.max_surplus=local.max_surplus;
            }
        }
        if(local.min_mass<out.min_mass ||
           (local.min_mass==out.min_mass && local.min_mass_body<out.min_mass_body)) {
            out.min_mass=local.min_mass;out.min_mass_body=local.min_mass_body;
            out.min_mass_comp=local.min_mass_comp;out.min_mass_k=local.min_mass_k;
            out.min_mass_surplus=local.min_mass_surplus;
        }
        if(local.max_component_value>out.max_component_value ||
           (local.max_component_value==out.max_component_value && local.max_comp_body<out.max_comp_body)) {
            out.max_component_value=local.max_component_value;out.max_comp_body=local.max_comp_body;
        }
    }
    require(out.bodies==14307150,"body universe");
    if(full_ledger) {
        Fnv1a64 ledger;
        std::map<u64,u64> k_hist;
        for(size_t i=0;i<bodies.size();i++) {
            const u64 surplus=45*ledger_mass[i]-4*D;
            const u64 theta=std::min<u64>(ledger_k[i],ledger_rho[i]);
            ledger.add(bodies[i]);ledger.add(ledger_k[i]);ledger.add(ledger_mass[i]);
            ledger.add(ledger_component[i]);ledger.add(surplus);ledger.add(ledger_rho[i]);ledger.add(theta);
            k_hist[ledger_k[i]]++;
            if(ledger_rho[i]>out.repair_max) {out.repair_max=ledger_rho[i];out.repair_max_count=1;}
            else if(ledger_rho[i]==out.repair_max) out.repair_max_count++;
            if(theta>out.theta_max) {out.theta_max=theta;out.theta_max_count=1;}
            else if(theta==out.theta_max) out.theta_max_count++;
            if(bodies[i]==out.max_body) out.max_body_repair=ledger_rho[i];
            if(ledger_mass[i]>out.mass_max ||
               (ledger_mass[i]==out.mass_max && bodies[i]<out.mass_max_body)) {
                out.mass_max=ledger_mass[i];out.mass_max_body=bodies[i];
            }
            if(int(ledger_component[i])<out.min_component_value ||
               (int(ledger_component[i])==out.min_component_value && bodies[i]<out.min_comp_body)) {
                out.min_component_value=int(ledger_component[i]);out.min_comp_body=bodies[i];
            }
        }
        out.fnv7=ledger.value;
        auto it=k_hist.rbegin();require(it!=k_hist.rend(),"k histogram");++it;
        require(it!=k_hist.rend(),"second k histogram");
        out.second_k=it->first;out.second_k_count=it->second;
    }
    if(verbose) std::cout<<"SCAN Q "<<q<<" EDGES "<<out.edges<<" BODIES "<<out.bodies
        <<" CHECKS "<<out.checks<<" COVERS "<<out.covers<<" BAD "<<out.bad
        <<" MAX_K "<<(out.max_k==std::numeric_limits<u64>::max()?"INF":std::to_string(out.max_k))
        <<" MAX_BODY "<<labels(out.max_body)<<" MASS "<<out.max_mass<<" COMPONENTS "<<out.max_comp
        <<" SURPLUS "<<out.max_surplus<<" MAX_K_COUNT "<<out.max_k_count
        <<" MIN_K "<<out.min_k<<" MIN_BODY "<<labels(out.min_body)
        <<" POSITIVE "<<out.positive<<" NONPOSITIVE "<<out.nonpositive
        <<" MIN_MASS "<<out.min_mass<<" MIN_MASS_BODY "<<labels(out.min_mass_body)
        <<" MIN_MASS_C "<<out.min_mass_comp<<" MIN_MASS_K "<<out.min_mass_k
        <<" MIN_MASS_SURPLUS "<<out.min_mass_surplus
        <<" MAX_C "<<out.max_component_value<<" MAX_C_BODY "<<labels(out.max_comp_body)
        <<" RIGHT_EQUALITIES "<<out.right_equalities
        <<" DIGEST_XOR "<<std::hex<<out.digest_xor<<" DIGEST_SUM "<<out.digest_sum;
    if(full_ledger) std::cout<<" FNV7 "<<out.fnv7<<std::dec
        <<" SECOND_K "<<out.second_k<<" SECOND_K_COUNT "<<out.second_k_count
        <<" REPAIR_MAX "<<out.repair_max<<" REPAIR_MAX_COUNT "<<out.repair_max_count
        <<" MAX_BODY_REPAIR "<<out.max_body_repair
        <<" THETA_MAX "<<out.theta_max<<" THETA_MAX_COUNT "<<out.theta_max_count
        <<" MASS_MAX "<<out.mass_max<<" MASS_MAX_BODY "<<labels(out.mass_max_body)
        <<" MIN_C "<<out.min_component_value<<" MIN_C_BODY "<<labels(out.min_comp_body);
    std::cout<<std::dec<<"\n";
    return out;
}

int main(int argc,char**argv) {
    auto [D,cells]=build_cells();
    auto rows=build_rows(D,cells);
    auto tables=build_body_tables(cells);
    std::vector<u32> bodies;bodies.reserve(14307150);
    u32 body=(u32(1)<<9)-1,limit=u32(1)<<30;
    while(body<limit) {
        bodies.push_back(body);
        u32 low=body&-body,ripple=body+low;
        if(ripple>=limit||ripple==0) break;
        body=ripple|(((ripple^body)>>2)/low);
    }
    require(bodies.size()==14307150,"body materialization");
    std::cout<<"READY D "<<D<<" CELLS "<<cells.size()<<" ROWS "<<rows.size()
             <<" MASS_COEFF "<<tables.mass.size()<<" COMPONENT_COEFF "<<tables.component.size()
             <<" MIN_A "<<rows.front().activation<<" MAX_A "<<rows.back().activation
             <<" THREADS "<<omp_get_max_threads()<<"\n";
    if(argc>1) {
        for(int i=1;i<argc;i++) scan_cutoff(std::stoull(argv[i]),D,tables,rows,bodies);
        return 0;
    }
    // Monotone composite predicate: no residual cover has direct threshold > Q.
    u64 lo=rows.front().activation-1,hi=17548;
    while(lo+1<hi) {
        u64 mid=(lo+hi)/2;
        auto s=scan_cutoff(mid,D,tables,rows,bodies);
        if(s.bad==0) hi=mid; else lo=mid;
    }
    auto below=scan_cutoff(hi-1,D,tables,rows,bodies);
    auto at=scan_cutoff(hi,D,tables,rows,bodies);
    std::cout<<"COMPOSITE_MIN_Q "<<hi<<" BELOW_BAD "<<below.bad<<" AT_BAD "<<at.bad<<"\n";
}
