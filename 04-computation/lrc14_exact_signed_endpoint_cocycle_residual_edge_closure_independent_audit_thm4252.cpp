// Independent exact audit for THM-4252.
//
// This path reads only the labelled prefix masks emitted by the primary,
// builds a fresh joint wall arrangement for all thirty pool speeds and the
// named pair, integrates every prefix repair directly, and recursively
// enumerates all nine-bodies.  It does not use the primitive H observable,
// fixed-cell atom coefficients, colex ranks, Gosper body enumeration, or a
// superset zeta transform.
//
// Every gate uses require(), which remains active under NDEBUG.  No C/C++
// assert, floating point, or sampling participates.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr std::array<int, 30> POOL = {
    8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
    120,126,132,143,145,168,170,176,190,193,240,252,264,286,290};
constexpr u64 EXPECTED_BODIES = UINT64_C(14307150);

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool ok, const std::string& message) {
    if (!ok) fail(message);
}

std::string dec(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x) { s.push_back(static_cast<char>('0' + x % 10)); x /= 10; }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

i64 checked_lcm(i64 a, i64 b) {
    const i128 x = static_cast<i128>(a / std::gcd(a,b)) * b;
    require(x <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(x);
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    i128 residue = static_cast<i128>(speed) * (left + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

struct Cell { i64 width; u32 failed; bool pair_safe; };

std::pair<i64,std::vector<Cell>> joint_cells(int q, int r) {
    i64 grid = 1;
    for (int s : POOL) grid = checked_lcm(grid, 14LL*s);
    grid = checked_lcm(grid, 14LL*q);
    grid = checked_lcm(grid, 14LL*r);
    std::vector<i64> walls = {0,grid};
    auto add = [&](int s) {
        const i64 unit = grid/(14LL*s);
        for (int i=0;i<s;++i) {
            walls.push_back((14LL*i+1)*unit);
            walls.push_back((14LL*i+13)*unit);
        }
    };
    for (int s:POOL) add(s);
    add(q); add(r);
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    std::vector<Cell> cells;
    cells.reserve(walls.size()-1);
    for (std::size_t i=1;i<walls.size();++i) {
        const i64 left=walls[i-1], right=walls[i];
        u32 failed=0;
        for (int j=0;j<30;++j) {
            if (!safe_midpoint(POOL[j],grid,left,right)) failed |= u32{1}<<j;
        }
        cells.push_back({right-left,failed,
            safe_midpoint(q,grid,left,right)&&safe_midpoint(r,grid,left,right)});
    }
    return {grid,std::move(cells)};
}

struct Parsed {
    int q = 0;
    int r = 0;
    std::vector<u32> deck;
    u64 declared_size = 0;
    u64 declared_fnv = 0;
    i128 primary_mass_num = 0;
    i128 primary_mass_den = 0;
};

i128 parse_i128(const std::string& text) {
    require(!text.empty(), "empty integer token");
    std::size_t index = 0;
    bool negative = false;
    if (text[index] == '-') { negative = true; ++index; }
    require(index < text.size(), "sign-only integer token");
    require(text.size() - index <= 38, "integer token exceeds safe parser range");
    i128 answer = 0;
    for (; index < text.size(); ++index) {
        require('0' <= text[index] && text[index] <= '9', "bad integer digit");
        answer = 10 * answer + (text[index] - '0');
    }
    return negative ? -answer : answer;
}

std::map<std::pair<int,int>, Parsed> parse_primary(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open primary output");
    std::map<std::pair<int,int>, Parsed> answer;
    Parsed current;
    bool inside = false;
    std::string line;
    while (std::getline(in,line)) {
        if (line.rfind("BEGIN_PAIR ", 0) == 0) {
            require(!inside, "nested pair block");
            const std::string pair = line.substr(std::string("BEGIN_PAIR ").size());
            const std::size_t comma = pair.find(',');
            require(comma != std::string::npos, "bad pair header");
            current = Parsed{};
            current.q = std::stoi(pair.substr(0, comma));
            current.r = std::stoi(pair.substr(comma + 1));
            inside = true;
        } else if (line.rfind("PREFIX_CERT ",0)==0) {
            require(inside, "prefix outside pair block");
            std::istringstream s(line);
            std::string prefix_label, size_label, fnv_label;
            require(static_cast<bool>(s >> prefix_label >> size_label >>
                                      current.declared_size >> fnv_label >>
                                      std::hex >> current.declared_fnv),
                    "truncated prefix row");
            require(prefix_label == "PREFIX_CERT" && size_label == "SIZE" &&
                        fnv_label == "FNV",
                    "bad prefix row labels");
        } else if (line.rfind("PREFIX_MASKS_HEX ",0)==0) {
            require(inside, "mask list outside pair block");
            std::string rest=line.substr(std::string("PREFIX_MASKS_HEX ").size());
            std::stringstream s(rest);
            std::string token;
            while (std::getline(s,token,',')) {
                std::size_t consumed = 0;
                const unsigned long value = std::stoul(token, &consumed, 16);
                require(consumed == token.size() && value < (1UL << 30),
                        "bad repair mask token");
                current.deck.push_back(static_cast<u32>(value));
            }
        } else if (line.rfind("DECISIVE_MASS_NUM ", 0) == 0) {
            require(inside, "mass row outside pair block");
            std::istringstream s(line);
            std::string label, numerator, den_label, denominator;
            require(static_cast<bool>(s >> label >> numerator >> den_label >> denominator),
                    "truncated decisive mass row");
            require(label == "DECISIVE_MASS_NUM" && den_label == "DEN",
                    "bad decisive mass row");
            current.primary_mass_num = parse_i128(numerator);
            current.primary_mass_den = parse_i128(denominator);
        } else if (line.rfind("END_PAIR ", 0) == 0) {
            require(inside, "end outside pair block");
            const std::string pair = line.substr(std::string("END_PAIR ").size());
            require(pair == std::to_string(current.q) + "," +
                                std::to_string(current.r),
                    "mismatched pair block terminator");
            require(current.declared_size == current.deck.size() &&
                        !current.deck.empty() && current.primary_mass_den > 0,
                    "incomplete parsed pair block");
            const auto key = std::make_pair(current.q, current.r);
            require(answer.emplace(key, current).second, "duplicate pair block");
            inside = false;
        }
    }
    require(!inside && answer.size() == 3, "primary pair universe changed");
    return answer;
}

struct Fnv { u64 x=UINT64_C(0xcbf29ce484222325); void add(u64 w) {
    for(int j=0;j<8;++j){x^=(w>>(8*j))&0xffu;x*=UINT64_C(0x100000001b3);}
}};

std::string labels(u32 mask) {
    std::string s;
    for(int i=0;i<30;++i) if(mask&(u32{1}<<i)) {
        if(!s.empty())s+=','; s+=std::to_string(POOL[i]);
    }
    return s;
}

struct BodyAudit { u64 count=0, checks=0, max_checks=0, failures=0; u32 worst=0; };

void enumerate_bodies_recursive(int next, int need, u32 mask,
                                const std::vector<u32>& deck, BodyAudit& a) {
    if (need==0) {
        ++a.count;
        u64 used=0;
        for(u32 repair:deck){++used;if((repair&mask)==0)break;}
        a.checks+=used;
        if(used==deck.size() && (deck.back()&mask)!=0) ++a.failures;
        if(used>a.max_checks || (used==a.max_checks && mask<a.worst)) {
            a.max_checks=used;a.worst=mask;
        }
        return;
    }
    for(int bit=next;bit<=30-need;++bit) {
        enumerate_bodies_recursive(bit+1,need-1,mask|(u32{1}<<bit),deck,a);
    }
}

struct Expected {
    int q;
    int r;
    i64 grid;
    u64 cells;
    i64 pair_ticks;
    u64 prefix_size;
    u64 prefix_fnv;
    i128 minimum_margin;
    u64 body_checks;
    u32 worst_body;
};

constexpr std::array<Expected, 3> EXPECTED = {{
    {466, 699, INT64_C(4250190144039840), UINT64_C(9453),
     INT64_C(3238240109744640), UINT64_C(481),
     UINT64_C(0xa5461a33cd4e1e8c), static_cast<i128>(UINT64_C(9699882064362)),
     UINT64_C(413551394), UINT32_C(0x0f722000)},
    {616, 769, INT64_C(14027451591273120), UINT64_C(9789),
     INT64_C(10305882801751680), UINT64_C(709),
     UINT64_C(0xf07e55a1f5c2b1b6), static_cast<i128>(UINT64_C(950979820896)),
     UINT64_C(418577253), UINT32_C(0x0ce07400)},
    {721, 769, INT64_C(1444827513901131360), UINT64_C(10097),
     INT64_C(1061505928580423040), UINT64_C(672),
     UINT64_C(0x88a16a55ce12cb5a),
     static_cast<i128>(UINT64_C(3535141146680682)),
     UINT64_C(414845236), UINT32_C(0x06cc3001)}
}};

i128 gcd128(i128 left, i128 right) {
    if (left < 0) left = -left;
    if (right < 0) right = -right;
    while (right != 0) {
        const i128 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

std::string fraction(i128 numerator, i128 denominator) {
    require(denominator > 0, "nonpositive fraction denominator");
    const i128 divisor = gcd128(numerator, denominator);
    return dec(numerator / divisor) + "/" + dec(denominator / divisor);
}

void audit_pair(const Expected& expected, const Parsed& parsed) {
        const int q = expected.q;
        const int r = expected.r;
        require(parsed.q == q && parsed.r == r &&
                    parsed.declared_size == expected.prefix_size &&
                    parsed.declared_fnv == expected.prefix_fnv,
                "parsed prefix declaration changed");
        Fnv hash;
        for(u32 repair:parsed.deck){require(std::popcount(repair)==8,"bad repair weight");hash.add(repair);}
        require(hash.x==parsed.declared_fnv,"prefix FNV mismatch");
        std::vector<u32> unique=parsed.deck;
        std::sort(unique.begin(),unique.end());
        require(std::adjacent_find(unique.begin(),unique.end())==unique.end(),"duplicate repair");

        auto [grid,cells]=joint_cells(q,r);
        i64 pair_ticks=0;
        for(const Cell& c:cells)if(c.pair_safe)pair_ticks+=c.width;
        require(grid == expected.grid && cells.size() == expected.cells &&
                    pair_ticks == expected.pair_ticks,
                "fresh joint wall geometry changed");
        i128 minimum_margin=0;
        u32 minimum_repair=0;
        i64 last_mass=0;
        for(std::size_t index=0;index<parsed.deck.size();++index) {
            const u32 repair=parsed.deck[index];
            i64 mass=0;
            for(const Cell& c:cells) {
                if(c.pair_safe && (c.failed&~repair)==0)mass+=c.width;
            }
            const i128 margin=static_cast<i128>(63)*mass-static_cast<i128>(4)*grid;
            require(margin>=0,"inactive repair in claimed prefix");
            if(index==0||margin<minimum_margin){minimum_margin=margin;minimum_repair=repair;}
            if (index + 1 == parsed.deck.size()) last_mass = mass;
        }
        require(minimum_margin == expected.minimum_margin,
                "minimum direct prefix margin changed");
        require(static_cast<i128>(last_mass) * parsed.primary_mass_den ==
                    parsed.primary_mass_num * grid,
                "direct final-repair mass disagrees with primary component mass");
        BodyAudit audit;
        enumerate_bodies_recursive(0,9,0,parsed.deck,audit);
        require(audit.count==EXPECTED_BODIES,"body universe mismatch");
        require(audit.failures==0,"prefix fails universal body coverage");
        require(audit.max_checks==parsed.deck.size(),"prefix is not order-minimal");
        require(audit.checks == expected.body_checks &&
                    audit.worst == expected.worst_body,
                "recursive body ledger changed");
        require((audit.worst&parsed.deck.back())==0,"last repair misses minimality witness");
        for(std::size_t i=0;i+1<parsed.deck.size();++i)
            require((audit.worst&parsed.deck[i])!=0,"earlier repair misses witness");

        std::cout<<"BEGIN_PAIR "<<q<<','<<r<<'\n';
        std::cout<<"PAIR "<<q<<','<<r<<" GRID "<<grid<<" JOINT_CELLS "<<cells.size()
                 <<" PAIR_TICKS "<<pair_ticks<<'/'<<grid<<'\n';
        std::cout<<"PREFIX SIZE "<<parsed.deck.size()<<" FNV "<<std::hex<<hash.x
                 <<std::dec<<" MIN_MARGIN_NUM "<<dec(minimum_margin)
                 <<" MIN_REPAIR {"<<labels(minimum_repair)<<"}\n";
        std::cout<<"BODY_RECURSION COUNT "<<audit.count<<" FAILURES "<<audit.failures
                 <<" CHECKS "<<audit.checks<<" MAX_CHECKS "<<audit.max_checks
                 <<" MINIMALITY_WITNESS {"<<labels(audit.worst)<<"}\n";
        std::cout<<"FINAL_REPAIR_DIRECT_MASS "<<last_mass<<'/'<<grid
                 <<" SCALED_MARGIN_63MU_MINUS4 "
                 <<fraction(static_cast<i128>(63)*last_mass-static_cast<i128>(4)*grid,
                             grid)
                 <<" CROSS_PRIMARY PASS\n";
        std::cout<<"CHECKS PASS fresh_joint_walls,direct_midpoint_mass,recursive_bodies,"
                    "prefix_order_minimality\n";
        std::cout<<"END_PAIR "<<q<<','<<r<<'\n';
}

} // namespace

int main(int argc,char**argv) {
    try {
        require(argc==2,"usage: independent primary.out");
        const auto parsed = parse_primary(argv[1]);
        std::cout << "THM4252_INDEPENDENT_JOINT_WALL_REPLAY_V1\n";
        std::cout << "PAIRS " << EXPECTED.size() << " BODIES_PER_PAIR "
                  << EXPECTED_BODIES << '\n';
        for (const Expected& expected : EXPECTED) {
            const auto found = parsed.find({expected.q, expected.r});
            require(found != parsed.end(), "named pair missing from primary output");
            audit_pair(expected, found->second);
        }
        std::cout << "GLOBAL_VERDICT PASS FRESH_JOINT_WALLS THREE_PREFIXES\n";
        return 0;
    } catch(const std::exception&e){
        std::cerr<<"THM4252_INDEPENDENT_ERROR "<<e.what()<<'\n';
        return 1;
    }
}
