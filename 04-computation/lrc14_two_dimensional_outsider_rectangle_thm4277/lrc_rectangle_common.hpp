#ifndef LRC14_RECTANGLE_COMMON_HPP
#define LRC14_RECTANGLE_COMMON_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace lrc_rectangle {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr std::array<int,30> POOL = {
    8,10,15,16,20,30,40,42,60,63,
    80,84,85,88,95,120,126,132,143,145,
    168,170,176,190,193,240,252,264,286,290};
constexpr i64 D = INT64_C(18241159416480);
constexpr u64 REPAIR_COUNT = UINT64_C(5852925);
constexpr u64 BODY_COUNT = UINT64_C(14307150);
constexpr u64 ORDER_SEED = UINT64_C(0x4245422842334245);
constexpr std::size_t CANDIDATES = 16384;
constexpr std::size_t HOSTILE_BUDGET = 8192;
constexpr int Q0 = 450;
constexpr int Q1 = 499;
constexpr int R0 = 600;
constexpr int R1 = 650;
constexpr std::size_t PAIR_COUNT = 2550;

[[noreturn]] inline void fail(const std::string& message) {
    throw std::runtime_error(message);
}

inline void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

inline std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(),result.end());
    return result;
}

inline std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::nouppercase << std::setw(16) << std::setfill('0')
        << value;
    return out.str();
}

inline i128 gcd128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) { const i128 r = a % b; a = b; b = r; }
    return a;
}

struct U256 {
    std::array<std::uint32_t,8> limb{};
};

inline U256 multiply_nonnegative_i128(i128 left,i128 right) {
    require(left>=0 && right>=0,"negative factor in U256 product");
    const u128 a=static_cast<u128>(left);
    const u128 b=static_cast<u128>(right);
    std::array<std::uint32_t,4> al{},bl{};
    for (unsigned i=0;i<4;++i) {
        al[i]=static_cast<std::uint32_t>(a>>(32*i));
        bl[i]=static_cast<std::uint32_t>(b>>(32*i));
    }
    U256 product;
    for (unsigned i=0;i<4;++i) {
        std::uint64_t carry=0;
        for (unsigned j=0;j<4;++j) {
            const unsigned k=i+j;
            const std::uint64_t current=
                static_cast<std::uint64_t>(al[i])*bl[j]+product.limb[k]+carry;
            product.limb[k]=static_cast<std::uint32_t>(current);
            carry=current>>32;
        }
        unsigned k=i+4;
        while (carry!=0) {
            require(k<8,"U256 carry overflow");
            const std::uint64_t current=static_cast<std::uint64_t>(product.limb[k])+carry;
            product.limb[k]=static_cast<std::uint32_t>(current);
            carry=current>>32;
            ++k;
        }
    }
    return product;
}

inline bool u256_less(const U256& left,const U256& right) {
    for (int i=7;i>=0;--i) {
        if (left.limb[static_cast<std::size_t>(i)]!=right.limb[static_cast<std::size_t>(i)])
            return left.limb[static_cast<std::size_t>(i)]<right.limb[static_cast<std::size_t>(i)];
    }
    return false;
}

inline bool product_less_exact(i128 a,i128 b,i128 c,i128 d) {
    return u256_less(multiply_nonnegative_i128(a,b),
                     multiply_nonnegative_i128(c,d));
}

inline void check_u256_products() {
    require(product_less_exact(3,7,5,5),"U256 small comparison failed");
    require(!product_less_exact(5,5,3,7),"U256 reverse comparison failed");
    require(!product_less_exact(9,11,9,11),"U256 equality comparison failed");
    const i128 a=static_cast<i128>(1)<<120;
    const i128 b=static_cast<i128>(1)<<110;
    const i128 c=static_cast<i128>(1)<<119;
    const i128 d=static_cast<i128>(1)<<110;
    require(!product_less_exact(a,b,c,d),"U256 high-limb comparison failed");
    require(product_less_exact(c,d,a,b),"U256 high-limb reverse failed");
    const U256 square=multiply_nonnegative_i128(a+b,a+b);
    const U256 reverse=multiply_nonnegative_i128(a+b,a+b);
    require(!u256_less(square,reverse) && !u256_less(reverse,square),
            "U256 deterministic square failed");
}

inline std::string fraction(i128 numerator, i128 denominator) {
    require(denominator != 0,"zero fraction denominator");
    if (denominator < 0) { numerator=-numerator; denominator=-denominator; }
    const i128 divisor = gcd128(numerator,denominator);
    return decimal(numerator/divisor) + "/" + decimal(denominator/divisor);
}

inline u64 splitmix64(u64 x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

class Fnv64 {
  public:
    void byte(std::uint8_t value) {
        state_ ^= value;
        state_ *= UINT64_C(0x100000001b3);
    }
    void word(u64 value) {
        for (unsigned byte_index=0; byte_index<8; ++byte_index)
            byte(static_cast<std::uint8_t>((value>>(8*byte_index))&0xffu));
    }
    u64 value() const { return state_; }
  private:
    u64 state_ = UINT64_C(0xcbf29ce484222325);
};

inline u32 next_mask(u32 mask) {
    const u32 low = mask & (~mask + 1u);
    const u32 ripple = mask + low;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low);
}

inline std::string labels(u32 mask) {
    std::string result;
    for (int bit=0; bit<30; ++bit) if (mask & (u32{1}<<bit)) {
        if (!result.empty()) result += ',';
        result += std::to_string(POOL[bit]);
    }
    return result;
}

inline i64 lcm64(i64 a, i64 b) {
    const i64 divisor = std::gcd(a,b);
    const i128 value = static_cast<i128>(a/divisor)*b;
    require(value<=std::numeric_limits<i64>::max(),"i64 lcm overflow");
    return static_cast<i64>(value);
}

struct Cell { i64 left=0; i64 right=0; u32 failed=0; };

inline bool pool_safe_midpoint(int speed, i64 left, i64 right) {
    const i64 modulus = 2*D;
    const i64 residue = static_cast<i64>(
        (static_cast<i128>(speed)*(left+right))%modulus);
    return static_cast<i128>(7)*residue>=D &&
           static_cast<i128>(7)*residue<=static_cast<i128>(13)*D;
}

inline std::vector<Cell> build_pool_cells() {
    i64 check=1;
    for (int speed:POOL) check=lcm64(check,14LL*speed);
    require(check==D,"pool denominator changed");
    std::vector<i64> walls={0,D};
    walls.reserve(2+2*std::accumulate(POOL.begin(),POOL.end(),0));
    for (int speed:POOL) {
        const i64 unit=D/(14LL*speed);
        for (int tooth=0; tooth<speed; ++tooth) {
            walls.push_back((14LL*tooth+1)*unit);
            walls.push_back((14LL*tooth+13)*unit);
        }
    }
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    require(walls.size()==7134,"pool wall census changed");
    std::vector<Cell> result;
    result.reserve(walls.size()-1);
    for (std::size_t i=0;i+1<walls.size();++i) {
        u32 failed_mask=0;
        for (int bit=0;bit<30;++bit)
            if (!pool_safe_midpoint(POOL[bit],walls[i],walls[i+1]))
                failed_mask |= u32{1}<<bit;
        result.push_back({walls[i],walls[i+1],failed_mask});
    }
    return result;
}

struct Geometry {
    u32 repair=0;
    i64 width=0;
    std::vector<std::pair<i64,i64>> components;
};

inline Geometry build_geometry(u32 repair,const std::vector<Cell>& cells) {
    require(std::popcount(repair)==8,"repair arity changed");
    Geometry result;
    result.repair=repair;
    std::size_t i=0;
    while (i<cells.size()) {
        if ((cells[i].failed & ~repair)!=0) { ++i; continue; }
        const i64 left=cells[i].left;
        i64 right=cells[i].right;
        ++i;
        while (i<cells.size() && (cells[i].failed & ~repair)==0) {
            right=cells[i].right;
            ++i;
        }
        result.components.push_back({left,right});
        result.width += right-left;
    }
    require(result.width>0 && !result.components.empty(),"empty repair geometry");
    return result;
}

struct OrderedRepair { u64 key=0; u32 mask=0; };

inline std::vector<std::pair<int,int>> rectangle_pairs() {
    std::vector<std::pair<int,int>> result;
    result.reserve(PAIR_COUNT);
    for (int q=Q0;q<=Q1;++q) for (int r=R0;r<=R1;++r) {
        require(q<r,"rectangle order failure");
        result.push_back({q,r});
    }
    require(result.size()==PAIR_COUNT,"pair census changed");
    return result;
}

struct BodyScan {
    u64 bodies=0;
    u64 failures=0;
    u64 checks=0;
    u64 max_checks=0;
    u32 first_failure=std::numeric_limits<u32>::max();
    u32 worst_body=0;
    u32 worst_repair=0;
};

inline std::vector<u32> all_bodies() {
    std::vector<u32> result;
    result.reserve(BODY_COUNT);
    u32 body=(u32{1}<<9)-1;
    const u32 limit=u32{1}<<30;
    while (body!=0 && body<limit) {
        result.push_back(body);
        const u32 next=next_mask(body);
        if (next<=body) break;
        body=next;
    }
    require(result.size()==BODY_COUNT,"body census changed");
    return result;
}

inline BodyScan parallel_body_scan(const std::vector<u32>& deck,
                                   const std::vector<u32>& bodies) {
    const unsigned threads=std::max(1u,std::min(16u,std::thread::hardware_concurrency()));
    std::vector<BodyScan> local(threads);
    std::vector<std::thread> workers;
    for (unsigned lane=0;lane<threads;++lane) workers.emplace_back([&,lane]() {
        BodyScan report;
        const std::size_t begin=bodies.size()*lane/threads;
        const std::size_t end=bodies.size()*(lane+1)/threads;
        for (std::size_t index=begin;index<end;++index) {
            const u32 body=bodies[index];
            u64 used=0;
            u32 witness=0;
            for (u32 repair:deck) {
                ++used;
                if ((repair&body)==0) { witness=repair; break; }
            }
            ++report.bodies;
            report.checks+=used;
            if (witness==0) {
                ++report.failures;
                report.first_failure=std::min(report.first_failure,body);
            }
            if (used>report.max_checks ||
                (used==report.max_checks && body<report.worst_body)) {
                report.max_checks=used;
                report.worst_body=body;
                report.worst_repair=witness;
            }
        }
        local[lane]=report;
    });
    for (auto& worker:workers) worker.join();
    BodyScan result;
    for (const auto& report:local) {
        result.bodies+=report.bodies;
        result.failures+=report.failures;
        result.checks+=report.checks;
        result.first_failure=std::min(result.first_failure,report.first_failure);
        if (report.max_checks>result.max_checks ||
            (report.max_checks==result.max_checks && report.worst_body<result.worst_body)) {
            result.max_checks=report.max_checks;
            result.worst_body=report.worst_body;
            result.worst_repair=report.worst_repair;
        }
    }
    return result;
}

} // namespace lrc_rectangle

#endif
