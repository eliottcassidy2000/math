// FINITE-EXACT independent native interval verification of the network head.
// All-parity nonadditive head: primitive 1<=a<b<c<=H, each nonzero mod three.
// The Cursor/literal engine below is unchanged from the inherited odd-head audit.
// Primary engine: all six literal sheet assignments, exact endpoints on
// denominator 42abc, and a simultaneous three-pointer contact scan.
// No carrier congruence, lattice-direction classifier, or raw roof formula
// is used to compute the projections.
//
// Definitions independently transcribed from the native sheet construction in
// 04-computation/lrc14_one_ray_overnight_hexagon_sep05.py.
// Reproduce: c++ -std=c++17 -O3 -DNDEBUG <this.cpp> -o /tmp/lrc14_literal_sep06
//           /tmp/lrc14_literal_sep06 535
// Fifth-wave native nonadditive head. No carrier engine or math producer imported.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <sstream>

using I = std::int64_t;
using V = std::array<I, 3>;
using Wide = __int128_t;

static I gates=0;
static void need(bool value, const std::string& label) {
    ++gates;
    if (!value) throw std::runtime_error(label);
}

static std::string fraction(I p, I q) {
    const I g = std::gcd(p, q);
    return std::to_string(p / g) + "/" + std::to_string(q / g);
}

static std::string triple(const V& w) {
    return std::to_string(w[0]) + "," + std::to_string(w[1]) + "," + std::to_string(w[2]);
}

struct Cursor {
    I left, right, width;
    I step, radius, denominator;
    int index, speed, color;
    bool valid;

    Cursor(int speed_, int color_, I multiplier, I denominator_)
        : step(42 * multiplier), radius(3 * multiplier), denominator(denominator_),
          index(0), speed(speed_), color(color_), valid(true) {
        // Sorted residues (3*owner-speed*color) mod(3*speed) are precisely
        // r0, r0+3, ..., r0+3*(speed-1), r0=(-speed*color) mod3.
        if (color == 0) {
            left = 0;
            right = radius;
        } else {
            const I residue = (3 - (speed * color) % 3) % 3;
            const I center = 14 * residue * multiplier;
            left = center - radius;
            right = center + radius;
        }
        width = right - left;
    }

    void advance() {
        ++index;
        if (color == 0) {
            if (index > speed) {
                valid = false;
                return;
            }
            if (index == speed) {
                left = denominator - radius;
                right = denominator;
            } else {
                left = index * step - radius;
                right = index * step + radius;
            }
        } else {
            if (index == speed) {
                valid = false;
                return;
            }
            left += step;
            right += step;
        }
        width = right - left;
    }
};

struct Native {
    I denominator;
    V projections{};
    I mass = 0;
    I contacts = 0;
};

static Native literal(const V& w) {
    Native result;
    result.denominator = 42 * w[0] * w[1] * w[2];
    constexpr int permutations[6][3] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2},
        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    };
    for (const auto& assignment : permutations) {
        Cursor x(static_cast<int>(w[0]), assignment[0], result.denominator / (42 * w[0]), result.denominator);
        Cursor y(static_cast<int>(w[1]), assignment[1], result.denominator / (42 * w[1]), result.denominator);
        Cursor z(static_cast<int>(w[2]), assignment[2], result.denominator / (42 * w[2]), result.denominator);
        while (x.valid && y.valid && z.valid) {
            const I left = std::max(x.left, std::max(y.left, z.left));
            const I right = std::min(x.right, std::min(y.right, z.right));
            if (left < right) {
                // All three intervals meet. For each omitted sheet, use the
                // complete intersection length of the other two intervals,
                // capped by the complete omitted interval length.
                const I yz = std::min(y.right, z.right) - std::max(y.left, z.left);
                const I xz = std::min(x.right, z.right) - std::max(x.left, z.left);
                const I xy = std::min(x.right, y.right) - std::max(x.left, y.left);
                result.projections[0] += std::min(yz, x.width);
                result.projections[1] += std::min(xz, y.width);
                result.projections[2] += std::min(xy, z.width);
                result.mass += right - left;
                ++result.contacts;
            }
            // At least one interval ends. Ties advance simultaneously and
            // zero-length contacts never contribute.
            if (x.right == right) x.advance();
            if (y.right == right) y.advance();
            if (z.right == right) z.advance();
        }
    }
    return result;
}

struct Leader {
    I numerator = 0, denominator = 1;
    V speeds{};
    V projections{};
    I contacts = 0;

    void update(I numerator_, const Native& row, const V& w) {
        if (Wide(numerator_) * denominator > Wide(numerator) * row.denominator) {
            numerator = numerator_;
            denominator = row.denominator;
            speeds = w;
            projections = row.projections;
            contacts = row.contacts;
        }
    }

    void print(const std::string& label) const {
        std::cout << label << " " << fraction(numerator, denominator)
                  << " w=" << triple(speeds) << " E="
                  << fraction(projections[0], denominator) << ","
                  << fraction(projections[1], denominator) << ","
                  << fraction(projections[2], denominator)
                  << " native_contacts=" << contacts << "\n";
    }
};

static void controls() {
    struct Control { V w; V p; V q; I mp, mq, contacts; };
    const Control cases[] = {
        {{1, 5, 7}, {8, 6, 4}, {245, 49, 35}, 8, 245, 2},
        {{1, 5, 11}, {6, 6, 6}, {77, 77, 77}, 6, 77, 2},
        {{1, 19, 23}, {12, 12, 12}, {437, 161, 161}, 12, 437, 4},
        {{17, 23, 25}, {106, 12, 2546}, {4025, 425, 68425}, 744, 68425, 4},
        {{19, 23, 29}, {156, 192, 3840}, {4669, 3857, 88711}, 154, 12673, 6},
        {{1, 1201, 1205}, {310116, 516, 516}, {10130435, 8435, 8435}, 310116, 10130435, 172},
        {{1, 599, 607}, {38940, 132, 132}, {2545151, 4249, 4249}, 38940, 2545151, 44},
        {{5, 1001, 1003}, {156680, 122, 848864}, {7028021, 5015, 35140105}, 767384, 35140105, 58}
    };
    for (const auto& c : cases) {
        const Native row = literal(c.w);
        for (int i = 0; i < 3; ++i)
            need(Wide(row.projections[i]) * c.q[i] == Wide(c.p[i]) * row.denominator,
                 "known literal/raw control mismatch " + triple(c.w));
        need(Wide(row.mass) * c.mq == Wide(c.mp) * row.denominator,
             "known physical mass mismatch " + triple(c.w));
        need(row.contacts == 3 * c.contacts, "native/raw contact multiplicity mismatch " + triple(c.w) + " native=" + std::to_string(row.contacts));
        std::cout << "CONTROL w=" << triple(c.w) << " E="
                  << fraction(row.projections[0], row.denominator) << ","
                  << fraction(row.projections[1], row.denominator) << ","
                  << fraction(row.projections[2], row.denominator)
                  << " mass=" << fraction(row.mass, row.denominator)
                  << " native_contacts=" << row.contacts << "\n";
    }
}



int main(int argc,char** argv) {
 try {
  const I H=argc>1?std::stoi(argv[1]):535;
  need(H>=11 && H<=535,"H bound11..535");
  controls();std::cout.flush();
  I universe=0,additive=0,rows=0,empty=0,contacts=0,fail=0;
  I outside_low_rows=0,outside_low_fail=0;
  Leader max_min,max_mass,max_outside;
  std::vector<V> min_equalities,mass_equalities;
  std::uint64_t digest=14695981039346656037ULL;
  auto hash=[&](I value) {auto x=std::uint64_t(value);for(int i=0;i<8;++i){digest^=x&255U;digest*=1099511628211ULL;x>>=8;}};
  for(I c=1;c<=H;++c) {
   if(c%3==0)continue;
   for(I b=1;b<c;++b) {
    if(b%3==0)continue;
    for(I a=1;a<b;++a) {
     if(a%3==0 || std::gcd(a,std::gcd(b,c))!=1)continue;
     ++universe;
     if(a+b==c) {++additive;continue;}
     const V w{a,b,c};const Native row=literal(w);++rows;
     const I minimum=*std::min_element(row.projections.begin(),row.projections.end());
     const I maximum=*std::max_element(row.projections.begin(),row.projections.end());
     max_min.update(minimum,row,w);max_mass.update(row.mass,row,w);
     if(Wide(140)*minimum>Wide(11)*row.denominator) {++fail;std::cout<<"FAILURE "<<triple(w)<<" min="<<fraction(minimum,row.denominator)<<"\n";}
     if(Wide(140)*minimum==Wide(11)*row.denominator)min_equalities.push_back(w);
     if(Wide(140)*row.mass==Wide(11)*row.denominator)mass_equalities.push_back(w);
     const bool norm4=(c==2*a+b || c==a+2*b || 2*b==a+c);
     const bool norm5=(c==2*(a+b) || 2*c==a+2*b || 2*c==2*a+b || 2*b==2*a+c);
     if(!norm4 && !norm5) {++outside_low_rows;max_outside.update(maximum,row,w);outside_low_fail+=Wide(140)*maximum>=Wide(11)*row.denominator;}
     empty+=row.contacts==0;contacts+=row.contacts;
     for(I x:w)hash(x);hash(row.denominator);for(I x:row.projections)hash(x);hash(row.mass);hash(row.contacts);
    }
   }
   if(c%50==1) {std::cerr<<"PROGRESS c="<<c<<" nonadditive="<<rows<<" failures="<<fail<<"\n";}
  }
  need(universe==rows+additive,"partition");need(fail==0,"selected nonadditive finite ceiling");
  const std::vector<V> expected=H>=20?std::vector<V>{{2,11,20}}:std::vector<V>{};
  need(min_equalities==expected && mass_equalities==expected,"both equality loci");
  std::cout<<"UNIVERSE H="<<H<<" primitive_unit="<<universe<<" additive="<<additive<<" nonadditive="<<rows<<" empty="<<empty<<" contacts="<<contacts<<"\n";
  max_min.print("MAX_NONADDITIVE_MIN");max_mass.print("MAX_NONADDITIVE_PHYSICAL");max_outside.print("MAX_ALL_E_OUTSIDE_NORM4_NORM5");
  std::cout<<"OUTSIDE_LOW rows="<<outside_low_rows<<" all_E_ge_target="<<outside_low_fail<<"\n";
  std::cout<<"EQUALITIES selected="<<min_equalities.size()<<" physical="<<mass_equalities.size()<<" both(2,11,20)whenpresent\n";
  std::cout<<"FNV1A64 "<<digest<<"\n";
  std::cout<<"FINITE-EXACT PASS; independent raw replay and all-height reduction still required.\n";
 }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}
 return 0;
}
