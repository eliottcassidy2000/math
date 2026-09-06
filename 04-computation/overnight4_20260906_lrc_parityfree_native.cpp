// FINITE-EXACT independent native interval verification of the network head.
// Adapted all-parity head: primitive 1<=a<b<c<=63, each nonzero modulo three.
// The Cursor/literal engine below is unchanged from the inherited odd-head audit.
// Primary engine: all six literal sheet assignments, exact endpoints on
// denominator 42abc, and a simultaneous three-pointer contact scan.
// No carrier congruence, lattice-direction classifier, or raw roof formula
// is used to compute the projections.
//
// Definitions independently transcribed from the native sheet construction in
// 04-computation/lrc14_one_ray_overnight_hexagon_sep05.py.
// Reproduce: c++ -std=c++17 -O3 -DNDEBUG <this.cpp> -o /tmp/lrc14_literal_sep06
//           /tmp/lrc14_literal_sep06 <raw-head.tsv>
// Required argument: standalone raw producer TSV. Every eligible row is generated
// independently here; all three E columns, physical mass, and contact count agree.

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
    need(argc==2,"supply the raw producer head63 TSV");
    std::ifstream input(argv[1]);need(bool(input),"cannot open input TSV");
    std::string line;std::getline(input,line);
    if(!line.empty() && line.back()=='\r')line.pop_back();
    need(line=="a\tb\tc\tdenominator\tE0_numerator\tE1_numerator\tE2_numerator\tmass_numerator\traw_carriers","TSV header");
    std::map<V,std::array<I,6>> raw;
    while(std::getline(input,line)) {
      if(line.empty())continue;
      std::istringstream stream(line);V w{};std::array<I,6> record{};
      need(bool(stream>>w[0]>>w[1]>>w[2]),"speed row parse");
      for(auto& x:record)need(bool(stream>>x),"record parse");
      std::string excess;need(!(stream>>excess),"extra row column");
      need(raw.emplace(w,record).second,"duplicate raw row");
    }
    need(raw.size()==10074,"independent H63 input universe cardinality");
    std::cout<<"FINITE-EXACT all-parity native six-sheet interval audit H=63\n";
    controls();
    I rows=0,old_failures=0,empty=0,total_contacts=0;
    Leader maximum_min,maximum_mass;
    std::vector<V> min_equalities,mass_equalities;
    for(I c=1;c<=63;++c) {
     if(c%3==0)continue;
     for(I b=1;b<c;++b) {
      if(b%3==0)continue;
      for(I a=1;a<b;++a) {
       if(a%3==0 || std::gcd(a,std::gcd(b,c))!=1)continue;
       const V w{a,b,c}; const Native row=literal(w);
       auto iter=raw.find(w);need(iter!=raw.end(),"missing eligible raw row "+triple(w));
       const auto& expected=iter->second;
       need(row.denominator==expected[0],"denominator "+triple(w));
       for(int i=0;i<3;++i)need(row.projections[i]==expected[i+1],"native E mismatch "+triple(w));
       need(row.mass==expected[4],"native mass mismatch "+triple(w));
       need(row.contacts==3*expected[5],"native contact multiplicity "+triple(w));
       const I minimum=*std::min_element(row.projections.begin(),row.projections.end());
       need(row.mass<=minimum,"physical versus network order");
       need(Wide(55)*minimum<=Wide(6)*row.denominator,"selected ceiling");
       if(Wide(55)*minimum==Wide(6)*row.denominator)min_equalities.push_back(w);
       if(Wide(55)*row.mass==Wide(6)*row.denominator)mass_equalities.push_back(w);
       old_failures+=Wide(77)*row.mass>Wide(6)*row.denominator;
       empty+=row.contacts==0;total_contacts+=row.contacts;++rows;
       maximum_min.update(minimum,row,w);maximum_mass.update(row.mass,row,w);
      }
     }
    }
    need(rows==10074 && rows==I(raw.size()),"full independently enumerated H63 universe");
    const std::vector<V> equality{{1,10,11}};
    need(min_equalities==equality && mass_equalities==equality,"both exact equality loci");
    need(old_failures==151,"old parity ceiling hostile count");
    const Native zero=literal({1,2,4});
    need(zero.contacts==0 && zero.mass==0,"empty-defect hostile included");
    const Native gap=literal({4,7,11});
    need(Wide(gap.mass)*2156==Wide(215)*gap.denominator,"physical strict-gap control");
    need(Wide(*std::min_element(gap.projections.begin(),gap.projections.end()))*2156==Wide(223)*gap.denominator,"network strict-gap control");
    const Native parity=literal({2,5,7});
    const V pp{22,6,1},pq{245,49,10};
    for(int i=0;i<3;++i)need(Wide(parity.projections[i])*pq[i]==Wide(pp[i])*parity.denominator,"mixed parity projection control");
    std::cout<<"UNIVERSE rows="<<rows<<" empty="<<empty<<" native_contacts="<<total_contacts<<" old_physical_ceiling_failures="<<old_failures<<"\n";
    maximum_min.print("MAX_SELECTED");maximum_mass.print("MAX_PHYSICAL");
    std::cout<<"EQUALITIES both=(1,10,11)\n";
    std::cout<<"PASS gates="<<gates<<"; every raw E, physical mass and 3x contact count matches native intervals.\n";
    std::cout<<"SCOPE: finite H63 certificate; unbounded theorem requires the separate analytic reduction.\n";
 } catch(const std::exception& e) {std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}
 return 0;
}
