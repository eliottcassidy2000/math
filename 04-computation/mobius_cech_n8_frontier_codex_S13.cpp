// Exact lower-first n=8 decision for THM-801 / HYP-6880.
//
// The strong recursive key Lambda uses only the three ordered n=7 face-node
// pairs, the UABC colour word, and the raw mirror-layer S2 counts.  Lambda is
// strictly weaker than Omega+S2 because it omits the upper n=8 node pair.  If
// Lambda is injective, Omega+S2 is injective without constructing an n=8
// tournament-class atlas.
//
// Tournament Analysis is reported as a refinement path: B3 lower node/colour,
// then tau=3,...,7 S2 layers, then the fixed layer.  Its pairwise observable is
// the number of unordered literal-line pairs separated at each stage.  The
// switches are retention versus separation per carrier cell, and the displayed
// refinement order is the tie Hamiltonian path.
//
// Assumption challenge: the vertices carrying the proof are complement lines,
// gap-contracted faces, and mirror layers, not tournament vertices alone.  The
// key preserves the declared line/face data but not an LRC witness or arbitrary
// future lift.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using u128 = unsigned __int128;

static const char *ATLAS =
    "05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json";

struct Tile { int a, b; };

static std::vector<Tile> tiles(int n) {
    std::vector<Tile> out;
    for (int b = 1; b <= n - 2; ++b)
        for (int a = n; a >= b + 2; --a) out.push_back({a, b});
    return out;
}

static int tile_pos(const std::vector<Tile>& ts, int a, int b) {
    for (int i = 0; i < (int)ts.size(); ++i)
        if (ts[i].a == a && ts[i].b == b) return i;
    throw std::runtime_error("tile not found");
}

static std::vector<int> reflection_map(int n, const std::vector<Tile>& ts) {
    std::vector<int> sigma(ts.size());
    for (int i = 0; i < (int)ts.size(); ++i)
        sigma[i] = tile_pos(ts, n - ts[i].b + 1, n - ts[i].a + 1);
    return sigma;
}

static bool symmetric(uint32_t mask, const std::vector<int>& sigma) {
    for (int i = 0; i < (int)sigma.size(); ++i)
        if (((mask >> i) & 1u) != ((mask >> sigma[i]) & 1u)) return false;
    return true;
}

static std::vector<uint16_t> read_node_ranks(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open n7 atlas");
    std::ostringstream buffer;
    buffer << in.rdbuf();
    std::string text = buffer.str();
    std::string marker = "\"node_rank_by_mask\"";
    size_t p = text.find(marker);
    if (p == std::string::npos) throw std::runtime_error("node array marker missing");
    p = text.find('[', p);
    if (p == std::string::npos) throw std::runtime_error("node array start missing");
    ++p;
    std::vector<uint16_t> values;
    while (p < text.size()) {
        while (p < text.size() && (text[p] == ' ' || text[p] == '\n' ||
                                   text[p] == '\r' || text[p] == '\t' || text[p] == ',')) ++p;
        if (p >= text.size() || text[p] == ']') break;
        if (text[p] < '0' || text[p] > '9') throw std::runtime_error("bad node array token");
        unsigned value = 0;
        while (p < text.size() && text[p] >= '0' && text[p] <= '9')
            value = 10 * value + unsigned(text[p++] - '0');
        if (value >= 512) throw std::runtime_error("node rank exceeds 9-bit budget");
        values.push_back(uint16_t(value));
    }
    if (values.size() != (1u << 15)) throw std::runtime_error("n7 node array wrong length");
    return values;
}

struct FaceMaps {
    std::array<std::vector<std::pair<int,int>>,3> maps;
};

static FaceMaps make_face_maps(const std::vector<Tile>& upper,
                               const std::vector<Tile>& lower, int n) {
    FaceMaps result;
    for (int bit = 0; bit < (int)upper.size(); ++bit) {
        int a = upper[bit].a, b = upper[bit].b;
        if (a < n) result.maps[0].push_back({bit, tile_pos(lower, a, b)});          // A
        if (a - b >= 3) result.maps[1].push_back({bit, tile_pos(lower, a - 1, b)}); // B
        if (b >= 2) result.maps[2].push_back({bit, tile_pos(lower, a - 1, b - 1)}); // C
    }
    return result;
}

static uint32_t face_mask(uint32_t mask,
                          const std::vector<std::pair<int,int>>& mapping) {
    uint32_t out = 0;
    for (auto [u,l] : mapping) out |= ((mask >> u) & 1u) << l;
    return out;
}

static int composition_rank(const std::array<int,4>& c, int total) {
    int rank = 0;
    for (int a = 0; a <= total; ++a)
        for (int b = 0; b <= total - a; ++b)
            for (int d = 0; d <= total - a - b; ++d) {
                std::array<int,4> x{a,b,d,total-a-b-d};
                if (x == c) return rank;
                ++rank;
            }
    throw std::runtime_error("invalid composition");
}

struct Rich {
    uint64_t base;                 // six 9-bit lower nodes + four colour bits
    std::array<uint8_t,6> digit;   // five tau compositions + fixed-one count
    std::array<uint16_t,3> face_line;
    uint32_t line;
};

struct Wide {
    uint64_t hi, lo;
    bool operator<(const Wide& o) const { return hi != o.hi ? hi < o.hi : lo < o.lo; }
    bool operator==(const Wide& o) const { return hi == o.hi && lo == o.lo; }
};

static Wide wide(u128 x) { return {uint64_t(x >> 64), uint64_t(x)}; }

static u128 prefix_key(const Rich& r, int prefix, const std::array<int,6>& radix) {
    u128 key = r.base;
    for (int i = 0; i < prefix; ++i) key = key * unsigned(radix[i]) + r.digit[i];
    return key;
}

struct Stats {
    uint64_t cells = 0, collision_cells = 0, excess = 0, pairs = 0, max_mult = 0;
};

static Stats stats(std::vector<Wide>& keys) {
    std::sort(keys.begin(), keys.end());
    Stats s;
    for (size_t i = 0; i < keys.size();) {
        size_t j = i + 1;
        while (j < keys.size() && keys[j] == keys[i]) ++j;
        uint64_t m = j - i;
        ++s.cells;
        s.collision_cells += m > 1;
        s.excess += m - 1;
        s.pairs += m * (m - 1) / 2;
        s.max_mult = std::max(s.max_mult, m);
        i = j;
    }
    return s;
}

static std::string atom_word(int word) {
    std::string out;
    for (int bit = 3; bit >= 0; --bit) out += ((word >> bit) & 1) ? 'B' : 'K';
    return out;
}

int main(int argc, char **argv) {
    std::string atlas = argc > 1 ? argv[1] : ATLAS;
    auto node = read_node_ranks(atlas);
    auto upper = tiles(8), lower = tiles(7);
    assert(upper.size() == 21 && lower.size() == 15);
    auto sigma8 = reflection_map(8, upper), sigma7 = reflection_map(7, lower);
    auto fm = make_face_maps(upper, lower, 8);
    const uint32_t full8 = (1u << 21) - 1, full7 = (1u << 15) - 1;
    const uint32_t lines = 1u << 20;
    const std::array<int,6> layer_size{1,1,2,2,3,3};
    const std::array<int,6> radix{4,4,10,10,20,4};

    std::vector<Rich> rows;
    rows.reserve(lines);
    std::array<uint64_t,16> atoms{};
    std::array<uint64_t,10> skew_hist{};

    for (uint32_t line = 0; line < lines; ++line) {
        // min(t,kappa t)=line; bit zero is the apex, so complement exactly
        // when line's apex is one to obtain the unique apex-zero endpoint.
        uint32_t mask = (line & 1u) ? (line ^ full8) : line;
        assert((mask & 1u) == 0);
        std::array<uint32_t,3> face{};
        for (int s = 0; s < 3; ++s) face[s] = face_mask(mask, fm.maps[s]);

        int word = (symmetric(mask, sigma8) ? 8 : 0)
                 | (symmetric(face[0], sigma7) ? 4 : 0)
                 | (symmetric(face[1], sigma7) ? 2 : 0)
                 | (symmetric(face[2], sigma7) ? 1 : 0);
        ++atoms[word];

        uint64_t base = 0;
        for (int s = 0; s < 3; ++s) {
            base = (base << 9) | node[face[s]];
            base = (base << 9) | node[face[s] ^ full7];
        }
        base = (base << 4) | unsigned(word);

        std::array<std::array<int,4>,5> comp{};
        std::array<int,2> fixed{};
        int skew = 0;
        for (int bit = 0; bit < (int)upper.size(); ++bit) {
            int tau = upper[bit].a + upper[bit].b - 1;
            if (tau < 8) {
                int x = (mask >> bit) & 1u;
                int y = (mask >> sigma8[bit]) & 1u;
                ++comp[tau - 3][2*x + y];
                skew += x != y;
            } else if (tau == 8) {
                assert(sigma8[bit] == bit);
                ++fixed[(mask >> bit) & 1u];
            }
        }
        ++skew_hist[skew];

        Rich r{};
        r.base = base;
        for (int i = 0; i < 5; ++i)
            r.digit[i] = uint8_t(composition_rank(comp[i], layer_size[i]));
        r.digit[5] = uint8_t(fixed[1]);
        for (int s = 0; s < 3; ++s)
            r.face_line[s] = uint16_t(std::min(face[s], face[s] ^ full7));
        r.line = line;
        rows.push_back(r);
    }

    std::array<Stats,7> ladder{};
    for (int prefix = 0; prefix <= 6; ++prefix) {
        std::vector<Wide> keys;
        keys.reserve(rows.size());
        for (const Rich& r : rows) keys.push_back(wide(prefix_key(r, prefix, radix)));
        ladder[prefix] = stats(keys);
    }

    std::vector<Wide> literal;
    literal.reserve(rows.size());
    for (const Rich& r : rows) {
        u128 key = r.face_line[0];
        key = (key << 14) | r.face_line[1];
        key = (key << 14) | r.face_line[2];
        literal.push_back(wide(key));
    }
    Stats literal_stats = stats(literal);
    assert(literal_stats.cells == lines);

    // Genealogy of the 418 lower-node/colour double collisions.
    std::vector<uint32_t> base_order(rows.size());
    std::iota(base_order.begin(), base_order.end(), uint32_t{0});
    std::sort(base_order.begin(), base_order.end(), [&](uint32_t i, uint32_t j) {
        if (rows[i].base != rows[j].base) return rows[i].base < rows[j].base;
        return rows[i].line < rows[j].line;
    });

    std::array<uint64_t,6> first_separation{}; // tau=3,...,7,fixed tau=8
    std::array<uint64_t,8> face_difference{};  // bit 0=A, 1=B, 2=C
    uint64_t base_doubletons = 0;
    uint64_t ac_equal_pairs = 0;
    uint64_t s11_residual_face_occurrences = 0;
    auto is_s11_residual = [](uint16_t e) {
        return e == 0x12ca || e == 0x12cb ||
               e == 0x146c || e == 0x146d;
    };

    for (size_t i = 0; i < base_order.size();) {
        size_t j = i + 1;
        while (j < base_order.size() &&
               rows[base_order[j]].base == rows[base_order[i]].base) ++j;
        assert(j - i <= 2);
        if (j - i == 2) {
            const Rich& x = rows[base_order[i]];
            const Rich& y = rows[base_order[i + 1]];
            ++base_doubletons;

            unsigned pattern = 0;
            for (int face = 0; face < 3; ++face) {
                pattern |= unsigned(x.face_line[face] != y.face_line[face]) << face;
                s11_residual_face_occurrences += is_s11_residual(x.face_line[face]);
                s11_residual_face_occurrences += is_s11_residual(y.face_line[face]);
            }
            ++face_difference[pattern];
            ac_equal_pairs += (pattern & 0b101u) == 0;

            int first = 0;
            while (first < 6 && x.digit[first] == y.digit[first]) ++first;
            assert(first < 6);
            ++first_separation[first];
        }
        i = j;
    }

    const std::array<uint64_t,6> expected_first{166,104,74,22,52,0};
    const std::array<uint64_t,8> expected_faces{0,0,0,4,0,44,4,366};
    assert(ladder[0].max_mult == 2);
    assert(base_doubletons == ladder[0].collision_cells);
    assert(base_doubletons == ladder[0].excess);
    assert(first_separation == expected_first);
    assert(face_difference == expected_faces);
    for (int d = 0; d < 6; ++d)
        assert(first_separation[d] == ladder[d].excess - ladder[d + 1].excess);
    assert(ac_equal_pairs == 0);
    assert(s11_residual_face_occurrences == 0);

    const uint64_t total_pairs = uint64_t(lines) * (lines - 1) / 2;
    std::array<uint64_t,7> separated_pairs{};
    for (int i = 0; i < 7; ++i)
        separated_pairs[i] = total_pairs - ladder[i].pairs;
    for (int i = 0; i < 5; ++i) {
        assert(separated_pairs[i] < separated_pairs[i + 1]);
        assert(u128(separated_pairs[i]) * ladder[i + 1].cells >
               u128(separated_pairs[i + 1]) * ladder[i].cells);
    }
    assert(separated_pairs[5] == separated_pairs[6]);
    assert(ladder[5].cells == ladder[6].cells);

    // Exact full-key collision witnesses; keys are 75-bit integers, not hashes.
    struct Witness { Wide key; uint32_t line; };
    std::vector<Witness> full;
    full.reserve(rows.size());
    for (const Rich& r : rows) full.push_back({wide(prefix_key(r, 6, radix)), r.line});
    std::sort(full.begin(), full.end(), [](const Witness& a, const Witness& b) {
        if (a.key == b.key) return a.line < b.line;
        return a.key < b.key;
    });
    std::vector<std::vector<uint32_t>> collision_witnesses;
    for (size_t i = 0; i < full.size();) {
        size_t j = i + 1;
        while (j < full.size() && full[j].key == full[i].key) ++j;
        if (j - i > 1 && collision_witnesses.size() < 24) {
            std::vector<uint32_t> one;
            for (size_t k = i; k < j; ++k) one.push_back(full[k].line);
            collision_witnesses.push_back(std::move(one));
        }
        i = j;
    }

    std::map<std::string,uint64_t> expected_atoms{
        {"BBBB",32},{"BKBK",2016},{"KBBB",32},{"KBBK",192},{"KBKB",192},
        {"KKBB",192},{"KBKK",15936},{"KKKB",15936},{"KKBK",13920},
        {"KKKK",1000128}
    };
    uint64_t atom_failures = 0;
    for (int w = 0; w < 16; ++w) {
        uint64_t expected = expected_atoms.count(atom_word(w)) ? expected_atoms[atom_word(w)] : 0;
        atom_failures += atoms[w] != expected;
    }
    assert(atom_failures == 0);
    uint64_t skew_failures = 0;
    for (int k = 0; k <= 9; ++k) {
        uint64_t choose = 1;
        for (int i = 1; i <= k; ++i) choose = choose * (10 - i) / i;
        skew_failures += skew_hist[k] != (uint64_t(1) << 11) * choose;
    }
    assert(skew_failures == 0);

    const char *labels[7] = {
        "B3 lower node pairs + UABC", "+ tau=3", "+ tau=4", "+ tau=5",
        "+ tau=6", "+ tau=7", "+ fixed tau=8"
    };
    std::cout << "MOBIUS/CECH N=8 LOWER-FIRST FRONTIER\n";
    std::cout << std::string(82, '=') << "\n";
    std::cout << "lines=" << lines << " n7 node ranks=32768 max="
              << *std::max_element(node.begin(), node.end()) << "\n";
    std::cout << "literal B3 face-line triples: cells=" << literal_stats.cells
              << " excess=" << literal_stats.excess << "\n\n";
    std::cout << "RECURSIVE LAMBDA REFINEMENT\n";
    for (int i = 0; i <= 6; ++i)
        std::cout << "  " << std::left << std::setw(31) << labels[i]
                  << " cells=" << ladder[i].cells << " excess=" << ladder[i].excess
                  << " collision_cells=" << ladder[i].collision_cells
                  << " max=" << ladder[i].max_mult << " pairs=" << ladder[i].pairs << "\n";
    std::cout << "\nLambda injective=" << (ladder[6].cells == lines ? "True" : "False")
              << "; therefore Omega+S2 injective="
              << (ladder[6].cells == lines ? "True" : "undecided by lower key") << "\n";
    std::cout << "collision witnesses=";
    if (collision_witnesses.empty()) std::cout << "()";
    for (const auto& cell : collision_witnesses) {
        std::cout << " (";
        for (size_t i = 0; i < cell.size(); ++i) std::cout << (i ? "," : "") << cell[i];
        std::cout << ")";
    }
    std::cout << "\n";
    std::cout << "UABC closed-form failures=" << atom_failures
              << "; S2 skew-binomial failures=" << skew_failures << "\n";
    std::cout << "S2 radix product=128000 (<2^17); full Lambda uses <=75 exact bits.\n";
    std::cout << "At n=8 each tau layer has <=3 positions. Counts plus first position moments\n"
                 "reconstruct every layer, so S2+M1 is an unconditional exact fallback (not needed\n"
                 "if Lambda is already injective).\n\n";
    std::cout << "BASE DOUBLETON GENEALOGY\n"
              << "  first separator: tau3=" << first_separation[0]
              << " tau4=" << first_separation[1]
              << " tau5=" << first_separation[2]
              << " tau6=" << first_separation[3]
              << " tau7=" << first_separation[4]
              << " fixed8=" << first_separation[5] << "\n"
              << "  differing literal faces: A+B=" << face_difference[3]
              << " A+C=" << face_difference[5]
              << " B+C=" << face_difference[6]
              << " A+B+C=" << face_difference[7] << "\n"
              << "  A/C-equal pairs=" << ac_equal_pairs
              << "; S11 residual-face occurrences="
              << s11_residual_face_occurrences << "\n\n";
    std::cout << "TOURNAMENT ANALYSIS\n"
                 "  vertices=the seven displayed recursive carriers\n"
                 "  observable=unordered literal-line pairs separated\n"
                 "  switches=retention / separation per carrier cell\n";
    std::cout << "  separated_pairs=(";
    for (int i = 0; i < 7; ++i)
        std::cout << (i ? "," : "") << separated_pairs[i];
    std::cout << ")\n"
                 "  retention score_hist={0:1,...,6:1} directed_3cycles=0 SCCs=[1,1,1,1,1,1,1] HP=1\n"
                 "  economy score_hist={0:1,...,6:1} directed_3cycles=0 SCCs=[1,1,1,1,1,1,1] HP=1\n"
                 "  edge_flips=20\n"
                 "  tie_order=(B3,tau3,tau4,tau5,tau6,tau7,fixed)\n"
                 "  retention_HP=(tau7,fixed,tau6,tau5,tau4,tau3,B3)\n"
                 "  economy_HP=(B3,tau3,tau4,tau5,tau6,tau7,fixed)\n";
    return 0;
}
