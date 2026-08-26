#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Standalone independent exact verifier for the THM-4207
// (P,{50,51}) depth-eight/depth-nine claim.
//
// Independence from the maintained primary closure implementation:
//   * geometry is an endpoint-event overlay, initialized from the cyclic
//     danger state; it never classifies a midpoint and never imports the
//     original source;
//   * deletion masses are accumulated atom -> fixed-cardinality supersets
//     into combinadic arrays, rather than deletion -> failure submasks;
//   * nine-bodies are recursively generated and E8 edges are ordered by exact
//     mass, rather than Gosper-generated and SplitMix ordered;
//   * the 35 E9 edges are frozen data, not greedily rediscovered.

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> P = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr int Q = 50;
constexpr int R = 51;
constexpr u32 FULL = (u32{1} << 30) - 1;
constexpr u64 FNV_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr u64 FNV_PRIME = UINT64_C(0x100000001b3);

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string out;
    while (value != 0) {
        out.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex;
    out.width(16);
    out.fill('0');
    out << value;
    return out.str();
}

class Ledger {
  public:
    void add(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value_ ^= static_cast<std::uint8_t>((word >> shift) & 0xffu);
            value_ *= FNV_PRIME;
        }
    }
    u64 value() const { return value_; }
  private:
    u64 value_ = FNV_OFFSET;
};

i64 exact_lcm(i64 a, i64 b) {
    const i64 divisor = std::gcd(a, b);
    const i128 value = static_cast<i128>(a / divisor) * b;
    require(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
}

std::array<std::array<u64, 10>, 31> make_binomial() {
    std::array<std::array<u64, 10>, 31> c{};
    for (int n = 0; n <= 30; ++n) {
        c[n][0] = 1;
        for (int k = 1; k <= 9; ++k) {
            c[n][k] = (n == 0 ? 0 : c[n - 1][k] + c[n - 1][k - 1]);
        }
    }
    return c;
}

const auto BINOM = make_binomial();

std::size_t colex_rank(u32 mask, int cardinality) {
    u64 rank = 0;
    int ordinal = 1;
    while (mask != 0) {
        const int vertex = std::countr_zero(mask);
        rank += BINOM[vertex][ordinal];
        ++ordinal;
        mask &= mask - 1;
    }
    require(ordinal == cardinality + 1, "combinadic cardinality mismatch");
    require(rank < BINOM[30][cardinality], "combinadic rank escaped range");
    return static_cast<std::size_t>(rank);
}

template <class Callback>
void choose_recursive(int next, int need, u32 mask, Callback& callback) {
    if (need == 0) {
        callback(mask);
        return;
    }
    for (int vertex = next; vertex <= 30 - need; ++vertex) {
        choose_recursive(vertex + 1, need - 1,
                         mask | (u32{1} << vertex), callback);
    }
}

template <class Callback>
void for_each_subset(int cardinality, Callback&& callback) {
    auto local = std::forward<Callback>(callback);
    choose_recursive(0, cardinality, 0, local);
}

struct Event {
    i64 position;
    std::uint8_t speed;
    bool becomes_safe;
};

struct Overlay {
    i64 denominator = 0;
    std::size_t event_count = 0;
    std::size_t event_positions = 0;
    std::size_t all_intervals = 0;
    std::size_t pair_safe_intervals = 0;
    i64 pair_safe_mass = 0;
    std::map<u32, i64> atom_mass;
    u64 ledger = FNV_OFFSET;
};

Overlay build_event_overlay() {
    std::array<int, 32> speeds{};
    std::copy(P.begin(), P.end(), speeds.begin());
    speeds[30] = Q;
    speeds[31] = R;
    i64 denominator = 1;
    for (int speed : speeds) denominator = exact_lcm(denominator, 14LL * speed);

    std::vector<Event> events;
    for (int index = 0; index < 32; ++index) {
        const int speed = speeds[index];
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            // Immediately right of (14k+1)/(14s), the danger tooth has ended.
            events.push_back({(14LL * tooth + 1) * unit,
                              static_cast<std::uint8_t>(index), true});
            // Immediately right of (14k+13)/(14s), the next danger tooth begins.
            events.push_back({(14LL * tooth + 13) * unit,
                              static_cast<std::uint8_t>(index), false});
        }
    }
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        if (a.position != b.position) return a.position < b.position;
        if (a.speed != b.speed) return a.speed < b.speed;
        return a.becomes_safe < b.becomes_safe;
    });

    // Every speed is dangerous immediately to the right of zero because its
    // zero-centered danger tooth wraps around the circle.
    std::array<bool, 32> safe{};
    Overlay out;
    out.denominator = denominator;
    out.event_count = events.size();
    i64 previous = 0;
    std::size_t index = 0;
    auto retain_interval = [&](i64 right) {
        require(right >= previous, "overlay order failure");
        if (right == previous) return;
        ++out.all_intervals;
        if (safe[30] && safe[31]) {
            u32 failed = 0;
            for (int vertex = 0; vertex < 30; ++vertex) {
                if (!safe[vertex]) failed |= u32{1} << vertex;
            }
            const i64 length = right - previous;
            out.atom_mass[failed] += length;
            out.pair_safe_mass += length;
            ++out.pair_safe_intervals;
        }
    };
    while (index < events.size()) {
        const i64 position = events[index].position;
        require(position > 0 && position < denominator,
                "unexpected endpoint-event position");
        retain_interval(position);
        ++out.event_positions;
        std::array<unsigned char, 32> touched{};
        while (index < events.size() && events[index].position == position) {
            const Event event = events[index++];
            require(touched[event.speed] == 0,
                    "one speed changed state twice at one wall");
            touched[event.speed] = 1;
            safe[event.speed] = event.becomes_safe;
        }
        previous = position;
    }
    retain_interval(denominator);

    Ledger digest;
    digest.add(static_cast<u64>(out.denominator));
    digest.add(out.event_count);
    digest.add(out.event_positions);
    digest.add(out.all_intervals);
    digest.add(out.pair_safe_intervals);
    digest.add(static_cast<u64>(out.pair_safe_mass));
    for (const auto& [mask, mass] : out.atom_mass) {
        digest.add(mask);
        digest.add(static_cast<u64>(mass));
    }
    out.ledger = digest.value();
    return out;
}

void add_supersets(u32 available, int need, u32 current, i64 length,
                   int cardinality, std::vector<i64>& masses) {
    if (need == 0) {
        masses[colex_rank(current, cardinality)] += length;
        return;
    }
    while (std::popcount(available) >= need) {
        const u32 bit = available & (~available + 1u);
        available ^= bit;
        add_supersets(available, need - 1, current | bit, length,
                      cardinality, masses);
    }
}

std::vector<i64> build_cardinality_masses(const Overlay& overlay,
                                          int cardinality) {
    std::vector<i64> masses(BINOM[30][cardinality], 0);
    for (const auto& [failed, length] : overlay.atom_mass) {
        const int used = std::popcount(failed);
        if (used > cardinality) continue;
        add_supersets(FULL ^ failed, cardinality - used, failed, length,
                      cardinality, masses);
    }
    return masses;
}

i64 direct_deletion_mass(u32 deletion, const Overlay& overlay) {
    i64 mass = 0;
    for (const auto& [failed, length] : overlay.atom_mass) {
        if ((failed & ~deletion) == 0) mass += length;
    }
    return mass;
}

std::string labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((mask & (u32{1} << vertex)) == 0) continue;
        if (!first) out << ',';
        first = false;
        out << P[vertex];
    }
    return out.str();
}

u32 label_mask(const std::array<int, 9>& values) {
    u32 mask = 0;
    for (int value : values) {
        const auto found = std::find(P.begin(), P.end(), value);
        require(found != P.end(), "frozen separator label absent from P");
        mask |= u32{1} << std::distance(P.begin(), found);
    }
    require(std::popcount(mask) == 9, "frozen separator is not a nine-set");
    return mask;
}

struct Edge8 {
    u32 mask;
    i64 mass;
};

struct LayerSummary {
    u64 edges = 0;
    u64 equalities = 0;
    u64 ledger = FNV_OFFSET;
};

LayerSummary summarize_layer(const std::vector<i64>& masses,
                             i64 denominator) {
    LayerSummary out;
    Ledger digest;
    for (std::size_t rank = 0; rank < masses.size(); ++rank) {
        const i128 delta = static_cast<i128>(63) * masses[rank] -
                           static_cast<i128>(4) * denominator;
        if (delta == 0) ++out.equalities;
        if (delta >= 0) {
            ++out.edges;
            digest.add(rank);
            digest.add(static_cast<u64>(masses[rank]));
        }
    }
    out.ledger = digest.value();
    return out;
}

struct BodyAudit {
    u64 candidates = 0;
    u64 checks = 0;
    u64 monotone_extensions = 0;
    int minimum_extension_choices = 30;
    std::vector<u32> covers;
    u64 cover_ledger = FNV_OFFSET;
};

BodyAudit audit_all_nine_bodies(std::vector<Edge8> edges,
                                const std::vector<i64>& mass8,
                                const std::vector<i64>& mass9,
                                i64 denominator) {
    std::sort(edges.begin(), edges.end(), [](const Edge8& a, const Edge8& b) {
        if (a.mass != b.mass) return a.mass > b.mass;
        return a.mask < b.mask;
    });
    BodyAudit out;
    Ledger cover_digest;
    auto inspect = [&](u32 body) {
        ++out.candidates;
        const Edge8* missed = nullptr;
        for (const Edge8& edge : edges) {
            ++out.checks;
            if ((edge.mask & body) == 0) {
                missed = &edge;
                break;
            }
        }
        if (missed == nullptr) {
            out.covers.push_back(body);
            cover_digest.add(body);
            return;
        }

        // R is an E8 repair disjoint from B. There are exactly 30-8-9=13
        // pool labels outside B union R. Adding any one gives a disjoint E9
        // repair, and deletion monotonicity can only increase safe mass.
        const u32 spare = FULL & ~(body | missed->mask);
        const int choices = std::popcount(spare);
        require(choices == 13, "monotone extension choice count changed");
        out.minimum_extension_choices =
            std::min(out.minimum_extension_choices, choices);
        const u32 extension = missed->mask | (spare & (~spare + 1u));
        const i64 extension_mass = mass9[colex_rank(extension, 9)];
        require(extension_mass >= missed->mass,
                "deletion monotonicity failed in reconstructed masses");
        require(static_cast<i128>(63) * extension_mass >=
                    static_cast<i128>(4) * denominator,
                "deterministic E8-to-E9 extension is not lawful");
        require(mass8[colex_rank(missed->mask, 8)] == missed->mass,
                "E8 mass lookup changed during body audit");
        ++out.monotone_extensions;
    };
    for_each_subset(9, inspect);
    out.cover_ledger = cover_digest.value();
    return out;
}

struct FrozenSeparator {
    std::array<int, 9> values;
    i64 expected_delta;
    u64 expected_resolves;
};

const std::array<FrozenSeparator, 35> FROZEN = {{
    {{{15,16,30,60,120,170,176,190,240}}, 11332443973410LL, 23},
    {{{15,30,60,63,85,126,170,176,252}}, 6908621619942LL, 52},
    {{{15,16,30,60,63,120,126,240,252}}, 4205208449610LL, 9},
    {{{20,40,60,85,145,170,176,190,290}}, 4031205116652LL, 22},
    {{{8,40,42,80,120,132,190,252,264}}, 1161813182922LL, 73},
    {{{10,42,80,85,88,132,170,190,264}}, 3885931764792LL, 25},
    {{{15,30,60,63,85,120,126,170,252}}, 5035944515112LL, 15},
    {{{15,30,60,63,85,126,170,176,193}}, 1448234172LL, 35},
    {{{8,16,42,63,84,132,193,264,290}}, 8989654757580LL, 19},
    {{{15,16,30,60,85,168,170,190,290}}, 477508755150LL, 9},
    {{{10,40,85,120,145,170,176,240,286}}, 938229028002LL, 7},
    {{{20,40,60,85,143,170,176,190,286}}, 16138831115862LL, 4},
    {{{8,10,16,40,85,132,170,240,264}}, 595918160460LL, 18},
    {{{8,10,15,16,20,85,88,170,190}}, 702445254840LL, 16},
    {{{15,16,30,60,143,168,190,286,290}}, 1593975121830LL, 3},
    {{{8,16,42,63,84,132,168,176,264}}, 33337500840LL, 9},
    {{{15,16,30,60,120,143,190,240,286}}, 20233480454460LL, 3},
    {{{10,15,63,85,126,145,170,190,290}}, 739981581192LL, 5},
    {{{15,16,30,60,85,88,120,170,240}}, 19230377882940LL, 2},
    {{{8,16,42,63,85,88,132,168,264}}, 4121159517720LL, 4},
    {{{15,30,60,120,143,170,176,190,193}}, 9119096042742LL, 4},
    {{{8,16,40,85,132,170,176,264,286}}, 5460196721490LL, 7},
    {{{42,63,84,132,168,176,193,264,290}}, 1968341209002LL, 3},
    {{{10,15,42,85,143,170,193,286,290}}, 25124458174452LL, 1},
    {{{10,42,85,88,132,170,193,264,290}}, 21520242284562LL, 2},
    {{{8,10,42,85,132,168,170,176,264}}, 6646777574472LL, 4},
    {{{20,40,60,80,85,120,170,176,190}}, 17453318379462LL, 24},
    {{{10,42,85,88,132,170,176,264,290}}, 22838605917372LL, 1},
    {{{8,16,42,63,80,132,190,264,290}}, 6548342684310LL, 9},
    {{{10,42,85,120,170,176,193,240,290}}, 21973677633684LL, 2},
    {{{8,15,30,60,85,120,143,170,176}}, 2997034762932LL, 4},
    {{{42,63,80,126,132,145,190,252,264}}, 13343187024342LL, 7},
    {{{8,16,42,63,84,132,145,264,286}}, 2711801931930LL, 6},
    {{{15,16,30,60,85,145,168,170,190}}, 16187216478930LL, 7},
    {{{8,16,20,40,60,80,120,145,176}}, 12675717979470LL, 2}
}};

u64 verify_frozen_certificate(const Overlay& overlay,
                              const std::vector<i64>& mass9,
                              const std::vector<u32>& covers) {
    std::vector<unsigned char> unresolved(covers.size(), 1);
    std::size_t remaining = covers.size();
    std::set<u32> unique;
    Ledger digest;
    for (std::size_t index = 0; index < FROZEN.size(); ++index) {
        const FrozenSeparator& row = FROZEN[index];
        const u32 edge = label_mask(row.values);
        require(unique.insert(edge).second, "duplicate frozen E9 separator");
        const i64 direct_mass = direct_deletion_mass(edge, overlay);
        const i64 ranked_mass = mass9[colex_rank(edge, 9)];
        require(direct_mass == ranked_mass,
                "direct and reverse-incidence E9 masses disagree");
        const i128 delta = static_cast<i128>(63) * direct_mass -
                           static_cast<i128>(4) * overlay.denominator;
        require(delta == row.expected_delta,
                "frozen E9 separator delta changed");
        require(delta >= 0, "frozen E9 separator is unlawful");

        u64 resolves = 0;
        for (std::size_t body = 0; body < covers.size(); ++body) {
            if (unresolved[body] != 0 && (covers[body] & edge) == 0) {
                unresolved[body] = 0;
                ++resolves;
            }
        }
        require(resolves == row.expected_resolves,
                "frozen separator resolution count changed");
        require(resolves <= remaining, "separator remaining count underflow");
        remaining -= resolves;
        digest.add(edge);
        digest.add(static_cast<u64>(direct_mass));
        digest.add(resolves);
        digest.add(remaining);
        std::cout << "FROZEN_E9 index=" << (index + 1)
                  << " repair={" << labels(edge) << "}"
                  << " delta63=" << decimal(delta)
                  << " resolves=" << resolves
                  << " remaining=" << remaining << '\n';
    }
    require(remaining == 0, "35-edge separator left an E8 transversal");
    require(std::all_of(unresolved.begin(), unresolved.end(),
                        [](unsigned char value) { return value == 0; }),
            "unresolved vector disagrees with remaining count");
    return digest.value();
}

}  // namespace

int main() {
    try {
        const Overlay overlay = build_event_overlay();
        require(overlay.denominator == INT64_C(91205797082400),
                "pair denominator changed");
        require(overlay.event_count == 7464 && overlay.event_positions == 7312 &&
                    overlay.all_intervals == 7313 &&
                    overlay.pair_safe_intervals == 5358 &&
                    overlay.atom_mass.size() == 2296,
                "independent overlay census changed");
        require(overlay.pair_safe_mass == INT64_C(67011990405360),
                "pair-safe total mass changed");
        std::cout << "OVERLAY pair=50,51 denominator=" << overlay.denominator
                  << " events=" << overlay.event_count
                  << " event_positions=" << overlay.event_positions
                  << " intervals=" << overlay.all_intervals
                  << " pair_safe_intervals=" << overlay.pair_safe_intervals
                  << " atom_masks=" << overlay.atom_mass.size()
                  << " pair_safe_mass=" << overlay.pair_safe_mass
                  << " ledger=" << hex64(overlay.ledger) << '\n';

        std::vector<i64> mass8 = build_cardinality_masses(overlay, 8);
        std::vector<i64> mass9 = build_cardinality_masses(overlay, 9);
        const LayerSummary layer8 = summarize_layer(mass8, overlay.denominator);
        const LayerSummary layer9 = summarize_layer(mass9, overlay.denominator);
        require(layer8.edges == 311544 && layer8.equalities == 0,
                "independent E8 layer disagrees with claim");
        require(layer9.edges == 3159764 && layer9.equalities == 0,
                "independent E9 layer disagrees with claim");

        std::vector<Edge8> edges8;
        edges8.reserve(layer8.edges);
        auto retain8 = [&](u32 edge) {
            const i64 mass = mass8[colex_rank(edge, 8)];
            if (static_cast<i128>(63) * mass >=
                static_cast<i128>(4) * overlay.denominator) {
                edges8.push_back({edge, mass});
            }
        };
        for_each_subset(8, retain8);
        require(edges8.size() == layer8.edges,
                "recursive E8 reconstruction count drift");
        std::cout << "LAYERS E8_edges=" << layer8.edges
                  << " E8_equalities=" << layer8.equalities
                  << " E8_ledger=" << hex64(layer8.ledger)
                  << " E9_edges=" << layer9.edges
                  << " E9_equalities=" << layer9.equalities
                  << " E9_ledger=" << hex64(layer9.ledger) << '\n';

        const BodyAudit bodies = audit_all_nine_bodies(
            std::move(edges8), mass8, mass9, overlay.denominator);
        require(bodies.candidates == UINT64_C(14307150),
                "nine-body universe changed");
        require(bodies.checks == UINT64_C(39296089046),
                "nine-body incidence count changed");
        require(bodies.covers.size() == 436,
                "independent E8 transversal count disagrees with claim");
        require(bodies.monotone_extensions == bodies.candidates - 436 &&
                    bodies.minimum_extension_choices == 13,
                "E8-to-E9 deletion-monotonicity audit changed");
        std::cout << "E8_TRANSVERSALS bodies=" << bodies.candidates
                  << " checks=" << bodies.checks
                  << " covers=" << bodies.covers.size()
                  << " cover_ledger=" << hex64(bodies.cover_ledger)
                  << " noncovers_monotone_extended="
                  << bodies.monotone_extensions
                  << " extension_choices_each="
                  << bodies.minimum_extension_choices << '\n';

        const u64 separator_ledger =
            verify_frozen_certificate(overlay, mass9, bodies.covers);
        std::cout << "INDEPENDENT_DEPTH9_ACCEPT"
                  << " E8_nine_transversals=" << bodies.covers.size()
                  << " frozen_E9_edges=" << FROZEN.size()
                  << " separator_ledger=" << hex64(separator_ledger)
                  << " other_nine_sets=E8_TO_E9_BY_DELETION_MONOTONICITY"
                  << " conclusion=TAU_E9_GT_9\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
