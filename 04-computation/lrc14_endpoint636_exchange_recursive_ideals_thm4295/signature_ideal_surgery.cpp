#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <bit>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>

namespace {

constexpr std::size_t kJointRowsProbe = 421;
constexpr std::size_t kBodyRowsProbe = 14307150;

struct ProbePair {
    int q = 0;
    int r = 0;
    friend bool operator<(const ProbePair& a, const ProbePair& b) {
        return std::tie(a.q, a.r) < std::tie(b.q, b.r);
    }
};

using ProbeSignature = std::array<u64, 7>;

std::vector<std::string> probe_split(const std::string& s, char delimiter) {
    std::vector<std::string> out;
    std::stringstream in(s);
    std::string field;
    while (std::getline(in, field, delimiter)) out.push_back(field);
    return out;
}

std::vector<unsigned> probe_indices(const ProbeSignature& s) {
    std::vector<unsigned> out;
    for (unsigned i = 0; i < kJointRowsProbe; ++i)
        if ((s[i / 64] >> (i % 64)) & 1) out.push_back(i);
    return out;
}

bool probe_subset(const ProbeSignature& a, const ProbeSignature& b) {
    for (std::size_t w = 0; w < a.size(); ++w)
        if ((a[w] & ~b[w]) != 0) return false;
    return true;
}

std::map<ProbePair, ProbeSignature> probe_read_signatures(const char* path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open signatures");
    std::string line;
    require(std::getline(in, line) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::map<ProbePair, ProbeSignature> out;
    while (std::getline(in, line)) {
        const auto f = probe_split(line, ',');
        require(f.size() == 10, "bad signature row");
        ProbeSignature s{};
        unsigned weight = 0;
        for (unsigned w = 0; w < 7; ++w) {
            s[w] = std::stoull(f[w + 3], nullptr, 16);
            weight += std::popcount(s[w]);
        }
        require(weight == std::stoul(f[2]), "signature weight mismatch");
        require(out.emplace(ProbePair{std::stoi(f[0]), std::stoi(f[1])}, s).second,
                "duplicate signature row");
    }
    require(out.size() == 24223, "signature universe changed");
    return out;
}

std::vector<ProbePair> probe_read_rows(const char* path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open residual");
    std::vector<ProbePair> out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto f = probe_split(line, ',');
        require(f.size() == 2, "bad residual row");
        out.push_back({std::stoi(f[0]), std::stoi(f[1])});
    }
    require(out.size() == 22647, "current residual changed");
    return out;
}

std::vector<u32> probe_read_joint(const char* path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open joint deck");
    std::vector<u32> out;
    std::string token;
    while (in >> token) out.push_back(static_cast<u32>(std::stoul(token, nullptr, 16)));
    require(out.size() == kJointRowsProbe, "joint deck changed");
    return out;
}

ProbeSignature probe_deleted_signature(const std::string& arg) {
    ProbeSignature out{};
    for (const std::string& token : probe_split(arg, ',')) {
        const unsigned i = std::stoul(token);
        require(i < kJointRowsProbe, "deleted index out of range");
        require(((out[i / 64] >> (i % 64)) & 1) == 0, "duplicate deleted index");
        out[i / 64] |= u64{1} << (i % 64);
    }
    require(!probe_indices(out).empty(), "empty deletion set");
    return out;
}

u64 probe_pair_fnv(const std::vector<ProbePair>& rows) {
    FnvLocal fnv;
    for (auto p : rows) { fnv.add(p.q); fnv.add(p.r); }
    return fnv.state;
}

u32 probe_next_same_popcount(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}

std::vector<u32> probe_lost_bodies(const std::vector<u32>& joint,
                                   const std::set<unsigned>& deleted,
                                   u64& checks, u64& fnv_out) {
    std::vector<u32> lost;
    FnvLocal fnv;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < kBodyRowsProbe; ++ordinal) {
        bool retained_hit = false;
        for (unsigned i = 0; i < joint.size(); ++i) {
            if (deleted.contains(i)) continue;
            ++checks;
            if ((joint[i] & body) == 0) { retained_hit = true; break; }
        }
        if (!retained_hit) {
            bool old_hit = false;
            for (unsigned i : deleted) old_hit |= (joint[i] & body) == 0;
            require(old_hit, "inherited joint deck did not cover body");
            lost.push_back(body);
            fnv.add(body);
        }
        if (ordinal + 1 < kBodyRowsProbe) body = probe_next_same_popcount(body);
    }
    require(body == UINT32_C(0x3fe00000), "body endpoint changed");
    fnv_out = fnv.state;
    return lost;
}

struct ProbeClass { u64 count = 0; u32 least = 0; };

struct ProbeCandidate64 {
    u64 pattern = 0;
    u32 least = 0;
};

bool probe_cover_search64(const std::vector<ProbeCandidate64>& candidates,
                          u64 full, u64 covered, int remaining,
                          std::unordered_map<u64, int>& memo,
                          std::vector<u32>& path) {
    if (covered == full) return true;
    if (remaining == 0) return false;
    const auto remembered = memo.find(covered);
    if (remembered != memo.end() && remembered->second >= remaining)
        return false;
    memo[covered] = remaining;
    const u64 uncovered = full & ~covered;
    unsigned max_gain = 0;
    for (const auto& candidate : candidates)
        max_gain = std::max(max_gain,
                            static_cast<unsigned>(std::popcount(
                                candidate.pattern & uncovered)));
    if (max_gain == 0 ||
        (std::popcount(uncovered) + max_gain - 1) / max_gain > remaining)
        return false;

    unsigned chosen_bit = 0;
    std::size_t chosen_options = std::numeric_limits<std::size_t>::max();
    for (unsigned bit = 0; bit < 64; ++bit) {
        const u64 bit_mask = u64{1} << bit;
        if ((uncovered & bit_mask) == 0) continue;
        std::size_t options = 0;
        for (const auto& candidate : candidates)
            options += (candidate.pattern & bit_mask) != 0;
        if (options < chosen_options) {
            chosen_options = options;
            chosen_bit = bit;
        }
    }
    const u64 chosen_mask = u64{1} << chosen_bit;
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < candidates.size(); ++i)
        if (candidates[i].pattern & chosen_mask) order.push_back(i);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        const unsigned gain_a = std::popcount(candidates[a].pattern & uncovered);
        const unsigned gain_b = std::popcount(candidates[b].pattern & uncovered);
        if (gain_a != gain_b) return gain_a > gain_b;
        return candidates[a].pattern > candidates[b].pattern;
    });
    for (std::size_t i : order) {
        const u64 joined = covered | candidates[i].pattern;
        if (joined == covered) continue;
        path.push_back(candidates[i].least);
        if (probe_cover_search64(candidates, full, joined, remaining - 1,
                                 memo, path)) return true;
        path.pop_back();
    }
    return false;
}

bool probe_is_full(const std::vector<u64>& p, std::size_t n) {
    for (std::size_t w = 0; w < p.size(); ++w) {
        const u64 expected = (w + 1 == p.size() && n % 64)
            ? ((u64{1} << (n % 64)) - 1) : ~u64{0};
        if (p[w] != expected) return false;
    }
    return true;
}

bool probe_union_full(const std::vector<u64>& a, const std::vector<u64>& b,
                      std::size_t n) {
    for (std::size_t w = 0; w < a.size(); ++w) {
        const u64 expected = (w + 1 == a.size() && n % 64)
            ? ((u64{1} << (n % 64)) - 1) : ~u64{0};
        if ((a[w] | b[w]) != expected) return false;
    }
    return true;
}

struct ProbeBodyScan { u64 checks = 0; u64 failures = 0; u64 fnv = 0; };

ProbeBodyScan probe_scan_body_cover(const std::vector<u32>& deck) {
    ProbeBodyScan out;
    FnvLocal fnv;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < kBodyRowsProbe; ++ordinal) {
        u64 prefix = 0;
        for (u32 mask : deck) {
            ++prefix;
            ++out.checks;
            if ((mask & body) == 0) break;
        }
        if (prefix == deck.size() && (deck.back() & body) != 0) {
            ++out.failures;
            prefix = 0;
        }
        fnv.add(body); fnv.add(prefix);
        if (ordinal + 1 < kBodyRowsProbe) body = probe_next_same_popcount(body);
    }
    out.fnv = fnv.state;
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 7, "usage: probe JOINT SIGNATURES CURRENT DELETED_CSV FIBRE_OUT WITNESS_OUT");
        init_choose8_local();
        const auto joint = probe_read_joint(argv[1]);
        const auto signatures = probe_read_signatures(argv[2]);
        const auto current = probe_read_rows(argv[3]);
        ProbeSignature deleted_sig{};
        const std::string deletion_arg = argv[4];
        if (!deletion_arg.empty() && deletion_arg.front() == '@') {
            const auto qr = probe_split(deletion_arg.substr(1), ',');
            require(qr.size() == 2, "target deletion must be @q,r");
            const ProbePair target{std::stoi(qr[0]), std::stoi(qr[1])};
            const auto found = signatures.find(target);
            require(found != signatures.end(), "target signature absent");
            deleted_sig = found->second;
        } else {
            deleted_sig = probe_deleted_signature(deletion_arg);
        }
        const auto deleted_vec = probe_indices(deleted_sig);
        const std::set<unsigned> deleted(deleted_vec.begin(), deleted_vec.end());

        std::vector<ProbePair> fibre;
        std::map<std::vector<unsigned>, u64> shapes;
        for (auto p : current) {
            const auto found = signatures.find(p);
            require(found != signatures.end(), "residual row missing signature");
            if (probe_subset(found->second, deleted_sig)) {
                fibre.push_back(p);
                ++shapes[probe_indices(found->second)];
            }
        }
        require(!fibre.empty(), "empty current signature ideal");
        std::ofstream fibre_out(argv[5]);
        for (auto p : fibre) fibre_out << p.q << ',' << p.r << '\n';
        fibre_out.close();

        u64 retained_checks = 0, lost_fnv = 0;
        const auto lost = probe_lost_bodies(joint, deleted, retained_checks, lost_fnv);
        require(!lost.empty(), "deletion exposed no bodies");

        if (std::getenv("PROBE_CENSUS_ONLY") != nullptr) {
            std::cout << "THM4295_SIGNATURE_IDEAL_CENSUS_V1\nDELETED "
                      << deleted.size() << " INDICES";
            for (unsigned i : deleted) std::cout << ' ' << i;
            std::cout << "\nFIBRE_ROWS " << fibre.size() << " FNV "
                      << std::hex << probe_pair_fnv(fibre) << std::dec
                      << " SHAPES " << shapes.size() << " TOP636";
            for (auto p : fibre)
                if (p.r == 636) std::cout << ' ' << p.q << ',' << p.r;
            std::cout << "\nLOST_BODIES " << lost.size() << " FNV "
                      << std::hex << lost_fnv << std::dec
                      << " RETAINED_CHECKS " << retained_checks << '\n';
            for (std::size_t i = 0;
                 i < std::min<std::size_t>(lost.size(), 12); ++i)
                std::cout << "LOST " << i << ' ' << std::hex
                          << std::setw(8) << std::setfill('0') << lost[i]
                          << std::dec << std::setfill(' ') << '\n';
            std::cout << "SCOPE FINITE_EXACT_FIXED_POOL_CENSUS_ONLY\n";
            return 0;
        }

        const auto cells = build_pool_cells();
        std::vector<unsigned char> common_active(EXPECTED_REPAIRS, 1);
        u64 activity_tests = 0;
        for (auto p : fibre) {
            const ActiveUniverse active = build_active_universe(cells, p.q, p.r);
            for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank)
                common_active[rank] &= active.active[rank];
            activity_tests += EXPECTED_REPAIRS;
        }

        u64 common_count = 0, response_incidences = 0;
        FnvLocal active_fnv, response_fnv;
        std::map<std::vector<u64>, ProbeClass> classes;
        const std::size_t words = (lost.size() + 63) / 64;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            if (!common_active[rank]) continue;
            const u32 mask = unrank_colex8(rank);
            ++common_count; active_fnv.add(mask);
            std::vector<u64> response(words);
            for (std::size_t i = 0; i < lost.size(); ++i)
                if ((mask & lost[i]) == 0) response[i / 64] |= u64{1} << (i % 64);
            for (u64 word : response) response_incidences += std::popcount(word);
            response_fnv.add(mask);
            for (u64 word : response) response_fnv.add(word);
            auto& cls = classes[response];
            ++cls.count;
            if (cls.count == 1 || mask < cls.least) cls.least = mask;
        }

        int exact_min = -1;
        std::vector<u32> witness;
        for (const auto& [pattern, cls] : classes) {
            if (probe_is_full(pattern, lost.size())) {
                exact_min = 1; witness = {cls.least}; break;
            }
        }
        if (exact_min < 0) {
            for (auto a = classes.begin(); a != classes.end() && exact_min < 0; ++a)
                for (auto b = a; b != classes.end(); ++b)
                    if (probe_union_full(a->first, b->first, lost.size())) {
                        exact_min = 2; witness = {a->second.least, b->second.least};
                        break;
                    }
        }
        std::size_t maximal_classes = 0;
        if (exact_min < 0 && lost.size() <= 64) {
            std::vector<ProbeCandidate64> candidates;
            for (const auto& [pattern, cls] : classes) {
                require(pattern.size() == 1, "64-bit quotient shape changed");
                if (pattern[0] != 0) candidates.push_back({pattern[0], cls.least});
            }
            std::sort(candidates.begin(), candidates.end(),
                      [](const auto& a, const auto& b) {
                          const unsigned pa = std::popcount(a.pattern);
                          const unsigned pb = std::popcount(b.pattern);
                          if (pa != pb) return pa > pb;
                          return a.pattern > b.pattern;
                      });
            std::vector<ProbeCandidate64> maximal;
            for (const auto& candidate : candidates) {
                bool dominated = false;
                for (const auto& kept : maximal)
                    if ((candidate.pattern & ~kept.pattern) == 0) {
                        dominated = true;
                        break;
                    }
                if (!dominated) maximal.push_back(candidate);
            }
            candidates = std::move(maximal);
            maximal_classes = candidates.size();
            const u64 full = lost.size() == 64 ? ~u64{0}
                                               : ((u64{1} << lost.size()) - 1);
            for (int depth = 3; depth <= 12 && exact_min < 0; ++depth) {
                std::unordered_map<u64, int> memo;
                std::vector<u32> path;
                if (probe_cover_search64(candidates, full, 0, depth, memo,
                                         path)) {
                    exact_min = depth;
                    witness = std::move(path);
                }
            }
        }

        ProbeBodyScan scan;
        u64 rebuilt_fnv = 0;
        if (exact_min > 0) {
            std::vector<u32> rebuilt;
            for (unsigned i = 0; i < joint.size(); ++i)
                if (!deleted.contains(i)) rebuilt.push_back(joint[i]);
            rebuilt.insert(rebuilt.end(), witness.begin(), witness.end());
            std::set<u32> distinct(rebuilt.begin(), rebuilt.end());
            require(distinct.size() == rebuilt.size(), "rebuilt deck has duplicate mask");
            FnvLocal fnv;
            for (u32 mask : rebuilt) fnv.add(mask);
            rebuilt_fnv = fnv.state;
            scan = probe_scan_body_cover(rebuilt);
            require(scan.failures == 0, "response witness failed full body scan");
        }
        std::ofstream witness_out(argv[6]);
        for (u32 mask : witness)
            witness_out << std::hex << std::setw(8) << std::setfill('0') << mask << '\n';
        witness_out.close();

        std::cout << "THM4295_SIGNATURE_IDEAL_SURGERY_V1\nDELETED " << deleted.size() << " INDICES";
        for (unsigned i : deleted) std::cout << ' ' << i;
        std::cout << "\nFIBRE_ROWS " << fibre.size() << " FNV " << std::hex
                  << probe_pair_fnv(fibre) << std::dec << " SHAPES " << shapes.size()
                  << " TOP636";
        for (auto p : fibre) if (p.r == 636) std::cout << ' ' << p.q << ',' << p.r;
        std::cout << "\nLOST_BODIES " << lost.size() << " FNV " << std::hex << lost_fnv
                  << std::dec << " RETAINED_CHECKS " << retained_checks << '\n';
        for (std::size_t i = 0; i < std::min<std::size_t>(lost.size(), 12); ++i)
            std::cout << "LOST " << i << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << lost[i] << std::dec << std::setfill(' ') << '\n';
        std::cout << "COMMON_ACTIVE " << common_count << " FNV " << std::hex
                  << active_fnv.state << " RESPONSE_FNV " << response_fnv.state
                  << std::dec << " ACTIVITY_TESTS " << activity_tests
                  << " RESPONSE_CLASSES " << classes.size()
                  << " MAXIMAL_CLASSES " << maximal_classes
                  << " RESPONSE_INCIDENCES " << response_incidences << '\n'
                  << "EXACT_MIN_UP_TO_12 " << exact_min << " WITNESS";
        for (u32 mask : witness) std::cout << ' ' << std::hex << std::setw(8)
                                           << std::setfill('0') << mask;
        std::cout << std::dec << std::setfill(' ') << "\nREBUILT_FNV " << std::hex
                  << rebuilt_fnv << " BODY_ROW_FNV " << scan.fnv << std::dec
                  << " BODY_CHECKS " << scan.checks << " FAILURES " << scan.failures
                  << "\nSCOPE FINITE_EXACT_FIXED_POOL_NO_LRC14\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "THM4295_SIGNATURE_IDEAL_ERROR " << e.what() << '\n';
        return 1;
    }
}
