// Exact THM-4278 literal-wall survey of active rank-eight repairs disjoint
// from a prescribed union of failed rank-nine bodies.

#define TOP_LITERAL_LIBRARY_ONLY
#include "lrc14_endpoint_520_688_minimum_one_atom_literal_wall_audit_thm4278.cpp"
#undef TOP_LITERAL_LIBRARY_ONLY

#include <queue>
#include <unordered_set>

namespace {

u32 parse_mask(const std::string& word) {
    std::size_t used = 0;
    const u64 wide = std::stoull(word, &used, 0);
    require(used == word.size() && wide < (UINT64_C(1) << 30),
            "bad forbidden mask");
    return static_cast<u32>(wide);
}

std::string labels_survey(u32 mask) {
    return labels(mask);
}

bool hits_all(u32 probe, const std::vector<u32>& family) {
    return std::all_of(family.begin(), family.end(),
                       [&](u32 mask) { return (mask & probe) != 0; });
}

void choose_transversals(const std::vector<int>& bits, int cursor, int need,
                         u32 mask, const std::vector<u32>& family,
                         std::vector<u32>& answer) {
    if (need == 0) {
        if (hits_all(mask, family)) answer.push_back(mask);
        return;
    }
    for (int i = cursor; i <= static_cast<int>(bits.size()) - need; ++i) {
        choose_transversals(bits, i + 1, need - 1,
                            mask | (u32{1} << bits[i]), family, answer);
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 4, "usage: family-survey q r forbidden-mask");
        const int q = std::stoi(argv[1]);
        const int r = std::stoi(argv[2]);
        const u32 forbidden = parse_mask(argv[3]);
        const LiteralAtoms atoms = build_literal_atoms(q, r);

        std::vector<int> allowed;
        for (int bit = 0; bit < 30; ++bit) {
            if ((forbidden & (u32{1} << bit)) == 0) allowed.push_back(bit);
        }
        require(allowed.size() >= 8, "fewer than eight allowed labels");
        std::vector<u32> candidates;
        enumerate_allowed(allowed, 0, 8, 0, candidates);

        std::vector<u32> active;
        u64 equalities = 0;
        i128 min_margin = 0;
        i128 max_margin = 0;
        u32 min_margin_mask = 0;
        u32 max_margin_mask = 0;
        for (u32 repair : candidates) {
            const i64 mass = literal_mass(atoms, repair);
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) * atoms.geometry.grid;
            if (margin < 0) continue;
            if (margin == 0) ++equalities;
            if (active.empty()) {
                min_margin = max_margin = margin;
                min_margin_mask = max_margin_mask = repair;
            }
            if (margin < min_margin) {
                min_margin = margin;
                min_margin_mask = repair;
            }
            if (margin > max_margin) {
                max_margin = margin;
                max_margin_mask = repair;
            }
            active.push_back(repair);
        }
        require(!active.empty(), "active common family is empty");
        std::sort(active.begin(), active.end());
        Fnv active_fnv;
        for (u32 mask : active) active_fnv.add(mask);

        u32 support = 0;
        u32 core = (u32{1} << 30) - 1;
        std::array<u64, 30> frequency{};
        for (u32 mask : active) {
            support |= mask;
            core &= mask;
            for (int bit = 0; bit < 30; ++bit) {
                if ((mask & (u32{1} << bit)) != 0) ++frequency[bit];
            }
        }
        const u32 allowed_mask = ((u32{1} << 30) - 1) & ~forbidden;
        const u32 unused_allowed = allowed_mask & ~support;

        std::unordered_set<u32> active_set(active.begin(), active.end());
        std::unordered_set<u32> seen;
        std::vector<u64> component_sizes;
        std::vector<u32> component_representatives;
        u64 twice_edges = 0;
        u64 min_degree = std::numeric_limits<u64>::max();
        u64 max_degree = 0;
        u32 min_degree_mask = 0;
        u32 max_degree_mask = 0;
        for (u32 mask : active) {
            u64 degree = 0;
            for (int out = 0; out < 30; ++out) {
                if ((mask & (u32{1} << out)) == 0) continue;
                for (int in : allowed) {
                    if ((mask & (u32{1} << in)) != 0) continue;
                    const u32 neighbor = mask ^ (u32{1} << out) ^
                                         (u32{1} << in);
                    if (active_set.contains(neighbor)) ++degree;
                }
            }
            twice_edges += degree;
            if (degree < min_degree) {
                min_degree = degree;
                min_degree_mask = mask;
            }
            if (degree > max_degree) {
                max_degree = degree;
                max_degree_mask = mask;
            }
            if (seen.contains(mask)) continue;
            u64 size = 0;
            std::queue<u32> queue;
            queue.push(mask);
            seen.insert(mask);
            while (!queue.empty()) {
                const u32 current = queue.front();
                queue.pop();
                ++size;
                for (int out = 0; out < 30; ++out) {
                    if ((current & (u32{1} << out)) == 0) continue;
                    for (int in : allowed) {
                        if ((current & (u32{1} << in)) != 0) continue;
                        const u32 neighbor = current ^ (u32{1} << out) ^
                                             (u32{1} << in);
                        if (active_set.contains(neighbor) &&
                            seen.insert(neighbor).second) {
                            queue.push(neighbor);
                        }
                    }
                }
            }
            component_sizes.push_back(size);
            component_representatives.push_back(mask);
        }
        require(twice_edges % 2 == 0, "Johnson edge parity failed");
        std::sort(component_sizes.begin(), component_sizes.end(),
                  std::greater<u64>());

        std::vector<int> support_bits;
        for (int bit = 0; bit < 30; ++bit) {
            if ((support & (u32{1} << bit)) != 0) support_bits.push_back(bit);
        }
        int transversal_number = 0;
        std::vector<u32> minimum_transversals;
        for (int size = 1; size <= 4; ++size) {
            choose_transversals(support_bits, 0, size, 0, active,
                                minimum_transversals);
            if (!minimum_transversals.empty()) {
                transversal_number = size;
                break;
            }
        }

        if (q == 520 && r == 688 && forbidden == UINT32_C(0x27587008)) {
            require(active.size() == 72 &&
                        active_fnv.state == UINT64_C(0xed1bfbaf6eaa06a3) &&
                        active.front() == UINT32_C(0x00048ec1),
                    "top-pair primary/literal control changed");
        }
        if (q == 542 && r == 732 && forbidden == UINT32_C(0x153a5400)) {
            require(active.size() == 2172 &&
                        active_fnv.state == UINT64_C(0x829ae906d6b54c9a) &&
                        active.front() == UINT32_C(0x0001aa89),
                    "542 primary/literal control changed");
        }

        std::cout << "LRC14_ONE_ATOM_COMMON_ACTIVE_FAMILY_SURVEY_THM4278\n";
        std::cout << "PAIR " << q << ',' << r << " FORBIDDEN " << std::hex
                  << forbidden << std::dec << " SIZE "
                  << std::popcount(forbidden) << " LABELS {"
                  << labels_survey(forbidden) << "}\n";
        std::cout << "GRID " << atoms.geometry.grid << " JOINT_CELLS "
                  << atoms.geometry.cells.size() << " LITERAL_ATOMS "
                  << atoms.mass.size() << " ALLOWED " << allowed.size()
                  << " CANDIDATES " << candidates.size() << " ACTIVE "
                  << active.size() << " EQUALITIES " << equalities
                  << " ACTIVE_FNV " << std::hex << active_fnv.state
                  << std::dec << '\n';
        std::cout << "LEAST " << std::hex << active.front() << std::dec
                  << " LABELS {" << labels_survey(active.front())
                  << "} MIN_MARGIN " << decimal(min_margin) << " MIN_MASK "
                  << std::hex << min_margin_mask << std::dec
                  << " MAX_MARGIN " << decimal(max_margin) << " MAX_MASK "
                  << std::hex << max_margin_mask << std::dec << '\n';
        std::cout << "SUPPORT SIZE " << std::popcount(support) << " LABELS {"
                  << labels_survey(support) << "} CORE SIZE "
                  << std::popcount(core) << " LABELS {" << labels_survey(core)
                  << "} UNUSED_ALLOWED SIZE " << std::popcount(unused_allowed)
                  << " LABELS {" << labels_survey(unused_allowed) << "}\n";
        std::cout << "JOHNSON COMPONENTS " << component_sizes.size()
                  << " SIZES";
        for (u64 size : component_sizes) std::cout << ' ' << size;
        std::cout << " EDGES " << twice_edges / 2 << " DEGREE_RANGE ["
                  << min_degree << ',' << max_degree << "]\n";
        std::cout << "DEGREE_EXTREMA MIN_MASK " << std::hex
                  << min_degree_mask << std::dec << " LABELS {"
                  << labels_survey(min_degree_mask) << "} MAX_MASK "
                  << std::hex << max_degree_mask << std::dec << " LABELS {"
                  << labels_survey(max_degree_mask) << "}\n";
        for (std::size_t i = 0; i < component_sizes.size(); ++i) {
            if (component_sizes[i] != 1) continue;
            const u32 isolated = component_representatives[i];
            i128 isolated_margin =
                static_cast<i128>(63) * literal_mass(atoms, isolated) -
                static_cast<i128>(4) * atoms.geometry.grid;
            bool first_neighbor = true;
            i128 best_neighbor_margin = 0;
            u32 best_neighbor = 0;
            for (int out = 0; out < 30; ++out) {
                if ((isolated & (u32{1} << out)) == 0) continue;
                for (int in : allowed) {
                    if ((isolated & (u32{1} << in)) != 0) continue;
                    const u32 neighbor = isolated ^ (u32{1} << out) ^
                                         (u32{1} << in);
                    const i128 margin =
                        static_cast<i128>(63) * literal_mass(atoms, neighbor) -
                        static_cast<i128>(4) * atoms.geometry.grid;
                    if (first_neighbor || margin > best_neighbor_margin) {
                        first_neighbor = false;
                        best_neighbor_margin = margin;
                        best_neighbor = neighbor;
                    }
                }
            }
            std::cout << "ISOLATED MASK " << std::hex << isolated << std::dec
                      << " LABELS {" << labels_survey(isolated)
                      << "} MARGIN " << decimal(isolated_margin)
                      << " BEST_ONE_SWAP " << std::hex << best_neighbor
                      << std::dec << " LABELS {"
                      << labels_survey(best_neighbor) << "} MARGIN "
                      << decimal(best_neighbor_margin) << '\n';
        }
        std::cout << "TRANSVERSAL_NUMBER_LE4 " << transversal_number
                  << " MINIMUM_COUNT " << minimum_transversals.size()
                  << " MINIMUM_LABELS";
        for (u32 mask : minimum_transversals) {
            std::cout << " {" << labels_survey(mask) << '}';
        }
        std::cout << '\n';
        std::cout << "FREQUENCIES";
        for (int bit : support_bits) {
            std::cout << ' ' << POOL[bit] << ':' << frequency[bit];
        }
        std::cout << '\n';
        std::cout << "STATUS FINITE_EXACT LITERAL_WALL PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
