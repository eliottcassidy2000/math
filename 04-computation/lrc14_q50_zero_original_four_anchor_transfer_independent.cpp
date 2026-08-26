#define main thm4179_dependency_main
#include "lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179_independent_audit.cpp"
#undef main

#include <functional>
#include <map>

namespace {

using Anchor4 = std::array<int, 4>;

constexpr std::array<int, 4> Q4_X = {85, 95, 145, 193};
constexpr std::array<int, 8> Q4_BODY_SPINE =
    {88, 95, 145, 168, 193, 240, 252, 286};
constexpr std::array<std::array<int, 2>, 5> Q4_SELECTOR_PAIRS = {{
    {8, 80}, {16, 85}, {16, 170}, {80, 85}, {80, 170}}};
constexpr std::array<u64, 5> Q4_BODY_PRESENTATIONS = {6, 4, 3, 8, 6};
constexpr std::array<u64, 5> Q4_BODY_WITNESSES = {243, 439, 736, 278, 514};

std::string q4_anchor_text(const Anchor4& anchor) {
    return std::to_string(anchor[0]) + "," + std::to_string(anchor[1]) +
           "," + std::to_string(anchor[2]) + "," + std::to_string(anchor[3]);
}

std::string q4_global_labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (int index = 0; index < 30; ++index) {
        if ((mask & (u32{1} << index)) == 0) continue;
        if (!first) out << ',';
        first = false;
        out << POOL[index];
    }
    return out.str();
}

bool q4_is_original(int label) {
    return label == 120 || label == 126 || label == 143;
}

bool q4_divisor_complete(const Anchor4& anchor) {
    for (int modulus = 2; modulus <= 14; ++modulus) {
        bool covered = false;
        for (int label : anchor) covered = covered || label % modulus == 0;
        if (!covered) return false;
    }
    return true;
}

int q4_anchor_gcd(const Anchor4& anchor) {
    int value = 0;
    for (int label : anchor) value = std::gcd(value, label);
    return value;
}

u32 q4_pool_mask(const Anchor4& anchor) {
    u32 mask = 0;
    for (int label : anchor) {
        const auto found = std::find(POOL.begin(), POOL.end(), label);
        require(found != POOL.end(), "four-anchor label is outside the pool");
        mask |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    require(std::popcount(mask) == 4, "four-anchor has repeated labels");
    return mask;
}

u32 q4_label_mask(std::initializer_list<int> labels) {
    u32 mask = 0;
    for (int label : labels) {
        const auto found = std::find(POOL.begin(), POOL.end(), label);
        require(found != POOL.end(), "declared label is outside the pool");
        mask |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    return mask;
}

template <std::size_t N>
u32 q4_array_mask(const std::array<int, N>& labels) {
    u32 mask = 0;
    for (int label : labels) {
        const auto found = std::find(POOL.begin(), POOL.end(), label);
        require(found != POOL.end(), "array label is outside the pool");
        mask |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    return mask;
}

std::vector<Anchor4> q4_enumerate_anchors() {
    std::vector<Anchor4> anchors;
    for (int i = 0; i < 30; ++i) {
        if (q4_is_original(POOL[i])) continue;
        for (int j = i + 1; j < 30; ++j) {
            if (q4_is_original(POOL[j])) continue;
            for (int k = j + 1; k < 30; ++k) {
                if (q4_is_original(POOL[k])) continue;
                for (int ell = k + 1; ell < 30; ++ell) {
                    if (q4_is_original(POOL[ell])) continue;
                    const Anchor4 anchor =
                        {POOL[i], POOL[j], POOL[k], POOL[ell]};
                    if (q4_anchor_gcd(anchor) == 1 &&
                        q4_divisor_complete(anchor)) {
                        anchors.push_back(anchor);
                    }
                }
            }
        }
    }
    return anchors;
}

struct Q4Layer {
    int arity = 0;
    u64 candidates = 0;
    u64 equalities = 0;
    std::vector<u32> edges;
    i64 minimum_edge_delta = std::numeric_limits<i64>::max();
    i64 maximum_nonedge_delta = std::numeric_limits<i64>::min();
    u64 edge_hash = 0;
};

Q4Layer q4_build_global_layer(const Geometry& geometry, int arity) {
    require(arity == 6 || arity == 7, "unexpected four-anchor deletion arity");
    std::unordered_map<u32, i64> failure_mass;
    failure_mass.reserve(geometry.atoms.size());
    for (const Atom& atom : geometry.atoms) {
        if (atom.failed_newcomer || std::popcount(atom.failed_pool) > arity) {
            continue;
        }
        failure_mass[atom.failed_pool] += atom.right - atom.left;
    }

    Q4Layer layer;
    layer.arity = arity;
    Fnv1a64 edge_hash;
    edge_hash.add_u64_le(static_cast<u64>(arity));
    const u32 limit = u32{1} << 30;
    u32 deletion = (u32{1} << arity) - 1;
    while (deletion != 0 && deletion < limit) {
        ++layer.candidates;
        i64 mass = 0;
        u32 subset = deletion;
        while (true) {
            const auto found = failure_mass.find(subset);
            if (found != failure_mass.end()) mass += found->second;
            if (subset == 0) break;
            subset = (subset - 1) & deletion;
        }
        const i64 delta = static_cast<i64>(
            static_cast<i128>(63) * mass -
            static_cast<i128>(4) * geometry.denominator);
        if (delta == 0) ++layer.equalities;
        if (delta >= 0) {
            layer.edges.push_back(deletion);
            layer.minimum_edge_delta =
                std::min(layer.minimum_edge_delta, delta);
            edge_hash.add_u64_le(deletion);
            edge_hash.add_i64_le(mass);
            edge_hash.add_i64_le(delta);
        } else {
            layer.maximum_nonedge_delta =
                std::max(layer.maximum_nonedge_delta, delta);
        }
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    require(layer.candidates == binomial(30, arity),
            "global layer candidate count mismatch");
    require(!layer.edges.empty(), "global repair layer is empty");
    layer.edge_hash = edge_hash.value();
    return layer;
}

std::vector<int> q4_optional_pool_indices(u32 anchor_mask) {
    std::vector<int> indices;
    for (int index = 0; index < 30; ++index) {
        if ((anchor_mask & (u32{1} << index)) == 0) indices.push_back(index);
    }
    require(indices.size() == 26, "four-anchor optional set is not size 26");
    return indices;
}

u32 q4_compress(u32 global_mask, const std::vector<int>& local_to_pool) {
    u32 local_mask = 0;
    for (std::size_t local = 0; local < local_to_pool.size(); ++local) {
        if ((global_mask & (u32{1} << local_to_pool[local])) != 0) {
            local_mask |= u32{1} << local;
        }
    }
    return local_mask;
}

u32 q4_expand(u32 local_mask, const std::vector<int>& local_to_pool) {
    u32 global_mask = 0;
    for (std::size_t local = 0; local < local_to_pool.size(); ++local) {
        if ((local_mask & (u32{1} << local)) != 0) {
            global_mask |= u32{1} << local_to_pool[local];
        }
    }
    return global_mask;
}

std::vector<u32> q4_active_local_edges(
    const std::vector<u32>& global_edges,
    u32 anchor_mask,
    const std::vector<int>& local_to_pool) {
    std::vector<u32> active;
    for (u32 edge : global_edges) {
        if ((edge & anchor_mask) == 0) active.push_back(q4_compress(edge, local_to_pool));
    }
    require(!active.empty(), "anchor repair hypergraph is empty");
    return active;
}

struct Q4CoverEnumeration {
    std::vector<u32> covers;
    u64 nodes = 0;
    std::array<u64, 7> level_sizes{};
};

Q4CoverEnumeration q4_enumerate_covers_through_six(
    const std::vector<u32>& edges) {
    std::array<std::unordered_set<u32>, 7> seen;
    seen[0].insert(0);
    Q4CoverEnumeration result;
    std::function<void(u32)> search = [&](u32 chosen) {
        ++result.nodes;
        u32 uncovered = 0;
        for (u32 edge : edges) {
            if ((edge & chosen) == 0) {
                uncovered = edge;
                break;
            }
        }
        if (uncovered == 0) {
            result.covers.push_back(chosen);
            return;
        }
        const int depth = std::popcount(chosen);
        if (depth == 6) return;
        for (u32 branch = uncovered; branch != 0; branch &= branch - 1) {
            const u32 bit = branch & (~branch + 1u);
            const u32 child = chosen | bit;
            if (seen[depth + 1].insert(child).second) search(child);
        }
    };
    search(0);
    for (int depth = 0; depth <= 6; ++depth) {
        result.level_sizes[depth] = seen[depth].size();
    }
    for (u32 cover : result.covers) {
        require(is_cover(cover, edges), "enumerated object is not a cover");
    }
    return result;
}

bool q4_is_blocked_anchor(const Anchor4& anchor) {
    if (anchor[2] != 252 || anchor[3] != 286) return false;
    const bool x_member =
        std::find(Q4_X.begin(), Q4_X.end(), anchor[1]) != Q4_X.end();
    const bool left = anchor[0] == 80 && x_member;
    const bool right = anchor[1] == 240 &&
        std::find(Q4_X.begin(), Q4_X.end(), anchor[0]) != Q4_X.end();
    return left || right;
}

int q4_expected_blockers(const Anchor4& anchor) {
    require(q4_is_blocked_anchor(anchor), "requested blocker count for positive row");
    if (anchor[0] == 80) return anchor[1] == 85 ? 1 : 3;
    return anchor[0] == 85 ? 2 : 5;
}

u32 q4_declared_d7_cover(const Anchor4& anchor) {
    require(q4_is_blocked_anchor(anchor), "requested d7 cover for positive row");
    u32 cover = q4_label_mask({8, 88, 168});
    const int x = anchor[0] == 80 ? anchor[1] : anchor[0];
    if (anchor[0] == 80) cover |= q4_label_mask({240});
    else cover |= q4_label_mask({80});
    for (int label : Q4_X) {
        if (label != x) cover |= q4_label_mask({label});
    }
    require(std::popcount(cover) == 7, "declared d7 cover is not size seven");
    return cover;
}

i64 q4_deletion_mass(const Geometry& geometry, u32 deletion) {
    i64 mass = 0;
    for (const Atom& atom : geometry.atoms) {
        if (atom.failed_newcomer) continue;
        if ((atom.failed_pool & ~deletion) == 0) mass += atom.right - atom.left;
    }
    return mass;
}

std::string q4_level_text(const std::array<u64, 7>& levels) {
    std::ostringstream out;
    for (int depth = 0; depth <= 6; ++depth) {
        if (depth != 0) out << ',';
        out << levels[depth];
    }
    return out.str();
}

void q4_audit_body_union(const std::vector<Anchor4>& anchors) {
    std::vector<int> zero_original_indices;
    for (int index = 0; index < 30; ++index) {
        if (!q4_is_original(POOL[index])) zero_original_indices.push_back(index);
    }
    require(zero_original_indices.size() == 27, "zero-original ground set changed");
    std::vector<u32> anchor_masks;
    for (const Anchor4& anchor : anchors) anchor_masks.push_back(q4_pool_mask(anchor));

    std::array<u64, 33> histogram{};
    u64 distinct = 0;
    u64 presentations = 0;
    const u32 limit = u32{1} << 27;
    u32 local_body = (u32{1} << 10) - 1;
    while (local_body != 0 && local_body < limit) {
        u32 global_body = 0;
        for (int local = 0; local < 27; ++local) {
            if ((local_body & (u32{1} << local)) != 0) {
                global_body |= u32{1} << zero_original_indices[local];
            }
        }
        int multiplicity = 0;
        for (u32 anchor : anchor_masks) {
            multiplicity += (global_body & anchor) == anchor;
        }
        if (multiplicity > 0) {
            ++distinct;
            presentations += static_cast<u64>(multiplicity);
            ++histogram[multiplicity];
        }
        const u32 next = next_combination(local_body);
        if (next <= local_body) break;
        local_body = next;
    }
    require(distinct == 1138494 && presentations == 3230304,
            "independent zero-original body union changed");
    const std::map<int, u64> expected = {
        {1, 262761}, {2, 332137}, {3, 204900}, {4, 178178},
        {5, 59560}, {6, 66399}, {7, 13188}, {8, 11229},
        {9, 6149}, {10, 2666}, {11, 658}, {12, 493},
        {13, 140}, {15, 36}};
    for (int multiplicity = 1; multiplicity <= 32; ++multiplicity) {
        const auto found = expected.find(multiplicity);
        const u64 wanted = found == expected.end() ? 0 : found->second;
        require(histogram[multiplicity] == wanted,
                "independent body multiplicity histogram changed");
    }
    std::cout << "ZERO_ORIGINAL_BODY_UNION presentations=" << presentations
              << " distinct=" << distinct << " histogram=";
    bool first = true;
    for (const auto& [multiplicity, count] : expected) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << multiplicity << ':' << count;
    }
    std::cout << "\n";
}

void q4_run_audit() {
    force_binary_stdout();
    const Geometry geometry = build_joint_geometry();
    require(geometry.denominator == 91205797082400 &&
                geometry.walls.size() == 7214 &&
                geometry.atoms.size() == 7213 &&
                geometry.atom_hash == UINT64_C(0x0ccd305ae47c79ea),
            "joint geometry control changed");
    const std::vector<Anchor4> anchors = q4_enumerate_anchors();
    require(anchors.size() == 32, "four-anchor count is not 32");
    require(std::all_of(anchors.begin(), anchors.end(), [](const Anchor4& anchor) {
                return std::find(anchor.begin(), anchor.end(), 286) != anchor.end();
            }), "286 is not forced in every four-anchor");

    std::cout << "AUDIT q50_zero_original_four_anchor_joint_wall_cpp20\n";
    std::cout << "JOINT_GEOMETRY denominator=" << geometry.denominator
              << " walls=" << geometry.walls.size()
              << " atoms=" << geometry.atoms.size()
              << " labelled_atom_fnv1a64_le=" << hex64(geometry.atom_hash)
              << "\n";
    std::cout << "FOUR_ANCHORS " << anchors.size() << " all_contain_286=yes\n";
    for (const Anchor4& anchor : anchors) {
        std::cout << "ANCHOR " << q4_anchor_text(anchor) << "\n";
    }

    const Q4Layer layer6 = q4_build_global_layer(geometry, 6);
    const Q4Layer layer7 = q4_build_global_layer(geometry, 7);
    require(layer6.edges.size() == 85324 && layer6.equalities == 0,
            "independent global d6 layer changed");
    require(layer7.edges.size() == 821737 && layer7.equalities == 0,
            "independent global d7 layer changed");
    for (const Q4Layer* layer : {&layer6, &layer7}) {
        std::cout << "GLOBAL_LAYER d=" << layer->arity
                  << " candidates=" << layer->candidates
                  << " edges=" << layer->edges.size()
                  << " equalities=" << layer->equalities
                  << " min_edge_delta=" << layer->minimum_edge_delta
                  << " max_nonedge_delta=" << layer->maximum_nonedge_delta
                  << " edge_fnv1a64_le=" << hex64(layer->edge_hash) << "\n";
    }

    struct BlockerPresentation {
        Anchor4 anchor;
        u32 anchor_mask;
        u32 blocker;
    };
    std::vector<BlockerPresentation> blocker_presentations;
    int positive_d6 = 0;
    int blocked_rows = 0;

    for (const Anchor4& anchor : anchors) {
        const u32 anchor_mask = q4_pool_mask(anchor);
        const std::vector<int> local_to_pool = q4_optional_pool_indices(anchor_mask);
        const std::vector<u32> active6 =
            q4_active_local_edges(layer6.edges, anchor_mask, local_to_pool);
        ExactCoverSearch budget_six(active6, 6);
        const SearchResult decision = budget_six.run();
        if (!q4_is_blocked_anchor(anchor)) {
            require(!decision.found, "unexpected independent d6 cover");
            ++positive_d6;
            std::cout << "D6_ROW anchor=" << q4_anchor_text(anchor)
                      << " edges=" << active6.size()
                      << " tau_gt_6=yes search_nodes=" << decision.nodes << "\n";
            continue;
        }

        require(decision.found && std::popcount(decision.cover) == 6,
                "declared blocked row did not have a six-cover");
        ++blocked_rows;
        const Q4CoverEnumeration enumeration =
            q4_enumerate_covers_through_six(active6);
        require(static_cast<int>(enumeration.covers.size()) ==
                    q4_expected_blockers(anchor),
                "independent d6 blocker count changed");
        require(std::all_of(enumeration.covers.begin(), enumeration.covers.end(),
                            [](u32 cover) { return std::popcount(cover) == 6; }),
                "independent cover below depth six");
        std::cout << "D6_ROW anchor=" << q4_anchor_text(anchor)
                  << " edges=" << active6.size()
                  << " tau=6 blockers=" << enumeration.covers.size()
                  << " enum_nodes=" << enumeration.nodes
                  << " enum_levels=" << q4_level_text(enumeration.level_sizes)
                  << "\n";
        for (u32 local_blocker : enumeration.covers) {
            const u32 blocker = q4_expand(local_blocker, local_to_pool);
            blocker_presentations.push_back({anchor, anchor_mask, blocker});
            std::cout << "D6_BLOCKER anchor=" << q4_anchor_text(anchor)
                      << " K=" << q4_global_labels(blocker) << "\n";
        }

        const std::vector<u32> active7 =
            q4_active_local_edges(layer7.edges, anchor_mask, local_to_pool);
        const u32 declared_cover = q4_declared_d7_cover(anchor);
        const u32 local_cover = q4_compress(declared_cover, local_to_pool);
        require(q4_expand(local_cover, local_to_pool) == declared_cover,
                "declared d7 cover used an anchor");
        require(is_cover(local_cover, active7), "declared d7 cover failed");
        std::cout << "D7_ROW anchor=" << q4_anchor_text(anchor)
                  << " edges=" << active7.size()
                  << " tau=7 cover=" << q4_global_labels(declared_cover)
                  << "\n";
    }
    require(positive_d6 == 24 && blocked_rows == 8 &&
                blocker_presentations.size() == 27,
            "independent 24/8/27 split changed");

    std::map<u32, std::vector<BlockerPresentation>> body_presentations;
    for (const BlockerPresentation& presentation : blocker_presentations) {
        const u32 body = presentation.anchor_mask | presentation.blocker;
        bool crossed = false;
        for (u32 repair : layer7.edges) {
            if ((repair & body) == 0) {
                crossed = true;
                break;
            }
        }
        require(crossed, "independent depth-seven blocker did not cross");
        body_presentations[body].push_back(presentation);
    }
    require(body_presentations.size() == 5,
            "independent blockers did not collapse to five body states");

    std::array<u32, 5> expected_bodies{};
    const u32 spine = q4_array_mask(Q4_BODY_SPINE);
    for (int index = 0; index < 5; ++index) {
        expected_bodies[index] = spine | q4_array_mask(Q4_SELECTOR_PAIRS[index]);
        const auto found = body_presentations.find(expected_bodies[index]);
        require(found != body_presentations.end(), "expected obstruction body absent");
        u64 witness_count = 0;
        u32 first_witness = 0;
        for (u32 repair : layer7.edges) {
            if ((repair & expected_bodies[index]) != 0) continue;
            if (first_witness == 0) first_witness = repair;
            ++witness_count;
        }
        require(found->second.size() == Q4_BODY_PRESENTATIONS[index] &&
                    witness_count == Q4_BODY_WITNESSES[index],
                "independent body incidence count changed");
        std::cout << "BLOCKER_BODY index=" << index
                  << " selector=" << Q4_SELECTOR_PAIRS[index][0] << ','
                  << Q4_SELECTOR_PAIRS[index][1]
                  << " labels=" << q4_global_labels(expected_bodies[index])
                  << " presentations=" << found->second.size()
                  << " d7_witnesses=" << witness_count
                  << " first_witness=" << q4_global_labels(first_witness)
                  << "\n";
    }

    std::map<int, u64> pattern_counts;
    std::map<int, u32> pattern_examples;
    for (u32 repair : layer7.edges) {
        int pattern = 0;
        for (int index = 0; index < 5; ++index) {
            if ((repair & expected_bodies[index]) == 0) pattern |= 1 << index;
        }
        if (pattern == 0) continue;
        ++pattern_counts[pattern];
        pattern_examples.emplace(pattern, repair);
    }
    const std::map<int, u64> expected_patterns = {
        {30, 27}, {21, 5}, {6, 284}, {10, 11}, {20, 91}, {24, 172},
        {1, 238}, {2, 117}, {4, 329}, {8, 68}, {16, 219}};
    require(pattern_counts == expected_patterns && pattern_counts.count(31) == 0,
            "independent witness pattern census changed");
    for (const auto& [pattern, count] : pattern_counts) {
        std::cout << "WITNESS_PATTERN mask=" << pattern
                  << " count=" << count
                  << " example=" << q4_global_labels(pattern_examples.at(pattern))
                  << "\n";
    }

    const std::array<u32, 2> bank = {
        q4_label_mask({8, 10, 15, 20, 143, 176, 290}),
        q4_label_mask({10, 15, 20, 85, 120, 143, 176})};
    int bank_union = 0;
    for (u32 repair : bank) {
        require(std::binary_search(layer7.edges.begin(), layer7.edges.end(), repair),
                "independent bank member is not a d7 repair");
        int pattern = 0;
        for (int index = 0; index < 5; ++index) {
            if ((repair & expected_bodies[index]) == 0) pattern |= 1 << index;
        }
        bank_union |= pattern;
        const i64 mass = q4_deletion_mass(geometry, repair);
        const i64 delta = static_cast<i64>(
            static_cast<i128>(63) * mass -
            static_cast<i128>(4) * geometry.denominator);
        require(delta > 0, "independent bank member is not strict");
        std::cout << "BANK_REPAIR labels=" << q4_global_labels(repair)
                  << " pattern=" << pattern << " mass_ticks=" << mass
                  << " delta=" << delta << "\n";
    }
    require(bank_union == 31, "independent two-repair bank is incomplete");
    std::cout << "COMMON_ONE_REPAIR_BANK no EXPLICIT_TWO_REPAIR_BANK yes\n";

    q4_audit_body_union(anchors);
    std::cout << "Q50_FOUR_DISJOINT_ORIGINAL_ANCHOR_SLICES 9758710\n";
    std::cout << "VERDICT PASS\n";
}

}  // namespace

int main() {
    try {
        q4_run_audit();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
