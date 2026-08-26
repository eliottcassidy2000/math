#define main thm4179_dependency_main
#include "lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179_independent_audit.cpp"
#undef main

#include <functional>
#include <map>

namespace {

constexpr std::array<int, 3> HIERARCHY_ORIGINAL = {120, 126, 143};

bool hierarchy_is_original(int label) {
    return std::find(HIERARCHY_ORIGINAL.begin(), HIERARCHY_ORIGINAL.end(),
                     label) != HIERARCHY_ORIGINAL.end();
}

u32 hierarchy_mask_of_labels(const std::vector<int>& labels) {
    u32 mask = 0;
    for (int label : labels) {
        const auto found = std::find(POOL.begin(), POOL.end(), label);
        require(found != POOL.end(), "hierarchy label is outside the pool");
        mask |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    require(std::popcount(mask) == static_cast<int>(labels.size()),
            "hierarchy label list has a duplicate");
    return mask;
}

u32 hierarchy_mask_of_labels(std::initializer_list<int> labels) {
    return hierarchy_mask_of_labels(std::vector<int>(labels));
}

std::string hierarchy_labels(u32 mask) {
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

bool hierarchy_good(u32 mask) {
    int common_gcd = 0;
    for (int index = 0; index < 30; ++index) {
        if ((mask & (u32{1} << index)) != 0) {
            common_gcd = std::gcd(common_gcd, POOL[index]);
        }
    }
    if (common_gcd != 1) return false;
    for (int modulus = 2; modulus <= 14; ++modulus) {
        bool covered = false;
        for (int index = 0; index < 30; ++index) {
            if ((mask & (u32{1} << index)) != 0 &&
                POOL[index] % modulus == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;
    }
    return true;
}

template <class Callback>
void hierarchy_for_each_mask(const std::vector<int>& indices,
                             int cardinality,
                             Callback&& callback) {
    require(cardinality >= 0 && cardinality <= static_cast<int>(indices.size()),
            "invalid hierarchy combination cardinality");
    if (cardinality == 0) {
        callback(u32{0});
        return;
    }
    const u32 limit = u32{1} << static_cast<int>(indices.size());
    u32 local = (u32{1} << cardinality) - 1;
    while (local != 0 && local < limit) {
        u32 global = 0;
        u32 bits = local;
        while (bits != 0) {
            const u32 bit = bits & (~bits + 1u);
            const int local_index = std::countr_zero(bit);
            global |= u32{1} << indices[local_index];
            bits ^= bit;
        }
        callback(global);
        const u32 next = next_combination(local);
        if (next <= local) break;
        local = next;
    }
}

struct HierarchyLayer {
    int arity = 0;
    u64 candidates = 0;
    u64 equalities = 0;
    std::vector<u32> edges;
    i64 minimum_edge_delta = std::numeric_limits<i64>::max();
    i64 maximum_nonedge_delta = std::numeric_limits<i64>::min();
    u64 edge_hash = 0;
};

HierarchyLayer hierarchy_build_layer(const Geometry& geometry, int arity) {
    require(arity >= 5 && arity <= 7, "unexpected hierarchy arity");
    std::unordered_map<u32, i64> atom_mass;
    atom_mass.reserve(geometry.atoms.size());
    for (const Atom& atom : geometry.atoms) {
        if (atom.failed_newcomer || std::popcount(atom.failed_pool) > arity) {
            continue;
        }
        atom_mass[atom.failed_pool] += atom.right - atom.left;
    }

    HierarchyLayer layer;
    layer.arity = arity;
    Fnv1a64 hash;
    hash.add_u64_le(static_cast<u64>(arity));
    const u32 limit = u32{1} << 30;
    u32 deletion = (u32{1} << arity) - 1;
    while (deletion != 0 && deletion < limit) {
        ++layer.candidates;
        i64 mass = 0;
        u32 subset = deletion;
        while (true) {
            const auto found = atom_mass.find(subset);
            if (found != atom_mass.end()) mass += found->second;
            if (subset == 0) break;
            subset = (subset - 1) & deletion;
        }
        const i64 delta = static_cast<i64>(
            static_cast<i128>(63) * mass -
            static_cast<i128>(4) * geometry.denominator);
        if (delta == 0) ++layer.equalities;
        if (delta >= 0) {
            layer.edges.push_back(deletion);
            layer.minimum_edge_delta = std::min(layer.minimum_edge_delta, delta);
            hash.add_u64_le(deletion);
            hash.add_i64_le(mass);
            hash.add_i64_le(delta);
        } else {
            layer.maximum_nonedge_delta =
                std::max(layer.maximum_nonedge_delta, delta);
        }
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    require(layer.candidates == binomial(30, arity),
            "hierarchy layer candidate count changed");
    layer.edge_hash = hash.value();
    return layer;
}

using AnchorLayers = std::array<std::vector<u32>, 7>;

bool hierarchy_lexicographic_mask_less(u32 left, u32 right) {
    while (left != 0 && right != 0) {
        const int left_index = std::countr_zero(left);
        const int right_index = std::countr_zero(right);
        if (left_index != right_index) return left_index < right_index;
        left &= left - 1;
        right &= right - 1;
    }
    return left == 0 && right != 0;
}

AnchorLayers hierarchy_enumerate_anchors(const std::vector<int>& ground) {
    AnchorLayers layers;
    std::vector<u32> previous;
    for (int size = 1; size <= 6; ++size) {
        hierarchy_for_each_mask(ground, size, [&](u32 mask) {
            if (!hierarchy_good(mask)) return;
            for (u32 old : previous) {
                if ((mask & old) == old) return;
            }
            layers[size].push_back(mask);
        });
        std::sort(layers[size].begin(), layers[size].end(),
                  hierarchy_lexicographic_mask_less);
        previous.insert(previous.end(), layers[size].begin(), layers[size].end());
    }
    const std::array<std::size_t, 7> expected = {0, 0, 0, 0, 32, 297, 24};
    const std::array<u64, 7> expected_hash = {
        0, UINT64_C(0xcbf29ce484222325), UINT64_C(0xcbf29ce484222325),
        UINT64_C(0xcbf29ce484222325), UINT64_C(0xf9f6fd7d15fa2ddb),
        UINT64_C(0x3841dfd22c768456), UINT64_C(0x11739190505473e1)};
    for (int size = 1; size <= 6; ++size) {
        require(layers[size].size() == expected[size],
                "hierarchy anchor-layer count changed");
        Fnv1a64 hash;
        for (u32 anchor : layers[size]) hash.add_u64_le(anchor);
        require(hash.value() == expected_hash[size],
                "hierarchy anchor-layer labels changed");
    }
    const u32 forced_286 = hierarchy_mask_of_labels({286});
    for (int size : {4, 5}) {
        for (u32 anchor : layers[size]) {
            require((anchor & forced_286) != 0,
                    "286 is not forced in a size-four/five anchor");
        }
    }
    const u32 size6_core = hierarchy_mask_of_labels({42, 63, 132, 286});
    const u32 left = hierarchy_mask_of_labels({8, 16, 88, 176});
    const u32 right = hierarchy_mask_of_labels({10, 20, 30, 170, 190, 290});
    for (u32 anchor : layers[6]) {
        require((anchor & size6_core) == size6_core &&
                    std::popcount(anchor & left) == 1 &&
                    std::popcount(anchor & right) == 1,
                "size-six product classification changed");
    }
    return layers;
}

std::vector<u32> hierarchy_active_edges(const std::vector<u32>& edges,
                                        u32 anchor) {
    std::vector<u32> active;
    active.reserve(edges.size());
    for (u32 edge : edges) {
        if ((edge & anchor) == 0) active.push_back(edge);
    }
    require(!active.empty(), "hierarchy active edge set is empty");
    return active;
}

bool hierarchy_is_cover(u32 chosen, const std::vector<u32>& edges) {
    for (u32 edge : edges) {
        if ((edge & chosen) == 0) return false;
    }
    return true;
}

std::vector<int> hierarchy_optional_vertices(u32 anchor) {
    std::vector<int> vertices;
    for (int index = 0; index < 30; ++index) {
        if ((anchor & (u32{1} << index)) == 0) vertices.push_back(index);
    }
    return vertices;
}

std::vector<u32> hierarchy_complete_covers(const std::vector<u32>& edges,
                                           u32 anchor,
                                           int budget) {
    const std::vector<int> vertices = hierarchy_optional_vertices(anchor);
    std::vector<u32> covers;
    for (int size = 0; size <= budget; ++size) {
        hierarchy_for_each_mask(vertices, size, [&](u32 candidate) {
            if (hierarchy_is_cover(candidate, edges)) covers.push_back(candidate);
        });
    }
    std::sort(covers.begin(), covers.end());
    require(std::adjacent_find(covers.begin(), covers.end()) == covers.end(),
            "duplicate hierarchy cover");
    return covers;
}

std::vector<u32> hierarchy_filter_covers(const std::vector<u32>& candidates,
                                         const std::vector<u32>& edges) {
    std::vector<u32> survivors;
    for (u32 candidate : candidates) {
        if (hierarchy_is_cover(candidate, edges)) survivors.push_back(candidate);
    }
    return survivors;
}

u64 hierarchy_update_digest(u64 digest,
                            u32 anchor,
                            u64 active_count,
                            const std::vector<u32>& covers,
                            const std::vector<u32>& survivors) {
    auto add = [&](u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            digest ^= static_cast<std::uint8_t>((word >> shift) & 0xffu);
            digest *= FNV1A64_PRIME;
        }
    };
    add(anchor);
    add(active_count);
    add(covers.size());
    for (u32 cover : covers) add(cover);
    add(survivors.size());
    for (u32 survivor : survivors) add(survivor);
    return digest;
}

struct HierarchyPresentation {
    u32 anchor;
    u32 cover;
};

void hierarchy_audit_size_five(const AnchorLayers& anchors,
                               const HierarchyLayer& layer5,
                               const HierarchyLayer& layer6,
                               const HierarchyLayer& layer7) {
    u64 digest = FNV1A64_OFFSET;
    std::map<int, u64> histogram;
    int empty_rows = 0;
    int blocked_rows = 0;
    std::vector<HierarchyPresentation> residual;
    for (u32 anchor : anchors[5]) {
        const std::vector<u32> active5 =
            hierarchy_active_edges(layer5.edges, anchor);
        const std::vector<u32> covers =
            hierarchy_complete_covers(active5, anchor, 5);
        if (covers.empty()) ++empty_rows;
        else ++blocked_rows;
        for (u32 cover : covers) ++histogram[std::popcount(cover)];
        const std::vector<u32> active6 =
            hierarchy_active_edges(layer6.edges, anchor);
        const std::vector<u32> survivors =
            hierarchy_filter_covers(covers, active6);
        for (u32 cover : survivors) residual.push_back({anchor, cover});
        digest = hierarchy_update_digest(
            digest, anchor, active5.size(), covers, survivors);
    }
    require(empty_rows == 41 && blocked_rows == 256,
            "independent size-five row split changed");
    require(histogram == std::map<int, u64>{{3, 33}, {4, 1229}, {5, 20244}},
            "independent complete size-five cover histogram changed");
    require(residual.size() == 15,
            "independent size-five E6 residual count changed");
    require(digest == UINT64_C(0x9dbf51ecd6f8c2d5),
            "independent size-five row ledger changed: " + hex64(digest));

    std::map<u32, std::vector<HierarchyPresentation>> bodies;
    for (const HierarchyPresentation& presentation : residual) {
        bodies[presentation.anchor | presentation.cover].push_back(presentation);
        std::cout << "SIZE5_E6_BLOCKER anchor="
                  << hierarchy_labels(presentation.anchor)
                  << " K=" << hierarchy_labels(presentation.cover) << "\n";
    }
    const u32 body16 = hierarchy_mask_of_labels(
        {16, 88, 95, 145, 168, 170, 193, 240, 252, 286});
    const u32 body80 = hierarchy_mask_of_labels(
        {80, 88, 95, 145, 168, 170, 193, 240, 252, 286});
    require(bodies.size() == 2 && bodies[body16].size() == 9 &&
                bodies[body80].size() == 6,
            "independent size-five residual bodies changed");
    const std::map<u32, std::pair<u64, u32>> expected = {
        {body16, {736, hierarchy_mask_of_labels({8, 10, 15, 20, 80, 85, 290})}},
        {body80, {514, hierarchy_mask_of_labels({8, 10, 15, 16, 20, 85, 120})}}};
    for (const auto& [body, body_presentations] : bodies) {
        u64 witness_count = 0;
        u32 first_witness = 0;
        for (u32 edge : layer7.edges) {
            if ((edge & body) != 0) continue;
            if (first_witness == 0) first_witness = edge;
            ++witness_count;
        }
        require(witness_count == expected.at(body).first &&
                    first_witness == expected.at(body).second,
                "independent size-five E7 incidence changed");
        std::cout << "SIZE5_RESIDUAL_BODY labels=" << hierarchy_labels(body)
                  << " presentations=" << body_presentations.size()
                  << " E7_witnesses=" << witness_count
                  << " first_witness=" << hierarchy_labels(first_witness)
                  << "\n";
    }
    std::cout << "SIZE5_SUMMARY no_E5_cover_rows=" << empty_rows
              << " E5_cover_rows=" << blocked_rows
              << " complete_E5_covers="
              << histogram[3] + histogram[4] + histogram[5]
              << " cover_size_hist=3:" << histogram[3]
              << ",4:" << histogram[4] << ",5:" << histogram[5]
              << " E6_blocker_presentations=" << residual.size()
              << " E6_blocker_bodies=" << bodies.size()
              << " row_ledger_fnv1a64_le=" << hex64(digest) << "\n";
}

void hierarchy_audit_size_six(const AnchorLayers& anchors,
                              const HierarchyLayer& layer5,
                              const HierarchyLayer& layer6) {
    u64 digest = FNV1A64_OFFSET;
    int empty_rows = 0;
    int blocked_rows = 0;
    std::vector<HierarchyPresentation> covers_all;
    for (u32 anchor : anchors[6]) {
        const std::vector<u32> active5 =
            hierarchy_active_edges(layer5.edges, anchor);
        const std::vector<u32> covers =
            hierarchy_complete_covers(active5, anchor, 4);
        if (covers.empty()) ++empty_rows;
        else ++blocked_rows;
        const std::vector<u32> active6 =
            hierarchy_active_edges(layer6.edges, anchor);
        const std::vector<u32> survivors =
            hierarchy_filter_covers(covers, active6);
        require(survivors.empty(), "independent size-six E6 blocker survived");
        for (u32 cover : covers) covers_all.push_back({anchor, cover});
        digest = hierarchy_update_digest(
            digest, anchor, active5.size(), covers, survivors);
    }
    require(empty_rows == 22 && blocked_rows == 2 && covers_all.size() == 4,
            "independent size-six row/cover split changed");
    require(digest == UINT64_C(0xa4f776dc340f0260),
            "independent size-six row ledger changed");
    const u32 exceptional88 =
        hierarchy_mask_of_labels({42, 63, 88, 132, 286, 290});
    const u32 exceptional176 =
        hierarchy_mask_of_labels({42, 63, 132, 176, 286, 290});
    int count88 = 0;
    int count176 = 0;
    for (const HierarchyPresentation& presentation : covers_all) {
        count88 += presentation.anchor == exceptional88;
        count176 += presentation.anchor == exceptional176;
        const u32 body = presentation.anchor | presentation.cover;
        u64 witness_count = 0;
        u32 first_witness = 0;
        for (u32 edge : layer6.edges) {
            if ((edge & body) != 0) continue;
            if (first_witness == 0) first_witness = edge;
            ++witness_count;
        }
        require(witness_count > 0,
                "independent size-six E5 cover did not cross E6");
        std::cout << "SIZE6_E5_BLOCKER anchor="
                  << hierarchy_labels(presentation.anchor)
                  << " K=" << hierarchy_labels(presentation.cover)
                  << " E6_witnesses=" << witness_count
                  << " first_witness=" << hierarchy_labels(first_witness)
                  << "\n";
    }
    require(count88 == 3 && count176 == 1,
            "independent size-six exceptional rows changed");
    std::cout << "SIZE6_SUMMARY no_E5_cover_rows=" << empty_rows
              << " E5_cover_rows=" << blocked_rows
              << " complete_E5_covers=" << covers_all.size()
              << " cover_size_hist=4:4 E6_blockers=0"
              << " row_ledger_fnv1a64_le=" << hex64(digest) << "\n";
}

void hierarchy_audit_body_universe(const std::vector<int>& ground,
                                   const AnchorLayers& anchors) {
    std::array<u64, 7> presentations{};
    presentations[4] = anchors[4].size() * binomial(23, 6);
    presentations[5] = anchors[5].size() * binomial(22, 5);
    presentations[6] = anchors[6].size() * binomial(21, 4);
    require(presentations[4] == 3230304 && presentations[5] == 7821198 &&
                presentations[6] == 143640,
            "independent hierarchy presentation counts changed");

    u64 valid = 0;
    std::array<u64, 7> first_anchor_size{};
    hierarchy_for_each_mask(ground, 10, [&](u32 body) {
        if (!hierarchy_good(body)) return;
        ++valid;
        int first = 0;
        for (int size : {4, 5, 6}) {
            for (u32 anchor : anchors[size]) {
                if ((body & anchor) == anchor) {
                    first = size;
                    break;
                }
            }
            if (first != 0) break;
        }
        require(first != 0,
                "valid zero-original ten-body has no size-four/five/six anchor");
        ++first_anchor_size[first];
    });
    require(valid == 1491665 && first_anchor_size[4] == 1138494 &&
                first_anchor_size[5] == 348712 &&
                first_anchor_size[6] == 4459,
            "independent zero-original body hierarchy changed");
    std::cout << "BODY_HIERARCHY ambient=" << binomial(27, 10)
              << " valid_primitive_divisor_complete=" << valid
              << " size4_union=" << first_anchor_size[4]
              << " after_size5=" << first_anchor_size[4] + first_anchor_size[5]
              << " after_size6=" << valid
              << " remainder_after_size4="
              << first_anchor_size[5] + first_anchor_size[6]
              << " remainder_after_size5=" << first_anchor_size[6]
              << " remainder_after_size6=0"
              << " presentations_by_anchor_size=4:" << presentations[4]
              << ",5:" << presentations[5] << ",6:" << presentations[6]
              << "\n";
}

void hierarchy_run_audit() {
    force_binary_stdout();
    const Geometry geometry = build_joint_geometry();
    require(geometry.denominator == 91205797082400 &&
                geometry.walls.size() == 7214 &&
                geometry.atoms.size() == 7213 &&
                geometry.atom_hash == UINT64_C(0x0ccd305ae47c79ea),
            "independent joint geometry control changed");

    std::vector<int> ground;
    for (int index = 0; index < 30; ++index) {
        if (!hierarchy_is_original(POOL[index])) ground.push_back(index);
    }
    require(ground.size() == 27, "zero-original hierarchy ground changed");
    const AnchorLayers anchors = hierarchy_enumerate_anchors(ground);

    std::cout << "AUDIT q50_zero_original_minimal_anchor_hierarchy_joint_wall_cpp20\n";
    std::cout << "JOINT_GEOMETRY denominator=" << geometry.denominator
              << " walls=" << geometry.walls.size()
              << " atoms=" << geometry.atoms.size()
              << " labelled_atom_fnv1a64_le=" << hex64(geometry.atom_hash)
              << "\n";
    for (int size = 1; size <= 6; ++size) {
        Fnv1a64 hash;
        for (u32 anchor : anchors[size]) hash.add_u64_le(anchor);
        std::cout << "MINIMAL_ANCHOR_LAYER size=" << size
                  << " count=" << anchors[size].size()
                  << " mask_fnv1a64_le=" << hex64(hash.value()) << "\n";
    }
    std::cout << "SIZE6_PRODUCT core=42,63,132,286"
              << " left=8,16,88,176 right=10,20,30,170,190,290\n";

    HierarchyLayer layer5 = hierarchy_build_layer(geometry, 5);
    HierarchyLayer layer6 = hierarchy_build_layer(geometry, 6);
    HierarchyLayer layer7 = hierarchy_build_layer(geometry, 7);
    for (HierarchyLayer* layer : {&layer5, &layer6, &layer7}) {
        std::sort(layer->edges.begin(), layer->edges.end(),
                  hierarchy_lexicographic_mask_less);
    }
    require(layer5.edges.size() == 3017 && layer5.equalities == 0 &&
                layer5.minimum_edge_delta == 9312918552LL &&
                layer5.maximum_nonedge_delta == -3452771070LL &&
                layer5.edge_hash == UINT64_C(0x6beebf7aa4701ef0),
            "independent global E5 layer changed");
    require(layer6.edges.size() == 85324 && layer6.equalities == 0 &&
                layer6.minimum_edge_delta == 17425422LL &&
                layer6.maximum_nonedge_delta == -131971518LL &&
                layer6.edge_hash == UINT64_C(0xe9951e7b4a3fba86),
            "independent global E6 layer changed");
    require(layer7.edges.size() == 821737 && layer7.equalities == 0 &&
                layer7.minimum_edge_delta == 17425422LL &&
                layer7.maximum_nonedge_delta == -20115648LL &&
                layer7.edge_hash == UINT64_C(0x808ca531b29c5e83),
            "independent global E7 layer changed");
    for (const HierarchyLayer* layer : {&layer5, &layer6, &layer7}) {
        std::cout << "GLOBAL_LAYER d=" << layer->arity
                  << " candidates=" << layer->candidates
                  << " edges=" << layer->edges.size()
                  << " equalities=" << layer->equalities
                  << " min_edge_delta=" << layer->minimum_edge_delta
                  << " max_nonedge_delta=" << layer->maximum_nonedge_delta
                  << " edge_fnv1a64_le=" << hex64(layer->edge_hash) << "\n";
    }

    hierarchy_audit_size_five(anchors, layer5, layer6, layer7);
    hierarchy_audit_size_six(anchors, layer5, layer6);
    hierarchy_audit_body_universe(ground, anchors);
    std::cout << "Q50_FOUR_DISJOINT_ORIGINAL_ANCHOR_SLICES 10111881\n";
    std::cout << "VERDICT PASS\n";
}

}  // namespace

int main() {
    try {
        hierarchy_run_audit();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
