#define JOINT421_LITERAL_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/verify_joint421_literal_r670.cpp"
#undef JOINT421_LITERAL_LIBRARY_ONLY

#include <fstream>
#include <iomanip>

namespace {

constexpr std::array<u32, 6> CONTINUATION_4281 = {{
    UINT32_C(0x003884c8), UINT32_C(0x00c4c124),
    UINT32_C(0x02a05206), UINT32_C(0x10e05240),
    UINT32_C(0x0004409f), UINT32_C(0x00c0c125)}};
constexpr std::array<std::size_t, 7> INACTIVE_JOINT_INDICES = {{
    75, 107, 139, 374, 394, 405, 417}};

void write_masks(const std::filesystem::path& path,
                 const std::vector<u32>& masks) {
    std::ofstream output(path);
    require(static_cast<bool>(output), "cannot create carrier output");
    for (u32 mask : masks)
        output << std::hex << std::setw(8) << std::setfill('0') << mask
               << '\n';
    require(output.good(), "failed writing carrier output");
}

bool rational_less_nonnegative(i128 a, i128 b, i128 c, i128 d) {
    require(a >= 0 && c >= 0 && b > 0 && d > 0,
            "invalid nonnegative rational comparison");
    bool reversed = false;
    while (true) {
        const i128 aq = a / b;
        const i128 cq = c / d;
        if (aq != cq) return reversed ? aq > cq : aq < cq;
        a %= b;
        c %= d;
        if (a == 0 || c == 0) {
            if (a == 0 && c == 0) return false;
            const bool raw_less = a == 0;
            return reversed ? !raw_less : raw_less;
        }
        std::swap(a, b);
        std::swap(c, d);
        reversed = !reversed;
    }
}

struct LiteralPartition {
    std::vector<u32> active;
    std::vector<u32> inactive;
    std::array<u64, 3> active_parts{};
    std::array<u64, 3> inactive_parts{};
    u64 equalities = 0;
    u32 weakest_active_mask = 0;
    i128 weakest_active_margin = 0;
    u32 closest_inactive_mask = 0;
    i128 closest_inactive_deficit = 0;
    bool weakest_set = false;
    bool closest_set = false;
};

LiteralPartition partition_literal(const Geometry& geometry,
                                   const std::vector<u32>& carrier) {
    LiteralPartition out;
    for (std::size_t index = 0; index < carrier.size(); ++index) {
        const u32 repair = carrier[index];
        i64 mass = 0;
        for (const Cell& cell : geometry.cells) {
            if (cell.pair_safe && (cell.failed_pool & ~repair) == 0)
                mass += cell.width;
        }
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * geometry.grid;
        const std::size_t part = index < 8524 ? 0 : (index < 8945 ? 1 : 2);
        out.equalities += margin == 0;
        if (margin >= 0) {
            out.active.push_back(repair);
            ++out.active_parts[part];
            if (!out.weakest_set || rational_less_nonnegative(
                    margin, geometry.grid, out.weakest_active_margin,
                    geometry.grid) ||
                (margin == out.weakest_active_margin &&
                 repair < out.weakest_active_mask)) {
                out.weakest_set = true;
                out.weakest_active_mask = repair;
                out.weakest_active_margin = margin;
            }
        } else {
            out.inactive.push_back(repair);
            ++out.inactive_parts[part];
            const i128 deficit = -margin;
            if (!out.closest_set || rational_less_nonnegative(
                    deficit, geometry.grid, out.closest_inactive_deficit,
                    geometry.grid) ||
                (deficit == out.closest_inactive_deficit &&
                 repair < out.closest_inactive_mask)) {
                out.closest_set = true;
                out.closest_inactive_mask = repair;
                out.closest_inactive_deficit = deficit;
            }
        }
    }
    require(out.weakest_set && out.closest_set,
            "literal carrier partition lacks a side");
    return out;
}

std::string hex8(u32 value) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 10,
                "usage: audit REPLAY_BAND PREFIX_416_704 PREFIX_416_700 "
                "PREFIX_520_700 PREFIX_384_694 PREFIX_520_688 JOINT_DECK "
                "CARRIER_OUT VULNERABLE_OUT");
        const std::filesystem::path band(argv[1]);
        const Parsed prefix704 = parse_probe(argv[2]);
        const Parsed prefix416 = parse_probe(argv[3]);
        const Parsed prefix520 = parse_probe(argv[4]);
        const Parsed prefix384 = parse_probe(argv[5]);
        const Parsed prefix688 = parse_probe(argv[6]);
        require(prefix704.q == 416 && prefix704.r == 704 &&
                    prefix704.deck.size() == 2608 &&
                    prefix704.declared_prefix_fnv ==
                        UINT64_C(0x18ff663a123e684e) &&
                    prefix416.q == 416 && prefix416.r == 700 &&
                    prefix416.deck.size() == 5894 &&
                    prefix416.declared_prefix_fnv ==
                        UINT64_C(0x701c0233c8f8abeb) &&
                    prefix520.q == 520 && prefix520.r == 700 &&
                    prefix520.deck.size() == 4557 &&
                    prefix520.declared_prefix_fnv ==
                        UINT64_C(0xd8466558bcd8ef9d) &&
                    prefix384.q == 384 && prefix384.r == 694 &&
                    prefix384.deck.size() == 5307 &&
                    prefix384.declared_prefix_fnv ==
                        UINT64_C(0x7b9a46c60e8514f8) &&
                    prefix688.q == 520 && prefix688.r == 688 &&
                    prefix688.deck.size() == 5398 &&
                    prefix688.declared_prefix_fnv ==
                        UINT64_C(0x6ab471da88c8e1d1),
                "canonical prefix identity changed");

        std::vector<u32> carrier;
        std::set<u32> seen;
        for (const auto& [q, r] : EXPECTED_PAIRS) {
            const Parsed parsed = parse_probe(
                band / (std::to_string(q) + "_" + std::to_string(r) +
                        ".out"));
            require(parsed.q == q && parsed.r == r,
                    "canonical band transcript mismatch");
            append_unseen_literal421(parsed.deck, seen, carrier);
        }
        require(carrier.size() == 4675 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0xce4e76ec11df057c),
                "base carrier changed");
        append_unseen_literal421(prefix704.deck, seen, carrier);
        require(carrier.size() == 4733 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0xa7b046289655c733),
                "round-one carrier changed");
        const std::set<u32> round1_seen = seen;
        std::vector<u32> novel416;
        std::vector<u32> novel520;
        for (u32 repair : prefix416.deck)
            if (!round1_seen.count(repair)) novel416.push_back(repair);
        for (u32 repair : prefix520.deck)
            if (!round1_seen.count(repair)) novel520.push_back(repair);
        append_unseen_literal421(novel416, seen, carrier);
        append_unseen_literal421(novel520, seen, carrier);
        require(carrier.size() == 7986 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0xbaef1d2f49444638),
                "round-two carrier changed");
        append_unseen_literal421(prefix384.deck, seen, carrier);
        require(carrier.size() == 8319 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0xe08b227730f6793c),
                "round-three carrier changed");
        append_unseen_literal421(prefix688.deck, seen, carrier);
        require(carrier.size() == 8518 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0x1603e3fe970f8428),
                "THM-4271 carrier changed");
        append_unseen_literal421(
            std::vector<u32>(THM4276_AUGMENTATION.begin(),
                             THM4276_AUGMENTATION.end()),
            seen, carrier);
        require(carrier.size() == 8524 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0x5ddb84a44f5d2ad7),
                "THM-4276 carrier changed");

        const std::vector<u32> joint = read_joint421(argv[7]);
        for (u32 mask : joint)
            require(!seen.count(mask), "joint mask overlaps base carrier");
        append_unseen_literal421(joint, seen, carrier);
        require(carrier.size() == 8945 &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0x3212efa05dd18c00),
                "joint carrier changed");
        for (u32 witness : CONTINUATION_4281)
            require(std::popcount(witness) == 8 && !seen.count(witness),
                    "continuation witness invalid/duplicate");
        append_unseen_literal421(
            std::vector<u32>(CONTINUATION_4281.begin(),
                             CONTINUATION_4281.end()),
            seen, carrier);
        require(carrier.size() == 8951 && seen.size() == carrier.size() &&
                    fnv_literal421(carrier) ==
                        UINT64_C(0x188f82ab9dd1695a),
                "final THM-4281 carrier changed");
        write_masks(argv[8], carrier);

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "labelled body universe changed");
        Fnv body_universe_ledger;
        for (u32 body : bodies) body_universe_ledger.add(body);

        const Geometry geometry = build_joint_geometry(256, 663);
        const LiteralPartition final_partition =
            partition_literal(geometry, carrier);
        require(final_partition.active.size() +
                    final_partition.inactive.size() == carrier.size(),
                "literal carrier partition lost masks");

        const std::vector<u32> active_joint = active_literal421(256, 663, joint);
        std::vector<std::size_t> inactive_joint_indices;
        std::vector<u32> inactive_joint;
        Fnv inactive_joint_ledger;
        for (std::size_t index = 0; index < joint.size(); ++index) {
            if (std::find(active_joint.begin(), active_joint.end(),
                          joint[index]) != active_joint.end())
                continue;
            inactive_joint_indices.push_back(index);
            inactive_joint.push_back(joint[index]);
            inactive_joint_ledger.add(index);
            inactive_joint_ledger.add(joint[index]);
        }
        require(inactive_joint_indices ==
                    std::vector<std::size_t>(INACTIVE_JOINT_INDICES.begin(),
                                             INACTIVE_JOINT_INDICES.end()),
                "literal inactive joint signature changed");

        const std::vector<u32> vulnerable =
            failures_literal421(bodies, active_joint);
        std::ofstream vulnerable_out(argv[9]);
        require(static_cast<bool>(vulnerable_out),
                "cannot create vulnerable output");
        vulnerable_out << "body,deleted_pattern\n";
        Fnv body_ledger;
        Fnv row_ledger;
        u64 deleted_incidences = 0;
        std::vector<u64> vulnerable_patterns;
        vulnerable_patterns.reserve(vulnerable.size());
        for (u32 body : vulnerable) {
            u64 pattern = 0;
            for (std::size_t index = 0; index < inactive_joint.size(); ++index)
                if ((body & inactive_joint[index]) == 0)
                    pattern |= u64{1} << index;
            require(pattern > 0 && pattern < 128,
                    "vulnerable body lacks deleted response");
            body_ledger.add(body);
            row_ledger.add(body);
            row_ledger.add(pattern);
            deleted_incidences += std::popcount(pattern);
            vulnerable_patterns.push_back(pattern);
            vulnerable_out << hex8(body) << ',' << std::hex << std::setw(2)
                           << std::setfill('0') << pattern << std::dec
                           << std::setfill(' ') << '\n';
        }
        require(vulnerable_out.good(), "failed writing vulnerable output");

        const BodyAudit closure =
            audit_bodies(bodies, final_partition.active);
        Fnv active_ledger;
        for (u32 mask : final_partition.active) active_ledger.add(mask);
        Fnv inactive_ledger;
        for (u32 mask : final_partition.inactive) inactive_ledger.add(mask);
        u64 response_incidences = 0;
        u64 response_failures = 0;
        u64 minimum_hits = std::numeric_limits<u64>::max();
        u64 maximum_hits = 0;
        u32 minimum_body = 0;
        Fnv response_ledger;
        for (std::size_t row = 0; row < vulnerable.size(); ++row) {
            const u32 body = vulnerable[row];
            std::vector<u32> hits;
            for (u32 mask : final_partition.active)
                if ((body & mask) == 0) hits.push_back(mask);
            response_failures += hits.empty();
            response_incidences += hits.size();
            if (hits.size() < minimum_hits) {
                minimum_hits = hits.size();
                minimum_body = body;
            }
            maximum_hits = std::max<u64>(maximum_hits, hits.size());
            response_ledger.add(body);
            response_ledger.add(vulnerable_patterns[row]);
            response_ledger.add(hits.size());
            for (u32 hit : hits) response_ledger.add(hit);
        }

        require(body_universe_ledger.state ==
                    UINT64_C(0x0d0aed8567285749) &&
                    geometry.grid == INT64_C(291858550663680) &&
                    geometry.cells.size() == 8969 &&
                    geometry.pair_ticks == INT64_C(214426971029040) &&
                    final_partition.active.size() == 5427 &&
                    active_ledger.state == UINT64_C(0xfe7da30cec531ca5) &&
                    final_partition.inactive.size() == 3524 &&
                    inactive_ledger.state == UINT64_C(0x47830f63011b6cca) &&
                    final_partition.equalities == 0 &&
                    final_partition.active_parts ==
                        std::array<u64, 3>{5008, 414, 5} &&
                    final_partition.inactive_parts ==
                        std::array<u64, 3>{3516, 7, 1} &&
                    final_partition.weakest_active_mask ==
                        UINT32_C(0x27048240) &&
                    final_partition.weakest_active_margin ==
                        static_cast<i128>(18878014764LL) &&
                    final_partition.closest_inactive_mask ==
                        UINT32_C(0x0121d204) &&
                    final_partition.closest_inactive_deficit ==
                        static_cast<i128>(6616222578LL) &&
                    active_joint.size() == 414 &&
                    inactive_joint.size() == 7 &&
                    inactive_joint_ledger.state ==
                        UINT64_C(0x38812f899b3636a8) &&
                    vulnerable.size() == 71 &&
                    body_ledger.state == UINT64_C(0x414d30a143d2e22d) &&
                    row_ledger.state == UINT64_C(0x287281de2900cca7) &&
                    deleted_incidences == 83 &&
                    response_incidences == 4575 && response_failures == 0 &&
                    minimum_hits == 7 &&
                    minimum_body == UINT32_C(0x2d306400) &&
                    maximum_hits == 228 &&
                    response_ledger.state == UINT64_C(0x454f40855bd2aac7) &&
                    closure.bodies == UINT64_C(14307150) &&
                    closure.failures == 0 &&
                    closure.checks == UINT64_C(499416624) &&
                    closure.max_checks == UINT64_C(5315) &&
                    closure.worst_body == UINT32_C(0x1d106401),
                "frozen detached carrier663 audit changed");

        std::cout << "THM4281_CARRIER663_DETACHED_RECONSTRUCTION_AUDIT_V1\n"
                  << "FINAL_CARRIER " << carrier.size() << " FNV " << std::hex
                  << fnv_literal421(carrier) << std::dec << '\n'
                  << "BODY_UNIVERSE " << bodies.size() << " FNV " << std::hex
                  << body_universe_ledger.state << std::dec << '\n'
                  << "LITERAL_GEOMETRY PAIR 256,663 GRID " << geometry.grid
                  << " CELLS " << geometry.cells.size() << " PAIR_TICKS "
                  << geometry.pair_ticks << '\n'
                  << "CARRIER_PARTITION ACTIVE "
                  << final_partition.active.size() << " ACTIVE_FNV "
                  << std::hex << active_ledger.state << " INACTIVE "
                  << std::dec << final_partition.inactive.size()
                  << " INACTIVE_FNV " << std::hex << inactive_ledger.state
                  << std::dec << " EQUALITIES " << final_partition.equalities
                  << " ACTIVE_PARTS " << final_partition.active_parts[0] << ','
                  << final_partition.active_parts[1] << ','
                  << final_partition.active_parts[2] << " INACTIVE_PARTS "
                  << final_partition.inactive_parts[0] << ','
                  << final_partition.inactive_parts[1] << ','
                  << final_partition.inactive_parts[2] << '\n'
                  << "WEAKEST_ACTIVE MASK " << std::hex
                  << final_partition.weakest_active_mask << std::dec
                  << " MARGIN_NUM "
                  << decimal(final_partition.weakest_active_margin)
                  << " DEN " << geometry.grid << " REDUCED "
                  << fraction(final_partition.weakest_active_margin,
                              geometry.grid) << '\n'
                  << "CLOSEST_INACTIVE MASK " << std::hex
                  << final_partition.closest_inactive_mask << std::dec
                  << " DEFICIT_NUM "
                  << decimal(final_partition.closest_inactive_deficit)
                  << " DEN " << geometry.grid << " REDUCED "
                  << fraction(final_partition.closest_inactive_deficit,
                              geometry.grid) << '\n'
                  << "JOINT_PARTITION ACTIVE " << active_joint.size()
                  << " INACTIVE " << inactive_joint.size()
                  << " INACTIVE_INDEX_MASK_FNV " << std::hex
                  << inactive_joint_ledger.state << std::dec << '\n'
                  << "VULNERABLE_BODIES " << vulnerable.size() << " BODY_FNV "
                  << std::hex << body_ledger.state << " ROW_FNV "
                  << row_ledger.state << std::dec << " DELETED_INCIDENCES "
                  << deleted_incidences << '\n'
                  << "VULNERABLE_ACTIVE_CARRIER_RESPONSES INCIDENCES "
                  << response_incidences << " FAILURES " << response_failures
                  << " MIN_HITS " << minimum_hits << " MIN_BODY " << std::hex
                  << minimum_body << std::dec << " MAX_HITS " << maximum_hits
                  << " FNV " << std::hex << response_ledger.state << std::dec
                  << '\n'
                  << "FULL_BODY_SCAN BODIES " << closure.bodies
                  << " FAILURES " << closure.failures << " CHECKS "
                  << closure.checks << " MAX_CHECKS " << closure.max_checks
                  << " WORST_BODY " << std::hex << closure.worst_body
                  << std::dec << '\n'
                  << "VERDICT "
                  << (closure.failures == 0 && response_failures == 0
                          ? "PASS DETACHED_LITERAL_FINAL_CARRIER_CLOSES_256_663"
                          : "FAIL")
                  << '\n';
        require(closure.bodies == EXPECTED_BODIES && closure.failures == 0 &&
                    response_failures == 0,
                "final carrier does not close (256,663)");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4281_CARRIER663_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
