// Independent exact certificate consumer for the support<=350 nonjoint
// deletion, its 84 residual witness hyperedges, and the 37-mask repair.

#define THRESHOLD_DELETION_REPLAY_MAIN threshold_replay_hidden_main_for_k350
#include "threshold_deletion_replay.cpp"
#undef THRESHOLD_DELETION_REPLAY_MAIN

#include <array>

namespace {

constexpr std::size_t kObligations350 = 84;

struct ObligationCert350 {
    PairAgent pair;
    u32 body = 0;
};

struct PatternCert350 {
    std::array<u64, 2> words{};
    bool operator<(const PatternCert350& other) const {
        return words < other.words;
    }
    bool operator==(const PatternCert350& other) const = default;
    void set(std::size_t bit) {
        words[bit / 64] |= UINT64_C(1) << (bit % 64);
    }
    bool contains(std::size_t bit) const {
        return (words[bit / 64] >> (bit % 64)) & UINT64_C(1);
    }
    bool empty() const { return words[0] == 0 && words[1] == 0; }
};

std::string hex16_cert350(u64 value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

std::vector<ObligationCert350> read_obligations_cert350(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open k350 failures");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "k350 failure header changed");
    std::vector<ObligationCert350> out;
    std::set<std::tuple<int, int, u32>> distinct;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank k350 failure row");
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token, extra;
        require(static_cast<bool>(fields >> q >> r >> token) &&
                    !(fields >> extra),
                "malformed k350 failure row");
        const u32 body = parse_mask_agent(token);
        require(q > 0 && q < r && std::popcount(body) == 9 &&
                    distinct.insert({q, r, body}).second,
                "invalid/duplicate k350 failure");
        out.push_back({PairAgent{q, r}, body});
    }
    require(input.eof() && out.size() == kObligations350,
            "k350 obligation count changed");
    return out;
}

std::vector<unsigned> read_packing_cert350(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open k350 packing");
    std::vector<unsigned> packing;
    std::set<unsigned> distinct;
    unsigned index = 0;
    while (input >> index) {
        require(index < kObligations350 && distinct.insert(index).second,
                "invalid/duplicate k350 packing index");
        packing.push_back(index);
    }
    require(input.eof() && packing.size() == 37,
            "k350 packing size changed");
    return packing;
}

u64 masks_fnv_cert350(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 18,
                "usage: k350_cert JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL "
                "OLD_UNION REPAIRS76 ADDITIONS4 DELETE73 ADDITIONS10 "
                "DELETE350 FAILURES84 RETAIN37 REPAIRED_DELETE5104 "
                "PACKING37 ATLAS_OUT EDGE_OUT");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> c596 = reconstruct_c596_thr(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9]);
        const std::vector<u32> additions10 =
            read_mixed617(argv[10], 10, kAddition10FnvThr);
        const std::vector<u32> deletion350 = read_threshold_ledger(argv[11]);
        const std::vector<ObligationCert350> obligations =
            read_obligations_cert350(argv[12]);
        const std::vector<u32> retention = read_threshold_ledger(argv[13]);
        const std::vector<u32> repaired_deletion =
            read_threshold_ledger(argv[14]);
        const std::vector<unsigned> packing =
            read_packing_cert350(argv[15]);
        require(deletion350.size() == 5141 && retention.size() == 37 &&
                    repaired_deletion.size() == 5104,
                "k350 deletion/retention counts changed");

        const std::set<PairKeyThr> old_union =
            read_pair_set_thr(argv[6], 1624, kOldUnionFnvThr);
        const std::vector<PairAgent> raw_band = read_band_agent(argv[5], 595);
        require(raw_band.size() == 394, "raw r>=595 band changed");
        std::vector<PairAgent> target;
        std::set<PairKeyThr> target_set;
        for (PairAgent pair : raw_band)
            if (pair.r >= 596 || !old_union.contains({pair.q, pair.r})) {
                target.push_back(pair);
                target_set.insert({pair.q, pair.r});
            }
        require(target.size() == 391, "target changed");
        std::vector<Geometry> target_geometries;
        for (PairAgent pair : target)
            target_geometries.push_back(build_geometry(pair.q, pair.r));
        for (const auto& obligation : obligations)
            require(target_set.contains(
                        {obligation.pair.q, obligation.pair.r}),
                    "k350 failure escaped target");

        // Independent complete support scan and exact threshold equality.
        std::vector<u32> expected_deletion;
        Fnv sign_ledger, support_ledger;
        u64 equalities = 0;
        for (u32 mask : c596) {
            u64 active_rows = 0;
            Fnv one_support;
            for (std::size_t row = 0; row < target.size(); ++row) {
                const Margin exact = margin(target_geometries[row], mask);
                const bool active = exact.ticks >= 0;
                sign_ledger.add(mask);
                sign_ledger.add(active);
                equalities += exact.ticks == 0;
                if (active) {
                    ++active_rows;
                    one_support.add(target[row].q);
                    one_support.add(target[row].r);
                }
            }
            support_ledger.add(mask);
            support_ledger.add(active_rows);
            support_ledger.add(one_support.state);
            if (!joint_set.contains(mask) && active_rows <= 350)
                expected_deletion.push_back(mask);
        }
        require(equalities == 0 && deletion350 == expected_deletion,
                "ledger is not exact protected-nonjoint support<=350 set");

        const std::unordered_set<u32> deletion_set(deletion350.begin(),
                                                   deletion350.end());
        const std::unordered_set<u32> retention_set(retention.begin(),
                                                    retention.end());
        for (u32 mask : retention)
            require(deletion_set.contains(mask), "retention escaped D350");
        std::vector<u32> expected_repaired_deletion;
        for (u32 mask : deletion350)
            if (!retention_set.contains(mask))
                expected_repaired_deletion.push_back(mask);
        require(repaired_deletion == expected_repaired_deletion,
                "repaired deletion is not D350 minus retention in D350 order");

        std::vector<u32> deleted_carrier;
        for (u32 mask : c596)
            if (!deletion_set.contains(mask)) deleted_carrier.push_back(mask);
        for (u32 mask : additions10) deleted_carrier.push_back(mask);
        std::vector<Geometry> obligation_geometries;
        for (const auto& obligation : obligations)
            obligation_geometries.push_back(
                build_geometry(obligation.pair.q, obligation.pair.r));

        // Verify that every imported failure really has zero survivor in the
        // D350 carrier, then reconstruct all residual responses from D350.
        for (std::size_t index = 0; index < obligations.size(); ++index)
            for (u32 mask : deleted_carrier)
                require(margin(obligation_geometries[index], mask).ticks < 0 ||
                            (mask & obligations[index].body) != 0,
                        "imported failure already has a D350 survivor");

        std::map<PatternCert350, std::vector<u32>> types;
        std::map<u32, PatternCert350> mask_patterns;
        Fnv responder_ledger;
        for (u32 mask : deletion350) {
            PatternCert350 pattern;
            for (std::size_t index = 0; index < obligations.size(); ++index)
                if ((mask & obligations[index].body) == 0 &&
                    margin(obligation_geometries[index], mask).ticks >= 0)
                    pattern.set(index);
            if (pattern.empty()) continue;
            responder_ledger.add(mask);
            types[pattern].push_back(mask);
            mask_patterns.emplace(mask, pattern);
        }

        std::ofstream atlas_out(argv[16]);
        std::ofstream edge_out(argv[17]);
        require(static_cast<bool>(atlas_out) && static_cast<bool>(edge_out),
                "cannot create independent k350 ledgers");
        atlas_out << "w0,w1,count,least_mask\n";
        Fnv type_ledger;
        u64 responders = 0;
        for (const auto& [pattern, masks] : types) {
            responders += masks.size();
            const u32 least = *std::min_element(masks.begin(), masks.end());
            type_ledger.add(pattern.words[0]);
            type_ledger.add(pattern.words[1]);
            type_ledger.add(masks.size());
            type_ledger.add(least);
            atlas_out << hex16_cert350(pattern.words[0]) << ','
                      << hex16_cert350(pattern.words[1]) << ',' << masks.size()
                      << ',' << hex8(least) << '\n';
        }
        require(responders == 70 && types.size() == 53,
                "k350 responder quotient changed");

        edge_out << "index,q,r,body_hex,witness_count,witness_fnv\n";
        Fnv obligation_ledger, edge_ledger;
        for (std::size_t index = 0; index < obligations.size(); ++index) {
            std::vector<u32> witnesses;
            for (const auto& [mask, pattern] : mask_patterns)
                if (pattern.contains(index)) witnesses.push_back(mask);
            require(!witnesses.empty(), "empty residual witness edge");
            Fnv witness_ledger;
            for (u32 mask : witnesses) witness_ledger.add(mask);
            const auto& obligation = obligations[index];
            obligation_ledger.add(obligation.pair.q);
            obligation_ledger.add(obligation.pair.r);
            obligation_ledger.add(obligation.body);
            edge_ledger.add(obligation.pair.q);
            edge_ledger.add(obligation.pair.r);
            edge_ledger.add(obligation.body);
            edge_ledger.add(witnesses.size());
            edge_ledger.add(witness_ledger.state);
            edge_out << index << ',' << obligation.pair.q << ','
                     << obligation.pair.r << ',' << hex8(obligation.body) << ','
                     << witnesses.size() << ',' << std::hex
                     << witness_ledger.state << std::dec << '\n';
        }
        require(atlas_out.good() && edge_out.good(),
                "independent k350 ledger write failed");

        // Upper certificate: the supplied 37 retained masks cover all 84
        // residual obligations.
        PatternCert350 cover_union;
        Fnv retention_response_ledger;
        for (u32 mask : retention) {
            const auto found = mask_patterns.find(mask);
            require(found != mask_patterns.end(),
                    "retained mask has empty residual response");
            cover_union.words[0] |= found->second.words[0];
            cover_union.words[1] |= found->second.words[1];
            retention_response_ledger.add(mask);
            retention_response_ledger.add(found->second.words[0]);
            retention_response_ledger.add(found->second.words[1]);
        }
        require(cover_union.words[0] == UINT64_MAX &&
                    cover_union.words[1] == (UINT64_C(1) << 20) - 1,
                "retention37 does not cover all 84 obligations");

        // Lower certificate: every possible D350 response hits at most one of
        // these 37 packed obligations.  Therefore any retention cover has at
        // least 37 masks.
        for (const auto& [pattern, masks] : types) {
            (void)masks;
            unsigned hits = 0;
            for (unsigned index : packing) hits += pattern.contains(index);
            require(hits <= 1, "packing collision in a response type");
        }
        Fnv packing_ledger;
        for (unsigned index : packing) packing_ledger.add(index);

        std::vector<u32> final_carrier;
        const std::unordered_set<u32> repaired_set(repaired_deletion.begin(),
                                                   repaired_deletion.end());
        for (u32 mask : c596)
            if (!repaired_set.contains(mask)) final_carrier.push_back(mask);
        for (u32 mask : additions10) final_carrier.push_back(mask);
        const u64 rank8 = std::count_if(
            final_carrier.begin(), final_carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });

        std::cout << "LRC14_K350_RETENTION_CERTIFICATE_V1\n"
                  << "TARGET_ROWS 391 C596 9019 C596_FNV " << std::hex
                  << masks_fnv_agent(c596) << std::dec << '\n'
                  << "SUPPORT_SIGN_CELLS " << c596.size() * target.size()
                  << " SIGN_FNV " << std::hex << sign_ledger.state << std::dec
                  << " EQUALITIES " << equalities << " SUPPORT_LEDGER_FNV "
                  << std::hex << support_ledger.state << std::dec << '\n'
                  << "DELETE350 " << deletion350.size() << " FNV " << std::hex
                  << masks_fnv_cert350(deletion350) << std::dec
                  << " EXACT_NONJOINT_THRESHOLD 1\n"
                  << "OBLIGATIONS " << obligations.size() << " FNV "
                  << std::hex << obligation_ledger.state << std::dec
                  << " EDGE_FNV " << std::hex << edge_ledger.state << std::dec
                  << '\n'
                  << "RESPONDERS " << responders << " FNV " << std::hex
                  << responder_ledger.state << std::dec << " TYPES "
                  << types.size() << " TYPE_FNV " << std::hex
                  << type_ledger.state << std::dec << '\n'
                  << "COVER 37 UNION " << std::hex << cover_union.words[0]
                  << ':' << cover_union.words[1] << std::dec
                  << " RESPONSE_LEDGER_FNV " << std::hex
                  << retention_response_ledger.state << std::dec << '\n'
                  << "PACKING 37 FNV " << std::hex << packing_ledger.state
                  << std::dec << " MAX_RESPONSE_HITS 1\n"
                  << "MIN_RETENTION 37 LOWER_PACKING 37 UPPER_COVER 37\n"
                  << "REPAIRED_DELETE " << repaired_deletion.size() << " FNV "
                  << std::hex << masks_fnv_cert350(repaired_deletion) << std::dec
                  << " FINAL_CARRIER " << final_carrier.size() << " FNV "
                  << std::hex << masks_fnv_agent(final_carrier) << std::dec
                  << " RANK8 " << rank8 << " RANK9 "
                  << final_carrier.size() - rank8 << " JOINT_RETAINED 421\n"
                  << "SCOPE EXACT_NONJOINT_SUPPORT350_RESIDUAL_HYPERGRAPH_"
                     "PACKING_COVER_CERTIFICATE_FIXED_391_ROWS_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "K350_CERTIFICATE_ERROR " << error.what() << '\n';
        return 1;
    }
}
