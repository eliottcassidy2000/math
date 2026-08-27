// CEGAR descent for the THM-4254 unlabelled repair-mask carrier.
//
// This is intentionally outside the repository.  It reads the 59 frozen
// primary transcripts directly, reconstructs their first-occurrence union,
// audits the first resistant row (416,704), and scans a semantic residual edge
// file with either the original carrier, the smallest one-witness enrichment,
// or all novel masks in the row's order-minimal full-deck prefix.

// The exact endpoint arithmetic comes from the already audited discovery
// implementation.  Union parsing, CEGAR selection, layer grouping, ledgers,
// and stopping logic are independent of the THM-4254 compact checker.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
#define CASCADE_LIBRARY_ONLY
#include "cascade_pair_exhaustive_primary.cpp"
#undef CASCADE_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

namespace {

struct EdgeCegar {
    int q = 0;
    int r = 0;
};

constexpr std::array<EdgeCegar, 59> FROZEN_BAND_CEGAR = {{
    {616,755},{616,756},{616,757},{616,758},{616,759},{616,760},{616,761},
    {616,762},{616,763},{616,764},{616,765},{616,766},{616,767},{616,768},
    {698,755},{698,757},
    {704,755},{704,757},{704,758},{704,759},{704,761},{704,762},{704,763},
    {704,764},{704,765},
    {721,755},{721,757},{721,758},{721,759},{721,761},{721,762},{721,763},
    {721,764},{721,765},{721,766},{721,767},{721,768},
    {726,755},{726,757},{726,758},{726,761},
    {732,755},{732,757},{732,761},{732,762},{732,763},
    {744,762},{744,763},{744,765},{744,766},{744,768},
    {750,762},{750,763},{750,765},{750,766},{750,768},
    {765,766},{765,768},{766,768}
}};

struct PrefixCegar {
    int q = 0;
    int r = 0;
    u64 stated_size = 0;
    u64 stated_fnv = 0;
    std::vector<u32> masks;
};

u64 parse_hex_cegar(const std::string& word) {
    std::size_t used = 0;
    const u64 value = std::stoull(word, &used, 16);
    require(used == word.size(), "bad hexadecimal token");
    return value;
}

PrefixCegar read_prefix_cegar(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open prefix transcript");
    PrefixCegar answer;
    bool saw_pair = false;
    bool saw_prefix = false;
    bool saw_masks = false;
    bool saw_verdict = false;
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("PAIR ", 0) == 0) {
            std::istringstream row(line.substr(5));
            char comma = 0;
            row >> answer.q >> comma >> answer.r;
            require(static_cast<bool>(row) && comma == ',', "bad pair row");
            saw_pair = true;
        } else if (line.rfind("PREFIX_CERT ", 0) == 0) {
            std::istringstream row(line);
            std::string tag, size_tag, fnv_tag, fnv_word;
            row >> tag >> size_tag >> answer.stated_size >> fnv_tag >> fnv_word;
            require(static_cast<bool>(row) && tag == "PREFIX_CERT" &&
                    size_tag == "SIZE" && fnv_tag == "FNV",
                    "bad prefix declaration");
            answer.stated_fnv = parse_hex_cegar(fnv_word);
            saw_prefix = true;
        } else if (line.rfind("PREFIX_MASKS_HEX ", 0) == 0) {
            std::string rest = line.substr(17);
            std::size_t begin = 0;
            while (true) {
                const std::size_t comma = rest.find(',', begin);
                const std::string word = rest.substr(
                    begin, comma == std::string::npos ? std::string::npos
                                                      : comma - begin);
                require(!word.empty(), "empty prefix mask");
                const u64 wide = parse_hex_cegar(word);
                require(wide < (UINT64_C(1) << 30), "mask leaves pool");
                answer.masks.push_back(static_cast<u32>(wide));
                if (comma == std::string::npos) break;
                begin = comma + 1;
            }
            saw_masks = true;
        } else if (line.rfind("VERDICT EVERY_BODY_CLOSED", 0) == 0) {
            saw_verdict = true;
        }
    }
    require(saw_pair && saw_prefix && saw_masks && saw_verdict,
            "incomplete closure transcript");
    require(answer.masks.size() == answer.stated_size,
            "prefix length mismatch");
    FnvLocal ledger;
    std::set<u32> distinct;
    for (u32 mask : answer.masks) {
        require(std::popcount(mask) == 8, "repair is not an eight-set");
        require(distinct.insert(mask).second, "prefix repeats a mask");
        ledger.add(mask);
    }
    require(ledger.state == answer.stated_fnv, "prefix FNV mismatch");
    return answer;
}

std::vector<EdgeCegar> read_edges_cegar(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open residual edge file");
    std::vector<EdgeCegar> edges;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        std::istringstream row(line);
        EdgeCegar edge;
        char comma = 0;
        row >> edge.q >> comma >> edge.r;
        require(static_cast<bool>(row) && comma == ',' && edge.q > 0 &&
                edge.q < edge.r, "bad residual edge row");
        edges.push_back(edge);
    }
    std::sort(edges.begin(), edges.end(), [](const EdgeCegar& a,
                                             const EdgeCegar& b) {
        return a.r != b.r ? a.r > b.r : a.q < b.q;
    });
    require(!edges.empty(), "empty residual edge file");
    return edges;
}

i128 atom_mass_cegar(const AtomData& atoms, u32 repair) {
    i128 mass = 0;
    for (const auto& [failure, value] : atoms.mass) {
        if ((failure & ~repair) == 0) mass += value;
    }
    return mass;
}

void fnv_add_i128_cegar(FnvLocal& ledger, i128 value) {
    const __uint128_t bits = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

struct CarrierScanCegar {
    u64 active = 0;
    ScanResult body;
    u32 best_disjoint_repair = 0;
    i128 best_disjoint_margin = 0;
    i128 denominator = 0;
};

CarrierScanCegar scan_carrier_cegar(const std::vector<Cell>& pool_cells,
                                    int q, int r,
                                    const std::vector<u32>& carrier) {
    const i64 g = std::gcd(q, r);
    const PrimitivePair primitive = build_primitive(q / g, r / g);
    const AtomData atoms = build_cocycle_atoms(pool_cells, primitive, g);
    const i128 denominator = static_cast<i128>(primitive.grid) * g * COMMON;
    std::vector<u32> active;
    active.reserve(carrier.size());
    std::vector<i128> margins;
    margins.reserve(carrier.size());
    for (u32 repair : carrier) {
        const i128 margin = static_cast<i128>(63) *
                                atom_mass_cegar(atoms, repair) -
                            static_cast<i128>(4) * denominator;
        margins.push_back(margin);
        if (margin >= 0) active.push_back(repair);
    }
    CarrierScanCegar answer;
    answer.active = active.size();
    answer.body = scan_bodies(active);
    answer.denominator = denominator;
    if (answer.body.failures != 0) {
        bool first = true;
        for (std::size_t index = 0; index < carrier.size(); ++index) {
            if ((carrier[index] & answer.body.first_failure) != 0) continue;
            if (first || margins[index] > answer.best_disjoint_margin ||
                (margins[index] == answer.best_disjoint_margin &&
                 carrier[index] < answer.best_disjoint_repair)) {
                first = false;
                answer.best_disjoint_margin = margins[index];
                answer.best_disjoint_repair = carrier[index];
            }
        }
        require(!first && answer.best_disjoint_margin < 0,
                "failure has active disjoint carrier repair");
    }
    return answer;
}

std::string mask_hex_cegar(u32 mask) {
    std::ostringstream out;
    out << std::hex << mask;
    return out.str();
}

}  // namespace

#ifndef CARRIER_CEGAR_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        require(argc == 8,
                "usage: cegar PACKET_ROOT FULL_416_704 RESIDUAL_CSV "
                "START_ENDPOINT MIN_ENDPOINT MODE base_compare(0|1)");
        const std::filesystem::path packet(argv[1]);
        const PrefixCegar row_prefix = read_prefix_cegar(argv[2]);
        require(row_prefix.q == 416 && row_prefix.r == 704 &&
                row_prefix.stated_size == 2608 &&
                row_prefix.stated_fnv == UINT64_C(0x18ff663a123e684e),
                "416,704 discovery prefix changed");
        const int start_endpoint = std::stoi(argv[4]);
        const int min_endpoint = std::stoi(argv[5]);
        const std::string mode(argv[6]);
        const bool compare_base = std::stoi(argv[7]) != 0;
        require(start_endpoint >= min_endpoint, "bad endpoint range");
        require(mode == "base" || mode == "one" || mode == "prefix",
                "unknown enrichment mode");

        std::vector<u32> base_carrier;
        std::set<u32> base_seen;
        FnvLocal incidence;
        u64 incidences = 0;
        for (const EdgeCegar& expected : FROZEN_BAND_CEGAR) {
            const auto path = packet / "replay_band" /
                (std::to_string(expected.q) + "_" +
                 std::to_string(expected.r) + ".out");
            const PrefixCegar prefix = read_prefix_cegar(path);
            require(prefix.q == expected.q && prefix.r == expected.r,
                    "band transcript/pair mismatch");
            incidence.add(static_cast<u64>(prefix.q));
            incidence.add(static_cast<u64>(prefix.r));
            incidence.add(prefix.masks.size());
            for (u32 mask : prefix.masks) {
                incidence.add(mask);
                ++incidences;
                if (base_seen.insert(mask).second) base_carrier.push_back(mask);
            }
        }
        FnvLocal base_ledger;
        for (u32 mask : base_carrier) base_ledger.add(mask);
        require(incidences == 56419 &&
                incidence.state == UINT64_C(0xacc8347addf27ac3) &&
                base_carrier.size() == 4675 &&
                base_ledger.state == UINT64_C(0xce4e76ec11df057c),
                "base carrier fingerprint changed");

        std::vector<u32> novel;
        std::set<u32> all_seen = base_seen;
        for (u32 mask : row_prefix.masks) {
            if (all_seen.insert(mask).second) novel.push_back(mask);
        }
        FnvLocal novel_ledger;
        for (u32 mask : novel) novel_ledger.add(mask);
        require(novel.size() == 58, "novel-prefix count changed");

        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "pool-cell count changed");
        const i64 row_g = std::gcd(row_prefix.q, row_prefix.r);
        const PrimitivePair row_primitive = build_primitive(
            row_prefix.q / row_g, row_prefix.r / row_g);
        const AtomData row_atoms = build_cocycle_atoms(
            pool_cells, row_primitive, row_g);
        const i128 row_denominator =
            static_cast<i128>(row_primitive.grid) * row_g * COMMON;
        const ScanResult full_prefix_scan = scan_bodies(row_prefix.masks);
        require(full_prefix_scan.failures == 0 &&
                full_prefix_scan.max_checks == row_prefix.masks.size(),
                "full row prefix lost closure/order minimality");
        const CarrierScanCegar base_row = scan_carrier_cegar(
            pool_cells, 416, 704, base_carrier);
        require(base_row.body.failures == 1 &&
                base_row.body.first_failure == UINT32_C(0x2f903000),
                "base 416,704 stopping witness changed");

        FnvLocal atom_component_ledger;
        u32 first_witness = 0;
        for (u32 mask : novel) {
            const i128 mass = atom_mass_cegar(row_atoms, mask);
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) * row_denominator;
            require(margin >= 0, "novel prefix mask is inactive at seed row");
            const ComponentReport component = component_report(
                pool_cells, row_primitive, row_g, mask);
            require(component.exact_mass_num == mass,
                    "atom/component control mismatch");
            atom_component_ledger.add(mask);
            fnv_add_i128_cegar(atom_component_ledger, mass);
            fnv_add_i128_cegar(atom_component_ledger, margin);
            if (first_witness == 0 &&
                (mask & base_row.body.first_failure) == 0) {
                first_witness = mask;
            }
        }
        require(first_witness != 0, "no novel counterexample repair");

        std::vector<u32> one_carrier = base_carrier;
        one_carrier.push_back(first_witness);
        const CarrierScanCegar one_row = scan_carrier_cegar(
            pool_cells, 416, 704, one_carrier);
        require(one_row.body.failures == 0,
                "one-witness enrichment does not close seed row");

        std::vector<u32> prefix_carrier = base_carrier;
        prefix_carrier.insert(prefix_carrier.end(), novel.begin(), novel.end());
        FnvLocal enriched_ledger;
        for (u32 mask : prefix_carrier) enriched_ledger.add(mask);
        const CarrierScanCegar prefix_row = scan_carrier_cegar(
            pool_cells, 416, 704, prefix_carrier);
        require(prefix_row.body.failures == 0,
                "prefix enrichment does not close seed row");

        const std::vector<u32>* carrier = &base_carrier;
        if (mode == "one") carrier = &one_carrier;
        if (mode == "prefix") carrier = &prefix_carrier;
        FnvLocal selected_ledger;
        for (u32 mask : *carrier) selected_ledger.add(mask);

        std::cout << "CARRIER_CEGAR_DESCENT_V1\n";
        std::cout << "BASE PREFIX_INCIDENCES " << incidences
                  << " INCIDENCE_FNV " << std::hex << incidence.state
                  << " UNION " << std::dec << base_carrier.size()
                  << " UNION_FNV " << std::hex << base_ledger.state << std::dec
                  << "\n";
        std::cout << "SEED PAIR 416,704 FULL_PREFIX "
                  << row_prefix.masks.size() << " FULL_PREFIX_FNV " << std::hex
                  << row_prefix.stated_fnv << std::dec
                  << " OVERLAP " << row_prefix.masks.size() - novel.size()
                  << " NOVEL " << novel.size() << " NOVEL_FNV " << std::hex
                  << novel_ledger.state << std::dec
                  << " FULL_PREFIX_MAX_CHECKS " << full_prefix_scan.max_checks
                  << "\n";
        std::cout << "SEED_BASE ACTIVE " << base_row.active
                  << " FAILURES " << base_row.body.failures
                  << " FIRST_FAILURE " << mask_hex_cegar(
                         base_row.body.first_failure)
                  << " FIRST_FAILURE_LABELS {" << labels(
                         base_row.body.first_failure) << "}"
                  << " BEST_DISJOINT_REPAIR " << mask_hex_cegar(
                         base_row.best_disjoint_repair)
                  << " BEST_DISJOINT_LABELS {" << labels(
                         base_row.best_disjoint_repair) << "}"
                  << " BEST_DISJOINT_MARGIN_NUM "
                  << decimal(base_row.best_disjoint_margin)
                  << " DEN " << decimal(base_row.denominator) << "\n";
        std::cout << "CEGAR_FIRST_WITNESS " << mask_hex_cegar(first_witness)
                  << " LABELS {" << labels(first_witness) << "}"
                  << " ONE_UNION " << one_carrier.size()
                  << " ONE_ROW_FAILURES " << one_row.body.failures
                  << " PREFIX_UNION " << prefix_carrier.size()
                  << " PREFIX_UNION_FNV " << std::hex << enriched_ledger.state
                  << " ATOM_COMPONENT_LEDGER " << atom_component_ledger.state
                  << std::dec << " PREFIX_ROW_FAILURES "
                  << prefix_row.body.failures << "\n";
        std::cout << "SELECTED MODE " << mode << " UNION " << carrier->size()
                  << " FNV " << std::hex << selected_ledger.state << std::dec
                  << " RANGE " << start_endpoint << ".." << min_endpoint
                  << " BASE_COMPARE " << compare_base << "\n";

        const std::vector<EdgeCegar> all_edges = read_edges_cegar(argv[3]);
        int current_layer = std::numeric_limits<int>::max();
        u64 layer_rows = 0;
        u64 layer_resistant = 0;
        u64 layer_base_resistant = 0;
        u64 layer_active_sum = 0;
        u64 layer_checks = 0;
        u64 layer_max_checks = 0;
        FnvLocal layer_ledger;
        std::vector<std::string> resistant_lines;
        bool found_layer = false;

        auto flush_layer = [&]() {
            if (current_layer == std::numeric_limits<int>::max()) return false;
            std::cout << "LAYER " << current_layer << " ROWS " << layer_rows
                      << " RESISTANT " << layer_resistant
                      << " BASE_RESISTANT " << layer_base_resistant
                      << " ACTIVE_SUM " << layer_active_sum
                      << " CHECKS " << layer_checks
                      << " MAX_CHECKS " << layer_max_checks
                      << " ROW_LEDGER " << std::hex << layer_ledger.state
                      << std::dec << "\n";
            for (const std::string& line : resistant_lines) {
                std::cout << line << '\n';
            }
            return layer_resistant != 0;
        };

        for (const EdgeCegar& edge : all_edges) {
            if (edge.r > start_endpoint || edge.r < min_endpoint) continue;
            if (current_layer != edge.r) {
                if (flush_layer()) {
                    found_layer = true;
                    break;
                }
                current_layer = edge.r;
                layer_rows = layer_resistant = layer_base_resistant = 0;
                layer_active_sum = layer_checks = layer_max_checks = 0;
                layer_ledger = FnvLocal{};
                resistant_lines.clear();
            }
            const CarrierScanCegar selected = scan_carrier_cegar(
                pool_cells, edge.q, edge.r, *carrier);
            u64 base_failures = 0;
            if (compare_base) {
                const CarrierScanCegar base = scan_carrier_cegar(
                    pool_cells, edge.q, edge.r, base_carrier);
                base_failures = base.body.failures;
                if (base_failures != 0) ++layer_base_resistant;
            }
            ++layer_rows;
            if (selected.body.failures != 0) ++layer_resistant;
            layer_active_sum += selected.active;
            layer_checks += selected.body.checks;
            layer_max_checks = std::max(layer_max_checks,
                                        selected.body.max_checks);
            layer_ledger.add(static_cast<u64>(edge.q));
            layer_ledger.add(static_cast<u64>(edge.r));
            layer_ledger.add(selected.active);
            layer_ledger.add(selected.body.failures);
            layer_ledger.add(selected.body.checks);
            layer_ledger.add(selected.body.max_checks);
            layer_ledger.add(selected.body.first_failure);
            layer_ledger.add(selected.best_disjoint_repair);
            fnv_add_i128_cegar(layer_ledger,
                              selected.best_disjoint_margin);
            fnv_add_i128_cegar(layer_ledger, selected.denominator);
            if (selected.body.failures != 0) {
                std::ostringstream row;
                row << "RESISTANT PAIR " << edge.q << ',' << edge.r
                    << " ACTIVE " << selected.active
                    << " FAILURES " << selected.body.failures
                    << " FIRST_FAILURE " << std::hex
                    << selected.body.first_failure << std::dec
                    << " FIRST_FAILURE_LABELS {"
                    << labels(selected.body.first_failure) << "}"
                    << " BEST_DISJOINT_REPAIR " << std::hex
                    << selected.best_disjoint_repair << std::dec
                    << " BEST_DISJOINT_LABELS {"
                    << labels(selected.best_disjoint_repair) << "}"
                    << " BEST_DISJOINT_MARGIN_NUM "
                    << decimal(selected.best_disjoint_margin)
                    << " DEN " << decimal(selected.denominator)
                    << " BASE_FAILURES " << base_failures;
                resistant_lines.push_back(row.str());
            }
        }
        if (!found_layer) found_layer = flush_layer();
        std::cout << "STOP FIRST_RESISTANT_LAYER "
                  << (found_layer ? current_layer : -1)
                  << " VERDICT " << (found_layer ? "BOUNDARY_FOUND"
                                                   : "RANGE_CLOSED")
                  << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CARRIER_CEGAR_ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif
