// Detached audit of the endpoint636 agent's exact 14-for-14 carrier exchange.
// The included maintained source imports no project code and reconstructs
// literal walls and all rank-nine bodies directly.

#define main thm4287_original_detached_main
#include "04-computation/lrc14_repaired_carrier_endpoint637_descent_thm4287/detached_endpoint637_literal_audit.cpp"
#undef main

#include <sstream>

namespace {

constexpr std::array<u32, 14> kExchangeDelete = {
    0x00003e1a, 0x000132a3, 0x00017464, 0x00033388,
    0x000a16c2, 0x000f8118, 0x00142a1a, 0x00154348,
    0x00184ba0, 0x001aa260, 0x00202c2b, 0x002066a4,
    0x002b018a, 0x0030c2a2};
constexpr std::array<u32, 14> kExchangeAdd = {
    0x18468880, 0x080e8281, 0x22081017, 0x08422a82,
    0x004cac40, 0x19c04044, 0x00c08ec0, 0x10443016,
    0x01609124, 0x10413209, 0x01611640, 0x00606449,
    0x0128d084, 0x08806449};
constexpr u64 kExchangeDeleteFnv = UINT64_C(0xa497f155f01aee9e);
constexpr u64 kExchangeAddFnv = UINT64_C(0x8c648463d5cede1b);
constexpr u64 kExchangeCarrierFnv = UINT64_C(0x8062ce6d5728da1f);
constexpr std::array<u64, 9> kExpectedActiveCounts = {
    5670, 7666, 7646, 3567, 3249, 7972, 8493, 6898, 6089};
constexpr std::array<u64, 9> kExpectedActiveFnv = {
    UINT64_C(0xc663178b4940d1f2), UINT64_C(0xea36324274e62ea6),
    UINT64_C(0x19b2fa64bf54036e), UINT64_C(0x16c10d3de13f45a4),
    UINT64_C(0x2d1b9c1611108c9d), UINT64_C(0x4557b5dad6866cd1),
    UINT64_C(0xec6edd802340bcaf), UINT64_C(0xda97157e50969930),
    UINT64_C(0x387b60e510c195de)};

struct ExchangeFailures {
    std::vector<u32> a;
    std::vector<u32> b;
};

ExchangeFailures read_exchange_failures_literal(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open exchange failure ledger");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "exchange failure header changed");
    ExchangeFailures out;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int q = 0;
        int r = 0;
        std::string token;
        row >> q >> r >> token;
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(row && used == token.size() && r == 636 &&
                    wide < (UINT64_C(1) << 30) &&
                    std::popcount(static_cast<u32>(wide)) == 9,
                "malformed exchange failure body");
        if (q == 100) out.a.push_back(static_cast<u32>(wide));
        else if (q == 256) out.b.push_back(static_cast<u32>(wide));
        else fail("unexpected exchange failure pair");
    }
    require(out.a.size() == 64 &&
                body_fnv(out.a) == UINT64_C(0xd611500ea833ff83) &&
                out.b.size() == 37 &&
                body_fnv(out.b) == UINT64_C(0xee7792a8a2fd51c9),
            "exchange failure families changed");
    return out;
}

std::vector<u32> build_original_carrier(
    const std::filesystem::path& inherited_path,
    const std::filesystem::path& additions_path,
    const std::filesystem::path& witness_path) {
    std::vector<u32> carrier =
        read_masks(inherited_path, kInheritedCount, kInheritedFnv);
    const std::vector<u32> additions =
        read_masks(additions_path, kAdditionCount, kAdditionsFnv);
    const std::vector<u32> witness =
        read_masks(witness_path, kWitnessCount, kWitnessFnv);
    std::set<u32> distinct(carrier.begin(), carrier.end());
    for (u32 mask : additions) {
        require(distinct.insert(mask).second, "addition overlaps carrier");
        carrier.push_back(mask);
    }
    require(distinct.insert(kPriorRepair).second, "prior repair overlaps");
    carrier.push_back(kPriorRepair);
    for (u32 mask : witness) {
        require(distinct.insert(mask).second, "witness overlaps carrier");
        carrier.push_back(mask);
    }
    require(carrier.size() == kRepairedCount &&
                mask_fnv(carrier) == kRepairedFnv,
            "original repaired carrier changed");
    return carrier;
}

std::vector<i128> literal_margins(Pair pair,
                                  const std::vector<u32>& masks) {
    const Geometry geometry = build_geometry(pair);
    const CarrierIndex index = build_index(masks);
    const std::vector<i64> mass = classify(geometry, masks, index);
    std::vector<i128> margins;
    margins.reserve(masks.size());
    for (i64 value : mass)
        margins.push_back(static_cast<i128>(63) * value -
                          static_cast<i128>(4) * geometry.grid);
    return margins;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 8,
                "usage: detached-exchange BASE8951 ADD45 WITNESS9 JOINT421 "
                "CURRENT22682 FAILURES THREADS");
        const std::vector<u32> original =
            build_original_carrier(argv[1], argv[2], argv[3]);
        const std::vector<u32> joint =
            read_masks(argv[4], 421, UINT64_C(0x20d63dd42fe8150e));
        const std::set<u32> original_set(original.begin(), original.end());
        const std::set<u32> joint_set(joint.begin(), joint.end());
        Fnv delete_ledger;
        Fnv add_ledger;
        std::set<u32> deleted;
        for (u32 mask : kExchangeDelete) {
            require(original_set.contains(mask) && !joint_set.contains(mask) &&
                        deleted.insert(mask).second,
                    "invalid exchange deletion");
            delete_ledger.add(mask);
        }
        require(delete_ledger.state == kExchangeDeleteFnv,
                "exchange deletion identity changed");
        std::set<u32> added;
        for (u32 mask : kExchangeAdd) {
            require(!original_set.contains(mask) && added.insert(mask).second,
                    "invalid exchange addition");
            add_ledger.add(mask);
        }
        require(add_ledger.state == kExchangeAddFnv,
                "exchange addition identity changed");
        std::vector<u32> replaced;
        for (u32 mask : original)
            if (!deleted.contains(mask)) replaced.push_back(mask);
        replaced.insert(replaced.end(), kExchangeAdd.begin(), kExchangeAdd.end());
        require(replaced.size() == 9006 &&
                    mask_fnv(replaced) == kExchangeCarrierFnv,
                "exchange carrier identity changed");

        const std::vector<Pair> band = read_boundary_band(argv[5]);
        require(band.size() == 9, "exchange boundary changed");
        std::cout << "INCOMING_SIGNAL_DETACHED_EXCHANGE_AUDIT_V1\n"
                  << "EXCHANGE DELETE 14 FNV " << std::hex
                  << delete_ledger.state << " ADD 14 FNV " << add_ledger.state
                  << " CARRIER 9006 FNV " << mask_fnv(replaced) << std::dec
                  << '\n';

        Fnv deletion_margin_ledger;
        std::vector<i128> deletion_min(kExchangeDelete.size());
        std::vector<i128> deletion_max(kExchangeDelete.size());
        bool first_layer = true;
        u64 deletion_equalities = 0;
        for (Pair pair : band) {
            const std::vector<u32> masks(kExchangeDelete.begin(),
                                         kExchangeDelete.end());
            const std::vector<i128> margins = literal_margins(pair, masks);
            for (std::size_t index = 0; index < margins.size(); ++index) {
                require(margins[index] < 0,
                        "deleted mask is not strictly boundary-inactive");
                deletion_equalities += margins[index] == 0;
                if (first_layer || margins[index] < deletion_min[index])
                    deletion_min[index] = margins[index];
                if (first_layer || margins[index] > deletion_max[index])
                    deletion_max[index] = margins[index];
                deletion_margin_ledger.add(pair.q);
                deletion_margin_ledger.add(pair.r);
                deletion_margin_ledger.add(kExchangeDelete[index]);
                add_i128(deletion_margin_ledger, margins[index]);
            }
            first_layer = false;
        }
        for (std::size_t index = 0; index < kExchangeDelete.size(); ++index)
            std::cout << "DELETE_MARGIN " << std::hex << std::setw(8)
                      << std::setfill('0') << kExchangeDelete[index] << std::dec
                      << std::setfill(' ') << " RANGE "
                      << decimal(deletion_min[index]) << ".."
                      << decimal(deletion_max[index]) << '\n';
        std::cout << "DELETE_STRICT_CELLS "
                  << band.size() * kExchangeDelete.size()
                  << " EQUALITIES " << deletion_equalities << " FNV "
                  << std::hex << deletion_margin_ledger.state << std::dec
                  << '\n';

        const ExchangeFailures failures =
            read_exchange_failures_literal(argv[6]);
        const std::vector<u32> additions(kExchangeAdd.begin(),
                                         kExchangeAdd.end());
        const std::vector<i128> margins_a =
            literal_margins({100, 636}, additions);
        const std::vector<i128> margins_b =
            literal_margins({256, 636}, additions);
        u64 cover_a = 0;
        u64 cover_b = 0;
        Fnv response_ledger;
        for (std::size_t index = 0; index < additions.size(); ++index) {
            u64 response_a = 0;
            u64 response_b = 0;
            if (margins_a[index] >= 0)
                for (std::size_t body = 0; body < failures.a.size(); ++body)
                    if ((additions[index] & failures.a[body]) == 0)
                        response_a |= UINT64_C(1) << body;
            if (margins_b[index] >= 0)
                for (std::size_t body = 0; body < failures.b.size(); ++body)
                    if ((additions[index] & failures.b[body]) == 0)
                        response_b |= UINT64_C(1) << body;
            cover_a |= response_a;
            cover_b |= response_b;
            response_ledger.add(additions[index]);
            response_ledger.add(response_a);
            response_ledger.add(response_b);
            add_i128(response_ledger, margins_a[index]);
            add_i128(response_ledger, margins_b[index]);
            std::cout << "ADD_RESPONSE " << std::hex << std::setw(8)
                      << std::setfill('0') << additions[index] << " LO "
                      << std::setw(16) << response_a << " HI " << std::setw(10)
                      << response_b << std::dec << std::setfill(' ')
                      << " MARGIN_A " << decimal(margins_a[index])
                      << " MARGIN_B " << decimal(margins_b[index]) << '\n';
        }
        require(cover_a == ~UINT64_C(0) &&
                    cover_b == (UINT64_C(1) << 37) - 1,
                "exchange additions do not cover 101 hostile bodies");
        std::cout << "ADD_COVER 64,37 RESPONSE_MARGIN_FNV " << std::hex
                  << response_ledger.state << std::dec << '\n';

        const CarrierIndex replaced_index = build_index(replaced);
        const unsigned requested_threads =
            static_cast<unsigned>(std::stoul(argv[7]));
        require(requested_threads >= 1 && requested_threads <= 8,
                "thread count outside 1..8");
        std::vector<PairAudit> audits(band.size());
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        for (unsigned thread = 0; thread < requested_threads; ++thread) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= band.size()) break;
                    audits[index] = audit_pair(band[index], replaced,
                                               replaced_index);
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        Fnv boundary_ledger;
        u64 total_bodies = 0;
        u64 total_checks = 0;
        u64 total_failures = 0;
        u64 total_equalities = 0;
        for (std::size_t index = 0; index < audits.size(); ++index) {
            const PairAudit& audit = audits[index];
            require(audit.bodies == kBodyCount && audit.failures.empty() &&
                        audit.equalities == 0 &&
                        audit.active_carrier == kExpectedActiveCounts[index] &&
                        audit.active_carrier_fnv == kExpectedActiveFnv[index],
                    "detached exchange boundary row changed");
            total_bodies += audit.bodies;
            total_checks += audit.body_checks;
            total_failures += audit.failures.size();
            total_equalities += audit.equalities;
            boundary_ledger.add(audit.pair.q);
            boundary_ledger.add(audit.pair.r);
            boundary_ledger.add(audit.active_carrier);
            boundary_ledger.add(audit.active_carrier_fnv);
            boundary_ledger.add(audit.body_checks);
            boundary_ledger.add(audit.failures.size());
            boundary_ledger.add(body_fnv(audit.failures));
            add_i128(boundary_ledger, audit.minimum_active_margin);
            add_i128(boundary_ledger, audit.maximum_inactive_margin);
            std::cout << "PAIR " << audit.pair.q << ',' << audit.pair.r
                      << " ACTIVE " << audit.active_carrier << " ACTIVE_FNV "
                      << std::hex << audit.active_carrier_fnv << std::dec
                      << " ACTIVE_MARGIN_MIN "
                      << decimal(audit.minimum_active_margin)
                      << " INACTIVE_MARGIN_MAX "
                      << decimal(audit.maximum_inactive_margin)
                      << " BODY_CHECKS " << audit.body_checks
                      << " FAILURES " << audit.failures.size() << '\n';
        }
        require(total_bodies == 9 * kBodyCount && total_failures == 0 &&
                    total_equalities == 0,
                "detached aggregate exchange closure changed");
        std::cout << "BOUNDARY ROWS 9 BODIES " << total_bodies
                  << " BODY_CHECKS " << total_checks << " FAILURES "
                  << total_failures << " EQUALITIES " << total_equalities
                  << " FNV " << std::hex << boundary_ledger.state << std::dec
                  << '\n'
                  << "SCOPE DETACHED_LITERAL_NINE_ROW_EXCHANGE_NO_PHYSICAL_"
                     "ENTRY_NO_LRC14\n"
                  << "VERDICT PASS STRICT_DELETIONS_RESPONSES_AND_EXHAUSTIVE_"
                     "BOUNDARY_CLOSURE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DETACHED_EXCHANGE_ERROR " << error.what() << '\n';
        return 1;
    }
}
