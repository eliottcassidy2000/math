// Direct proof replay for the 422-mask signature-surgery deck obtained from
// the THM-4281 421 deck by deleting its seven (256,663)-inactive masks and
// appending an explicit eight-mask cover witness.  This replay makes no
// cardinality-minimality claim.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace {

constexpr std::array<std::size_t, 7> DELETED_663 =
    {75, 107, 139, 374, 394, 405, 417};
constexpr std::array<u32, 7> DELETED_MASKS_663 = {
    UINT32_C(0x21948006), UINT32_C(0x128c8900),
    UINT32_C(0x10259240), UINT32_C(0x20027016),
    UINT32_C(0x10a41016), UINT32_C(0x2051d200),
    UINT32_C(0x01188550)};
constexpr std::array<u32, 8> REPLACEMENTS_663 = {
    UINT32_C(0x0110a550), UINT32_C(0x04871108),
    UINT32_C(0x10241207), UINT32_C(0x1042d088),
    UINT32_C(0x12848902), UINT32_C(0x21141284),
    UINT32_C(0x30249140), UINT32_C(0x31202206)};
constexpr u64 EXPECTED_JOINT_FNV = UINT64_C(0x20d63dd42fe8150e);

u64 mask_fnv663(const std::vector<u32>& deck) {
    FnvLocal ledger;
    for (u32 mask : deck) ledger.add(mask);
    return ledger.state;
}

i128 atom_mass663(const AtomData& atoms, u32 repair) {
    i128 mass = 0;
    for (const auto& [failure, value] : atoms.mass)
        if ((failure & ~repair) == 0) mass += value;
    return mass;
}

std::vector<u32> read_joint(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open joint deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "bad joint mask token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "joint mask arity/distinctness changed");
        deck.push_back(mask);
    }
    require(deck.size() == 421 && mask_fnv663(deck) == EXPECTED_JOINT_FNV,
            "joint deck identity changed");
    for (std::size_t slot = 0; slot < DELETED_663.size(); ++slot)
        require(deck[DELETED_663[slot]] == DELETED_MASKS_663[slot],
                "deleted mask identity changed");
    return deck;
}

std::vector<u32> read_obligations663(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open obligation CSV");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "body,response",
            "obligation header changed");
    std::vector<u32> bodies;
    FnvLocal body_ledger;
    FnvLocal obligation_ledger;
    while (std::getline(input, line)) {
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos, "malformed obligation row");
        const u32 body = static_cast<u32>(
            std::stoul(line.substr(0, comma), nullptr, 16));
        const unsigned response = static_cast<unsigned>(
            std::stoul(line.substr(comma + 1), nullptr, 16));
        require(std::popcount(body) == 9 && response > 0 && response < 128 &&
                    (bodies.empty() || bodies.back() < body),
                "invalid obligation row");
        bodies.push_back(body);
        body_ledger.add(body);
        obligation_ledger.add(body);
        obligation_ledger.add(response);
    }
    require(bodies.size() == 71 &&
                body_ledger.state == UINT64_C(0x414d30a143d2e22d) &&
                obligation_ledger.state == UINT64_C(0x287281de2900cca7),
            "obligation universe changed");
    return bodies;
}

struct Signature663 {
    int q = 0;
    int r = 0;
    u64 inactive_count = 0;
    std::array<u64, 7> words{};
};

std::vector<Signature663> read_signatures663(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature CSV");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::vector<Signature663> rows;
    FnvLocal signature_ledger;
    FnvLocal row_ledger;
    while (std::getline(input, line)) {
        std::vector<std::string> fields;
        std::stringstream stream(line);
        std::string field;
        while (std::getline(stream, field, ',')) fields.push_back(field);
        require(fields.size() == 10, "malformed signature row");
        Signature663 row;
        row.q = std::stoi(fields[0]);
        row.r = std::stoi(fields[1]);
        row.inactive_count = std::stoull(fields[2]);
        for (std::size_t word = 0; word < row.words.size(); ++word) {
            row.words[word] = std::stoull(fields[word + 3], nullptr, 16);
            signature_ledger.add(row.words[word]);
        }
        require(row.q > 0 && row.q < row.r &&
                    (rows.empty() ||
                     std::pair{rows.back().q, rows.back().r} <
                         std::pair{row.q, row.r}),
                "signature edges are not ordered/distinct");
        u64 count = 0;
        for (u64 word : row.words) count += std::popcount(word);
        require(count == row.inactive_count &&
                    (row.words.back() >> (421 % 64)) == 0,
                "signature weight/padding changed");
        row_ledger.add(row.q);
        row_ledger.add(row.r);
        row_ledger.add(row.inactive_count);
        for (u64 word : row.words) row_ledger.add(word);
        rows.push_back(row);
    }
    require(rows.size() == 24223 &&
                signature_ledger.state == UINT64_C(0x1f991a30499f0ae0) &&
                row_ledger.state == UINT64_C(0x3652b5590b330704),
            "full signature ledger changed");
    return rows;
}

std::array<u64, 7> deleted_signature_words() {
    std::array<u64, 7> words{};
    for (std::size_t index : DELETED_663)
        words[index / 64] |= UINT64_C(1) << (index % 64);
    return words;
}

struct KeyControl {
    int q = 0;
    int r = 0;
    u64 active = 0;
    u32 first_inactive = 0;
    i128 minimum_margin = 0;
    bool seen = false;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 6,
                "usage: surgery-deck JOINT_DECK OBLIGATIONS SIGNATURES "
                "SURGERY_DECK_OUT COMMON_OUT");
        const std::vector<u32> joint = read_joint(argv[1]);
        const std::vector<u32> obligations = read_obligations663(argv[2]);
        const std::vector<Signature663> signatures =
            read_signatures663(argv[3]);
        std::array<bool, 421> deleted{};
        for (std::size_t index : DELETED_663) deleted[index] = true;
        std::vector<u32> surgery;
        surgery.reserve(422);
        for (std::size_t index = 0; index < joint.size(); ++index)
            if (!deleted[index]) surgery.push_back(joint[index]);
        surgery.insert(surgery.end(), REPLACEMENTS_663.begin(),
                       REPLACEMENTS_663.end());
        require(surgery.size() == 422 &&
                    std::set<u32>(surgery.begin(), surgery.end()).size() ==
                        surgery.size(),
                "surgery deck arity/distinctness changed");
        const u64 surgery_fnv = mask_fnv663(surgery);
        std::ofstream deck_out(argv[4]);
        require(static_cast<bool>(deck_out), "cannot create surgery deck");
        for (u32 mask : surgery)
            deck_out << std::hex << std::setw(8) << std::setfill('0') << mask
                     << std::dec << std::setfill(' ') << '\n';
        require(deck_out.good(), "failed writing surgery deck");

        u64 obligation_incidences = 0;
        u64 minimum_obligation_response = UINT64_MAX;
        u64 maximum_obligation_response = 0;
        FnvLocal obligation_response_ledger;
        for (u32 body : obligations) {
            u64 response = 0;
            for (u32 mask : REPLACEMENTS_663) response += (body & mask) == 0;
            require(response > 0, "replacement cover misses an obligation");
            obligation_incidences += response;
            minimum_obligation_response =
                std::min(minimum_obligation_response, response);
            maximum_obligation_response =
                std::max(maximum_obligation_response, response);
            obligation_response_ledger.add(body);
            obligation_response_ledger.add(response);
        }

        u64 bodies = 0;
        u64 body_checks = 0;
        u64 body_failures = 0;
        FnvLocal failure_ledger;
        u32 body = (u32{1} << 9) - 1;
        while (body < (u32{1} << 30)) {
            ++bodies;
            bool covered = false;
            for (u32 mask : surgery) {
                ++body_checks;
                if ((body & mask) == 0) {
                    covered = true;
                    break;
                }
            }
            if (!covered) {
                ++body_failures;
                failure_ledger.add(body);
            }
            const u32 next = next_combination(body);
            if (next <= body) break;
            body = next;
        }
        require(bodies == UINT64_C(14307150) && body_failures == 0,
                "surgery deck body coverage failed");

        const std::array<u64, 7> deleted_words = deleted_signature_words();
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::ofstream common_out(argv[5]);
        require(static_cast<bool>(common_out), "cannot create common CSV");
        u64 candidates = 0;
        u64 common = 0;
        u64 activation_checks = 0;
        u64 equalities = 0;
        FnvLocal candidate_ledger;
        FnvLocal common_ledger;
        std::array<KeyControl, 5> controls = {
            KeyControl{256, 663}, KeyControl{256, 662},
            KeyControl{256, 287}, KeyControl{256, 558},
            KeyControl{256, 575}};
        for (const Signature663& row : signatures) {
            bool subset = true;
            for (std::size_t word = 0; word < row.words.size(); ++word)
                subset &= (row.words[word] & ~deleted_words[word]) == 0;
            if (!subset) continue;
            ++candidates;
            candidate_ledger.add(row.q);
            candidate_ledger.add(row.r);
            const i64 g = std::gcd(row.q, row.r);
            const PrimitivePair primitive = build_primitive(row.q / g,
                                                             row.r / g);
            const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
            const i128 denominator =
                static_cast<i128>(primitive.grid) * g * COMMON;
            u64 active = 0;
            u32 first_inactive = 0;
            i128 minimum_margin = 0;
            bool first = true;
            for (u32 mask : REPLACEMENTS_663) {
                ++activation_checks;
                const i128 margin = static_cast<i128>(63) *
                                        atom_mass663(atoms, mask) -
                                    static_cast<i128>(4) * denominator;
                equalities += margin == 0;
                if (first || margin < minimum_margin) {
                    minimum_margin = margin;
                    first = false;
                }
                if (margin < 0) {
                    if (first_inactive == 0) first_inactive = mask;
                } else {
                    ++active;
                }
            }
            if (first_inactive == 0) {
                ++common;
                common_ledger.add(row.q);
                common_ledger.add(row.r);
                common_out << row.q << ',' << row.r << '\n';
            }
            for (KeyControl& control : controls) {
                if (control.q == row.q && control.r == row.r) {
                    control.active = active;
                    control.first_inactive = first_inactive;
                    control.minimum_margin = minimum_margin;
                    control.seen = true;
                }
            }
        }
        require(common_out.good(), "failed writing common CSV");
        require(candidates == 395, "signature-subset census changed");
        for (const KeyControl& control : controls)
            require(control.seen, "key signature-subset control absent");

        std::cout << "THM4281_SIGNATURE_SURGERY_256_663_DECK_REPLAY_V1\n"
                  << "ORIGINAL 421 FNV " << std::hex << EXPECTED_JOINT_FNV
                  << std::dec << " DELETE 7 APPEND 8 SURGERY "
                  << surgery.size() << " FNV " << std::hex << surgery_fnv
                  << std::dec << '\n'
                  << "OBLIGATIONS " << obligations.size()
                  << " REPLACEMENT_INCIDENCES " << obligation_incidences
                  << " MIN_RESPONSE " << minimum_obligation_response
                  << " MAX_RESPONSE " << maximum_obligation_response
                  << " RESPONSE_FNV " << std::hex
                  << obligation_response_ledger.state << std::dec << '\n'
                  << "BODY_SCAN " << bodies << " CHECKS " << body_checks
                  << " FAILURES " << body_failures << " FAILURE_FNV "
                  << std::hex << failure_ledger.state << std::dec << '\n'
                  << "SIGNATURE_SUBSET_CANDIDATES " << candidates
                  << " CANDIDATE_FNV " << std::hex << candidate_ledger.state
                  << std::dec << " REPLACEMENT_TESTS " << activation_checks
                  << " EQUALITIES " << equalities << '\n'
                  << "SURGERY_COMMON " << common << " COMMON_FNV " << std::hex
                  << common_ledger.state << std::dec << '\n';
        for (const KeyControl& control : controls)
            std::cout << "CONTROL " << control.q << ',' << control.r
                      << " ACTIVE_REPLACEMENTS " << control.active
                      << " FIRST_INACTIVE " << std::hex
                      << control.first_inactive << std::dec << " MIN_MARGIN "
                      << decimal(control.minimum_margin) << '\n';
        std::cout << "VERDICT PASS DIRECT_BODY_AND_COMMON_ACTIVITY_REPLAY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SURGERY_DECK_ERROR " << error.what() << '\n';
        return 1;
    }
}
