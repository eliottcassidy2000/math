#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Pair = std::pair<int, int>;
using PairSet = std::set<Pair>;

namespace {

constexpr std::uint64_t FNV_OFFSET = 0xcbf29ce484222325ULL;
constexpr std::uint64_t FNV_PRIME = 0x100000001b3ULL;

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string &message) {
    if (!condition) fail(message);
}

std::vector<std::string> split_csv(const std::string &line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const auto comma = line.find(',', begin);
        fields.push_back(line.substr(begin, comma == std::string::npos
                                                ? comma
                                                : comma - begin));
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    return fields;
}

std::vector<Pair> read_pairs(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    require(input.good(), "cannot open pair ledger: " + path);
    std::vector<Pair> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line == "q,r") continue;
        const auto fields = split_csv(line);
        require(fields.size() == 2, "bad pair row in " + path + ": " + line);
        rows.emplace_back(std::stoi(fields[0]), std::stoi(fields[1]));
    }
    require(std::is_sorted(rows.begin(), rows.end()), "unsorted pair ledger: " + path);
    require(std::adjacent_find(rows.begin(), rows.end()) == rows.end(),
            "duplicate pair ledger: " + path);
    return rows;
}

PairSet as_set(const std::vector<Pair> &rows) {
    return PairSet(rows.begin(), rows.end());
}

std::vector<Pair> as_vector(const PairSet &rows) {
    return std::vector<Pair>(rows.begin(), rows.end());
}

PairSet intersection(const PairSet &left, const PairSet &right) {
    PairSet out;
    std::set_intersection(left.begin(), left.end(), right.begin(), right.end(),
                          std::inserter(out, out.end()));
    return out;
}

PairSet difference(const PairSet &left, const PairSet &right) {
    PairSet out;
    std::set_difference(left.begin(), left.end(), right.begin(), right.end(),
                        std::inserter(out, out.end()));
    return out;
}

PairSet unite(const PairSet &left, const PairSet &right) {
    PairSet out;
    std::set_union(left.begin(), left.end(), right.begin(), right.end(),
                   std::inserter(out, out.end()));
    return out;
}

bool subset(const PairSet &small, const PairSet &large) {
    return std::includes(large.begin(), large.end(), small.begin(), small.end());
}

std::uint64_t fnv_pairs(const std::vector<Pair> &rows) {
    std::uint64_t state = FNV_OFFSET;
    for (const auto &[q, r] : rows) {
        for (const std::uint64_t value : {static_cast<std::uint64_t>(q),
                                          static_cast<std::uint64_t>(r)}) {
            for (unsigned shift = 0; shift < 64; shift += 8) {
                state ^= (value >> shift) & 0xffULL;
                state *= FNV_PRIME;
            }
        }
    }
    return state;
}

std::string fnv_hex(const PairSet &rows) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16)
        << fnv_pairs(as_vector(rows));
    return out.str();
}

void write_pairs(const std::string &path, const PairSet &rows) {
    std::ofstream output(path, std::ios::binary);
    require(output.good(), "cannot create output: " + path);
    for (const auto &[q, r] : rows) output << q << ',' << r << '\n';
}

bool files_equal(const std::string &left_path, const std::string &right_path) {
    std::ifstream left(left_path, std::ios::binary);
    std::ifstream right(right_path, std::ios::binary);
    require(left.good() && right.good(), "cannot open files for byte comparison");
    std::istreambuf_iterator<char> left_it(left), right_it(right), end;
    while (left_it != end && right_it != end) {
        if (*left_it != *right_it) return false;
        ++left_it;
        ++right_it;
    }
    return left_it == end && right_it == end;
}

std::string show(const PairSet &rows) {
    std::ostringstream out;
    bool first = true;
    for (const auto &[q, r] : rows) {
        if (!first) out << ' ';
        first = false;
        out << '(' << q << ',' << r << ')';
    }
    return first ? "EMPTY" : out.str();
}

PairSet rows_at_max_endpoint(const PairSet &rows) {
    require(!rows.empty(), "empty row set has no maximum endpoint");
    int max_r = 0;
    for (const auto &[q, r] : rows) max_r = std::max(max_r, r);
    PairSet out;
    for (const auto &pair : rows) {
        if (pair.second == max_r) out.insert(pair);
    }
    return out;
}

PairSet derive_j19(const std::string &path, const PairSet &current) {
    std::ifstream input(path, std::ios::binary);
    require(input.good(), "cannot open signatures: " + path);
    PairSet result;
    std::string line;
    require(static_cast<bool>(std::getline(input, line)), "empty signatures ledger");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    require(line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "unexpected signatures header");
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_csv(line);
        require(fields.size() == 10, "bad signatures row: " + line);
        const Pair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        if (!current.contains(pair) || std::stoi(fields[2]) != 1) continue;
        std::array<std::uint64_t, 7> words{};
        for (int i = 0; i < 7; ++i) words[i] = std::stoull(fields[3 + i], nullptr, 16);
        bool exactly_19 = words[0] == (1ULL << 19);
        for (int i = 1; i < 7; ++i) exactly_19 = exactly_19 && words[i] == 0;
        if (exactly_19) result.insert(pair);
    }
    return result;
}

struct ScanAudit {
    PairSet successes;
    PairSet failures;
    PairSet all;
    std::vector<Pair> order;
    std::uint64_t body_failures = 0;
    std::uint64_t completed_rows = 0;
    std::uint64_t declared_failures = 0;
    bool stopped = false;
    bool saw_summary = false;
    bool saw_pass = false;
};

ScanAudit parse_descent_scan(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    require(input.good(), "cannot open descent scan: " + path);
    ScanAudit audit;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.starts_with("PAIR ")) {
            std::istringstream fields(line);
            std::vector<std::string> tokens;
            for (std::string token; fields >> token;) tokens.push_back(token);
            require(tokens.size() >= 4 && tokens[0] == "PAIR",
                    "bad PAIR line in descent scan: " + line);
            const auto comma = tokens[1].find(',');
            require(comma != std::string::npos,
                    "bad pair token in descent scan: " + tokens[1]);
            const Pair pair{std::stoi(tokens[1].substr(0, comma)),
                            std::stoi(tokens[1].substr(comma + 1))};
            const auto failure_it = std::find(tokens.begin(), tokens.end(), "FAILURES");
            require(failure_it != tokens.end() && failure_it + 1 != tokens.end(),
                    "PAIR line lacks failure count: " + line);
            const std::uint64_t failure_count = std::stoull(*(failure_it + 1));
            require(audit.all.insert(pair).second, "duplicate pair in descent scan");
            audit.order.push_back(pair);
            audit.body_failures += failure_count;
            (failure_count == 0 ? audit.successes : audit.failures).insert(pair);
            continue;
        }
        if (line.starts_with("COMPLETED_ROWS ")) {
            require(!audit.saw_summary, "duplicate descent completion summary");
            std::istringstream fields(line);
            std::string completed_label, failure_label, stopped_label;
            int stopped = -1;
            fields >> completed_label >> audit.completed_rows >> failure_label
                   >> audit.declared_failures >> stopped_label >> stopped;
            require(fields && completed_label == "COMPLETED_ROWS" &&
                        failure_label == "TOTAL_FAILURES" && stopped_label == "STOPPED" &&
                        (stopped == 0 || stopped == 1),
                    "bad descent completion summary: " + line);
            audit.stopped = stopped != 0;
            audit.saw_summary = true;
            continue;
        }
        if (line.starts_with("VERDICT PASS ")) audit.saw_pass = true;
    }
    require(audit.saw_summary && audit.saw_pass,
            "descent scan is partial or lacks PASS verdict");
    require(audit.completed_rows == audit.all.size(),
            "COMPLETED_ROWS differs from parsed PAIR rows");
    require(audit.declared_failures == audit.body_failures,
            "TOTAL_FAILURES differs from parsed body failures");
    require(std::is_sorted(audit.order.begin(), audit.order.end(),
                           [](Pair left, Pair right) {
                               return left.second != right.second
                                          ? left.second > right.second
                                          : left.first < right.first;
                           }),
            "descent PAIR rows are not in endpoint/q order");
    return audit;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        require(argc == 12,
                "usage: audit CURRENT SINGLETON SIGNATURES ENDPOINT J19 "
                "EXPECTED_UNION EXPECTED_FINAL DESCENT_SCAN OUT_ENDPOINT "
                "OUT_UNION OUT_FINAL");

        const auto current_vec = read_pairs(argv[1]);
        const auto singleton_vec = read_pairs(argv[2]);
        const PairSet current = as_set(current_vec);
        const PairSet singleton = as_set(singleton_vec);
        const PairSet endpoint_ledger = as_set(read_pairs(argv[4]));
        const PairSet j19_ledger = as_set(read_pairs(argv[5]));
        const PairSet expected_union = as_set(read_pairs(argv[6]));
        const PairSet expected_final = as_set(read_pairs(argv[7]));

        require(current.size() == 22647, "current cardinality mismatch");
        require(fnv_pairs(current_vec) == 0xdf5374d4aca67677ULL,
                "current FNV mismatch");
        require(singleton.size() == 1219, "singleton cardinality mismatch");
        require(fnv_pairs(singleton_vec) == 0x3706723fd3574334ULL,
                "singleton FNV mismatch");
        require(subset(singleton, current), "singleton union is not within current residual");

        const PairSet j19_derived = derive_j19(argv[3], current);
        require(j19_derived == j19_ledger,
                "independently derived J19 ideal differs from frozen ledger");
        require(j19_ledger.size() == 36, "J19 cardinality mismatch");
        require(fnv_hex(j19_ledger) == "5c8af37cf2f002e7", "J19 FNV mismatch");
        require(subset(j19_ledger, current), "J19 is not within current residual");

        require(endpoint_ledger.size() == 70, "endpoint cardinality mismatch");
        require(fnv_hex(endpoint_ledger) == "213e8d6fcbd0f1cd",
                "endpoint FNV mismatch");
        require(subset(endpoint_ledger, current), "endpoint ledger is not within current residual");
        const ScanAudit scan = parse_descent_scan(argv[8]);
        require(scan.successes == endpoint_ledger,
                "raw descent scan success rows differ from frozen endpoint ledger");
        require(scan.failures == PairSet{{100, 626}, {256, 626}},
                "raw descent scan failure rows differ from expected two-row frontier");
        require(scan.completed_rows == 72 && scan.body_failures == 1005 && scan.stopped,
                "raw descent completion totals changed");
        PairSet current_prefix;
        for (const auto &pair : current)
            if (pair.second >= 626) current_prefix.insert(pair);
        require(scan.all == current_prefix,
                "raw descent rows are not the complete current r>=626 prefix");

        const PairSet s_j = intersection(singleton, j19_ledger);
        const PairSet s_e = intersection(singleton, endpoint_ledger);
        const PairSet j_e = intersection(j19_ledger, endpoint_ledger);
        const PairSet all_three = intersection(s_j, endpoint_ledger);
        require(s_j.empty(), "singleton/J19 overlap is nonempty");
        require(s_e == PairSet{{410, 626}, {506, 626}},
                "singleton/endpoint overlap mismatch");
        require(j_e == PairSet{{338, 628}, {338, 636}},
                "J19/endpoint overlap mismatch");
        require(all_three.empty(), "unexpected triple overlap");

        const PairSet proof_union = unite(unite(singleton, j19_ledger), endpoint_ledger);
        const PairSet final_rows = difference(current, proof_union);
        require(proof_union.size() == 1321, "typed union cardinality mismatch");
        require(fnv_hex(proof_union) == "69bcc5c5fc86ac8e", "typed union FNV mismatch");
        require(final_rows.size() == 21326, "final residual cardinality mismatch");
        require(fnv_hex(final_rows) == "a05a372bf57ea064", "final residual FNV mismatch");
        require(proof_union == expected_union, "typed union differs from root frozen ledger");
        require(final_rows == expected_final, "final residual differs from root frozen ledger");
        require(intersection(proof_union, final_rows).empty(), "partition overlap");
        require(unite(proof_union, final_rows) == current, "partition does not recover current");

        const PairSet final_top = rows_at_max_endpoint(final_rows);
        require(final_top == PairSet{{100, 626}, {256, 626}},
                "final top-row set mismatch");

        write_pairs(argv[9], endpoint_ledger);
        write_pairs(argv[10], proof_union);
        write_pairs(argv[11], final_rows);
        require(files_equal(argv[4], argv[9]),
                "derived endpoint output is not byte-identical to packet ledger");
        require(files_equal(argv[6], argv[10]),
                "derived union output is not byte-identical to packet ledger");
        require(files_equal(argv[7], argv[11]),
                "derived residual output is not byte-identical to packet ledger");

        std::cout << "INDEPENDENT_CPP_TYPED_PROOF_GRAPH_AUDIT_V2\n";
        std::cout << "CURRENT rows=" << current.size() << " fnv=" << fnv_hex(current)
                  << " max=" << rows_at_max_endpoint(current).begin()->second
                  << " top=" << show(rows_at_max_endpoint(current)) << "\n";
        std::cout << "SINGLETON type=110_separate_421_mask_common_decks rows="
                  << singleton.size() << " fnv=" << fnv_hex(singleton) << "\n";
        std::cout << "J19 type=one_422_mask_two_replacement_common_deck rows="
                  << j19_ledger.size() << " fnv=" << fnv_hex(j19_ledger) << "\n";
        std::cout << "ENDPOINT type=one_9019_mask_mixed_rank_exchange_carrier rows="
                  << endpoint_ledger.size() << " fnv=" << fnv_hex(endpoint_ledger) << "\n";
        std::cout << "INTERSECTION singleton_j19=" << s_j.size()
                  << " singleton_endpoint=" << s_e.size() << " rows=" << show(s_e)
                  << " j19_endpoint=" << j_e.size() << " rows=" << show(j_e)
                  << " triple=" << all_three.size() << "\n";
        std::cout << "MEMBERSHIP_PROFILE singleton_only=1217 j19_only=34 "
                     "endpoint_only=66 singleton_and_endpoint=2 "
                     "j19_and_endpoint=2\n";
        std::cout << "TYPED_UNION rows=" << proof_union.size()
                  << " fnv=" << fnv_hex(proof_union)
                  << " max=" << rows_at_max_endpoint(proof_union).begin()->second
                  << " top=" << show(rows_at_max_endpoint(proof_union)) << "\n";
        std::cout << "FINAL rows=" << final_rows.size()
                  << " fnv=" << fnv_hex(final_rows)
                  << " max=" << final_top.begin()->second
                  << " top=" << show(final_top) << "\n";
        std::cout << "LEDGER_EQUALITY endpoint_bytes=PASS union_bytes=PASS "
                     "final_bytes=PASS j19_derivation=PASS raw_descent_scan=PASS "
                     "complete_endpoint_prefix=PASS partition=PASS\n";
        std::cout << "TYPE_GUARD no_single_common_deck_claim_for_1321 "
                     "j19_not_421 singleton_margins_use_literal_corrected_ledger\n";
        std::cout << "VERDICT PASS FINITE_EXACT_TYPED_PROOF_GRAPH_UNION LRC14_OPEN\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "VERDICT FAIL " << error.what() << '\n';
        return 1;
    }
}
