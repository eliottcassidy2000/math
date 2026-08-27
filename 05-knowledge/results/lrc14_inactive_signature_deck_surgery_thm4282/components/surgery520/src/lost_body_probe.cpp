#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr std::size_t kDeckSize = 421;
constexpr std::size_t kBodyCount = 14307150;
constexpr u64 kDeckFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr std::array<std::size_t, 5> kDeleted{57, 107, 222, 275, 345};

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (value >> (8 * byte)) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};

std::vector<u32> read_deck(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "bad mask token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "bad or duplicate mask");
        deck.push_back(mask);
    }
    Fnv ledger;
    for (u32 mask : deck) ledger.add(mask);
    require(deck.size() == kDeckSize && ledger.state == kDeckFnv,
            "deck identity changed");
    return deck;
}

u32 next_same_popcount(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: lost-body-probe DECK OBLIGATIONS_CSV");
        const std::vector<u32> deck = read_deck(argv[1]);
        std::array<bool, kDeckSize> deleted{};
        for (std::size_t index : kDeleted) deleted[index] = true;

        std::ofstream output(argv[2]);
        require(static_cast<bool>(output), "cannot create obligations CSV");
        output << "ordinal,body_hex,deleted_response\n";
        std::array<u64, 32> response_hist{};
        Fnv obligation_ledger;
        u64 retained_checks = 0;
        u64 obligations = 0;
        u32 body = (UINT32_C(1) << 9) - 1;
        for (std::size_t ordinal = 0; ordinal < kBodyCount; ++ordinal) {
            bool retained_hit = false;
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (deleted[index]) continue;
                ++retained_checks;
                if ((body & deck[index]) == 0) {
                    retained_hit = true;
                    break;
                }
            }
            if (!retained_hit) {
                unsigned response = 0;
                for (unsigned local = 0; local < kDeleted.size(); ++local) {
                    if ((body & deck[kDeleted[local]]) == 0)
                        response |= 1u << local;
                }
                require(response != 0, "original deck failed body coverage");
                ++response_hist[response];
                obligation_ledger.add(body);
                obligation_ledger.add(response);
                output << obligations++ << ',' << std::hex << std::setw(8)
                       << std::setfill('0') << body << std::dec
                       << std::setfill(' ') << ',' << response << '\n';
            }
            if (ordinal + 1 < kBodyCount) body = next_same_popcount(body);
        }
        require(body == UINT32_C(0x3fe00000), "body endpoint changed");
        require(output.good(), "failed writing obligations CSV");
        std::cout << "THM4281_SIGNATURE_SURGERY_LOST_BODY_PROBE_V1\n"
                  << "TARGET 520,663 DELETED_INDICES 57,107,222,275,345\n"
                  << "BODIES " << kBodyCount << " RETAINED_CHECKS "
                  << retained_checks << " OBLIGATIONS " << obligations
                  << " OBLIGATION_FNV " << std::hex << obligation_ledger.state
                  << std::dec << '\n';
        for (unsigned response = 1; response < response_hist.size(); ++response)
            if (response_hist[response])
                std::cout << "DELETED_RESPONSE " << std::hex << response
                          << std::dec << " COUNT " << response_hist[response]
                          << '\n';
        std::cout << "VERDICT PASS EXACT_ALL_NINE_BODIES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LOST_BODY_PROBE_ERROR " << error.what() << '\n';
        return 1;
    }
}
