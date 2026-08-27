#include <algorithm>
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
constexpr std::size_t TARGET_INDEX = 367;
constexpr u32 TARGET_MASK = UINT32_C(0x02188125);
constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0x20d63dd42fe8150e);
struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (value >> (8 * byte)) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
u32 next_same_popcount(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}
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
                "bad deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "bad deck mask");
        deck.push_back(mask);
    }
    Fnv ledger;
    for (u32 mask : deck) ledger.add(mask);
    require(deck.size() == 421 && ledger.state == EXPECTED_DECK_FNV &&
                deck[TARGET_INDEX] == TARGET_MASK,
            "deck identity changed");
    return deck;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: private-bodies JOINT_DECK");
        const std::vector<u32> deck = read_deck(argv[1]);
        std::vector<u32> private_bodies;
        u64 bodies = 0;
        u64 target_disjoint = 0;
        u64 other_checks = 0;
        u32 body = (UINT32_C(1) << 9) - 1;
        while (body < (UINT32_C(1) << 30)) {
            ++bodies;
            if ((body & TARGET_MASK) == 0) {
                ++target_disjoint;
                bool other = false;
                for (std::size_t index = 0; index < deck.size(); ++index) {
                    if (index == TARGET_INDEX) continue;
                    ++other_checks;
                    if ((body & deck[index]) == 0) {
                        other = true;
                        break;
                    }
                }
                if (!other) private_bodies.push_back(body);
            }
            const u32 next = next_same_popcount(body);
            if (next <= body) break;
            body = next;
        }
        require(bodies == UINT64_C(14307150) &&
                    target_disjoint == UINT64_C(497420),
                "body enumeration changed");
        Fnv ledger;
        for (u32 private_body : private_bodies) ledger.add(private_body);
        std::cout << "PRIVATE_BODIES_367_V1\n"
                  << "TARGET_INDEX " << TARGET_INDEX << " TARGET_MASK "
                  << std::hex << TARGET_MASK << std::dec << '\n'
                  << "BODIES " << bodies << " TARGET_DISJOINT "
                  << target_disjoint << " OTHER_CHECKS " << other_checks
                  << " PRIVATE " << private_bodies.size() << " FNV "
                  << std::hex << ledger.state << std::dec << '\n'
                  << "PRIVATE_MASKS_HEX ";
        for (std::size_t index = 0; index < private_bodies.size(); ++index) {
            if (index) std::cout << ',';
            std::cout << std::hex << std::setw(8) << std::setfill('0')
                      << private_bodies[index] << std::dec << std::setfill(' ');
        }
        std::cout << "\nVERDICT PASS EXACT_PRIVATE_BODY_FIBRE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PRIVATE_BODIES_ERROR " << error.what() << '\n';
        return 1;
    }
}
