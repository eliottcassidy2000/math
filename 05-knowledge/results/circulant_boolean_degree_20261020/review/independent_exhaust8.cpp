// Independent enumeration: ordinary binary coefficient order, split into two
// lookup tables. No Gray-code traversal and no producer source imports.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

using Wide = __uint128_t;

int main(int argc, char **argv) {
    if (argc != 2) return 1;
    std::ifstream input(argv[1]);
    unsigned dimension, source_types;
    input >> dimension >> source_types;
    if (!input || dimension != 29 || source_types != 16) return 2;
    std::vector<Wide> basis(dimension);
    for (unsigned k = 0; k < dimension; ++k) {
        uint64_t low, high;
        input >> low >> high;
        if (!input) return 3;
        basis[k] = Wide(low) | (Wide(high) << 64);
    }
    std::array<int, 256> target_class{};
    for (int &entry : target_class) {
        input >> entry;
        if (!input || entry < -1 || entry > 15) return 4;
    }
    constexpr unsigned split = 15;
    std::vector<Wide> bottom(1 << split), top(1 << (29 - split));
    for (unsigned c = 0; c < bottom.size(); ++c)
        for (unsigned bit = 0; bit < split; ++bit)
            if ((c >> bit) & 1) bottom[c] ^= basis[bit];
    for (unsigned c = 0; c < top.size(); ++c)
        for (unsigned bit = 0; bit < 29 - split; ++bit)
            if ((c >> bit) & 1) top[c] ^= basis[split + bit];
    std::array<uint64_t, 17> histogram{};
    uint64_t visited = 0;
    uint64_t first_bijection = 0;
    bool found_bijection = false;
    for (unsigned high = 0; high < top.size(); ++high) {
        for (unsigned low = 0; low < bottom.size(); ++low) {
            Wide packed = top[high] ^ bottom[low];
            unsigned image_classes = 0;
            bool allowed = true;
            ++visited;
            for (unsigned source = 0; source < source_types; ++source) {
                const int target = target_class[unsigned(packed & 255)];
                if (target < 0) { allowed = false; break; }
                image_classes |= 1u << target;
                packed >>= 8;
            }
            if (allowed) {
                const unsigned coverage = __builtin_popcount(image_classes);
                ++histogram[coverage];
                if (coverage == 16 && !found_bijection) {
                    found_bijection = true;
                    first_bijection = (uint64_t(high) << split) | low;
                }
            }
        }
    }
    std::cout << "{\"maps\":" << visited
              << ",\"first_coefficient_mask\":";
    if (found_bijection) std::cout << first_bijection;
    else std::cout << "null";
    std::cout << ",\"target_valued_type_coverage\":{";
    bool comma = false;
    for (unsigned count = 0; count <= 16; ++count) {
        if (!histogram[count]) continue;
        if (comma) std::cout << ',';
        comma = true;
        std::cout << '"' << count << "\":" << histogram[count];
    }
    std::cout << "}}\n";
    return histogram[16] ? 5 : 0;
}
