// Exact sparse one-apex inverse for the full upper perfect-number neighbor.
#include <cassert>
#include <cstdio>
#include <cstdint>

using Integer = __int128_t;

constexpr Integer power(int exponent) { return Integer(1) << exponent; }

// All previous marks occur strictly above the proposed new smallest mark.
Integer increment(Integer address, int previous, int position) {
    return power(position) - (address >> (position + 1)) - (position + 1 == previous);
}

bool find(int remaining, int previous, Integer address, Integer height, Integer target) {
    if (remaining == 1) {
        if (previous < 2 || height + increment(address, previous, 1) > target ||
            height + increment(address, previous, previous - 1) < target) {
            return false;
        }
        int low = 1;
        int high = previous - 1;
        while (low < high) {
            int middle = (low + high) / 2;
            if (height + increment(address, previous, middle) < target) low = middle + 1;
            else high = middle;
        }
        return height + increment(address, previous, low) == target;
    }
    for (int position = remaining; position < previous; ++position) {
        Integer next_height = height + increment(address, previous, position);
        Integer next_address = address + power(position);
        // All subsequent marks contribute less than 2^position in positive terms.
        if (next_height + power(position) - 2 < target) continue;
        // Their negative geometric weights have total multiplier at most 1/2.
        if (next_height - (next_address + power(position)) / 2 - remaining > target) continue;
        if (find(remaining - 1, position, next_address, next_height, target)) return true;
    }
    return false;
}

bool find_top_mark(int exponent, int mark_count) {
    const int top = 2 * exponent - 1;
    const Integer target = power(top) - power(exponent - 1) + 1;
    if (mark_count == 1) return power(top) + 1 == target;
    return find(mark_count - 1, top, power(top), power(top) + 1, target);
}

int main() {
    // A positive even-exponent control has marks 15,11,3.
    assert(find_top_mark(8, 3));
    for (int marks = 1; marks <= 8; ++marks) {
        if (find_top_mark(61, marks)) return 1;
        std::printf("exponent 61, %d marks: no upper hit\n", marks);
        std::fflush(stdout);
    }
}
