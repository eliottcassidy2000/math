#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

struct Rat {
  std::int64_t n;
  std::int64_t d;
  Rat(std::int64_t numerator = 0, std::int64_t denominator = 1) {
    if (denominator == 0) throw std::runtime_error("zero denominator");
    if (denominator < 0) numerator = -numerator, denominator = -denominator;
    const std::int64_t g = std::gcd(numerator < 0 ? -numerator : numerator,
                                    denominator);
    n = numerator / g;
    d = denominator / g;
  }
};

static Rat operator-(Rat a, Rat b) { return Rat(a.n * b.d - b.n * a.d, a.d * b.d); }
static Rat operator*(std::int64_t a, Rat b) { return Rat(a * b.n, b.d); }
static bool operator==(Rat a, Rat b) { return a.n == b.n && a.d == b.d; }
static bool operator<(Rat a, Rat b) { return a.n * b.d < b.n * a.d; }
static bool operator>=(Rat a, Rat b) { return !(a < b); }

static std::ostream& operator<<(std::ostream& stream, Rat value) {
  stream << value.n;
  if (value.d != 1) stream << '/' << value.d;
  return stream;
}

static void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

static Rat gap(std::int64_t speed, Rat phase) {
  std::int64_t residue = (speed * phase.n) % phase.d;
  if (residue < 0) residue += phase.d;
  return Rat(std::min(residue, phase.d - residue), phase.d);
}

static Rat clearance(const std::vector<int>& speeds, Rat phase) {
  Rat answer(1, 2);
  for (int speed : speeds) answer = std::min(answer, gap(speed, phase));
  return answer;
}

static bool width_gate(int first, int last) {
  return std::int64_t{27} * (13 * first - last) >=
         std::int64_t{4} * first * last;
}

static std::uint64_t choose(int n, int k) {
  if (k > n - k) k = n - k;
  std::uint64_t answer = 1;
  for (int i = 1; i <= k; ++i) answer = answer * (n - k + i) / i;
  return answer;
}

static std::uint64_t fnv64(const std::string& text) {
  std::uint64_t value = 14695981039346656037ULL;
  for (unsigned char byte : text) {
    value ^= byte;
    value *= 1099511628211ULL;
  }
  return value;
}

int main() {
  try {
    const Rat delta(1, 14), gate(2, 189);
    require(3 * gate == Rat(2, 63), "scale-three gate");
    require(Rat(2, 7 * 27) == gate, "q=27 gate");

    Rat residual_min(1, 2);
    for (int m = 3; m <= 500; ++m) {
      for (int r = 1; r <= 25; r += 2) {
        residual_min = std::min(residual_min, Rat(1, 2) - Rat(r, 28 * m));
      }
    }
    require(residual_min == Rat(17, 84) && residual_min >= delta,
            "residual endpoint clock");

    std::vector<int> valid_46;
    for (int first = 3; first <= 1000; ++first) {
      if (width_gate(first, first + 45)) valid_46.push_back(first);
      require(!width_gate(first, first + 46), "47-label maximality");
    }
    const std::vector<int> expected_starts{14,15,16,17,18,19,20,21,22};
    require(valid_46 == expected_starts, "46-label starts");
    require(140 * 140 - 16 * 1242 == -272, "47-label discriminant");
    require(144 * 144 - 16 * 1215 == 1296, "46-label discriminant");

    std::vector<int> block;
    for (int h = 15; h <= 60; ++h) block.push_back(h);
    const Rat left(1, 210), right(13, 840), block_width = right - left;
    require(block_width == Rat(3, 280), "block width");
    require(block_width - gate == Rat(1, 7560), "block surplus");
    require(clearance(block, left) >= delta && clearance(block, right) >= delta,
            "closed first-window endpoints");

    std::vector<int> doubled;
    for (int h : block) doubled.push_back(2 * h);
    const Rat endpoint_clock(211, 420);
    Rat endpoint_tail_min(1, 2);
    for (int r = 1; r <= 25; r += 2) {
      endpoint_tail_min = std::min(endpoint_tail_min, gap(r, endpoint_clock));
    }
    require(endpoint_tail_min == Rat(37, 84), "A=15 endpoint tail minimum");
    require(clearance(doubled, endpoint_clock) >= delta, "doubled body endpoint");

    require(gap(1, Rat(1, 420)) == Rat(1, 420), "hostile lower lift");
    require(gap(211, Rat(211, 420)) == Rat(1, 420), "hostile upper lift");
    require((211LL * 211LL) % 420 == 1, "hostile congruence");

    const Rat rescue_body(1, 105), rescue(53, 105);
    require(2 * rescue - Rat(1) == rescue_body, "rescue lift");
    std::vector<int> superbody = doubled;
    superbody.push_back(1);
    superbody.push_back(211);
    const Rat rescue_clearance = clearance(superbody, rescue);
    require(rescue_clearance == Rat(1, 7), "rescue clearance");
    require(gap(1, rescue) == Rat(52, 105), "rescue first tail");
    require(gap(211, rescue) == Rat(52, 105), "rescue second tail");

    const std::uint64_t families = choose(46, 11);
    require(families == 13340783196ULL, "family count");

    std::uint64_t all_width_families = 0;
    int first_minimum = -1, last_minimum = -1, maximum_label = 0;
    for (int minimum = 3; minimum <= 1000; ++minimum) {
      const int last_cap = 351 * minimum / (4 * minimum + 27);
      if (last_cap - minimum < 10) continue;
      if (first_minimum < 0) first_minimum = minimum;
      last_minimum = minimum;
      maximum_label = std::max(maximum_label, last_cap);
      all_width_families += choose(last_cap - minimum, 10);
    }
    require(first_minimum == 3 && last_minimum == 70,
            "all-width minimum range");
    require(maximum_label == 80, "all-width maximum label");
    require(all_width_families == 60301609751ULL,
            "all-width eleven-body count");
    const std::string ledger =
        "gate=2/189;scale3=2/63;q27=2/189;"
        "residual=17/84;valid46=14..22;max47disc=-272;"
        "block15_60=3/280;surplus=1/7560;"
        "clock=211/420;clock_tail=37/84;"
        "hostile=1/420,1/420;rescue=1/7;"
        "families=13340783196;all11=60301609751;mrange=3..70;maxlabel=80";
    const std::uint64_t semantic = fnv64(ledger);
    require(semantic == 0x56b59e7ee90c92f6ULL, "semantic FNV");

    std::cout << "THM4148_FIRST_WINDOW_WIDTH_ODD_TAIL_TRANSFER_INDEPENDENT_20260825\n";
    std::cout << "status=PASS;scope=independent integer-rational controls\n";
    std::cout << "width_gate=2/189;scale_gates=(2/63,2/189)\n";
    std::cout << "residual_clock_min=" << residual_min << '\n';
    std::cout << "consecutive_46_starts=(14,15,16,17,18,19,20,21,22)\n";
    std::cout << "consecutive_47_discriminant=-272;consecutive_46_discriminant=1296\n";
    std::cout << "block_15_60_interval=(" << left << ',' << right << ','
              << block_width << ")\n";
    std::cout << "block_15_60_surplus=" << block_width - gate << '\n';
    std::cout << "residual_clock=(211/420," << endpoint_tail_min << ")\n";
    std::cout << "moving_base_hostile=((1/420,1),(211/420,211),1/420)\n";
    std::cout << "moving_base_rescue=(1/105,53/105," << rescue_clearance << ")\n";
    std::cout << "eleven_subsets=" << families << '\n';
    std::cout << "all_width_body_minima=" << first_minimum << ".."
              << last_minimum << ";all_width_max_label=" << maximum_label
              << '\n';
    std::cout << "all_width_eleven_subsets=" << all_width_families << '\n';
    std::cout << "semantic_fnv64=" << std::hex << std::setw(16)
              << std::setfill('0') << semantic << std::dec << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "status=FAIL;error=" << error.what() << '\n';
    return 1;
  }
}
