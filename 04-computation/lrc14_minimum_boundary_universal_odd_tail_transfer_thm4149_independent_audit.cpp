#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
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

static Rat operator+(Rat a, Rat b) {
  return Rat(a.n * b.d + b.n * a.d, a.d * b.d);
}
static Rat operator-(Rat a, Rat b) {
  return Rat(a.n * b.d - b.n * a.d, a.d * b.d);
}
static Rat operator*(std::int64_t a, Rat b) { return Rat(a * b.n, b.d); }
static Rat operator/(Rat a, std::int64_t b) { return Rat(a.n, a.d * b); }
static bool operator==(Rat a, Rat b) { return a.n == b.n && a.d == b.d; }
static bool operator<(Rat a, Rat b) { return a.n * b.d < b.n * a.d; }
static bool operator>=(Rat a, Rat b) { return !(a < b); }

static void require(bool condition, const std::string& label) {
  if (!condition) throw std::runtime_error(label);
}

static Rat mod_one(Rat value) {
  std::int64_t residue = value.n % value.d;
  if (residue < 0) residue += value.d;
  return Rat(residue, value.d);
}

static Rat gap(int speed, Rat phase) {
  const Rat residue = mod_one(speed * phase);
  return std::min(residue, Rat(1) - residue);
}

static Rat clearance(const std::vector<int>& speeds, Rat phase) {
  Rat answer(1, 2);
  for (int speed : speeds) answer = std::min(answer, gap(speed, phase));
  return answer;
}

static Rat physical(Rat body_phase, int sheet) {
  require(sheet == 0 || sheet == 1, "sheet");
  return (body_phase + Rat(sheet)) / 2;
}

static bool both_sheets_bad(int p, int q, Rat body_phase) {
  const Rat delta(1, 14);
  for (int sheet = 0; sheet <= 1; ++sheet) {
    const Rat x = physical(body_phase, sheet);
    if (!(std::min(gap(p, x), gap(q, x)) < delta)) return false;
  }
  return true;
}

static std::vector<Rat> cross_walls(int p, int q) {
  const Rat delta(1, 14);
  std::vector<Rat> walls{Rat(0), Rat(1)};
  for (int twice_shift : {0, 1}) {
    const Rat shift(twice_shift, 2);
    for (int speed : {p, q}) {
      for (int integer = 0; integer < speed; ++integer) {
        for (int sign : {-1, 1}) {
          const Rat phase_wall = (Rat(integer) + sign * delta) / speed;
          walls.push_back(mod_one(2 * (phase_wall - shift)));
        }
      }
    }
  }
  std::sort(walls.begin(), walls.end());
  walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
  return walls;
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
    const Rat delta(1, 14);
    require(27 * (13 * 1 - 11) >= 4 * 1 * 11, "minimum-one cap positive");
    require(27 * (13 * 1 - 12) < 4 * 1 * 12, "minimum-one cap sharp");
    require(27 * (13 * 2 - 20) >= 4 * 2 * 20, "minimum-two cap positive");
    require(27 * (13 * 2 - 21) < 4 * 2 * 21, "minimum-two cap sharp");

    std::vector<std::pair<int, int>> pairs;
    for (int q = 3; q <= 25; q += 2) {
      for (int p = 1; p < q; p += 2) {
        if (std::gcd(p, q) == 1) pairs.emplace_back(p, q);
      }
    }
    require(pairs.size() == 68, "primitive residual census");

    std::vector<int> full_one, doubled_one;
    for (int h = 1; h <= 11; ++h) {
      full_one.push_back(h);
      doubled_one.push_back(2 * h);
    }
    const std::vector<Rat> one_y{Rat(1, 14), Rat(6, 77), Rat(3, 14)};
    const std::vector<int> one_sheet{0, 1, 0};
    const std::vector<Rat> expected_body{Rat(1, 14), Rat(6, 77), Rat(1, 14)};
    for (int i = 0; i < 3; ++i) {
      require(clearance(full_one, one_y[i]) == expected_body[i],
              "minimum-one body clock");
    }

    std::vector<int> one_count(3, 0);
    std::vector<Rat> one_min(3, Rat(1, 2));
    for (const auto& [p, q] : pairs) {
      int category = p >= 3 ? 0 : (q != 13 ? 1 : 2);
      const Rat x = physical(one_y[category], one_sheet[category]);
      const Rat tail = std::min(gap(p, x), gap(q, x));
      std::vector<int> row = doubled_one;
      row.push_back(p);
      row.push_back(q);
      require(tail >= delta && clearance(row, x) >= delta,
              "minimum-one residual clock");
      ++one_count[category];
      one_min[category] = std::min(one_min[category], tail);
    }
    require(one_count == std::vector<int>({56, 11, 1}),
            "minimum-one partition");
    require(one_min == std::vector<Rat>({Rat(3, 28), Rat(1, 14), Rat(3, 28)}),
            "minimum-one tail minima");

    const Rat first_left(1, 14), first_right(13, 154);
    const auto walls = cross_walls(1, 13);
    auto upper = std::upper_bound(walls.begin(), walls.end(), first_left);
    require(upper != walls.begin() && upper != walls.end(), "hostile wall location");
    const Rat wall_left = *(upper - 1), wall_right = *upper;
    require(wall_left == Rat(6, 91) && wall_right == Rat(8, 91),
            "hostile component walls");
    require(wall_left < first_left && first_right < wall_right,
            "first window strictly trapped");
    require(both_sheets_bad(1, 13, (first_left + first_right) / 2),
            "hostile cell sign");
    require(!both_sheets_bad(1, 13, wall_left) &&
            !both_sheets_bad(1, 13, wall_right), "hostile open endpoints");

    const Rat isolated(3, 14), epsilon(1, 10000);
    require(clearance(full_one, isolated) == delta, "isolated clock safe");
    require(gap(5, isolated) == delta && gap(9, isolated) == delta,
            "isolated owners");
    require(clearance(full_one, isolated - epsilon) < delta &&
            clearance(full_one, isolated + epsilon) < delta,
            "isolated hostile neighbours");

    std::vector<int> full_two, doubled_two;
    for (int h = 2; h <= 20; ++h) {
      full_two.push_back(h);
      doubled_two.push_back(2 * h);
    }
    require(clearance(full_two, Rat(1, 28)) == delta &&
            clearance(full_two, Rat(13, 280)) == delta,
            "minimum-two common endpoints");
    const std::vector<Rat> two_y{Rat(1, 28), Rat(1, 28), Rat(13, 280)};
    const std::vector<int> two_sheet{1, 0, 1};
    std::vector<int> two_count(3, 0);
    std::vector<Rat> two_min(3, Rat(1, 2));
    for (const auto& [p, q] : pairs) {
      int category = q <= 23 ? 0 : (p >= 5 ? 1 : 2);
      const Rat x = physical(two_y[category], two_sheet[category]);
      const Rat tail = std::min(gap(p, x), gap(q, x));
      std::vector<int> row = doubled_two;
      row.push_back(p);
      row.push_back(q);
      require(tail >= delta && clearance(row, x) >= delta,
              "minimum-two residual clock");
      ++two_count[category];
      two_min[category] = std::min(two_min[category], tail);
    }
    require(two_count == std::vector<int>({58, 8, 2}),
            "minimum-two partition");
    require(gap(25, physical(Rat(1, 28), 1)) == Rat(3, 56),
            "minimum-two upper hostile");
    require(gap(1, physical(Rat(1, 28), 0)) == Rat(1, 56) &&
            gap(3, physical(Rat(1, 28), 0)) == Rat(3, 56),
            "minimum-two lower hostiles");
    require(two_min[2] == Rat(9, 112), "minimum-two endpoint rescue");

    std::uint64_t all_eleven = 0;
    int first_minimum = -1, last_minimum = -1, maximum_label = 0;
    std::vector<std::uint64_t> first_counts;
    for (int minimum = 1; minimum <= 1000; ++minimum) {
      const int last_cap = 351 * minimum / (4 * minimum + 27);
      if (last_cap - minimum < 10) continue;
      if (first_minimum < 0) first_minimum = minimum;
      last_minimum = minimum;
      maximum_label = std::max(maximum_label, last_cap);
      const std::uint64_t count = choose(last_cap - minimum, 10);
      if (minimum <= 2) first_counts.push_back(count);
      all_eleven += count;
    }
    require(first_minimum == 1 && last_minimum == 70 && maximum_label == 80,
            "all-width range");
    require(first_counts == std::vector<std::uint64_t>({1, 43758}),
            "boundary body counts");
    require(all_eleven == 60301653510ULL, "all-width body count");

    const std::string ledger =
        "gate=2/189;m1cap=11;m2cap=20;pairs=68;"
        "m1partition=56,11,1;m1clocks=1/14:0,6/77:1,3/14:0;"
        "m1bodygaps=1/14,6/77,1/14;m1hostile=6/91..8/91;m1isolated=3/14;"
        "m2partition=58,8,2;m2clocks=1/28:1,1/28:0,13/280:1;"
        "m2hostile=3/56;m2rescue=9/112;"
        "all11=60301653510;added=43759;mrange=1..70;maxlabel=80";
    const std::uint64_t semantic = fnv64(ledger);
    require(semantic == 0x6784b7e0b01a759dULL, "semantic FNV");

    std::cout << "THM4149_MINIMUM_BOUNDARY_ODD_TAIL_TRANSFER_INDEPENDENT_20260825\n";
    std::cout << "status=PASS;scope=independent integer-rational audit\n";
    std::cout << "width_gate=2/189;minimum_caps=(1:11,2:20)\n";
    std::cout << "primitive_residual_pairs=" << pairs.size() << '\n';
    std::cout << "minimum_one_partition=(56,11,1)\n";
    std::cout << "minimum_one_clocks=((1/14,0),(6/77,1),(3/14,0))\n";
    std::cout << "minimum_one_body_gaps=(1/14,6/77,1/14)\n";
    std::cout << "minimum_one_first_window_hostile=(1,13;6/91,8/91)\n";
    std::cout << "minimum_one_isolated_clock=3/14\n";
    std::cout << "minimum_two_partition=(58,8,2)\n";
    std::cout << "minimum_two_clocks=((1/28,1),(1/28,0),(13/280,1))\n";
    std::cout << "minimum_two_hostile=3/56;minimum_two_rescue=9/112\n";
    std::cout << "added_eleven_bodies=43759\n";
    std::cout << "all_width_body_minima=" << first_minimum << ".." << last_minimum
              << ";all_width_max_label=" << maximum_label << '\n';
    std::cout << "all_width_eleven_bodies=" << all_eleven << '\n';
    std::cout << "semantic_fnv64=" << std::hex << std::setw(16)
              << std::setfill('0') << semantic << std::dec << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "status=FAIL;error=" << error.what() << '\n';
    return 1;
  }
}
