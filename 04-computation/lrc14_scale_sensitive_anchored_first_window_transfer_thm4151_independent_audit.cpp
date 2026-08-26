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
static Rat operator/(Rat a, std::int64_t b) { return Rat(a.n, a.d * b); }
static bool operator==(Rat a, Rat b) { return a.n == b.n && a.d == b.d; }
static bool operator<(Rat a, Rat b) { return a.n * b.d < b.n * a.d; }
static bool operator<=(Rat a, Rat b) { return !(b < a); }

static std::ostream& operator<<(std::ostream& out, Rat value) {
  out << value.n;
  if (value.d != 1) out << '/' << value.d;
  return out;
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
  Rat result(1, 2);
  for (int speed : speeds) result = std::min(result, gap(speed, phase));
  return result;
}

static int last_cap(int first) { return (156 * first + 13) / 16; }

static Rat width(int first, int last) {
  return Rat(13, 14 * last) - Rat(1, 14 * first);
}

static Rat anchored_gate(int first) {
  return Rat(4 * first - 1, 14LL * first * (12 * first + 1));
}

static std::uint64_t choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  if (k > n - k) k = n - k;
  std::uint64_t answer = 1;
  for (int j = 1; j <= k; ++j) {
    answer *= n - k + j;
    answer /= j;
  }
  return answer;
}

static std::uint64_t family_count(int bound) {
  std::uint64_t answer = 0;
  for (int first = 1; first <= bound; ++first) {
    const int available = std::min(bound, last_cap(first)) - first;
    if (available >= 10) answer += choose(available, 10);
  }
  return answer;
}

static int translation_minimum(int span) {
  return std::max(1, (16 * span - 13 + 139) / 140);
}

static bool window_has_safe_lift(int first, int last, int a, int b) {
  const Rat left(1, 14 * first), right(13, 14 * last);
  std::vector<Rat> walls{left, right};
  for (int speed : {a, b}) {
    for (int integer = -2; integer <= speed + 2; ++integer) {
      for (int sign : {-1, 1}) {
        const Rat lower_wall(2 * (14 * integer + sign), 14 * speed);
        for (Rat wall : {lower_wall, lower_wall - Rat(1)}) {
          if (left <= wall && wall <= right) walls.push_back(wall);
        }
      }
    }
  }
  std::sort(walls.begin(), walls.end());
  walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
  std::vector<Rat> probes = walls;
  for (std::size_t i = 1; i < walls.size(); ++i) {
    probes.push_back((walls[i - 1] + walls[i]) / 2);
  }
  for (Rat y : probes) {
    for (Rat lift : {y / 2, (Rat(1) + y) / 2}) {
      if (Rat(1, 14) <= std::min(gap(a, lift), gap(b, lift))) return true;
    }
  }
  return false;
}

static std::vector<Rat> physical_walls(const std::vector<int>& speeds) {
  std::vector<Rat> walls{Rat(0), Rat(1)};
  for (int speed : speeds) {
    for (int integer = 0; integer <= speed; ++integer) {
      for (int sign : {-1, 1}) {
        const Rat wall(14 * integer + sign, 14 * speed);
        if (Rat(0) <= wall && wall <= Rat(1)) walls.push_back(wall);
      }
    }
  }
  std::sort(walls.begin(), walls.end());
  walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
  return walls;
}

static std::vector<std::pair<Rat, Rat>> safe_components(
    const std::vector<int>& speeds) {
  const Rat delta(1, 14);
  const std::vector<Rat> walls = physical_walls(speeds);
  std::vector<std::pair<Rat, Rat>> pieces;
  for (std::size_t i = 1; i < walls.size(); ++i) {
    const Rat left = walls[i - 1], right = walls[i];
    const Rat middle = (left + right) / 2;
    if (clearance(speeds, middle) < delta) continue;
    require(delta <= clearance(speeds, left) &&
                delta <= clearance(speeds, right),
            "closed safe component endpoints");
    if (!pieces.empty() && pieces.back().second == left) {
      pieces.back().second = right;
    } else {
      pieces.emplace_back(left, right);
    }
  }
  for (Rat wall : walls) {
    if (clearance(speeds, wall) < delta) continue;
    bool contained = false;
    for (const auto& piece : pieces) {
      if (piece.first <= wall && wall <= piece.second) contained = true;
    }
    if (!contained) pieces.emplace_back(wall, wall);
  }
  std::sort(pieces.begin(), pieces.end());
  return pieces;
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
    for (int first = 1; first <= 1000; ++first) {
      const int cap = last_cap(first);
      require(16 * cap <= 156 * first + 13 &&
                  156 * first + 13 < 16 * (cap + 1),
              "affine cap");
      require(anchored_gate(first) <= width(first, cap), "cap admitted");
      require(width(first, cap + 1) < anchored_gate(first),
              "first cap failure");
    }

    Rat endpoint_floor(1, 2);
    for (int first = 1; first <= 500; ++first) {
      const Rat y0(1, 14 * first);
      const Rat upper_lift = (Rat(1) + y0) / 2;
      for (int speed = 1; speed <= 12 * first; speed += 2) {
        const Rat value = gap(speed, upper_lift);
        require(delta < value, "endpoint tail floor");
        endpoint_floor = std::min(endpoint_floor, value);
      }
      const int hostile = 12 * first + 1;
      const Rat carrier_left(6, 7 * hostile);
      const Rat carrier_right(8, 7 * hostile);
      require(y0 - carrier_left == Rat(1, 14LL * first * hostile),
              "odd wall slack");
      require(carrier_right - y0 == anchored_gate(first),
              "remaining carrier room");
      const Rat y1(13, 14 * (last_cap(first) + 1));
      require(carrier_left < y0 && y0 <= y1 && y1 < carrier_right &&
                  carrier_right < Rat(1, 7),
              "first-window hostile containment");
    }

    const Rat hostile_right(13, 154);
    require(gap(1, hostile_right / 2) == Rat(13, 308),
            "hostile tail one");
    require(gap(13, (Rat(1) + hostile_right) / 2) == Rat(15, 308),
            "hostile tail thirteen");
    std::vector<int> rescue_row;
    for (int h = 1; h <= 11; ++h) rescue_row.push_back(2 * h);
    rescue_row.push_back(1);
    rescue_row.push_back(13);
    const Rat rescue_clearance = clearance(rescue_row, Rat(5, 24));
    require(rescue_clearance == Rat(1, 12), "hostile rescue");

    for (int span = 0; span <= 2000; ++span) {
      const int first = translation_minimum(span);
      require(first + span <= last_cap(first), "translation admitted");
      if (first > 1) {
        require(first - 1 + span > last_cap(first - 1),
                "translation minimal");
      }
    }
    require(translation_minimum(46) == 6 && last_cap(6) == 59,
            "47-label block translate");

    int primitive_pairs = 0, direct_rows = 0;
    for (int q = 3; q <= 51; q += 2) {
      for (int p = 1; p < q; p += 2) {
        if (std::gcd(p, q) != 1) continue;
        ++primitive_pairs;
        for (int first = 1; first <= 5; ++first) {
          require(window_has_safe_lift(first, last_cap(first), p, q),
                  "direct primitive wall-cell row");
          ++direct_rows;
        }
      }
    }
    require(primitive_pairs == 271 && direct_rows == 1355,
            "direct primitive row census");
    for (int first = 1; first <= 25; ++first) {
      require(!window_has_safe_lift(first, last_cap(first) + 1, 1,
                                    12 * first + 1),
              "direct carrier-sharp first failure");
    }

    std::vector<int> actual_hostile;
    for (int h = 2; h <= 21; ++h) actual_hostile.push_back(2 * h);
    const std::vector<std::pair<Rat, Rat>> expected_hostile_components{
        {Rat(1, 56), Rat(13, 588)},
        {Rat(281, 588), Rat(27, 56)},
        {Rat(29, 56), Rat(307, 588)},
        {Rat(575, 588), Rat(55, 56)},
    };
    require(safe_components(actual_hostile) == expected_hostile_components,
            "actual hostile body-safe components");
    actual_hostile.push_back(1);
    actual_hostile.push_back(25);
    const std::vector<Rat> hostile_walls = physical_walls(actual_hostile);
    require(hostile_walls.size() == 946, "actual hostile wall count");
    Rat hostile_wall_max(0), hostile_midpoint_max(0);
    for (std::size_t i = 0; i < hostile_walls.size(); ++i) {
      hostile_wall_max =
          std::max(hostile_wall_max, clearance(actual_hostile, hostile_walls[i]));
      if (i) {
        hostile_midpoint_max = std::max(
            hostile_midpoint_max,
            clearance(actual_hostile, (hostile_walls[i - 1] + hostile_walls[i]) / 2));
      }
    }
    require(hostile_wall_max == Rat(12, 175) &&
                hostile_midpoint_max == Rat(1, 16),
            "actual hostile wall/midpoint maxima");
    require(safe_components(actual_hostile).empty(),
            "actual cap hostile is unsafe");

    const std::vector<std::pair<int, std::uint64_t>> counts{
        {20, 75582ULL},
        {40, 792864735ULL},
        {80, 3548681310136ULL},
        {120, 392890789426111ULL},
        {160, 10500809430042208ULL},
        {200, 131378242150108190ULL},
    };
    for (const auto& row : counts) {
      require(family_count(row.first) == row.second, "finite family count");
    }
    require(family_count(20) + 1 == 75583, "combined minimum-one body");

    std::uint64_t density_num = 1, density_den = 1;
    for (int j = 0; j < 10; ++j) {
      density_num *= 35;
      density_den *= 39;
    }
    require(density_num == 2758547353515625ULL &&
                density_den == 8140406085191601ULL,
            "density fraction");

    const std::string ledger =
        "gate=(4m-1)/(14m(12m+1));cap=floor((156m+13)/16);"
        "endpoint=odd<=12m;hostile_b=12m+1;slack=1/(14mb);"
        "mechanism_hostile=(m,M,a,b)=(1,11,1,13);rescue=5/24,1/12;"
        "translate=ceil((16D-13)/140);span46_min=6;cap6=59;block6_52=47;"
        "direct_primitive_rows=1355;direct_first_failures=25;"
        "actual_cap_hostile=2{2..21}+{1,25}:unsafe;"
        "hostile_wall_probes=946;wallmax=12/175;midmax=1/16;"
        "anchored_counts20,40,80,120,160,200="
        "75582,792864735,3548681310136,392890789426111,"
        "10500809430042208,131378242150108190;"
        "density=(35/39)^10=2758547353515625/8140406085191601;"
        "combined_add_one_for_N>=11";
    const std::uint64_t semantic = fnv64(ledger);

    std::cout << "THM4151_SCALE_SENSITIVE_ANCHORED_FIRST_WINDOW_TRANSFER_"
                 "INDEPENDENT_20260825\n";
    std::cout << "status=PASS;scope=independent integer-rational controls\n";
    std::cout << "gate=(4m-1)/(14m(12m+1));integral_gate=16M<=156m+13\n";
    std::cout << "endpoint_test_floor_through_m500=" << endpoint_floor << '\n';
    std::cout << "sharp_carrier=b=12m+1;left_slack=1/(14mb);"
                 "remaining_room=gate\n";
    std::cout << "mechanism_hostile=(m=1,M=11,tails=(1,13),"
                 "window=(1/14,13/154))\n";
    std::cout << "hostile_rescue=(phase=5/24,clearance=" << rescue_clearance
              << ")\n";
    std::cout << "translation_span_gate=140m>=16D-13\n";
    std::cout << "span46_first_translate=6;cap_at_6=59;"
                 "block_6_52_labels=47;max_block_at_6=54\n";
    std::cout << "direct_wall_cells=1355_primitive_passes;"
                 "25_first_failures_confirmed\n";
    std::cout << "actual_cap_hostile=2{2..21}+{1,25};safe_components=0\n";
    std::cout << "actual_cap_hostile_probes=946_walls;wallmax=12/175;"
                 "midmax=1/16\n";
    std::cout << "anchored_family_counts=";
    for (std::size_t i = 0; i < counts.size(); ++i) {
      if (i) std::cout << ',';
      std::cout << counts[i].first << ':' << counts[i].second;
    }
    std::cout << '\n';
    std::cout << "combined_family_rule=anchored+1_for_N>=11;"
                 "combined_count_at_20=75583\n";
    std::cout << "asymptotic_density=" << density_num << '/' << density_den
              << '\n';
    std::cout << "semantic_fnv64=" << std::hex << std::setw(16)
              << std::setfill('0') << semantic << std::dec << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "status=FAIL;error=" << error.what() << '\n';
    return 1;
  }
}
