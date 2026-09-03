#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr int kBodySize = 11;
constexpr int kFirstClock = 15;
constexpr int kLastClock = 43;
constexpr std::uint16_t kFullDivisorMask = (std::uint16_t{1} << 13) - 1;

struct TailMask {
  std::uint64_t danger_zero{};
  std::uint64_t danger_one{};
};

struct PacketScore {
  std::uint16_t deficit{};
  std::uint16_t covering_pairs{};
};

struct Fitness {
  int certified_clocks{};
  int total_deficit{};
  int min_covering_pairs{};
  std::array<PacketScore, kLastClock + 1> clock{};
};

class XorShift64 {
 public:
  explicit XorShift64(std::uint64_t seed) : state_(seed == 0 ? 1 : seed) {}

  std::uint64_t next() {
    state_ ^= state_ << 13;
    state_ ^= state_ >> 7;
    state_ ^= state_ << 17;
    return state_;
  }

  int uniform(const int bound) {
    return static_cast<int>(next() % static_cast<std::uint64_t>(bound));
  }

  double unit() {
    return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0);
  }

 private:
  std::uint64_t state_;
};

int abs_residue(const int value, const int modulus) {
  const int residue = ((value % modulus) + modulus) % modulus;
  return std::min(residue, modulus - residue);
}

std::uint16_t divisor_mask(const int speed) {
  std::uint16_t mask = 0;
  for (int divisor = 2; divisor <= 14; ++divisor) {
    if (speed % divisor == 0) {
      mask |= std::uint16_t{1} << (divisor - 2);
    }
  }
  return mask;
}

std::uint64_t safe_mask(const int speed, const int clock) {
  std::uint64_t mask = 0;
  for (int residue = 0; residue < clock; ++residue) {
    if (14 * abs_residue(speed * residue, clock) >= clock) {
      mask |= std::uint64_t{1} << residue;
    }
  }
  return mask;
}

std::vector<TailMask> build_tail_masks(const int clock) {
  std::vector<TailMask> tails;
  tails.reserve(clock);
  for (int tail = 1; tail < 2 * clock; tail += 2) {
    TailMask mask;
    for (int residue = 0; residue < clock; ++residue) {
      if (14 * abs_residue(tail * residue, 2 * clock) < 2 * clock) {
        mask.danger_zero |= std::uint64_t{1} << residue;
      }
      if (14 * abs_residue(tail * (residue + clock), 2 * clock) <
          2 * clock) {
        mask.danger_one |= std::uint64_t{1} << residue;
      }
    }
    tails.push_back(mask);
  }
  return tails;
}

PacketScore score_packet(const std::uint64_t packet,
                         const std::vector<TailMask>& tails) {
  int minimum_deficit = std::numeric_limits<int>::max();
  int covering_pairs = 0;
  for (std::size_t i = 0; i < tails.size(); ++i) {
    for (std::size_t j = i; j < tails.size(); ++j) {
      const std::uint64_t uncovered_zero =
          packet & ~(tails[i].danger_zero | tails[j].danger_zero);
      const std::uint64_t uncovered_one =
          packet & ~(tails[i].danger_one | tails[j].danger_one);
      const int deficit = std::popcount(uncovered_zero) +
                          std::popcount(uncovered_one);
      minimum_deficit = std::min(minimum_deficit, deficit);
      covering_pairs += deficit == 0 ? 1 : 0;
    }
  }
  return {static_cast<std::uint16_t>(minimum_deficit),
          static_cast<std::uint16_t>(covering_pairs)};
}

std::string format_body(const std::array<int, kBodySize>& body) {
  std::string result;
  for (std::size_t i = 0; i < body.size(); ++i) {
    if (i != 0) {
      result += ',';
    }
    result += std::to_string(body[i]);
  }
  return result;
}

bool divisor_complete(const std::array<int, kBodySize>& body,
                      const std::vector<std::uint16_t>& divisor_masks) {
  std::uint16_t mask = 0;
  for (const int speed : body) {
    mask |= divisor_masks[speed];
  }
  return mask == kFullDivisorMask;
}

bool less_fitness(const Fitness& left, const Fitness& right) {
  if (left.certified_clocks != right.certified_clocks) {
    return left.certified_clocks < right.certified_clocks;
  }
  if (left.total_deficit != right.total_deficit) {
    return left.total_deficit < right.total_deficit;
  }
  return left.min_covering_pairs > right.min_covering_pairs;
}

int energy(const Fitness& fitness) {
  return 1000 * fitness.certified_clocks + fitness.total_deficit;
}

void print_candidate(const char* label,
                     const std::array<int, kBodySize>& body,
                     const Fitness& fitness) {
  std::cout << label << " body=" << format_body(body)
            << " certified_clocks=" << fitness.certified_clocks
            << " total_deficit=" << fitness.total_deficit
            << " min_covering_pairs=" << fitness.min_covering_pairs << '\n';
  std::cout << label << " clock_profile=";
  bool first = true;
  for (int clock = kFirstClock; clock <= kLastClock; ++clock) {
    if (!first) {
      std::cout << ',';
    }
    first = false;
    std::cout << clock << ':' << fitness.clock[clock].deficit << '/'
              << fitness.clock[clock].covering_pairs;
  }
  std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  int maximum = 500;
  std::uint64_t iterations = 2'000'000;
  std::uint64_t seed = 6124341;
  if (argc >= 2) {
    maximum = std::stoi(argv[1]);
  }
  if (argc >= 3) {
    iterations = std::stoull(argv[2]);
  }
  if (argc >= 4) {
    seed = std::stoull(argv[3]);
  }
  if (argc > 4 || maximum < 14) {
    std::cerr << "usage: " << argv[0]
              << " [maximum>=14] [iterations] [seed]\n";
    return 2;
  }

  std::vector<std::uint16_t> divisor_masks(maximum + 1);
  std::array<std::vector<std::uint64_t>, kLastClock + 1> safe_masks;
  std::array<std::vector<TailMask>, kLastClock + 1> tails;
  std::array<std::unordered_map<std::uint64_t, PacketScore>, kLastClock + 1>
      caches;
  for (int speed = 1; speed <= maximum; ++speed) {
    divisor_masks[speed] = divisor_mask(speed);
  }
  for (int clock = kFirstClock; clock <= kLastClock; ++clock) {
    safe_masks[clock].resize(maximum + 1);
    for (int speed = 1; speed <= maximum; ++speed) {
      safe_masks[clock][speed] = safe_mask(speed, clock);
    }
    tails[clock] = build_tail_masks(clock);
  }

  auto evaluate = [&](const std::array<int, kBodySize>& body) {
    Fitness result;
    result.min_covering_pairs = std::numeric_limits<int>::max();
    for (int clock = kFirstClock; clock <= kLastClock; ++clock) {
      std::uint64_t packet = (std::uint64_t{1} << clock) - 1;
      for (const int speed : body) {
        packet &= safe_masks[clock][speed];
      }
      auto [iterator, inserted] = caches[clock].try_emplace(packet);
      if (inserted) {
        iterator->second = score_packet(packet, tails[clock]);
      }
      result.clock[clock] = iterator->second;
      if (iterator->second.deficit != 0) {
        ++result.certified_clocks;
        result.total_deficit += iterator->second.deficit;
      } else {
        result.min_covering_pairs =
            std::min(result.min_covering_pairs,
                     static_cast<int>(iterator->second.covering_pairs));
      }
    }
    if (result.min_covering_pairs == std::numeric_limits<int>::max()) {
      result.min_covering_pairs = 0;
    }
    return result;
  };

  std::array<int, kBodySize> current = maximum >= 408
      ? std::array<int, kBodySize>{
            160, 180, 185, 186, 190, 208, 210, 297, 300, 396, 408}
      : std::array<int, kBodySize>{1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14};
  if (!divisor_complete(current, divisor_masks)) {
    std::cerr << "inherited starting body is not divisor-complete\n";
    return 1;
  }
  Fitness current_fitness = evaluate(current);
  if (maximum >= 408 && current_fitness.certified_clocks != 4) {
    std::cerr << "inherited starting profile changed\n";
    print_candidate("START_MISMATCH", current, current_fitness);
    return 1;
  }
  std::array<int, kBodySize> best = current;
  Fitness best_fitness = current_fitness;
  print_candidate("START", current, current_fitness);

  XorShift64 random(seed);
  std::uint64_t accepted = 0;
  std::uint64_t divisor_rejects = 0;
  std::uint64_t duplicate_rejects = 0;
  for (std::uint64_t iteration = 1; iteration <= iterations; ++iteration) {
    std::array<int, kBodySize> proposal = current;
    const int index = random.uniform(kBodySize);
    proposal[index] = 1 + random.uniform(maximum);
    std::sort(proposal.begin(), proposal.end());
    if (std::adjacent_find(proposal.begin(), proposal.end()) != proposal.end()) {
      ++duplicate_rejects;
      continue;
    }
    if (!divisor_complete(proposal, divisor_masks)) {
      ++divisor_rejects;
      continue;
    }
    const Fitness proposal_fitness = evaluate(proposal);
    const int delta = energy(proposal_fitness) - energy(current_fitness);
    const double progress = static_cast<double>(iteration) /
                            static_cast<double>(iterations);
    const double temperature = 40.0 * (1.0 - progress) + 0.05;
    const bool accept = delta <= 0 || random.unit() < std::exp(-delta / temperature);
    if (accept) {
      current = proposal;
      current_fitness = proposal_fitness;
      ++accepted;
    }
    if (less_fitness(proposal_fitness, best_fitness)) {
      best = proposal;
      best_fitness = proposal_fitness;
      print_candidate("IMPROVEMENT", best, best_fitness);
      if (best_fitness.certified_clocks == 0) {
        break;
      }
    }
    if (iteration % 100'000 == 0 && best_fitness.certified_clocks > 0) {
      // Periodically restart at the best point rather than letting a hot
      // trajectory forget the exact obstruction already found.
      current = best;
      current_fitness = best_fitness;
    }
  }

  print_candidate("BEST", best, best_fitness);
  std::uint64_t cache_entries = 0;
  for (int clock = kFirstClock; clock <= kLastClock; ++clock) {
    cache_entries += caches[clock].size();
  }
  std::cout << "maximum=" << maximum << " iterations=" << iterations
            << " seed=" << seed << " accepted=" << accepted
            << " divisor_rejects=" << divisor_rejects
            << " duplicate_rejects=" << duplicate_rejects
            << " cache_entries=" << cache_entries << '\n';
  std::cout << "RESULT="
            << (best_fitness.certified_clocks == 0 ? "OWNER_ESCAPE"
                                                    : "NO_ESCAPE_FOUND")
            << '\n';
}
