#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr int kFirstMaximum = 31;
constexpr int kLastMaximum = 40;
constexpr int kChoose = 10;
constexpr int kFirstClock = 15;
constexpr int kLastClock = 43;
constexpr std::uint16_t kFullDivisorMask = (std::uint16_t{1} << 13) - 1;

struct TailWord {
  std::uint64_t eligible{};
  std::uint64_t owner{};
  std::uint64_t danger_zero{};
  std::uint64_t danger_one{};
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

std::vector<TailWord> tail_words(const int clock) {
  std::vector<TailWord> words;
  words.reserve(clock);
  for (int tail = 1; tail < 2 * clock; tail += 2) {
    TailWord word;
    for (int residue = 0; residue < clock; ++residue) {
      if (7 * abs_residue(tail * residue, clock) < clock) {
        word.eligible |= std::uint64_t{1} << residue;
        const int nearest = (2 * tail * residue + clock) / (2 * clock);
        if (nearest % 2 != 0) {
          word.owner |= std::uint64_t{1} << residue;
        }
      }
      if (14 * abs_residue(tail * residue, 2 * clock) < 2 * clock) {
        word.danger_zero |= std::uint64_t{1} << residue;
      }
      if (14 * abs_residue(tail * (residue + clock), 2 * clock) <
          2 * clock) {
        word.danger_one |= std::uint64_t{1} << residue;
      }
    }
    words.push_back(word);
  }
  return words;
}

bool has_owner_word_certificate(const std::uint64_t packet,
                                const std::vector<TailWord>& tails) {
  if (packet == 0) {
    return false;
  }
  std::unordered_set<std::uint64_t> words;
  for (const TailWord& tail : tails) {
    if ((packet & ~tail.eligible) == 0) {
      words.insert(packet & tail.owner);
    }
  }
  for (const std::uint64_t word : words) {
    if (words.contains(packet ^ word)) {
      return false;
    }
  }
  return true;
}

bool has_literal_lift_certificate(const std::uint64_t packet,
                                  const std::vector<TailWord>& tails) {
  if (packet == 0) {
    return false;
  }
  for (std::size_t i = 0; i < tails.size(); ++i) {
    for (std::size_t j = i; j < tails.size(); ++j) {
      const std::uint64_t zero_cover =
          (tails[i].danger_zero | tails[j].danger_zero) & packet;
      const std::uint64_t one_cover =
          (tails[i].danger_one | tails[j].danger_one) & packet;
      if (zero_cover == packet && one_cover == packet) {
        return false;
      }
    }
  }
  return true;
}

std::string format_core(const std::array<int, 11>& core) {
  std::string result;
  for (std::size_t i = 0; i < core.size(); ++i) {
    if (i != 0) {
      result += ',';
    }
    result += std::to_string(core[i]);
  }
  return result;
}

std::uint64_t binomial(int n, int k) {
  if (k < 0 || k > n) {
    return 0;
  }
  k = std::min(k, n - k);
  std::uint64_t value = 1;
  for (int j = 1; j <= k; ++j) {
    value = value * static_cast<std::uint64_t>(n - k + j) /
            static_cast<std::uint64_t>(j);
  }
  return value;
}

}  // namespace

int main(const int argc, char** argv) {
  bool reverse_clocks = false;
  if (argc == 2 && std::string(argv[1]) == "--reverse-clocks") {
    reverse_clocks = true;
  } else if (argc != 1) {
    std::cerr << "usage: " << argv[0] << " [--reverse-clocks]\n";
    return 2;
  }

  std::array<std::uint16_t, kLastMaximum + 1> divisor_masks{};
  for (int speed = 1; speed <= kLastMaximum; ++speed) {
    divisor_masks[speed] = divisor_mask(speed);
  }

  std::array<std::array<std::uint64_t, kLastMaximum + 1>, kLastClock + 1>
      safe_masks{};
  std::array<std::vector<TailWord>, kLastClock + 1> tails{};
  std::array<std::unordered_map<std::uint64_t, bool>, kLastClock + 1> caches{};
  for (int clock = kFirstClock; clock <= kLastClock; ++clock) {
    for (int speed = 1; speed <= kLastMaximum; ++speed) {
      safe_masks[clock][speed] = safe_mask(speed, clock);
    }
    tails[clock] = tail_words(clock);
  }

  std::cout << "LRC14 ARBITRARY ENTRY OWNER-WORD CLOSURE, EXACT MAXIMA 31..40\n";
  std::cout << "clock_order=" << (reverse_clocks ? "descending" : "ascending")
            << " clocks=" << kFirstClock << ".." << kLastClock << '\n';

  std::uint64_t grand_searched = 0;
  std::uint64_t grand_divisor_complete = 0;
  std::uint64_t grand_uncertified = 0;
  for (int maximum = kFirstMaximum; maximum <= kLastMaximum; ++maximum) {
    std::uint64_t searched = 0;
    std::uint64_t divisor_complete = 0;
    std::uint64_t primitive = 0;
    std::uint64_t nonprimitive = 0;
    std::map<int, std::uint64_t> selected_clock_histogram;
    std::vector<std::array<int, 11>> uncertified;
    std::vector<std::array<int, 11>> last_clock_boundary;

    std::array<int, kChoose> choice{};
    std::iota(choice.begin(), choice.end(), 1);
    while (true) {
      ++searched;
      std::array<int, 11> core{};
      std::copy(choice.begin(), choice.end(), core.begin());
      core.back() = maximum;

      std::uint16_t cover = divisor_masks[maximum];
      int content = maximum;
      for (const int speed : choice) {
        cover |= divisor_masks[speed];
        content = std::gcd(content, speed);
      }
      if (cover == kFullDivisorMask) {
        ++divisor_complete;
        if (content == 1) {
          ++primitive;
        } else {
          ++nonprimitive;
        }

        int selected_clock = 0;
        for (int offset = 0; offset <= kLastClock - kFirstClock; ++offset) {
          const int clock = reverse_clocks ? kLastClock - offset
                                           : kFirstClock + offset;
          std::uint64_t packet = (std::uint64_t{1} << clock) - 1;
          for (const int speed : core) {
            packet &= safe_masks[clock][speed];
          }
          auto [iterator, inserted] = caches[clock].try_emplace(packet, false);
          if (inserted) {
            const bool owner_certificate =
                has_owner_word_certificate(packet, tails[clock]);
            const bool literal_certificate =
                has_literal_lift_certificate(packet, tails[clock]);
            if (owner_certificate != literal_certificate) {
              std::cerr << "owner/literal mismatch clock=" << clock
                        << " packet=" << packet << '\n';
              return 1;
            }
            iterator->second = owner_certificate;
          }
          if (iterator->second) {
            selected_clock = clock;
            ++selected_clock_histogram[clock];
            break;
          }
        }
        if (selected_clock == 0) {
          uncertified.push_back(core);
        } else if (!reverse_clocks && selected_clock == kLastClock) {
          last_clock_boundary.push_back(core);
        }
      }

      int index = kChoose - 1;
      while (index >= 0 &&
             choice[index] == (maximum - 1) - (kChoose - 1 - index)) {
        --index;
      }
      if (index < 0) {
        break;
      }
      ++choice[index];
      for (int j = index + 1; j < kChoose; ++j) {
        choice[j] = choice[j - 1] + 1;
      }
    }

    grand_searched += searched;
    grand_divisor_complete += divisor_complete;
    grand_uncertified += uncertified.size();
    std::uint64_t histogram_total = 0;
    for (const auto& [clock, count] : selected_clock_histogram) {
      static_cast<void>(clock);
      histogram_total += count;
    }
    if (searched != binomial(maximum - 1, kChoose) ||
        divisor_complete != primitive + nonprimitive ||
        histogram_total + uncertified.size() != divisor_complete) {
      std::cerr << "layer accounting failure maximum=" << maximum << '\n';
      return 1;
    }
    std::cout << "maximum=" << maximum << " searched=" << searched
              << " divisor_complete=" << divisor_complete
              << " primitive=" << primitive
              << " nonprimitive=" << nonprimitive
              << " uncertified=" << uncertified.size() << '\n';
    std::cout << "selected_clock_histogram=";
    bool first = true;
    for (const auto& [clock, count] : selected_clock_histogram) {
      if (!first) {
        std::cout << ',';
      }
      first = false;
      std::cout << clock << ':' << count;
    }
    std::cout << '\n';
    constexpr std::size_t kDisplayLimit = 40;
    for (std::size_t i = 0; i < std::min(kDisplayLimit, uncertified.size()); ++i) {
      std::cout << "uncertified_core=" << format_core(uncertified[i]) << '\n';
    }
    for (const auto& core : last_clock_boundary) {
      std::cout << "last_clock_boundary_core=" << format_core(core) << '\n';
    }
  }

  const std::array<int, 11> ap_core = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  std::uint64_t ap_packet = (std::uint64_t{1} << 12) - 1;
  for (const int speed : ap_core) {
    ap_packet &= safe_mask(speed, 12);
  }
  const std::vector<TailWord> ap_tails = tail_words(12);
  const bool ap_owner_certificate =
      has_owner_word_certificate(ap_packet, ap_tails);
  const bool ap_literal_certificate =
      has_literal_lift_certificate(ap_packet, ap_tails);
  int ap_clearance_numerator = 24;
  for (const int speed : ap_core) {
    ap_clearance_numerator =
        std::min(ap_clearance_numerator, abs_residue(2 * speed * 5, 24));
  }
  ap_clearance_numerator =
      std::min(ap_clearance_numerator, abs_residue(1 * 5, 24));
  ap_clearance_numerator =
      std::min(ap_clearance_numerator, abs_residue(9 * 5, 24));
  if (!ap_owner_certificate || !ap_literal_certificate ||
      ap_clearance_numerator != 2) {
    std::cerr << "non-divisor-complete AP control failed\n";
    return 1;
  }
  std::cout << "non_divisor_complete_control=H_1_to_11 missing=12,13,14"
               " universal_owner_clock=12 tails=1,9 witness=5/24"
               " clearance=1/12\n";
  std::cout << "grand_searched=" << grand_searched
            << " grand_divisor_complete=" << grand_divisor_complete
            << " grand_uncertified=" << grand_uncertified << '\n';
  std::uint64_t packet_audits = 0;
  for (int clock = kFirstClock; clock <= kLastClock; ++clock) {
    packet_audits += caches[clock].size();
  }
  std::cout << "owner_literal_packet_audits=" << packet_audits << '\n';
  std::cout << "RESULT=" << (grand_uncertified == 0 ? "PASS" : "SURVIVORS")
            << '\n';
}
