// Exact sparse degree-21 Macaulay ranks for the full-support-five FC chart.
//
// Build and run, for example:
//
//   g++ -O3 -DNDEBUG -std=c++20 \
//     04-computation/degree_four_support5_degree21_sparse_rank_scout_20260803.cpp \
//     -o /tmp/degree_four_support5_rank && \
//   /tmp/degree_four_support5_rank
//
// The reducer is deterministic.  A row is held in a reusable dense accumulator,
// while echelon pivots are stored sparsely.  Fixed-degree lexicographic order
// makes the first generator triangular and keeps the memory footprint far below
// that of the 13972 x 12650 dense matrix.  All arithmetic is exact in F_p.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Exp = std::array<std::uint8_t, 5>;

struct Entry {
  std::uint16_t column;
  std::uint8_t value;
};

struct SparseRow {
  std::vector<Entry> entries;
};

struct Polynomial {
  int degree = 0;
  std::vector<std::pair<Exp, std::uint8_t>> terms;
};

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

std::uint64_t choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t answer = 1;
  for (int i = 1; i <= k; ++i) {
    answer = answer * static_cast<std::uint64_t>(n - k + i) /
             static_cast<std::uint64_t>(i);
  }
  return answer;
}

int mod_pow(int base, int exponent, int prime) {
  int answer = 1;
  while (exponent) {
    if (exponent & 1) answer = static_cast<int>(1LL * answer * base % prime);
    base = static_cast<int>(1LL * base * base % prime);
    exponent >>= 1;
  }
  return answer;
}

int mod_inverse(int value, int prime) {
  require(value % prime != 0, "attempted inversion of zero");
  return mod_pow(value, prime - 2, prime);
}

bool is_prime(int value) {
  if (value < 2) return false;
  for (int d = 2; d * d <= value; ++d) {
    if (value % d == 0) return false;
  }
  return true;
}

void compositions_desc_rec(int total, int parts, int position, Exp& current,
                           std::vector<Exp>& out) {
  if (parts == 1) {
    current[position] = static_cast<std::uint8_t>(total);
    out.push_back(current);
    return;
  }
  for (int first = total; first >= 0; --first) {
    current[position] = static_cast<std::uint8_t>(first);
    compositions_desc_rec(total - first, parts - 1, position + 1, current, out);
  }
}

std::vector<Exp> compositions_desc(int total, int parts) {
  require(parts >= 1 && parts <= 5, "composition arity");
  Exp current{};
  std::vector<Exp> out;
  out.reserve(static_cast<std::size_t>(choose(total + parts - 1, parts - 1)));
  compositions_desc_rec(total, parts, 0, current, out);
  return out;
}

std::uint32_t encode_exp(const Exp& exponent) {
  // Every exponent here is at most 21, so five base-22 digits are collision-free.
  std::uint32_t code = 0;
  for (int i = 0; i < 5; ++i) code = code * 22 + exponent[i];
  return code;
}

struct ModularMoments {
  int prime;
  std::array<int, 87> factorial{};
  std::array<int, 87> inverse_factorial{};

  explicit ModularMoments(int p) : prime(p) {
    require(is_prime(prime), "modulus must be prime");
    require(prime > 86, "denominator guard requires p > 86");
    factorial[0] = 1;
    for (int i = 1; i <= 86; ++i) factorial[i] = 1LL * factorial[i - 1] * i % prime;
    inverse_factorial[86] = mod_inverse(factorial[86], prime);
    for (int i = 86; i >= 1; --i) {
      inverse_factorial[i - 1] = 1LL * inverse_factorial[i] * i % prime;
    }
  }

  int multinomial(int total, const Exp& exponents, int parts) const {
    int answer = factorial[total];
    for (int i = 0; i < parts; ++i) {
      answer = 1LL * answer * inverse_factorial[exponents[i]] % prime;
    }
    return answer;
  }

  int coefficient_c(int a, int b) const {
    int answer = 0;
    for (int k = 0; k <= std::min(a, b); ++k) {
      if ((a - k) % 3 != 0 || (b - k) % 3 != 0) continue;
      int i = (a - k) / 3;
      int j = (b - k) / 3;
      int total = i + j + k;
      int term = factorial[total];
      term = 1LL * term * inverse_factorial[i] % prime;
      term = 1LL * term * inverse_factorial[j] % prime;
      term = 1LL * term * inverse_factorial[k] % prime;
      term = 1LL * term * mod_pow(3, k, prime) % prime;
      answer += term;
      if (answer >= prime) answer -= prime;
    }
    return answer;
  }

  int mu(int a, int b) const {
    require(a + b + 2 <= 86, "simplex moment degree guard");
    int answer = 2 % prime;
    answer = 1LL * answer * factorial[a] % prime;
    answer = 1LL * answer * factorial[b] % prime;
    answer = 1LL * answer * coefficient_c(a, b) % prime;
    answer = 1LL * answer * inverse_factorial[a + b + 2] % prime;
    return answer;
  }
};

constexpr std::array<std::array<int, 2>, 5> kBasis{{
    {{0, 1}}, {{2, 0}}, {{1, 2}}, {{3, 1}}, {{0, 4}},
}};

Polynomial moment_polynomial(int degree, int variable_count,
                             const std::array<int, 5>& retained,
                             const ModularMoments& moments) {
  Polynomial polynomial;
  polynomial.degree = degree;
  for (const Exp& counts : compositions_desc(degree, variable_count)) {
    int z_power = 0;
    int zbar_power = 0;
    for (int i = 0; i < variable_count; ++i) {
      z_power += counts[i] * kBasis[retained[i]][0];
      zbar_power += counts[i] * kBasis[retained[i]][1];
    }
    int coefficient = moments.multinomial(degree, counts, variable_count);
    coefficient = 1LL * coefficient * moments.mu(z_power, zbar_power) % moments.prime;
    if (coefficient != 0) {
      polynomial.terms.emplace_back(counts, static_cast<std::uint8_t>(coefficient));
    }
  }
  require(!polynomial.terms.empty(), "moment polynomial vanished modulo p");
  return polynomial;
}

struct RankRecord {
  int rows = 0;
  int columns = 0;
  int rank = 0;
  std::size_t pivot_nonzeros = 0;
  std::size_t peak_pivot_nonzeros = 0;
  int selected_minor_determinant = 0;
  std::vector<int> selected_rows;
  std::vector<int> selected_columns;
};

class SparseReducer {
 public:
  SparseReducer(int target_degree, int variable_count, int prime)
      : target_degree_(target_degree),
        variable_count_(variable_count),
        prime_(prime),
        columns_(compositions_desc(target_degree, variable_count)),
        pivots_(columns_.size()),
        accumulator_(columns_.size(), 0) {
    require(prime < 256, "one-byte sparse coefficients require p < 256");
    require(columns_.size() < std::numeric_limits<std::uint16_t>::max(),
            "two-byte column indices");
    column_index_.reserve(columns_.size() * 2);
    for (std::size_t i = 0; i < columns_.size(); ++i) {
      column_index_.emplace(encode_exp(columns_[i]), static_cast<int>(i));
    }
  }

  int columns() const { return static_cast<int>(columns_.size()); }
  int rank() const { return rank_; }
  int rows() const { return rows_; }
  std::size_t pivot_nonzeros() const { return pivot_nonzeros_; }
  std::size_t peak_pivot_nonzeros() const { return peak_pivot_nonzeros_; }
  int selected_minor_determinant() const { return selected_minor_determinant_; }
  const std::vector<int>& selected_rows() const { return selected_rows_; }
  const std::vector<int>& selected_columns() const { return selected_columns_; }

  void add_generator(const Polynomial& polynomial) {
    require(polynomial.degree <= target_degree_, "generator above target degree");
    auto shifts = compositions_desc(target_degree_ - polynomial.degree,
                                    variable_count_);
    for (const Exp& shift : shifts) {
      std::fill(accumulator_.begin(), accumulator_.end(), 0);
      for (const auto& [term, coefficient] : polynomial.terms) {
        Exp monomial{};
        for (int i = 0; i < variable_count_; ++i) {
          monomial[i] = static_cast<std::uint8_t>(term[i] + shift[i]);
        }
        auto found = column_index_.find(encode_exp(monomial));
        require(found != column_index_.end(), "Macaulay output monomial");
        require(accumulator_[found->second] == 0,
                "distinct terms remain distinct under a fixed shift");
        accumulator_[found->second] = coefficient;
      }
      reduce_accumulator(rows_);
      ++rows_;
    }
  }

 private:
  void reduce_accumulator(int source_row) {
    int lead = 0;
    while (true) {
      while (lead < columns() && accumulator_[lead] == 0) ++lead;
      if (lead == columns()) return;
      if (pivots_[lead].entries.empty()) break;
      int factor = accumulator_[lead];
      for (const Entry& entry : pivots_[lead].entries) {
        int value = accumulator_[entry.column] - factor * entry.value;
        value %= prime_;
        if (value < 0) value += prime_;
        accumulator_[entry.column] = static_cast<std::uint8_t>(value);
      }
      require(accumulator_[lead] == 0, "normalized pivot cancellation");
      ++lead;
    }

    const int raw_pivot = accumulator_[lead];
    int inverse = mod_inverse(raw_pivot, prime_);
    SparseRow row;
    // Reserve a modest amount; geometric vector growth remains bounded by the
    // final persistent pivot storage and avoids a dense allocation per row.
    row.entries.reserve(64);
    for (int column = lead; column < columns(); ++column) {
      if (accumulator_[column] == 0) continue;
      int value = accumulator_[column] * inverse % prime_;
      row.entries.push_back(Entry{static_cast<std::uint16_t>(column),
                                  static_cast<std::uint8_t>(value)});
    }
    require(!row.entries.empty() && row.entries.front().column == lead &&
                row.entries.front().value == 1,
            "new normalized sparse pivot");
    pivot_nonzeros_ += row.entries.size();
    peak_pivot_nonzeros_ = std::max(peak_pivot_nonzeros_, row.entries.size());
    pivots_[lead] = std::move(row);
    selected_minor_determinant_ =
        1LL * selected_minor_determinant_ * raw_pivot % prime_;
    selected_rows_.push_back(source_row);
    selected_columns_.push_back(lead);
    ++rank_;
  }

  int target_degree_;
  int variable_count_;
  int prime_;
  std::vector<Exp> columns_;
  std::unordered_map<std::uint32_t, int> column_index_;
  std::vector<SparseRow> pivots_;
  std::vector<std::uint8_t> accumulator_;
  int rows_ = 0;
  int rank_ = 0;
  std::size_t pivot_nonzeros_ = 0;
  std::size_t peak_pivot_nonzeros_ = 0;
  int selected_minor_determinant_ = 1;
  std::vector<int> selected_rows_;
  std::vector<int> selected_columns_;
};

RankRecord rank_moment_system(int prime, int variable_count,
                              const std::array<int, 5>& retained,
                              const std::vector<int>& degrees,
                              int target_degree, bool progress) {
  ModularMoments moments(prime);
  SparseReducer reducer(target_degree, variable_count, prime);
  for (int degree : degrees) {
    Polynomial generator =
        moment_polynomial(degree, variable_count, retained, moments);
    const int rank_before = reducer.rank();
    const int rows_before = reducer.rows();
    reducer.add_generator(generator);
    if (progress) {
      std::cout << "  generator=M" << std::setw(2) << degree
                << " rows_added=" << std::setw(5) << reducer.rows() - rows_before
                << " rank_added=" << std::setw(5) << reducer.rank() - rank_before
                << " rank=" << std::setw(5) << reducer.rank()
                << " pivot_nnz=" << reducer.pivot_nonzeros()
                << " peak_row_nnz=" << reducer.peak_pivot_nonzeros() << '\n';
      std::cout.flush();
    }
  }
  return RankRecord{reducer.rows(), reducer.columns(), reducer.rank(),
                    reducer.pivot_nonzeros(), reducer.peak_pivot_nonzeros(),
                    reducer.selected_minor_determinant(), reducer.selected_rows(),
                    reducer.selected_columns()};
}

std::string compress_ordered_ranges(const std::vector<int>& indices) {
  require(!indices.empty(), "nonempty certificate index list");
  std::string answer;
  std::size_t begin = 0;
  while (begin < indices.size()) {
    std::size_t end = begin;
    while (end + 1 < indices.size() && indices[end + 1] == indices[end] + 1) ++end;
    if (!answer.empty()) answer += ',';
    answer += std::to_string(indices[begin]);
    if (end != begin) answer += '-' + std::to_string(indices[end]);
    begin = end + 1;
  }
  return answer;
}

int sign_to_sorted_order(const std::vector<int>& sequence) {
  std::vector<int> sorted = sequence;
  std::sort(sorted.begin(), sorted.end());
  require(std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end(),
          "permutation sign needs distinct indices");
  std::vector<int> permutation(sequence.size());
  for (std::size_t i = 0; i < sequence.size(); ++i) {
    auto found = std::lower_bound(sorted.begin(), sorted.end(), sequence[i]);
    require(found != sorted.end() && *found == sequence[i],
            "column occurs in its sorted order");
    permutation[i] = static_cast<int>(found - sorted.begin());
  }
  std::vector<bool> seen(sequence.size(), false);
  int cycles = 0;
  for (std::size_t i = 0; i < permutation.size(); ++i) {
    if (seen[i]) continue;
    ++cycles;
    int cursor = static_cast<int>(i);
    while (!seen[cursor]) {
      seen[cursor] = true;
      cursor = permutation[cursor];
    }
  }
  return ((static_cast<int>(sequence.size()) - cycles) & 1) ? -1 : 1;
}

void print_record(const std::string& name, int prime, const RankRecord& record,
                  bool print_certificate) {
  std::cout << name << " p=" << prime << " rows=" << record.rows
            << " columns=" << record.columns << " rank=" << record.rank
            << " nullity=" << record.columns - record.rank
            << " pivot_nnz=" << record.pivot_nonzeros
            << " peak_row_nnz=" << record.peak_pivot_nonzeros << '\n';
  if (print_certificate) {
    const int sign = sign_to_sorted_order(record.selected_columns);
    const int sorted_determinant =
        sign == 1 ? record.selected_minor_determinant
                  : prime - record.selected_minor_determinant;
    std::vector<int> sorted_columns = record.selected_columns;
    std::sort(sorted_columns.begin(), sorted_columns.end());
    std::cout << name << "_minor_determinant_pivot_order_mod_" << prime << '='
              << record.selected_minor_determinant << '\n';
    std::cout << name << "_pivot_to_natural_column_sign=" << sign << '\n';
    std::cout << name << "_minor_determinant_natural_columns_mod_" << prime << '='
              << sorted_determinant << '\n';
    std::cout << name << "_selected_rows="
              << compress_ordered_ranges(record.selected_rows) << '\n';
    std::cout << name << "_selected_columns_natural_order="
              << compress_ordered_ranges(sorted_columns) << '\n';
  }
}

void combinatorial_shape_controls() {
  const std::vector<int> degrees{3, 6, 9, 12, 15, 18, 21};
  int rows5 = 0;
  int rows4 = 0;
  for (int degree : degrees) {
    rows5 += choose(21 - degree + 4, 4);
    rows4 += choose(21 - degree + 3, 3);
  }
  require(rows5 == 13972, "five-variable degree-21 row count");
  require(choose(25, 4) == 12650, "five-variable degree-21 column count");
  require(rows4 == 2926, "four-variable degree-21 row count");
  require(choose(24, 3) == 2024, "four-variable degree-21 column count");

  auto formal_coefficient = [](const std::vector<int>& generator_degrees,
                               int target_degree) {
    std::int64_t answer = 0;
    const int count = static_cast<int>(generator_degrees.size());
    for (int mask = 0; mask < (1 << count); ++mask) {
      int shift = 0;
      int parity = 0;
      for (int i = 0; i < count; ++i) {
        if ((mask >> i) & 1) {
          shift += generator_degrees[i];
          parity ^= 1;
        }
      }
      if (shift <= target_degree) {
        const std::int64_t term = choose(target_degree - shift + 4, 4);
        answer += parity ? -term : term;
      }
    }
    return answer;
  };
  require(formal_coefficient({3, 6, 9, 12, 15}, 21) == 1705,
          "first-five complete-intersection coefficient");
  require(formal_coefficient(degrees, 21) == 1670,
          "seven-form formal degree-21 coefficient");
  require(formal_coefficient(degrees, 28) == 39,
          "seven-form formal degree-28 coefficient");
  require(formal_coefficient(degrees, 29) == -354,
          "seven-form formal degree-29 coefficient");
}

void denominator_hostile_control() {
  bool rejected = false;
  try {
    ModularMoments hostile(83);
    (void)hostile;
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  require(rejected, "p<=86 hostile denominator guard must reject");
}

}  // namespace

int main(int argc, char** argv) {
  try {
    combinatorial_shape_controls();
    denominator_hostile_control();
    const std::vector<int> all_degrees{3, 6, 9, 12, 15, 18, 21};
    const std::array<int, 5> all_variables{{0, 1, 2, 3, 4}};
    const std::array<int, 5> delete_e{{0, 1, 2, 3, 0}};

    std::cout << "DEGREE-FOUR CYCLIC EIGENSPACE: SUPPORT-FIVE SPARSE RANK\n";
    std::cout << "status=FINITE-EXACT deterministic sparse row reduction\n";
    std::cout << "hostile_prime_83_guard=REJECTED\n";
    std::cout << "ordering=fixed-degree descending lex on A,B,C,D,E\n";
    std::cout << "storage=dense reusable accumulator + sparse normalized pivots\n\n";

    std::vector<int> primes{101, 103};
    if (argc == 2) primes = {std::stoi(argv[1])};
    for (int prime : primes) {
      require(prime > 86 && prime < 256 && is_prime(prime),
              "argument must be a prime 86 < p < 256");
      std::cout << "PRIME " << prime << " guard=" << prime << ">86: PASS\n";

      // Independent control against THM-3321's already frozen delete-E rank.
      std::cout << "SUPPORT-FOUR CONTROL (delete E)\n";
      RankRecord control =
          rank_moment_system(prime, 4, delete_e, all_degrees, 21, false);
      print_record("delete_E", prime, control, false);
      require(control.rows == 2926 && control.columns == 2024 &&
                  control.rank == 2024,
              "support-four full-rank control");

      // Omit M15,M18,M21: a known-deficient hostile projective system.
      std::cout << "UNDER-GENERATED HOSTILE CONTROL (M3,M6,M9,M12 only)\n";
      RankRecord hostile = rank_moment_system(
          prime, 4, delete_e, std::vector<int>{3, 6, 9, 12}, 21, false);
      print_record("delete_E_short", prime, hostile, false);
      require(hostile.rank < hostile.columns, "hostile control must be deficient");

      std::cout << "SUPPORT-FIVE DEGREE-21 MAP\n";
      RankRecord full =
          rank_moment_system(prime, 5, all_variables, all_degrees, 21, true);
      // The first guarded prime carries the explicit maximal-minor row/column
      // certificate.  The second prime is an independent exact rank control
      // without duplicating that long deterministic index list.
      print_record("support_five", prime, full, prime == 101);
      require(full.rows == 13972 && full.columns == 12650,
              "support-five matrix shape");
      require(static_cast<int>(full.selected_rows.size()) == full.rank &&
                  static_cast<int>(full.selected_columns.size()) == full.rank,
              "minor certificate dimensions equal rank");
      require(std::is_sorted(full.selected_rows.begin(), full.selected_rows.end()) &&
                  std::adjacent_find(full.selected_rows.begin(),
                                     full.selected_rows.end()) ==
                      full.selected_rows.end(),
              "selected source rows are strictly increasing");
      {
        std::vector<int> sorted_columns = full.selected_columns;
        std::sort(sorted_columns.begin(), sorted_columns.end());
        require(std::adjacent_find(sorted_columns.begin(), sorted_columns.end()) ==
                    sorted_columns.end(),
                "selected pivot columns are distinct");
      }
      require(full.selected_minor_determinant != 0,
              "selected maximal minor determinant is nonzero");
      require(full.rank <= 10980,
              "complete-intersection ceiling rank<=10980 at degree 21");
      std::cout << "first_five_complete_intersection_lower_bound=1705\n";
      std::cout << "M18_max_new_rank=34 (35 cubic multipliers minus M3 overlap)\n";
      std::cout << "M21_max_new_rank=1\n";
      std::cout << "universal_degree21_nullity_lower_bound=1670\n";
      std::cout << "universal_degree21_rank_ceiling=10980\n";
      std::cout << "ceiling_attained=" << (full.rank == 10980 ? "YES" : "NO")
                << '\n';
      std::cout << "degree21_projective_exclusion=IMPOSSIBLE_FROM_THIS_MAP\n";
      std::cout << "formal_product_coefficient_degree28=39\n";
      std::cout << "formal_product_coefficient_degree29=-354\n";
      std::cout << "degree29=FIRST_NOT_RULED_OUT_BY_FORMAL_COUNT_ONLY\n";
      if (prime != primes.back()) std::cout << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
