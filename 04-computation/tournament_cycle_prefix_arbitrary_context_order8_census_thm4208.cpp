#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using i128 = __int128_t;

struct Row {
  int order;
  std::string label;
  long long h;
  std::array<long long, 11> left;
  std::array<long long, 11> right;
};

static std::string show(i128 value) {
  if (value == 0) return "0";
  const bool negative = value < 0;
  if (negative) value = -value;
  std::string answer;
  while (value != 0) {
    answer.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  if (negative) answer.push_back('-');
  std::reverse(answer.begin(), answer.end());
  return answer;
}

static i128 gcd128(i128 left, i128 right) {
  if (left < 0) left = -left;
  if (right < 0) right = -right;
  while (right != 0) {
    const i128 next = left % right;
    left = right;
    right = next;
  }
  return left;
}

static bool label_less(const Row &left_b, const Row &left_c,
                       const Row &right_b, const Row &right_c) {
  return std::pair(left_b.label, left_c.label) <
         std::pair(right_b.label, right_c.label);
}

struct Statistics {
  long long total = 0;
  long long negative = 0;
  long long failures = 0;
  long long equality = 0;
  bool initialized = false;
  i128 minimum_slack = 0;
  int minimum_slack_i = -1;
  int minimum_slack_j = -1;
  i128 minimum_numerator = 0;
  i128 minimum_denominator = 1;
  int minimum_ratio_i = -1;
  int minimum_ratio_j = -1;
  bool positive_initialized = false;
  i128 minimum_positive_slack = 0;
  int minimum_positive_slack_i = -1;
  int minimum_positive_slack_j = -1;
  i128 minimum_positive_numerator = 0;
  i128 minimum_positive_denominator = 1;
  int minimum_positive_ratio_i = -1;
  int minimum_positive_ratio_j = -1;
};

static void update(Statistics &stats, const std::vector<Row> &rows, int i,
                   int j, i128 value, i128 scale) {
  const i128 slack = value - i128(10764) * scale;
  ++stats.total;
  stats.negative += value < 0;
  stats.failures += slack < 0;
  stats.equality += slack == 0;
  if (!stats.initialized || slack < stats.minimum_slack ||
      (slack == stats.minimum_slack &&
       label_less(rows[i], rows[j], rows[stats.minimum_slack_i],
                  rows[stats.minimum_slack_j]))) {
    stats.minimum_slack = slack;
    stats.minimum_slack_i = i;
    stats.minimum_slack_j = j;
  }
  if (!stats.initialized ||
      value * stats.minimum_denominator < stats.minimum_numerator * scale ||
      (value * stats.minimum_denominator == stats.minimum_numerator * scale &&
       label_less(rows[i], rows[j], rows[stats.minimum_ratio_i],
                  rows[stats.minimum_ratio_j]))) {
    stats.minimum_numerator = value;
    stats.minimum_denominator = scale;
    stats.minimum_ratio_i = i;
    stats.minimum_ratio_j = j;
  }
  stats.initialized = true;
  if (slack <= 0) return;
  if (!stats.positive_initialized || slack < stats.minimum_positive_slack ||
      (slack == stats.minimum_positive_slack &&
       label_less(rows[i], rows[j], rows[stats.minimum_positive_slack_i],
                  rows[stats.minimum_positive_slack_j]))) {
    stats.minimum_positive_slack = slack;
    stats.minimum_positive_slack_i = i;
    stats.minimum_positive_slack_j = j;
  }
  if (!stats.positive_initialized ||
      value * stats.minimum_positive_denominator <
          stats.minimum_positive_numerator * scale ||
      (value * stats.minimum_positive_denominator ==
           stats.minimum_positive_numerator * scale &&
       label_less(rows[i], rows[j], rows[stats.minimum_positive_ratio_i],
                  rows[stats.minimum_positive_ratio_j]))) {
    stats.minimum_positive_numerator = value;
    stats.minimum_positive_denominator = scale;
    stats.minimum_positive_ratio_i = i;
    stats.minimum_positive_ratio_j = j;
  }
  stats.positive_initialized = true;
}

static void print_ratio(const std::string &name, i128 numerator,
                        i128 denominator, const Row &left, const Row &right) {
  const i128 divisor = gcd128(numerator, denominator);
  std::cout << name << ' ' << show(numerator / divisor) << " denominator "
            << show(denominator / divisor) << " B " << left.label << " C "
            << right.label << " H " << left.h << ' ' << right.h << '\n';
}

static void print_statistics(const std::string &prefix, const Statistics &stats,
                             const std::vector<Row> &rows) {
  std::cout << prefix << "ordered_pairs " << stats.total << '\n';
  std::cout << prefix << "negative " << stats.negative << '\n';
  std::cout << prefix << "bound_failures " << stats.failures << '\n';
  std::cout << prefix << "bound_equality " << stats.equality << '\n';
  std::cout << prefix << "minimum_slack " << show(stats.minimum_slack)
            << " B " << rows[stats.minimum_slack_i].label << " C "
            << rows[stats.minimum_slack_j].label << '\n';
  print_ratio(prefix + "minimum_ratio_numerator", stats.minimum_numerator,
              stats.minimum_denominator, rows[stats.minimum_ratio_i],
              rows[stats.minimum_ratio_j]);
  if (stats.positive_initialized) {
    std::cout << prefix << "minimum_positive_slack "
              << show(stats.minimum_positive_slack) << " B "
              << rows[stats.minimum_positive_slack_i].label << " C "
              << rows[stats.minimum_positive_slack_j].label << '\n';
    print_ratio(prefix + "minimum_positive_ratio_numerator",
                stats.minimum_positive_numerator,
                stats.minimum_positive_denominator,
                rows[stats.minimum_positive_ratio_i],
                rows[stats.minimum_positive_ratio_j]);
  }
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  int count;
  if (!(std::cin >> count)) return 2;
  std::vector<Row> rows(count);
  for (int k = 0; k < count; ++k) {
    int index;
    std::cin >> index >> rows[k].order >> rows[k].label >> rows[k].h;
    if (!std::cin || index != k) return 3;
    for (auto &value : rows[k].left) std::cin >> value;
    for (auto &value : rows[k].right) std::cin >> value;
    if (!std::cin) return 4;
  }

  std::array<long long, 11> left_max{};
  std::array<long long, 11> right_max{};
  long long h_max = 0;
  for (const Row &row : rows) {
    h_max = std::max(h_max, row.h);
    for (int k = 0; k < 11; ++k) {
      left_max[k] = std::max(left_max[k], std::llabs(row.left[k]));
      right_max[k] = std::max(right_max[k], std::llabs(row.right[k]));
    }
  }
  i128 dot_bound = 0;
  for (int k = 0; k < 11; ++k) {
    dot_bound += i128(left_max[k]) * right_max[k];
  }
  const i128 h4_bound = i128(h_max) * h_max * h_max * h_max;
  const i128 comparison_bound = dot_bound * h4_bound;
  if (comparison_bound >= (i128(1) << 126)) return 5;

  Statistics all;
  Statistics order8;
  for (int i = 0; i < count; ++i) {
    const i128 left_h2 = i128(rows[i].h) * rows[i].h;
    for (int j = 0; j < count; ++j) {
      i128 value = 0;
      for (int k = 0; k < 11; ++k) {
        value += i128(rows[i].left[k]) * rows[j].right[k];
      }
      const i128 scale = left_h2 * rows[j].h * rows[j].h;
      update(all, rows, i, j, value, scale);
      if (rows[i].order == 8 && rows[j].order == 8) {
        update(order8, rows, i, j, value, scale);
      }
    }
  }

  std::cout << "TOURNAMENT_CYCLE_PREFIX_ORDER8_CENSUS_EXACT_128\n";
  std::cout << "classes " << count << '\n';
  std::cout << "max_H " << h_max << '\n';
  std::cout << "absolute_dot_bound " << show(dot_bound) << '\n';
  std::cout << "ratio_cross_product_bound " << show(comparison_bound) << '\n';
  print_statistics("", all, rows);
  print_statistics("order8_sector_", order8, rows);
  std::cout << "PASS\n";
  return all.negative == 0 && all.failures == 0 && all.equality == 1 &&
                 order8.negative == 0 && order8.failures == 0
             ? 0
             : 1;
}
