// Dependency-free exact referee for THM-2514.
//
// The program works only with finite integer matrices and finite fields.  It
// checks the critical 13 x 7 -> K_14 chart, ordinary-degree reconstruction
// from any six cut phases, the exact affine covariance generators, the
// Fourier multiplier (including the alpha=0 Omega channel), and the sharp
// two-marginal loss boundary.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr int P = 13;
constexpr int Q = 7;
constexpr int OMEGA = 13;
constexpr int NV = 14;
constexpr int DIM = P * (Q - 1);
constexpr int MOD = 547;  // 547-1=2*3*7*13.

void require(bool condition, const string& message) {
  if (!condition) throw runtime_error(message);
}

int modn(long long x, int n) {
  x %= n;
  if (x < 0) x += n;
  return static_cast<int>(x);
}

int mod13(long long x) { return modn(x, P); }
int mod7(long long x) { return modn(x, Q); }
int modq(long long x) { return modn(x, MOD); }

int power_mod(int a, int e, int modulus = MOD) {
  long long answer = 1;
  long long base = modn(a, modulus);
  while (e > 0) {
    if (e & 1) answer = answer * base % modulus;
    base = base * base % modulus;
    e >>= 1;
  }
  return static_cast<int>(answer);
}

int inverse_mod(int a, int modulus) {
  require(modn(a, modulus) != 0, "attempted to invert zero");
  return power_mod(modn(a, modulus), modulus - 2, modulus);
}

int primitive_root_547() {
  const array<int, 4> prime_divisors{2, 3, 7, 13};
  for (int g = 2; g < MOD; ++g) {
    bool primitive = true;
    for (int ell : prime_divisors) {
      if (power_mod(g, (MOD - 1) / ell) == 1) primitive = false;
    }
    if (primitive) return g;
  }
  throw runtime_error("no primitive root modulo 547");
}

struct Edge {
  int u = 0;
  int v = 0;
  bool operator<(const Edge& other) const {
    return pair<int, int>{u, v} < pair<int, int>{other.u, other.v};
  }
  bool operator==(const Edge& other) const {
    return u == other.u && v == other.v;
  }
};

Edge make_edge(int u, int v) {
  if (u > v) swap(u, v);
  require(u != v, "K14 edge acquired a loop");
  return {u, v};
}

Edge chart_edge(int h, int r, int tau, int a, int c) {
  const int s = mod7(a * r + c);
  if (s == 0) return make_edge(OMEGA, h);
  return make_edge(mod13(h + tau * s), mod13(h - tau * s));
}

pair<int, int> ordered_pair(int h, int r, int tau, int a, int c) {
  const int s = mod7(a * r + c);
  return {mod13(h + tau * s), mod13(h - tau * s)};
}

int incidence(const Edge& edge, int vertex) {
  return static_cast<int>(edge.u == vertex || edge.v == vertex);
}

pair<int, int> basis_cell(int column) {
  return {column / (Q - 1), column % (Q - 1)};
}

int basis_value(int column, int h, int r) {
  const auto [h0, r0] = basis_cell(column);
  if (h != h0) return 0;
  if (r == r0) return 1;
  if (r == Q - 1) return -1;
  return 0;
}

array<int, NV> degree_basis(int column, int tau, int a, int c) {
  array<int, NV> answer{};
  const auto [h, r] = basis_cell(column);
  const Edge positive = chart_edge(h, r, tau, a, c);
  const Edge negative = chart_edge(h, Q - 1, tau, a, c);
  ++answer[positive.u];
  ++answer[positive.v];
  --answer[negative.u];
  --answer[negative.v];
  return answer;
}

array<int, P> radon_basis(int column, int tau, int a, int c) {
  array<int, P> answer{};
  const auto [h, r] = basis_cell(column);
  const int s_pos = mod7(a * r + c);
  const int s_neg = mod7(a * (Q - 1) + c);
  ++answer[mod13(h + tau * s_pos)];
  --answer[mod13(h + tau * s_neg)];
  return answer;
}

int rank_mod(vector<vector<int>> matrix) {
  if (matrix.empty()) return 0;
  const int rows = static_cast<int>(matrix.size());
  const int columns = static_cast<int>(matrix.front().size());
  for (auto& row : matrix) {
    require(static_cast<int>(row.size()) == columns, "ragged matrix");
    for (int& value : row) value = modq(value);
  }
  int rank = 0;
  for (int column = 0; column < columns && rank < rows; ++column) {
    int pivot = rank;
    while (pivot < rows && matrix[pivot][column] == 0) ++pivot;
    if (pivot == rows) continue;
    swap(matrix[rank], matrix[pivot]);
    const int inverse = inverse_mod(matrix[rank][column], MOD);
    for (int j = column; j < columns; ++j) {
      matrix[rank][j] = modq(1LL * matrix[rank][j] * inverse);
    }
    for (int i = rank + 1; i < rows; ++i) {
      const int factor = matrix[i][column];
      if (factor == 0) continue;
      for (int j = column; j < columns; ++j) {
        matrix[i][j] =
            modq(matrix[i][j] - 1LL * factor * matrix[rank][j]);
      }
    }
    ++rank;
  }
  return rank;
}

vector<vector<int>> degree_phase_matrix(int tau, int a,
                                        const vector<int>& phases) {
  vector<vector<int>> matrix;
  matrix.reserve(phases.size() * NV);
  for (int c : phases) {
    for (int vertex = 0; vertex < NV; ++vertex) {
      vector<int> row(DIM);
      for (int column = 0; column < DIM; ++column) {
        row[column] = degree_basis(column, tau, a, c)[vertex];
      }
      matrix.push_back(std::move(row));
    }
  }
  return matrix;
}

vector<vector<int>> marginal_matrix(int tau) {
  vector<vector<int>> matrix;
  matrix.reserve(2 * P);
  for (int slope : {tau, mod13(-tau)}) {
    for (int v = 0; v < P; ++v) {
      vector<int> row(DIM);
      for (int column = 0; column < DIM; ++column) {
        row[column] = radon_basis(column, slope, 1, 0)[v];
      }
      matrix.push_back(std::move(row));
    }
  }
  return matrix;
}

vector<vector<int>> symmetric_slope_matrix() {
  vector<vector<int>> matrix;
  matrix.reserve(6 * P);
  for (int tau = 1; tau <= 6; ++tau) {
    for (int v = 0; v < P; ++v) {
      vector<int> row(DIM);
      for (int column = 0; column < DIM; ++column) {
        const auto plus = radon_basis(column, tau, 1, 0);
        const auto minus = radon_basis(column, mod13(-tau), 1, 0);
        row[column] = plus[v] + minus[v];
      }
      matrix.push_back(std::move(row));
    }
  }
  return matrix;
}

using Defect = array<array<int, Q>, P>;

Defect marginal_hostile(int tau) {
  auto T = [](int h) { return mod13(h) == 0 ? 12 : -1; };
  Defect defect{};
  for (int h = 0; h < P; ++h) {
    const int packet = T(h) + T(h + tau) + T(h - tau);
    defect[h][0] = -T(h);
    defect[h][1] = packet;
    defect[h][2] = -packet;
    defect[h][3] = T(h);
  }
  return defect;
}

array<long long, P> radon(const Defect& defect, int tau, int a = 1,
                          int c = 0) {
  array<long long, P> answer{};
  for (int h = 0; h < P; ++h) {
    for (int r = 0; r < Q; ++r) {
      const int s = mod7(a * r + c);
      answer[mod13(h + tau * s)] += defect[h][r];
    }
  }
  return answer;
}

int zeta = 0;
int xi = 0;

int zeta_power(int exponent) { return power_mod(zeta, mod13(exponent)); }
int xi_power(int exponent) { return power_mod(xi, mod7(exponent)); }

int kernel_K(int lambda, int beta) {
  int answer = 0;
  for (int s = 0; s < Q; ++s) {
    answer = modq(answer + 1LL * zeta_power(-lambda * s) *
                              xi_power(-beta * s));
  }
  return answer;
}

int kernel_L(int lambda, int beta) {
  return modq(kernel_K(lambda, beta) + kernel_K(-lambda, beta) - 1);
}

int basis_transform(int column, int alpha, int beta) {
  const auto [h, r] = basis_cell(column);
  return modq(1LL * zeta_power(-alpha * h) *
              modq(xi_power(-beta * r) -
                   xi_power(-beta * (Q - 1))));
}

struct Gauge {
  int A;
  int H;
  int B;
  int C;
};

Edge affine_vertex_image(Edge edge, int A, int H) {
  auto move = [&](int vertex) {
    return vertex == OMEGA ? OMEGA : mod13(A * vertex + H);
  };
  return make_edge(move(edge.u), move(edge.v));
}

bool is_prime(int n) {
  if (n < 2) return false;
  for (int d = 2; d * d <= n; ++d) {
    if (n % d == 0) return false;
  }
  return true;
}

}  // namespace

int main() {
  try {
    const int primitive = primitive_root_547();
    zeta = power_mod(primitive, (MOD - 1) / P);
    xi = power_mod(primitive, (MOD - 1) / Q);
    require(power_mod(zeta, P) == 1 && zeta != 1,
            "bad thirteenth root");
    require(power_mod(xi, Q) == 1 && xi != 1, "bad seventh root");

    set<Edge> complete_edges;
    for (int u = 0; u < NV; ++u) {
      for (int v = u + 1; v < NV; ++v) complete_edges.insert({u, v});
    }
    require(complete_edges.size() == 91, "K14 edge count drifted");

    uint64_t chart_bijections = 0;
    uint64_t matching_rows = 0;
    uint64_t column_types = 0;
    for (int tau = 1; tau < P; ++tau) {
      set<int> differences;
      for (int s = 1; s < Q; ++s) differences.insert(mod13(-2 * tau * s));
      set<int> expected;
      for (int odd : {1, 3, 5, 7, 9, 11}) {
        expected.insert(mod13(tau * odd));
      }
      require(differences == expected, "circulant tournament orbit drifted");

      for (int a = 1; a < Q; ++a) {
        for (int c = 0; c < Q; ++c) {
          set<Edge> image;
          for (int h = 0; h < P; ++h) {
            array<int, NV> degrees{};
            for (int r = 0; r < Q; ++r) {
              const Edge edge = chart_edge(h, r, tau, a, c);
              image.insert(edge);
              ++degrees[edge.u];
              ++degrees[edge.v];
            }
            require(all_of(degrees.begin(), degrees.end(),
                           [](int degree) { return degree == 1; }),
                    "row left the cyclic one-factorization");
            ++matching_rows;
          }
          require(image == complete_edges, "chart is not a K14 bijection");
          ++chart_bijections;

          for (int r = 0; r < Q; ++r) {
            array<int, NV> degrees{};
            vector<vector<int>> adjacency(NV);
            for (int h = 0; h < P; ++h) {
              const Edge edge = chart_edge(h, r, tau, a, c);
              ++degrees[edge.u];
              ++degrees[edge.v];
              adjacency[edge.u].push_back(edge.v);
              adjacency[edge.v].push_back(edge.u);
            }
            const int s = mod7(a * r + c);
            if (s == 0) {
              require(degrees[OMEGA] == 13, "hub column is not a star");
              for (int v = 0; v < P; ++v) {
                require(degrees[v] == 1, "star leaf degree drifted");
              }
            } else {
              require(degrees[OMEGA] == 0, "cycle met Omega");
              for (int v = 0; v < P; ++v) {
                require(degrees[v] == 2, "Hamilton column is not 2-regular");
              }
              array<bool, NV> seen{};
              queue<int> frontier;
              frontier.push(0);
              seen[0] = true;
              while (!frontier.empty()) {
                const int v = frontier.front();
                frontier.pop();
                for (int w : adjacency[v]) {
                  if (!seen[w]) {
                    seen[w] = true;
                    frontier.push(w);
                  }
                }
              }
              for (int v = 0; v < P; ++v) {
                require(seen[v], "2-regular column is disconnected");
              }
            }
            ++column_types;
          }
        }
      }
    }

    const array<Gauge, 4> generators{{
        {1, 1, 1, 0},  // F_13 translation
        {2, 0, 1, 0},  // F_13 multiplier (order 12)
        {1, 0, 1, 1},  // F_7 translation
        {1, 0, 3, 0},  // F_7 multiplier (order 6)
    }};
    uint64_t covariance_rows = 0;
    for (const Gauge& gauge : generators) {
      const int Ainv = inverse_mod(gauge.A, P);
      for (int tau = 1; tau < P; ++tau) {
        for (int a = 1; a < Q; ++a) {
          for (int c = 0; c < Q; ++c) {
            const int old_tau = mod13(Ainv * tau);
            const int old_a = mod7(a * gauge.B);
            const int old_c = mod7(a * gauge.C + c);
            for (int h0 = 0; h0 < P; ++h0) {
              for (int r0 = 0; r0 < Q; ++r0) {
                const int h = mod13(gauge.A * h0 + gauge.H);
                const int r = mod7(gauge.B * r0 + gauge.C);
                const Edge old_edge =
                    chart_edge(h0, r0, old_tau, old_a, old_c);
                const Edge transported =
                    affine_vertex_image(old_edge, gauge.A, gauge.H);
                const Edge current = chart_edge(h, r, tau, a, c);
                require(transported == current,
                        "affine edge covariance failed");
                ++covariance_rows;
              }
            }
          }
        }
      }
    }

    uint64_t degree_formula_rows = 0;
    uint64_t phase_sum_rows = 0;
    for (int tau = 1; tau < P; ++tau) {
      for (int a = 1; a < Q; ++a) {
        const int ainv = inverse_mod(a, Q);
        for (int c = 0; c < Q; ++c) {
          const int rc = mod7(-ainv * c);
          for (int column = 0; column < DIM; ++column) {
            const auto degree = degree_basis(column, tau, a, c);
            const auto plus = radon_basis(column, tau, a, c);
            const auto minus = radon_basis(column, mod13(-tau), a, c);
            for (int x = 0; x < P; ++x) {
              const int expected_degree =
                  plus[x] + minus[x] - basis_value(column, x, rc);
              require(degree[x] == expected_degree,
                      "ordinary-degree/Radon identity failed");
              ++degree_formula_rows;
            }
            int omega_expected = 0;
            for (int h = 0; h < P; ++h) {
              omega_expected += basis_value(column, h, rc);
            }
            require(degree[OMEGA] == omega_expected,
                    "Omega column-total identity failed");
            ++degree_formula_rows;
          }
        }
        for (int column = 0; column < DIM; ++column) {
          for (int vertex = 0; vertex < NV; ++vertex) {
            int total = 0;
            for (int c = 0; c < Q; ++c) {
              total += degree_basis(column, tau, a, c)[vertex];
            }
            require(total == 0, "seven cut degrees do not sum to zero");
            ++phase_sum_rows;
          }
        }
      }
    }

    int minimum_six_phase_rank = DIM;
    int maximum_six_phase_rank = 0;
    int minimum_single_phase_rank = DIM;
    int maximum_single_phase_rank = 0;
    uint64_t six_phase_banks = 0;
    for (int tau = 1; tau < P; ++tau) {
      for (int a = 1; a < Q; ++a) {
        for (int c = 0; c < Q; ++c) {
          const int single_rank = rank_mod(degree_phase_matrix(tau, a, {c}));
          minimum_single_phase_rank = min(minimum_single_phase_rank, single_rank);
          maximum_single_phase_rank = max(maximum_single_phase_rank, single_rank);
        }
        for (int omitted = 0; omitted < Q; ++omitted) {
          vector<int> phases;
          for (int c = 0; c < Q; ++c) {
            if (c != omitted) phases.push_back(c);
          }
          const int rank = rank_mod(degree_phase_matrix(tau, a, phases));
          minimum_six_phase_rank = min(minimum_six_phase_rank, rank);
          maximum_six_phase_rank = max(maximum_six_phase_rank, rank);
          require(rank == DIM, "six cut phases failed to reconstruct");
          ++six_phase_banks;
        }
      }
    }
    require(minimum_single_phase_rank == 13 &&
                maximum_single_phase_rank == 13,
            "single-phase rank drifted");

    uint64_t multiplier_checks = 0;
    for (int lambda = 1; lambda < P; ++lambda) {
      require(kernel_L(lambda, 0) == 0,
              "zero cut character multiplier should vanish");
      ++multiplier_checks;
      for (int beta = 1; beta < Q; ++beta) {
        require(kernel_L(lambda, beta) != 0,
                "primitive ordinary-degree multiplier vanished");
        ++multiplier_checks;
      }
    }

    uint64_t fourier_checks = 0;
    uint64_t omega_checks = 0;
    for (int tau = 1; tau < P; ++tau) {
      for (int a = 1; a < Q; ++a) {
        for (int column = 0; column < DIM; ++column) {
          for (int alpha = 1; alpha < P; ++alpha) {
            for (int beta = 1; beta < Q; ++beta) {
              int lhs = 0;
              for (int c = 0; c < Q; ++c) {
                const auto degree = degree_basis(column, tau, a, c);
                int horizontal = 0;
                for (int x = 0; x < P; ++x) {
                  horizontal =
                      modq(horizontal + 1LL * degree[x] *
                                            zeta_power(-alpha * x));
                }
                lhs = modq(lhs + 1LL * horizontal * xi_power(-beta * c));
              }
              const int rhs =
                  modq(1LL * kernel_L(alpha * tau, beta) *
                       basis_transform(column, alpha, -beta * a));
              require(lhs == rhs, "finite-degree Fourier factorization failed");
              ++fourier_checks;
            }
          }
          for (int beta = 1; beta < Q; ++beta) {
            int lhs = 0;
            for (int c = 0; c < Q; ++c) {
              lhs = modq(lhs + 1LL * degree_basis(column, tau, a, c)[OMEGA] *
                                   xi_power(-beta * c));
            }
            const int rhs = basis_transform(column, 0, -beta * a);
            require(lhs == rhs, "Omega failed to recover alpha=0 mode");
            ++omega_checks;
          }
        }
      }
    }

    int marginal_rank_minimum = DIM;
    int marginal_rank_maximum = 0;
    long long hostile_l1 = -1;
    for (int tau = 1; tau < P; ++tau) {
      const int rank = rank_mod(marginal_matrix(tau));
      marginal_rank_minimum = min(marginal_rank_minimum, rank);
      marginal_rank_maximum = max(marginal_rank_maximum, rank);
      require(rank == 24, "antipodal marginal rank drifted");

      const Defect hostile = marginal_hostile(tau);
      long long current_l1 = 0;
      bool nonvertical = false;
      for (int h = 0; h < P; ++h) {
        int row_sum = 0;
        for (int r = 0; r < Q; ++r) {
          row_sum += hostile[h][r];
          current_l1 += llabs(hostile[h][r]);
          if (h > 0 && hostile[h][r] != hostile[0][r]) nonvertical = true;
        }
        require(row_sum == 0, "marginal hostile is not row-zero");
      }
      require(nonvertical, "marginal hostile became vertical");
      const auto plus = radon(hostile, tau);
      const auto minus = radon(hostile, mod13(-tau));
      require(all_of(plus.begin(), plus.end(), [](long long x) { return x == 0; }) &&
                  all_of(minus.begin(), minus.end(), [](long long x) { return x == 0; }),
              "integral marginal hostile did not vanish");
      if (hostile_l1 < 0) hostile_l1 = current_l1;
      require(current_l1 == hostile_l1, "hostile L1 depends on slope");
    }
    const int symmetric_rank = rank_mod(symmetric_slope_matrix());
    require(symmetric_rank == 72, "six converse-invariant sums rank drifted");

    uint64_t general_boundary_rows = 0;
    uint64_t critical_bijections = 0;
    for (int p = 3; p <= 31; p += 2) {
      if (!is_prime(p)) continue;
      for (int q = 2; q < p; ++q) {
        set<pair<int, int>> unordered;
        for (int h = 0; h < p; ++h) {
          for (int r = 0; r < q; ++r) {
            int x = modn(h + r, p);
            int y = modn(h - r, p);
            if (x > y) swap(x, y);
            unordered.insert({x, y});
          }
        }
        const bool injective = static_cast<int>(unordered.size()) == p * q;
        require(injective == (q <= (p + 1) / 2),
                "general antipodal injectivity boundary failed");
        if (q == (p + 1) / 2) {
          require(static_cast<int>(unordered.size()) == p * (p + 1) / 2,
                  "critical chart is not onto Sym^2");
          ++critical_bijections;
        }

        const int antipodal_pairs = (p - 1) / 2;
        const int bad_count = q - 2;
        set<int> bad;
        for (int k = 1; k <= antipodal_pairs &&
                        static_cast<int>(bad.size()) < bad_count;
             ++k) {
          bad.insert(k);
        }
        for (int k = 1; k <= antipodal_pairs &&
                        static_cast<int>(bad.size()) < bad_count;
             ++k) {
          bad.insert(p - k);
        }
        int full_good_pairs = 0;
        for (int k = 1; k <= antipodal_pairs; ++k) {
          if (!bad.count(k) && !bad.count(p - k)) ++full_good_pairs;
        }
        const int expected_full = max(0, (p - 2 * q + 3) / 2);
        require(full_good_pairs == expected_full,
                "sharp antipodal good-pair mask failed");
        ++general_boundary_rows;
      }
    }

    cout << "THM-2514 CYCLIC K14 ORDINARY-DEGREE RECONSTRUCTION REFEREE\n";
    cout << "field=F_547;primitive_root=" << primitive << ";zeta13=" << zeta
         << ";xi7=" << xi << "\n";
    cout << "chart_bijections=" << chart_bijections
         << ";matching_rows=" << matching_rows
         << ";column_type_checks=" << column_types << "\n";
    cout << "carrier=E(K14)=91;rows=13_perfect_matchings;columns=K1,13_plus_6_C13;"
            "oriented_connection=tau*{1,3,5,7,9,11}\n";
    cout << "affine_generator_covariance_rows=" << covariance_rows << "\n";
    cout << "degree_formula_rows=" << degree_formula_rows
         << ";seven_phase_sum_zero_rows=" << phase_sum_rows << "\n";
    cout << "single_phase_rank=" << minimum_single_phase_rank << ".."
         << maximum_single_phase_rank << ";six_phase_banks=" << six_phase_banks
         << ";six_phase_rank=" << minimum_six_phase_rank << ".."
         << maximum_six_phase_rank << ";domain_dim=78\n";
    cout << "multiplier_checks=" << multiplier_checks
         << ";primitive_L_nonzero=72;beta0_L_zero=12\n";
    cout << "finite_fourier_factorization_checks=" << fourier_checks
         << ";Omega_alpha0_recovery_checks=" << omega_checks << "\n";
    cout << "antipodal_marginal_rank=" << marginal_rank_minimum << ".."
         << marginal_rank_maximum << ";kernel_dim=" << DIM - marginal_rank_minimum
         << ";integral_hostiles=12;hostile_L1=" << hostile_l1 << "\n";
    cout << "six_converse_invariant_pair_sums_rank=" << symmetric_rank
         << ";vertical_kernel_dim=" << DIM - symmetric_rank << "\n";
    cout << "general_p_le_31_boundary_rows=" << general_boundary_rows
         << ";critical_Sym2_bijections=" << critical_bijections
         << ";injective_iff_q<=(p+1)/2;min_full_good_pairs=max(0,(p-2q+3)/2)\n";
    cout << "sharp_carry_boundary=star_degree_(13,1^13)_vs_cycle_degree_(2^13,0);"
            "no_vertex_automorphism\n";
    cout << "NONCLAIMS=signed_static_root_chart_not_runner_labels_not_owner_time_deep_ancestry;"
            "no_live_row_exclusion\n";
    cout << "ALL EXACT CHECKS PASSED\n";
  } catch (const exception& error) {
    cerr << "FAIL: " << error.what() << "\n";
    return 1;
  }
}
