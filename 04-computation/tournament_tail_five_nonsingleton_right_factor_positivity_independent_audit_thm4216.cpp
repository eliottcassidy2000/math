#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Clean-room literal referee for THM-4216.
//
// This program imports no tournament transfer or response implementation.  It
// rebuilds every audited tournament from labelled adjacency, counts directed
// paths by subset DP, constructs exposed-word capacities and rooted U/V states,
// and evaluates G_+ and R_+ literally.  The only theorem-level formula encoded
// after that reconstruction is THM-4208's displayed 11-jet response identity.

namespace {

using Count = std::int64_t;
using Out = std::vector<std::uint64_t>;

void need(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

struct Data {
  Out out;
  Count hamilton{};
  Count mass{};
  Count gate{};
  std::vector<std::vector<Count>> capacity;
  std::vector<std::array<Count, 2>> starts;
  std::vector<std::array<Count, 2>> ends;
  std::vector<Count> hamilton_starts;
};

std::size_t index(std::uint64_t mask, std::size_t vertex,
                  std::size_t order) {
  return static_cast<std::size_t>(mask) * order + vertex;
}

int bit_count(std::uint64_t value) {
  return __builtin_popcountll(value);
}

Data build_data(const Out& out) {
  const std::size_t order = out.size();
  need(order >= 1 && order < 63, "literal order range");
  const std::uint64_t size = std::uint64_t{1} << order;
  const std::uint64_t full = size - 1;

  std::vector<Count> start(static_cast<std::size_t>(size) * order, 0);
  std::vector<Count> finish(static_cast<std::size_t>(size) * order, 0);
  for (std::size_t vertex = 0; vertex < order; ++vertex) {
    const std::uint64_t singleton = std::uint64_t{1} << vertex;
    start[index(singleton, vertex, order)] = 1;
    finish[index(singleton, vertex, order)] = 1;
  }

  for (std::uint64_t mask = 1; mask < size; ++mask) {
    if ((mask & (mask - 1)) == 0) {
      continue;
    }
    for (std::size_t vertex = 0; vertex < order; ++vertex) {
      const std::uint64_t vertex_bit = std::uint64_t{1} << vertex;
      if ((mask & vertex_bit) == 0) {
        continue;
      }
      const std::uint64_t rest = mask ^ vertex_bit;
      for (std::size_t other = 0; other < order; ++other) {
        const std::uint64_t other_bit = std::uint64_t{1} << other;
        if ((rest & other_bit) == 0) {
          continue;
        }
        if ((out[vertex] & other_bit) != 0) {
          start[index(mask, vertex, order)] +=
              start[index(rest, other, order)];
        }
        if ((out[other] & vertex_bit) != 0) {
          finish[index(mask, vertex, order)] +=
              finish[index(rest, other, order)];
        }
      }
    }
  }

  std::vector<Count> subset_hamilton(static_cast<std::size_t>(size), 0);
  subset_hamilton[0] = 1;
  for (std::uint64_t mask = 1; mask < size; ++mask) {
    for (std::size_t vertex = 0; vertex < order; ++vertex) {
      subset_hamilton[static_cast<std::size_t>(mask)] +=
          finish[index(mask, vertex, order)];
    }
  }

  std::vector<std::vector<Count>> exposed(
      order, std::vector<Count>(order, 0));
  for (std::uint64_t left_mask = 1; left_mask < full; ++left_mask) {
    const std::uint64_t right_mask = full ^ left_mask;
    for (std::size_t left = 0; left < order; ++left) {
      if ((left_mask & (std::uint64_t{1} << left)) == 0) {
        continue;
      }
      const Count left_count = finish[index(left_mask, left, order)];
      if (left_count == 0) {
        continue;
      }
      for (std::size_t right = 0; right < order; ++right) {
        if ((right_mask & (std::uint64_t{1} << right)) != 0) {
          exposed[left][right] +=
              left_count * start[index(right_mask, right, order)];
        }
      }
    }
  }

  std::vector<std::vector<Count>> capacity(
      order, std::vector<Count>(order, 0));
  Count mass = 0;
  for (std::size_t left = 0; left < order; ++left) {
    for (std::size_t right = left + 1; right < order; ++right) {
      const Count value = exposed[left][right] + exposed[right][left];
      capacity[left][right] = value;
      capacity[right][left] = value;
      mass += value;
    }
  }

  std::vector<std::array<Count, 2>> starts(order, {0, 0});
  std::vector<std::array<Count, 2>> ends(order, {0, 0});
  for (std::uint64_t mask = 1; mask < size; ++mask) {
    const std::size_t parity = static_cast<std::size_t>(bit_count(mask) & 1);
    const Count complement_h =
        subset_hamilton[static_cast<std::size_t>(full ^ mask)];
    for (std::size_t vertex = 0; vertex < order; ++vertex) {
      starts[vertex][parity] +=
          start[index(mask, vertex, order)] * complement_h;
      ends[vertex][parity] +=
          finish[index(mask, vertex, order)] * complement_h;
    }
  }

  std::vector<Count> degree(order, 0);
  std::vector<Count> current(order, 0);
  Count edge_squares = 0;
  for (std::size_t left = 0; left < order; ++left) {
    for (std::size_t right = left + 1; right < order; ++right) {
      const Count value = capacity[left][right];
      degree[left] += value;
      degree[right] += value;
      edge_squares += value * value;
      if ((out[left] & (std::uint64_t{1} << right)) != 0) {
        current[left] += value;
        current[right] -= value;
      } else {
        current[left] -= value;
        current[right] += value;
      }
    }
  }
  Count degree_squares = 0;
  Count signed_current = 0;
  for (std::size_t vertex = 0; vertex < order; ++vertex) {
    degree_squares += degree[vertex] * degree[vertex];
    signed_current += degree[vertex] * current[vertex];
  }
  const Count disjoint_numerator =
      mass * mass + edge_squares - degree_squares;
  need(disjoint_numerator % 2 == 0, "integral disjoint-edge gate");
  const Count gate = disjoint_numerator / 2 + 2 * signed_current;

  std::vector<Count> hamilton_starts(order, 0);
  Count hamilton = 0;
  for (std::size_t vertex = 0; vertex < order; ++vertex) {
    hamilton_starts[vertex] = start[index(full, vertex, order)];
    hamilton += finish[index(full, vertex, order)];
  }
  need(hamilton == subset_hamilton[static_cast<std::size_t>(full)],
       "Hamilton total");

  return Data{out, hamilton, mass, gate, capacity, starts, ends,
              hamilton_starts};
}

Out transitive(std::size_t order) {
  Out out(order, 0);
  for (std::size_t left = 0; left < order; ++left) {
    for (std::size_t right = left + 1; right < order; ++right) {
      out[left] |= std::uint64_t{1} << right;
    }
  }
  return out;
}

Out cycle_three() {
  Out out(3, 0);
  out[0] |= std::uint64_t{1} << 1;
  out[1] |= std::uint64_t{1} << 2;
  out[2] |= std::uint64_t{1} << 0;
  return out;
}

Out ordinal_out(const Out& left, const Out& right) {
  need(!left.empty() && !right.empty(), "nonempty ordinal factors");
  const std::size_t left_order = left.size();
  const std::size_t right_order = right.size();
  need(left_order + right_order < 63, "ordinal order range");
  Out out(left_order + right_order, 0);
  const std::uint64_t right_mask =
      ((std::uint64_t{1} << right_order) - 1) << left_order;
  for (std::size_t vertex = 0; vertex < left_order; ++vertex) {
    out[vertex] = left[vertex] | right_mask;
  }
  for (std::size_t vertex = 0; vertex < right_order; ++vertex) {
    out[left_order + vertex] = right[vertex] << left_order;
  }
  return out;
}

Out q_tail(std::size_t tail) {
  if (tail == 0) {
    return cycle_three();
  }
  return ordinal_out(cycle_three(), transitive(tail));
}

Out labelled(std::size_t order, std::uint64_t code) {
  Out out(order, 0);
  std::size_t cursor = 0;
  for (std::size_t left = 0; left < order; ++left) {
    for (std::size_t right = left + 1; right < order; ++right) {
      if ((code & (std::uint64_t{1} << cursor)) != 0) {
        out[left] |= std::uint64_t{1} << right;
      } else {
        out[right] |= std::uint64_t{1} << left;
      }
      ++cursor;
    }
  }
  return out;
}

std::string label(const Out& out) {
  std::string answer;
  for (std::size_t left = 0; left < out.size(); ++left) {
    for (std::size_t right = left + 1; right < out.size(); ++right) {
      answer.push_back(
          (out[left] & (std::uint64_t{1} << right)) != 0 ? '1' : '0');
    }
  }
  return answer;
}

bool has_sink(const Out& out) {
  for (std::uint64_t row : out) {
    if (row == 0) {
      return true;
    }
  }
  return false;
}

struct DirectedMasses {
  std::vector<Count> degree;
  std::vector<Count> outgoing;
  std::vector<Count> incoming;
};

DirectedMasses directed_masses(const Data& data) {
  const std::size_t order = data.out.size();
  DirectedMasses result{std::vector<Count>(order, 0),
                        std::vector<Count>(order, 0),
                        std::vector<Count>(order, 0)};
  for (std::size_t left = 0; left < order; ++left) {
    for (std::size_t right = left + 1; right < order; ++right) {
      const Count value = data.capacity[left][right];
      result.degree[left] += value;
      result.degree[right] += value;
      if ((data.out[left] & (std::uint64_t{1} << right)) != 0) {
        result.outgoing[left] += value;
        result.incoming[right] += value;
      } else {
        result.incoming[left] += value;
        result.outgoing[right] += value;
      }
    }
  }
  return result;
}

using Jet = std::array<Count, 11>;

Jet left_jet(const Data& data) {
  const DirectedMasses masses = directed_masses(data);
  need(data.mass % 2 == 0, "even capacity mass");
  const Count s0 = data.mass / 2;
  const Count s1 = s0 + data.hamilton;
  Count q00 = 0;
  Count q01 = 0;
  Count q11 = 0;
  Count l0 = 0;
  Count l1 = 0;
  for (std::size_t vertex = 0; vertex < data.out.size(); ++vertex) {
    const Count linear = data.mass - masses.degree[vertex] +
                         4 * masses.outgoing[vertex];
    q00 += data.starts[vertex][0] * data.starts[vertex][0];
    q01 += data.starts[vertex][0] * data.starts[vertex][1];
    q11 += data.starts[vertex][1] * data.starts[vertex][1];
    l0 += data.starts[vertex][0] * linear;
    l1 += data.starts[vertex][1] * linear;
  }
  return Jet{
      data.hamilton * data.mass,
      2 * l0,
      2 * l1,
      2 * data.hamilton * s0,
      2 * data.hamilton * s1,
      2 * s0 * s0 + 6 * q00,
      4 * s0 * s1 + 12 * q01,
      2 * s1 * s1 + 6 * q11,
      2 * q00 - 10 * s0 * s0,
      4 * q01 - 20 * s0 * s1,
      2 * q11 - 10 * s1 * s1,
  };
}

Jet right_jet(const Data& data) {
  const DirectedMasses masses = directed_masses(data);
  need(data.mass % 2 == 0, "even capacity mass");
  const Count s0 = data.mass / 2;
  const Count s1 = s0 + data.hamilton;
  Count q00 = 0;
  Count q01 = 0;
  Count q11 = 0;
  Count l0 = 0;
  Count l1 = 0;
  for (std::size_t vertex = 0; vertex < data.out.size(); ++vertex) {
    const Count linear = data.mass - masses.degree[vertex] -
                         4 * masses.incoming[vertex];
    q00 += data.ends[vertex][0] * data.ends[vertex][0];
    q01 += data.ends[vertex][0] * data.ends[vertex][1];
    q11 += data.ends[vertex][1] * data.ends[vertex][1];
    l0 += data.ends[vertex][0] * linear;
    l1 += data.ends[vertex][1] * linear;
  }
  return Jet{
      data.hamilton * data.mass,
      data.hamilton * s0,
      data.hamilton * s1,
      l0,
      l1,
      s0 * s0,
      s0 * s1,
      s1 * s1,
      q00,
      q01,
      q11,
  };
}

Count dot(const Jet& left, const Jet& right) {
  Count answer = 0;
  for (std::size_t index_value = 0; index_value < left.size(); ++index_value) {
    answer += left[index_value] * right[index_value];
  }
  return answer;
}

Count literal_remainder(const Data& left, const Data& right) {
  const Data child = build_data(ordinal_out(left.out, right.out));
  return child.gate - right.hamilton * right.hamilton * left.gate -
         left.hamilton * left.hamilton * right.gate;
}

struct Decomposition {
  Count exact{};
  Count floor{};
  Count debt{};
  Count tau{};
};

Decomposition decompose(const Data& data) {
  const DirectedMasses masses = directed_masses(data);
  const Count h = data.hamilton;
  const Count w = data.mass;
  Count moment = 0;
  Count mixed = 0;
  Count tau = 0;
  Count avoidance = 0;
  for (std::size_t vertex = 0; vertex < data.out.size(); ++vertex) {
    const Count v = data.ends[vertex][0] + data.ends[vertex][1];
    const Count t = data.ends[vertex][1] - data.ends[vertex][0];
    need(t >= 0, "nonnegative Hamilton-start coordinate");
    need(t == data.hamilton_starts[vertex],
         "V chirality equals literal Hamilton starts");
    need(v - t == masses.incoming[vertex], "incoming fan identity");
    need(w - masses.degree[vertex] >= 0, "edge-avoidance mass");
    moment += v * v;
    mixed += v * t;
    tau += t * t;
    avoidance += (1143 * v + 9 * t) * (w - masses.degree[vertex]);
  }
  need(tau > 0, "strict endpoint-square boundary");
  const Count endpoint_defect =
      (w + h) * (w + 3 * h) - 3 * moment;
  const Count mixed_defect = (w + h) * h - mixed;
  need(endpoint_defect >= 0, "endpoint energy");
  need(mixed_defect >= 0, "mixed endpoint bound");
  const Count exact =
      353088 * h * h + 472302 * h * w + 118656 * w * w + avoidance -
      352116 * moment - 1170 * mixed + 18 * tau;
  const Count floor = 6 * (-33 * h * h + 274 * h * w + 214 * w * w);
  const Count debt = avoidance + 117372 * endpoint_defect +
                     1170 * mixed_defect + 18 * tau;
  need(exact - floor == debt, "exact four-debt certificate");
  need(debt > 0, "strict floor gap from tau");
  return Decomposition{exact, floor, debt, tau};
}

void fnv_update(std::uint64_t& digest, const std::string& text) {
  constexpr std::uint64_t prime = 1099511628211ULL;
  for (unsigned char byte : text) {
    digest ^= static_cast<std::uint64_t>(byte);
    digest *= prime;
  }
}

std::string join(const std::vector<Count>& values) {
  std::ostringstream stream;
  for (std::size_t position = 0; position < values.size(); ++position) {
    if (position != 0) {
      stream << ',';
    }
    stream << values[position];
  }
  return stream.str();
}

std::string join(const Jet& values) {
  return join(std::vector<Count>(values.begin(), values.end()));
}

}  // namespace

int main() {
  try {
    const Data p1 = build_data(transitive(1));
    const Data p2 = build_data(transitive(2));
    const Data q5 = build_data(q_tail(5));
    const Jet q5_lambda = left_jet(q5);
    const Jet expected_lambda = {
        1134, 232128, 233244, 1134, 1152, 117504,
        237276, 119844, -341856, -695052, -353268,
    };
    need(q5.hamilton == 3 && q5.mass == 378, "literal Q5 H/W");
    need(q5_lambda == expected_lambda, "literal Q5 response jet");

    need(q5_lambda[1] % 2 == 0 && q5_lambda[2] % 2 == 0 &&
             q5_lambda[6] % 2 == 0,
         "integral unary hW coefficients");
    need((q5_lambda[5] + q5_lambda[6] + q5_lambda[7]) % 4 == 0,
         "integral unary W-square coefficient");
    const Count coefficient_h2 = q5_lambda[2] + q5_lambda[7];
    const Count coefficient_hw =
        q5_lambda[0] + q5_lambda[1] / 2 + q5_lambda[2] / 2 +
        q5_lambda[6] / 2 + q5_lambda[7];
    const Count coefficient_w2 =
        (q5_lambda[5] + q5_lambda[6] + q5_lambda[7]) / 4;
    need(coefficient_h2 == 353088 && coefficient_hw == 472302 &&
             coefficient_w2 == 118656,
         "exact Q5 scalar coefficients");

    const Count q00_coefficient = q5_lambda[8];
    const Count q01_coefficient = q5_lambda[9];
    const Count q11_coefficient = q5_lambda[10];
    need((q00_coefficient + q01_coefficient + q11_coefficient) % 4 == 0,
         "quadratic moment divisibility");
    const Count moment_coefficient =
        (q00_coefficient + q01_coefficient + q11_coefficient) / 4;
    const Count mixed_coefficient =
        (-2 * q00_coefficient + 2 * q11_coefficient) / 4;
    const Count tau_coefficient =
        (q00_coefficient - q01_coefficient + q11_coefficient) / 4;
    need(moment_coefficient == -347544 && mixed_coefficient == -5706 &&
             tau_coefficient == -18,
         "V-to-endpoint quadratic coordinate change");
    need(moment_coefficient - 4 * 1143 == -352116 &&
             mixed_coefficient + 4 * 1134 == -1170 &&
             tau_coefficient + 4 * 9 == 18,
         "incoming-fan to edge-avoidance rewrite");

    std::vector<Data> contexts;
    for (std::size_t order = 1; order <= 4; ++order) {
      const std::size_t pair_count = order * (order - 1) / 2;
      const std::uint64_t count = std::uint64_t{1} << pair_count;
      for (std::uint64_t code = 0; code < count; ++code) {
        contexts.push_back(build_data(labelled(order, code)));
      }
    }
    need(contexts.size() == 75, "labelled contexts through order four");

    std::uint64_t semantic_digest = 14695981039346656037ULL;
    Count minimum_floor_gap = std::numeric_limits<Count>::max();
    Count minimum_nonsingleton_gap = std::numeric_limits<Count>::max();
    std::size_t nonsingleton_rows = 0;
    std::size_t all_adjacency_rows = 0;
    std::size_t strict_tau_rows = 0;
    for (const Data& context : contexts) {
      const Count direct = literal_remainder(q5, context);
      const Count response = dot(q5_lambda, right_jet(context));
      const Decomposition certificate = decompose(context);
      need(direct == response, "literal response-jet identity");
      need(direct == certificate.exact, "literal Q5 decomposition");
      need(direct > certificate.floor, "strict certified floor");
      minimum_floor_gap =
          std::min(minimum_floor_gap, direct - certificate.floor);
      ++strict_tau_rows;

      const Count marked_adjacencies =
          static_cast<Count>(context.out.size() - 1) * context.hamilton;
      need(context.mass >= marked_adjacencies,
           "all Hamilton-path adjacencies inject into exposed words");
      ++all_adjacency_rows;
      if (context.out.size() == 1) {
        need(direct == -180, "singleton Q5 boundary");
      } else {
        need(context.mass >= context.hamilton,
             "non-singleton marked-adjacency bound");
        need(direct > 2730 * context.hamilton * context.hamilton,
             "strict non-singleton Q5 floor");
        minimum_nonsingleton_gap =
            std::min(minimum_nonsingleton_gap,
                     direct - 2730 * context.hamilton * context.hamilton);
        ++nonsingleton_rows;
      }
      fnv_update(semantic_digest,
                 label(context.out) + "|" + std::to_string(direct) + "|" +
                     std::to_string(certificate.floor) + "|" +
                     std::to_string(certificate.debt) + "\n");
    }

    std::vector<Count> p1_values;
    std::vector<Count> p2_values;
    for (std::size_t tail = 0; tail <= 6; ++tail) {
      const Data q = build_data(q_tail(tail));
      p1_values.push_back(literal_remainder(q, p1));
      if (tail <= 5) {
        p2_values.push_back(literal_remainder(q, p2));
      }
    }
    const std::vector<Count> expected_p1 = {
        -72, -216, -468, -900, -1332, -180, 10764,
    };
    const std::vector<Count> expected_p2 = {
        -288, -684, -1368, -2232, -1512, 10584,
    };
    need(p1_values == expected_p1, "sharp P1 boundary table");
    need(p2_values == expected_p2, "sharp P2 boundary table");

    std::size_t later_tail_rows = 0;
    for (std::size_t tail : {std::size_t{6}, std::size_t{7}}) {
      const Data q = build_data(q_tail(tail));
      const Jet lambda = left_jet(q);
      for (const Data& context : contexts) {
        const Count value = dot(lambda, right_jet(context));
        need(value >= 10764 * context.hamilton * context.hamilton,
             "later-tail propagated lower-bound control");
        ++later_tail_rows;
      }
    }

    const Data two_cycle = build_data(ordinal_out(q5.out, q5.out));
    const Jet two_cycle_lambda = left_jet(two_cycle);
    std::size_t no_sink_controls = 0;
    for (const Data& context : contexts) {
      if (context.out.size() >= 3 && !has_sink(context.out)) {
        need(dot(two_cycle_lambda, right_jet(context)) > 0,
             "final-tail-five two-cycle no-sink control");
        ++no_sink_controls;
      }
    }
    need(no_sink_controls == 34, "labelled no-sink controls");

    std::ostringstream digest_stream;
    digest_stream << std::hex << std::setfill('0') << std::setw(16)
                  << semantic_digest;
    std::cout
        << "theorem=THM-4216\n"
        << "audit=INDEPENDENT_ACCEPT\n"
        << "construction=literal_adjacency_subset_path_dp_no_imported_transfer\n"
        << "q5_hamilton=" << q5.hamilton << "\n"
        << "q5_capacity_mass=" << q5.mass << "\n"
        << "q5_response_lambda=" << join(q5_lambda) << "\n"
        << "q5_unary_coefficients_h2_hw_w2_L0_L1_M00_M01_M11="
        << coefficient_h2 << ',' << coefficient_hw << ',' << coefficient_w2
        << ',' << q5_lambda[3] << ',' << q5_lambda[4] << ','
        << q5_lambda[8] << ',' << q5_lambda[9] << ',' << q5_lambda[10]
        << "\n"
        << "quadratic_change_m_p_tau=" << moment_coefficient << ','
        << mixed_coefficient << ',' << tau_coefficient << "\n"
        << "fan_rewrite_m_p_tau=-352116,-1170,18\n"
        << "floor_debt=E+117372*DeltaV+1170*((W+H)*H-p)+18*tau\n"
        << "literal_contexts_orders_1_to_4=" << contexts.size() << "\n"
        << "literal_q5_remainder_rows=" << contexts.size() << "\n"
        << "all_adjacency_injection_rows=" << all_adjacency_rows << "\n"
        << "strict_tau_rows=" << strict_tau_rows << "\n"
        << "minimum_strict_floor_gap=" << minimum_floor_gap << "\n"
        << "nonsingleton_positive_rows=" << nonsingleton_rows << "\n"
        << "minimum_strict_gap_over_2730H2=" << minimum_nonsingleton_gap
        << "\n"
        << "p2_hostiles_n0_to_n4="
        << join(std::vector<Count>(p2_values.begin(), p2_values.begin() + 5))
        << "\n"
        << "p1_hostiles_n0_to_n5="
        << join(std::vector<Count>(p1_values.begin(), p1_values.begin() + 6))
        << "\n"
        << "later_tail_jet_rows=" << later_tail_rows << "\n"
        << "two_cycle_hamilton=" << two_cycle.hamilton << "\n"
        << "two_cycle_capacity_mass=" << two_cycle.mass << "\n"
        << "final_tail_five_no_sink_controls=" << no_sink_controls << "\n"
        << "semantic_fnv1a64=" << digest_stream.str() << "\n";
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
