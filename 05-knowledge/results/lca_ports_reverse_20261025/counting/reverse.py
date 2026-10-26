"""Ordinary A000568: structural DAG, reverse potential, normalized residues.

The structural graph never depends on arithmetic weights or modular support.
Only the ordinary (zero clone-excess) transfer is implemented here.
"""
from dataclasses import dataclass


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def valid_order(n):
    if type(n) is not int or n < 0:
        raise ValueError('order must be a nonnegative integer')


def baseline(n):
    return n // 2 if n >= 3 else 0


def lower_bound(defect):
    return max(((value + 2*i - 1) // (2*i) for i, value in enumerate(defect) if i), default=0)


def factorial_twos(k):
    return k - k.bit_count()


def odd_factorials(n, modulus=None):
    values = [1]
    for j in range(1, n + 1):
        value = values[-1] * (j // (j & -j))
        values.append(value if modulus is None else value % modulus)
    return values


@dataclass(frozen=True)
class State:
    stage: int
    remaining: int
    defect: tuple


@dataclass(frozen=True)
class Edge:
    target: int
    multiplicity: int
    cost: int


class OrdinaryDAG:
    def __init__(self, n):
        valid_order(n)
        from math import gcd
        self.n = n
        top = n if n % 2 else max(0, n - 1)
        self.states = [State(top, n, (0,) * ((top + 1) // 2))]
        self.edges = [[]]
        self.layers = [[0]]
        self.indegree = [0]
        if n:
            for ell in range(top, 0, -2):
                next_ids = {}
                next_layer = []
                for source in self.layers[-1]:
                    state = self.states[source]
                    need(state.stage == ell, 'stage is part of the state')
                    s, defect = state.remaining, state.defect
                    need(len(defect) == (min(ell, s) + 1) // 2, 'exact affordable-coordinate mask')
                    b = lower_bound(defect)
                    selected = defect[(ell - 1) // 2] if ell <= s else 0
                    choices = (s,) if ell == 1 else range(s // ell + 1)
                    for k in choices:
                        rest = s - ell * k
                        dn = tuple(defect[(d - 1) // 2] + k * (d - gcd(ell, d))
                                   for d in range(1, min(ell - 2, rest) + 1, 2))
                        bp = lower_bound(dn)
                        intrinsic = k * (ell - 1) // 2 + ell * k * (k - 1) // 2
                        cost = intrinsic - factorial_twos(k) + k * (ell*b - selected) + rest * (k+b-bp)
                        need(cost >= 0 and bp <= b+k, 'nonnegative original integral cost')
                        key = State(max(0, ell - 2), rest, dn)
                        if key not in next_ids:
                            target = len(self.states)
                            next_ids[key] = target
                            next_layer.append(target)
                            self.states.append(key)
                            self.edges.append([])
                            self.indegree.append(0)
                        target = next_ids[key]
                        self.edges[source].append(Edge(target, k, cost))
                        self.indegree[target] += 1
                self.layers.append(next_layer)
        need(len(self.layers[-1]) == 1, 'one final terminal')
        self.terminal = self.layers[-1][0]
        need(self.states[self.terminal] == State(0, 0, ()), 'terminal state')

        # All states were generated from the initial vertex. Every positive
        # remainder has an all-singleton continuation, so all reach terminal.
        self.potential = [None] * len(self.states)
        self.potential[self.terminal] = 0
        self.reverse_relaxations = 0
        for source in range(len(self.states) - 2, -1, -1):
            outgoing = self.edges[source]
            need(bool(outgoing), 'reachable state has a continuation')
            need(all(edge.target > source and self.potential[edge.target] is not None
                     for edge in outgoing), 'topological reverse pass')
            self.potential[source] = min(edge.cost + self.potential[edge.target] for edge in outgoing)
            self.reverse_relaxations += len(outgoing)
        self.B = self.potential[0]
        need(self.B == baseline(n), 'reverse minimum agrees with proved termwise baseline')
        need(all(self.reduced_cost(i, edge) >= 0 for i, outgoing in enumerate(self.edges)
                 for edge in outgoing), 'all reduced costs are nonnegative')
        self._structure_cache = None

    def reduced_cost(self, source, edge):
        return edge.cost + self.potential[edge.target] - self.potential[source]

    def structure(self):
        if self._structure_cache is not None:
            return dict(self._structure_cache)
        edge_count = sum(map(len, self.edges))
        repeated = len(self.states) - len({(v.remaining, v.defect) for v in self.states})
        self._structure_cache = {'vertices_including_terminal': len(self.states), 'edges': edge_count,
                'merged_vertices': sum(d > 1 for d in self.indegree),
                'extra_incoming_edges': sum(max(0, d - 1) for d in self.indegree),
                'defect_coordinates_stored': sum(len(v.defect) for v in self.states),
                'stage_distinctions_lost_if_omitted': repeated,
                'reverse_relaxations': self.reverse_relaxations,
                'max_original_cost': max((e.cost for row in self.edges for e in row), default=0),
                'max_reduced_cost': max((self.reduced_cost(i, e) for i, row in enumerate(self.edges) for e in row), default=0),
                'baseline': self.B,
                'peak_weight_slots_two_layers': max((len(a) + len(b) for a,b in zip(self.layers,self.layers[1:])), default=1)}
        return dict(self._structure_cache)

    def residue(self, K, *, normalized=True, prune_edges=True, skip_zero_sources=True):
        """Normalized a(n)/2^B mod2^K, or absolute a(n) mod2^K."""
        if type(K) is not int or K < 1:
            raise ValueError('precision must be a positive integer')
        modulus = 1 << K
        odd = odd_factorials(self.n, modulus)
        transfer_cache = {}
        values = {0: 1}
        stats = {'precision_bits': K, 'edge_scans': 0, 'cost_pruned_edges': 0,
                 'zero_source_edges': 0, 'modular_updates': 0,
                 'multiply_operand_bit_sum': 0, 'max_residue_bits': 1,
                 'max_product_bits': 0, 'max_shifted_product_bits': 0,
                 'max_modular_sum_bits': 0, 'layer_residue_bit_sum': 1,
                 'odd_factorial_setup_products': self.n}
        for layer, next_layer in zip(self.layers, self.layers[1:]):
            destination = dict.fromkeys(next_layer, 0)
            for source in layer:
                value = values[source]
                state = self.states[source]
                u = self.n - state.remaining
                for edge in self.edges[source]:
                    stats['edge_scans'] += 1
                    shift = self.reduced_cost(source, edge) if normalized else edge.cost
                    if prune_edges and shift >= K:
                        stats['cost_pruned_edges'] += 1
                        continue
                    if skip_zero_sources and value == 0:
                        stats['zero_source_edges'] += 1
                        continue
                    key = (u, state.stage, edge.multiplicity)
                    if key not in transfer_cache:
                        k = edge.multiplicity
                        up = u + state.stage*k
                        denominator = odd[u] * pow(state.stage,k,modulus) * odd[k] % modulus
                        need(denominator % 2 == 1, 'odd modular denominator')
                        transfer_cache[key] = odd[up] * pow(denominator,-1,modulus) % modulus
                        need(transfer_cache[key] % 2 == 1, 'odd modular transfer')
                    factor = transfer_cache[key]
                    stats['multiply_operand_bit_sum'] += value.bit_length() + factor.bit_length()
                    product = value*factor
                    stats['max_product_bits'] = max(stats['max_product_bits'], product.bit_length())
                    if shift >= K:
                        contribution = 0
                    else:
                        shifted = product << shift
                        stats['max_shifted_product_bits'] = max(stats['max_shifted_product_bits'], shifted.bit_length())
                        contribution = shifted % modulus
                    target = edge.target
                    summed = destination[target] + contribution
                    stats['max_modular_sum_bits'] = max(stats['max_modular_sum_bits'], summed.bit_length())
                    destination[target] = summed % modulus
                    stats['modular_updates'] += 1
            values = destination
            stats['max_residue_bits'] = max(stats['max_residue_bits'], max((x.bit_length() for x in values.values()), default=0))
            stats['layer_residue_bit_sum'] += sum(x.bit_length() for x in values.values())
        result = values[self.terminal] * pow(odd[self.n], -1, modulus) % modulus
        stats.update({'transfer_cache_entries': len(transfer_cache),
                      'odd_inverse_calls': len(transfer_cache)+1,
                      'scalar_slots_peak': self.structure()['peak_weight_slots_two_layers']})
        need(stats['max_residue_bits'] <= K, 'residues stay within requested precision')
        need(stats['edge_scans'] == stats['modular_updates'] + stats['cost_pruned_edges'] + stats['zero_source_edges'],
             'all edge work accounted for')
        return result, stats

    def exact(self, *, audit_original=False):
        """Optional full-integer replay of the same production DAG."""
        odd = odd_factorials(self.n)
        cache = {}
        Z = [0] * len(self.states)
        Y = [0] * len(self.states) if audit_original else None
        Z[0] = 1
        if audit_original:
            Y[0] = 1
        coordinate_checks = 0
        for source,state in enumerate(self.states):
            if audit_original:
                power = self.potential[source] - self.B
                if power >= 0:
                    need(Y[source] << power == Z[source], 'exact reachable-lattice identity')
                else:
                    quotient, remainder = divmod(Y[source], 1 << -power)
                    need(remainder == 0 and quotient == Z[source], 'pathwise gauge divides exactly before merging')
                coordinate_checks += 1
            u = self.n-state.remaining
            for edge in self.edges[source]:
                key = (u,state.stage,edge.multiplicity)
                if key not in cache:
                    k = edge.multiplicity
                    denominator = odd[u] * state.stage**k * odd[k]
                    factor, remainder = divmod(odd[u+state.stage*k], denominator)
                    need(remainder == 0 and factor > 0 and factor % 2 == 1, 'positive odd exact transfer')
                    cache[key] = factor
                factor = cache[key]
                Z[edge.target] += factor * Z[source] << self.reduced_cost(source,edge)
                if audit_original:
                    Y[edge.target] += factor * Y[source] << edge.cost
        quotient, remainder = divmod(Z[self.terminal], odd[self.n])
        need(remainder == 0, 'exact final odd normalization')
        if audit_original:
            need(Y[self.terminal] == Z[self.terminal] << self.B, 'terminal baseline restored')
        return quotient << self.B, {'normalized_integer': quotient,
                'coordinate_checks': coordinate_checks,
                'max_Z_bits': max(map(int.bit_length,Z)),
                'max_Y_bits': max(map(int.bit_length,Y)) if audit_original else None}

    def prefix(self, multiplicities, *, stop_stage=0):
        """Audit utility: follow specified cycle choices, retaining the structural merge."""
        source = 0
        cost = 0
        while self.states[source].stage != stop_stage:
            state = self.states[source]
            need(state.stage > stop_stage, 'requested prefix stage exists')
            k = multiplicities.get(state.stage, 0)
            choices = [edge for edge in self.edges[source] if edge.multiplicity == k]
            need(len(choices) == 1, 'one structural edge per cycle choice')
            edge = choices[0]
            cost += edge.cost
            source = edge.target
        return source, cost
