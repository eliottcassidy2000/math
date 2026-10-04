"""Compile a shared checked frontier into coarse all-height descent families.

The initial graph contains only the previously saved 802 checked edges. New
queries are counted; family parameters never become unproved root seeds.
"""
from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path
import json

import adaptive_boundary_selector_20261004 as scalar
import adaptive_observation_union_20261004 as observation
from paley_geometry_observations_20261004 import residue as certificate_residue


def need(ok, message):
    if not ok:
        raise ValueError(message)


def oddpart(n):
    need(type(n) is int and n > 0, 'positive exact oddpart argument')
    return n >> scalar.v2(n)


@dataclass(frozen=True)
class Family:
    word: tuple
    P: int
    Q: int
    B: int
    residue: int
    nominal: int
    first_parameter: int

    def source(self, parameter):
        need(type(parameter) is int and parameter >= self.first_parameter, 'guarded family parameter')
        return self.residue+self.Q*parameter

    def endpoint(self, parameter):
        self.source(parameter)
        return oddpart(self.nominal+self.P*parameter)

    def extra_valuation(self, parameter):
        self.source(parameter)
        return scalar.v2(self.nominal+self.P*parameter)


def coarse_replay(source, word):
    """Last displayed exponent is a lower bound; all preceding ones are exact."""
    current, actual = source, []
    for position, lower in enumerate(word):
        current, exponent = observation.literal_step(current)
        need(exponent == lower if position+1 < len(word) else exponent >= lower,
             'exact head and guarded terminal division')
        actual.append(exponent)
    return current, tuple(actual)


def compile_family(word):
    word = tuple(word)
    need(bool(word) and all(type(a) is int and a > 0 for a in word), 'positive valuation word')
    P, Q, B = scalar.compose(word)
    if P >= Q:
        return None  # No unrefined nominal-contraction conclusion.
    residue = -B*pow(P, -1, Q) % Q
    need(residue > 0 and residue % 2, 'least positive odd coarse residue')
    nominal = (P*residue+B)//Q
    first = max(0, (nominal-residue)//(Q-P)+1)
    result = Family(word, P, Q, B, residue, nominal, first)
    n = result.source(first)
    endpoint, actual = coarse_replay(n, word)
    need(endpoint == result.endpoint(first) < n and nominal+P*first < n,
         'checked least covered boundary and strict nominal payment')
    need(actual[-1] == word[-1]+result.extra_valuation(first), 'retained extra final division')
    return result


def transport(head, tail):
    """Compose known exact head with a guarded tail; independently repay source."""
    p, q, b = scalar.compose(head)
    composed = compile_family(tuple(head)+tail.word)
    if composed is None:
        return None
    need((composed.P, composed.Q, composed.B) ==
         (tail.P*p, tail.Q*q, tail.P*b+tail.B*q), 'affine prefix/tail composition')
    pulled = (q*tail.residue-b)*pow(p, -1, q*tail.Q) % (q*tail.Q)
    need(pulled == composed.residue, 'coarse guard pullback')
    x0 = (p*composed.residue+b)//q
    offset = (x0-tail.residue)//tail.Q
    need(x0-tail.residue == tail.Q*offset, 'integer tail parameter offset')
    first = max(composed.first_parameter,
                max(0, (tail.first_parameter-offset+p-1)//p))
    composed = replace(composed, first_parameter=first)
    need(offset+p*first >= tail.first_parameter, 'transport retains the tail parameter threshold')
    return composed, offset, p


def valuation_class(family, extra):
    need(type(extra) is int and extra >= 0, 'nonnegative exact extra valuation')
    modulus = 1 << (extra+1)
    residue = ((1 << extra)-family.nominal)*pow(family.P, -1, modulus) % modulus
    return residue, modulus


def crt_coprime(a, m, b, n):
    return (a+m*((b-a)*pow(m, -1, n) % n)) % (m*n), m*n


def phase_and_valuation_parameter(family, source_phase, extra, phase_modulus=247):
    need(type(phase_modulus) is int and phase_modulus > 0 and phase_modulus % 2,
         'positive odd phase modulus')
    need(type(source_phase) is int and 0 <= source_phase < phase_modulus, 'canonical requested phase')
    binary, power = valuation_class(family, extra)
    phase = (source_phase-family.residue)*pow(family.Q, -1, phase_modulus) % phase_modulus
    residue, period = crt_coprime(binary, power, phase, phase_modulus)
    if residue < family.first_parameter:
        residue += ((family.first_parameter-residue+period-1)//period)*period
    return residue, period


def completion_exponents(family, hub):
    """All sufficiently large extra exponents for one supplied ternary-unit hub."""
    scalar.source_guard(hub)
    need(hub % 3 != 0, 'a nonempty odd route cannot end at a ternary multiple')
    exponent = next(a for a in (0, 1) if pow(2, a, 3)*hub % 3 == family.nominal % 3)
    modulus, period = 3, 2
    for _ in range(1, len(family.word)):
        modulus *= 3
        candidates = [exponent+d*period for d in range(3)
                      if pow(2, exponent+d*period, modulus)*hub % modulus == family.nominal % modulus]
        need(len(candidates) == 1, 'unique ternary digit of the completion exponent')
        exponent = candidates[0]
        period *= 3
    need(modulus == family.P, 'word length supplies the ternary denominator')
    threshold = family.nominal+family.P*family.first_parameter
    minimum = max(0, threshold.bit_length()-hub.bit_length())
    if hub << minimum < threshold:
        minimum += 1
    if hub == 1:
        minimum = max(minimum, 3-family.word[-1])  # Exclude a terminal root self-return.
    start = max(0, (minimum-exponent+period-1)//period)
    return exponent+period*start, period


def completed_certificate(family, hub, hub_certificate, parameter, hub_bit_cap=10000):
    """Construct an infinite completed subfamily; never expand its huge source."""
    need(type(parameter) is int and parameter >= 0, 'nonnegative completed-family parameter')
    need(type(hub_bit_cap) is int and hub_bit_cap >= 1, 'explicit supplied-hub expansion cap')
    codec = observation.codec
    codec.audit_certificate(hub_certificate)
    need(codec.expand(hub_certificate, bit_cap=hub_bit_cap) == hub, 'supplied hub source identity')
    first, period = completion_exponents(family, hub)
    extra = first+period*parameter
    word = family.word[:-1]+(family.word[-1]+extra,)
    cert = hub_certificate
    for exponent in reversed(word):
        numerator_mod9 = (pow(2, exponent, 9)*codec.mod3(cert, 2)-1) % 9
        need(numerator_mod9 % 3 == 0, 'completed inverse edge is integral')
        row = numerator_mod9//3
        least = codec.kappa(codec.mod3(cert, 2), row)
        need(exponent >= least and (exponent-least) % 6 == 0, 'completed inverse row/block')
        cert = codec.extend(cert, row, (exponent-least)//6)
    need(codec.ranks(cert) == (len(word)+codec.ranks(hub_certificate)[0],
                              sum(word)+len(word)+codec.ranks(hub_certificate)[1]),
         'completed first-hit clock includes the supplied suffix')
    return cert, extra


def completed_parameter_residue(family, hub, extra, modulus):
    scalar.source_guard(hub)
    need(type(extra) is int and extra >= 0, 'nonnegative exact extra exponent')
    need(type(modulus) is int and modulus >= 1, 'positive modular read precision')
    lifted = family.P*modulus
    numerator = (pow(2, extra, lifted)*hub-family.nominal) % lifted
    need(numerator % family.P == 0, 'exact ternary division in parameter reader')
    return numerator//family.P


def initial_graph():
    root = Path(__file__).resolve().parents[2]
    output = (root/'05-knowledge/results/adaptive_observation_union_20261004.out').read_text()
    line = next(line for line in output.splitlines() if line.startswith('INITIAL_GRAPH_JSON '))
    encoded = line.split(' ', 1)[1]
    need(sha256(encoded.encode()).hexdigest() ==
         '6cc583f405c674fd50b593ea4916f954c5f189b77873abb70ba9853e4b5620e8', 'frozen initial graph identity')
    data = json.loads(encoded)
    graph = observation.ObservationGraph()
    for source, target, exponent in data['edges']:
        graph.add_edge(source, target, exponent, ('initial802', source))
    return graph


def observe_to_ground(graph, source, grounded, query_cap, owner):
    current, word, queries, reads, seen = source, [], 0, 0, set()
    while current not in grounded:
        if current in seen or len(word) >= 128:
            return None, tuple(word), queries, reads
        seen.add(current)
        if current not in graph.edges:
            if queries >= query_cap:
                return None, tuple(word), queries, reads
            target, exponent = scalar.step(current)
            need(graph.add_edge(current, target, exponent, ('new-query', owner, queries)), 'fresh exact query')
            queries += 1
        else:
            reads += 1
        current, exponent = graph.edges[current]
        word.append(exponent)
    return current, tuple(word), queries, reads


def learn_segments(initial, requested, strategy, total_cap=128):
    need(strategy in ('frontier', 'original'), 'declared discovery strategy')
    need(type(total_cap) is int and total_cap >= 0, 'nonnegative exact query cap')
    graph, learned, queries, reads = initial.copy(), [], 0, 0
    while queries < total_cap:
        grounded = observation.least_closure(graph.edges)
        missing = [n for n in requested if n not in grounded]
        if not missing:
            break
        if strategy == 'frontier':
            demand, cyclic = observation.missing_frontier(graph, requested, grounded)
            if not demand:
                break
            source = min(demand, key=lambda n: (-len(demand[n]), n))
            beneficiaries = tuple(demand[source])
        else:
            source = min(missing)
            beneficiaries = (source,)
        endpoint, word, qcount, rcount = observe_to_ground(graph, source, grounded, total_cap-queries, strategy)
        queries += qcount
        reads += rcount
        if endpoint is None:
            break
        family = compile_family(word)
        if family is not None:
            need(source % family.Q == family.residue, 'observed source in learned family')
            parameter = (source-family.residue)//family.Q
            if parameter < family.first_parameter:
                family = None  # Keep the grounded point proof; do not promote this instance.
            else:
                need(family.endpoint(parameter) == endpoint, 'observed family instance reaches grounded join')
        learned.append((source, endpoint, family, beneficiaries, qcount, rcount))
    return graph, tuple(learned), queries, reads


def family_summary(family):
    return dict(length=len(family.word), valuation=sum(family.word), P=family.P, Q=family.Q,
                B=family.B, residue=family.residue, nominal=family.nominal,
                first_parameter=family.first_parameter)


def main():
    requested = tuple(range(1, 1024, 2))
    initial = initial_graph()
    print('FRONTIER FAMILY COMPILER: rooted finite discovery plus proved guarded all-height rules')
    print('Initial802 checked edges; requested512 odds1..1023; only ROOT1; query cap128, segment-word cap128')
    edge_graph, edge_trace = observation.frontier_expansion(initial, requested, cap=128)
    frontier_graph, families, fq, fr = learn_segments(initial, requested, 'frontier')
    original_graph, originals, oq, ore = learn_segments(initial, requested, 'original')
    need(all(row[2] is not None for row in families+originals),
         'every completed segment in this declared finite test has a paid family')
    need(edge_graph.edges == frontier_graph.edges == original_graph.edges, 'same final checked transition set')
    print('Strategy comparison:edge-only', len(edge_trace), 'new queries;',
          'frontier', fq, 'new queries,', fr, 'cached prefix reads, learned lengths', [len(x[2].word) for x in families],
          '; original', oq, 'new queries,', ore, 'cached prefix reads, learned lengths', [len(x[2].word) for x in originals])
    print('These reads are discovery-path reads only; graph-closure and verification work are separate.')
    transported = []
    for source, endpoint, family, beneficiaries, qcount, rcount in families:
        print('Frontier family', source, 'join', endpoint, 'beneficiaries', beneficiaries,
              json.dumps(family_summary(family), sort_keys=True))
        for owner in beneficiaries:
            path, head, cyclic = observation.route_walk(initial.edges, owner)
            need(not cyclic and path[-1] == source, 'checked observed head reaches this frontier')
            pulled, offset, multiplier = transport(head, family)
            need(owner == pulled.residue, 'declared original is least transported boundary')
            transported.append(pulled)
            print('  pulled family', owner, 'head length', len(head), 'tail parameter', [offset, multiplier],
                  json.dumps(family_summary(pulled), sort_keys=True))
    all_families = [row[2] for row in families]+transported
    instances = phases = 0
    for family in all_families:
        for k in (family.first_parameter, family.first_parameter+1, 2, 7, 19, 247, 1024):
            k = max(k, family.first_parameter)
            n = family.source(k)
            endpoint, actual = coarse_replay(n, family.word)
            need(endpoint == family.endpoint(k) < n, 'independent all-height family sample')
            need(actual[-1] == family.word[-1]+family.extra_valuation(k), 'exact extra valuation')
            instances += 1
        for h in range(16):
            residue, period = valuation_class(family, h)
            k = residue if residue >= family.first_parameter else residue+period
            need(family.extra_valuation(k) == h and family.extra_valuation(k+period) == h,
                 'exact dyadic parameter refinement')
            for source_phase in range(247):
                k, period = phase_and_valuation_parameter(family, source_phase, h)
                n = family.source(k)
                endpoint, actual = coarse_replay(n, family.word)
                need(n % 247 == source_phase and family.extra_valuation(k) == h,
                     'phase and true terminal clock jointly realized')
                predicted = (family.P*n+family.B)*pow(2, -(sum(family.word)+h), 247) % 247
                need(predicted == endpoint % 247 and sum(actual) == sum(family.word)+h,
                     'incoming joined mod36 clock uses actual extra division')
                phases += 1
    print('Family replay instances:', instances, '; complete phase/extra-clock controls:', phases)

    shared = families[0][2]
    n0, n1 = shared.source(0), shared.source(247)
    need(n0 % 247 == n1 % 247 and shared.extra_valuation(0) == 0 and shared.extra_valuation(247) == 1,
         'same source phase with different true terminal clock')
    wrong = (shared.P*n1+shared.B)*pow(2, -sum(shared.word), 247) % 247
    need(wrong != shared.endpoint(247) % 247, 'nominal clock predicts the wrong actual phase')
    print('Clock hostile:k0/k247 share source phase', n0 % 247, 'but extra halvings0/1;',
          'nominal phase', wrong, 'actual', shared.endpoint(247) % 247)
    tail = compile_family((2,))
    need(tail.residue == 1 and tail.first_parameter == 1, 'root boundary remains separately typed')
    need(transport((1,), tail) is None and scalar.step(41)[0] == 31 > 27,
         'a frontier descent does not automatically repay original27')
    guarded, offset, multiplier = transport((4,), tail)
    need(guarded.residue == 5 and guarded.first_parameter == 1 and offset == 0 and multiplier == 3,
         'pullback must retain tail threshold even if a boundary has another root proof')
    need(coarse_replay(3, (1, 2))[0] == 1 and coarse_replay(27, (1, 2))[0] == 31,
         'rejected unrefined nominal family still contains paid individual subfamilies')
    print('Boundary hostile:27->41->31; tail descends but pullback unrefined nominal family is not paid.')
    zero, no_rules, queries, _ = learn_segments(initial, requested, 'frontier', total_cap=0)
    need(zero.edges == initial.edges and not no_rules and queries == 0, 'zero cap supplies no new evidence')
    partial, no_rules, queries, _ = learn_segments(initial, requested, 'frontier', total_cap=1)
    need(len(partial.edges) == 803 and not no_rules and queries == 1,
         'unfinished observation survives without promotion to a family or certificate')

    ranks = observation.least_closure(frontier_graph.edges)
    certificates = observation.graph_certificates(frontier_graph, ranks)
    need(all(n in certificates for n in requested), 'declared finite requested universe completed')
    for n in requested:
        need(certificates[n] == observation.codec.encode_source(n), 'post-closure independent first-hit audit')
    print('Grounded final graph:', len(frontier_graph.edges), 'edges;', len(certificates),
          'certified vertices; requested', len(requested), 'all complete.')
    # Completion uses only hubs already grounded in the original802-edge graph.
    old_ranks = observation.least_closure(initial.edges)
    supplied_hubs = observation.graph_certificates(initial, old_ranks)
    completions, expanded, maximum_extra = 0, 0, 0
    for family in all_families:
        for hub in (1, 5, 157, 797):
            residues19 = []
            for parameter in range(3):
                cert, extra = completed_certificate(family, hub, supplied_hubs[hub], parameter)
                maximum_extra = max(maximum_extra, extra)
                for modulus in (81, 4096, 361, 247, 223, 233):
                    k = completed_parameter_residue(family, hub, extra, modulus)
                    expected = (family.residue+family.Q*k) % modulus
                    need(certificate_residue(cert, modulus) == expected, 'independent completed-source modular reader')
                residues19.append(certificate_residue(cert, 19))
                if observation.codec.bit_bounds(cert)[1] <= 10000:
                    n = observation.codec.expand(cert)
                    k = ((hub << extra)-family.nominal)//family.P
                    need(k >= family.first_parameter and n == family.source(k), 'expanded completed family membership')
                    need(coarse_replay(n, family.word)[0] == hub, 'expanded supplied-hub return')
                    expanded += 1
                completions += 1
            if len(family.word) >= 3:
                need(len(set(residues19)) == 1, 'completed ternary exponent progression fixes first19 digit')
    try:
        completed_certificate(shared, 3, supplied_hubs[3], 0)
    except ValueError:
        pass
    else:
        raise ValueError('ternary-multiple hub accepted')
    print('Completed subfamilies:', completions, 'symbolic AST controls;', expanded,
          'expanded controls; largest extra valuation', maximum_extra,
          '; zero new odd-edge queries for completion.')
    returns223 = [((1, 1, 2, 1, 2, 1, 4), (1703, 4096, 910, 2187)),
                  ((1, 1, 1, 1, 2, 1, 1, 5), (6815, 8192, 5459, 6561))]
    return_controls = 0
    for word, expected in returns223:
        family = compile_family(word)
        need((family.residue, family.Q, family.nominal, family.P) == expected and family.first_parameter == 0,
             'shared compiler interface for prime-derived return tails')
        need(family.B % 223 == 0, 'retained223 return predicate')
        for stop in range(1, len(word)):
            p, q, _ = scalar.compose(word[:stop])
            need(p > q, 'every proper prefix has strict expanding coefficient')
        for k in range(16):
            source, current = family.source(k), family.source(k)
            for index in range(len(word)):
                current, exponent = observation.literal_step(current)
                need(current > source if index+1 < len(word) else current < source,
                     'unconditional dyadic first descent, without223 filter')
            need(current == family.endpoint(k), 'returned shared family endpoint')
            return_controls += 1
        k0 = -family.residue*pow(family.Q, -1, 223) % 223
        for lift in range(9):
            k = k0+223*lift
            need(family.source(k) % 223 == family.endpoint(k) % 223 == 0,
                 '223 filter selects modular return within the whole descent cylinder')
            return_controls += 1
    print('Prime-derived shared Family interface:', return_controls,
          'controls;223 filter retained for modular return, unnecessary for whole-cylinder first descent.')
    print('No family parameter is an oracle seed; general family conclusion is descent, not unconditional home.')


if __name__ == '__main__':
    main()
