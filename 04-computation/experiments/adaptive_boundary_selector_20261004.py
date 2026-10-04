"""Bounded, source-preserving Collatz proof-obligation selection.

REDUCED means a smaller home obligation, not a completed route. PENDING is an
explicit bounded-search result. A supplied child certificate can discharge a
reduction and compile it to the inherited first-hit inverse-ray codec.
"""
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import json

import entry_20260927_recursive as inherited
import mod19_recursive_observers_20261004 as observer


def need(ok, message):
    if not ok:
        raise ValueError(message)


def source_guard(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'positive odd integer source')


def v2(n):
    need(type(n) is int and n > 0, 'positive integer valuation argument')
    return (n & -n).bit_length()-1


def step(n):
    source_guard(n)
    a = v2(3*n+1)
    return (3*n+1) >> a, a


def compose(word):
    p, q, b = 1, 1, 0
    for a in word:
        need(type(a) is int and a > 0, 'positive integer valuation letter')
        p, q, b = 3*p, q << a, 3*b+q
    return p, q, b


def replay(n, word):
    path = [n]
    for a in word:
        need(n != 1, 'no root self-return in a first-hit prefix')
        n, actual = step(n)
        need(actual == a, 'exact valuation guard')
        path.append(n)
    return tuple(path)


@dataclass(frozen=True)
class Macro:
    kind: str
    source: int
    word: tuple
    endpoint: int


def macro(kind, source, word):
    word = tuple(word)
    path = replay(source, word)
    p, q, b = compose(word)
    need(p*source+b == q*path[-1], 'independent affine/literal replay')
    return Macro(kind, source, word, path[-1])


def literal_macro(kind, source, count):
    n, word = source, []
    for _ in range(count):
        if n == 1:
            break
        n, a = step(n)
        word.append(a)
    return macro(kind, source, word)


def siblings(n):
    """All proper odd ancestors under S(a)=4a+1, in decreasing order."""
    source_guard(n)
    result = []
    while n % 8 == 5:
        n = (n-1)//4
        result.append(n)
    return tuple(result)


def repeated112(n, maximum):
    """Exact word (1,1,2)^m; guard independent of every mod-19 phase."""
    source_guard(n)
    need(type(maximum) is int and maximum >= 0, 'nonnegative repetition cap')
    m = min(maximum, (v2(11*n+19)-1)//4)
    if m < 1:
        return None
    result = macro('repeat112', n, (1, 1, 2)*m)
    p, q = 27**m, 16**m
    need(11*q*result.endpoint == p*(11*n+19)-19*q,
         'rational-anchor affine identity')
    return result


@dataclass(frozen=True)
class Selection:
    source: int
    status: str
    kind: str
    endpoint: int
    dependencies: tuple
    macros: tuple
    phase19: tuple

    @property
    def word(self):
        return tuple(a for record in self.macros for a in record.word)


def audit_selection(result):
    source_guard(result.source)
    need(result.status in ('ROOT', 'REDUCED', 'PENDING'), 'declared proof status')
    n = result.source
    for record in result.macros:
        need(record.source == n, 'macro preserves actual current source')
        n = replay(n, record.word)[-1]
        need(n == record.endpoint, 'stored macro endpoint')
    need(n == result.endpoint, 'enclosing endpoint')
    p, q, b = compose(result.word)
    need(p*result.source+b == q*n, 'immutable-source affine transfer')
    if result.status == 'ROOT':
        need(result.source == 1 and not result.word and result.dependencies == (), 'root convention')
    elif result.status == 'PENDING':
        need(result.dependencies == () and result.kind == 'budget', 'no pending conclusion')
    else:
        need(bool(result.dependencies) and all(type(a) is int and 0 < a < result.source
                                             and a % 2 for a in result.dependencies),
             'strictly smaller positive odd proof dependencies')
        if result.kind == 'actual':
            need(result.dependencies == (n,), 'actual endpoint dependency')
        else:
            need(result.kind == 'common_future' and all(a == n or step(a)[0] == n for a in result.dependencies),
                 'actual route and each auxiliary route meet at the retained join')
    return result


def finish(source, status, kind, current, dependencies, records, phases):
    if status == 'REDUCED' and kind == 'common_future' and current < source:
        dependencies = tuple(dict.fromkeys((*dependencies, current)))
    return audit_selection(Selection(source, status, kind, current, tuple(dependencies),
                                    tuple(records), tuple(phases)))


@lru_cache(None)
def bank():
    return inherited.bank()


@dataclass(frozen=True)
class LearnedDescent:
    word: tuple
    residue: int
    modulus: int
    endpoint_at_residue: int


def discover_descent(source, cap=128):
    """Bounded discovery, followed by an all-cylinder proof check."""
    source_guard(source)
    need(source > 1 and type(cap) is int and cap >= 1, 'nonroot source and positive discovery cap')
    n, word = source, []
    for _ in range(cap):
        n, a = step(n)
        word.append(a)
        if n < source:
            break
    else:
        return None
    p, q, b = compose(word)
    need(q > p, 'positive carry forces contraction at an actual descent')
    residue = (q-b)*pow(p, -1, 2*q) % (2*q)
    endpoint = (p*residue+b)//q
    # A failed least-boundary check remains an undispatched obligation.
    if endpoint >= residue:
        return None
    current = residue
    for a in word:
        if current == 1:
            return None  # First-hit prefix typing is retained at the boundary.
        current, actual = step(current)
        need(actual == a, 'discovered exact cylinder')
    need(current == endpoint, 'least-boundary literal and affine agreement')
    need(p*residue+b == q*endpoint and residue % 2, 'exact cylinder at its least positive source')
    return LearnedDescent(tuple(word), residue, 2*q, endpoint)


def validate_learned(row):
    need(isinstance(row, LearnedDescent), 'typed learned cylinder')
    p, q, b = compose(row.word)
    need(q > p and row.modulus == 2*q and 0 < row.residue < 2*q and row.residue % 2,
         'contracting exact cylinder')
    need(replay(row.residue, row.word)[-1] == row.endpoint_at_residue < row.residue,
         'checked least-boundary certificate')
    need(p*row.residue+b == q*row.endpoint_at_residue, 'learned cylinder affine identity')


def select(source, budget=8, policy='adaptive', word_cap=128, phase_depth=2, learned=()):
    """Bounded macro grammar; residue phases annotate, never authorize a rule."""
    source_guard(source)
    need(type(budget) is int and budget >= 0, 'nonnegative integer macro budget')
    need(type(word_cap) is int and word_cap >= 3, 'integer per-macro word cap at least three')
    need(type(phase_depth) is int and phase_depth >= 1, 'positive mod-19 precision')
    need(policy in ('literal', 'inherited', 'adaptive'), 'known bounded policy')
    for row in learned:
        validate_learned(row)
    if source == 1:
        return finish(source, 'ROOT', 'root', 1, (), (), ())
    n, records, phases = source, [], []
    for _ in range(budget):
        phases.append((n, observer.chart112_observer(n, phase_depth)))
        matched = next((row for row in learned if len(row.word) <= word_cap
                        and n % row.modulus == row.residue), None)
        if matched is not None:
            chosen = macro('learned_cylinder', n, matched.word)
            need(chosen.endpoint < n, 'learned cylinder own-source descent')
        elif policy == 'literal':
            chosen = literal_macro('step', n, 1)
        elif (3*n+1) & (3*n) == 0:
            chosen = literal_macro('to_one', n, 1)
        else:
            matches = [(row['J'], q) for q, row in bank().items()
                       if row['J']+1 <= word_cap and n % (1 << row['K']) == row['residue']]
            if matches:
                length, q = min(matches)
                chosen = literal_macro('core:'+str(q), n, length+1)
                need(chosen.endpoint < n, 'inherited core own-source conclusion')
            else:
                if policy == 'adaptive':
                    import virtual_contraction_ladders_20261004 as virtual
                    recognized = virtual.recognize(n)
                    if recognized is not None:
                        chart, child = recognized
                        if child < source and chart.j+1 <= word_cap:
                            chosen = macro('virtual:'+str(chart.s), n,
                                           (1,)*chart.j+(2*chart.s+1,))
                            need(step(child)[0] == chosen.endpoint, 'virtual common future')
                            records.append(chosen)
                            return finish(source, 'REDUCED', 'common_future', chosen.endpoint,
                                          (child,), records, phases)
                    children = tuple(a for a in siblings(n) if a < source)
                    if children:
                        chosen = literal_macro('sibling_join', n, 1)
                        records.append(chosen)
                        return finish(source, 'REDUCED', 'common_future', chosen.endpoint,
                                      children, records, phases)
                    chosen = repeated112(n, word_cap//3)
                else:
                    chosen = None
                if chosen is None:
                    if n % 4 == 1:
                        chosen = literal_macro('step', n, 1)
                    elif v2(n+5) >= 4:
                        m = min((v2(n+5)-1)//3, word_cap//2)
                        chosen = macro('repeat12', n, (1, 2)*m)
                    else:
                        length = min(v2(n+1)-1, word_cap)
                        need(length >= 1, 'complete fallback branch')
                        chosen = macro('repeat1', n, (1,)*length)
        records.append(chosen)
        n = chosen.endpoint
        if n < source:
            return finish(source, 'REDUCED', 'actual', n, (n,), records, phases)
    return finish(source, 'PENDING', 'budget', n, (), records, phases)


def discharge(result, child_certificate, bit_cap=10000):
    """Compile a discharged smaller obligation to a genuine first-hit home AST."""
    audit_selection(result)
    need(type(bit_cap) is int and bit_cap > 0, 'positive explicit certificate expansion cap')
    codec = observer.ray_codec()
    codec.audit_certificate(child_certificate)
    if result.status == 'ROOT':
        need(child_certificate == codec.ROOT, 'root supplied certificate')
        return codec.ROOT
    need(result.status == 'REDUCED', 'pending search is not a certificate')
    child = codec.expand(child_certificate, bit_cap=bit_cap)
    need(child in result.dependencies, 'certificate discharges an advertised obligation')
    if result.kind == 'actual' or child == result.endpoint:
        cert = child_certificate
    else:
        cert = codec.ROOT if child == 1 else child_certificate.parent
        need(codec.expand(cert, bit_cap=bit_cap) == result.endpoint, 'retained common-future suffix')
    suffix_rank = codec.ranks(cert)
    path = replay(result.source, result.word)
    for i in reversed(range(len(result.word))):
        source, target, a = path[i], path[i+1], result.word[i]
        least = codec.kappa(target % 9, source % 3)
        need(a >= least and (a-least) % 6 == 0, 'inverse row/block decoding')
        cert = codec.extend(cert, source % 3, (a-least)//6)
    need(codec.expand(cert, bit_cap=bit_cap) == result.source,
         'compiled original source')
    need(codec.ranks(cert) == (len(result.word)+suffix_rank[0],
                              sum(a+1 for a in result.word)+suffix_rank[1]),
         'exact compiled odd and ordinary first-hit ranks')
    codec.audit_certificate(cert)
    return cert


def virtual_selection(source, word_cap=128):
    """Expose one guarded virtual rule, independently of policy precedence."""
    source_guard(source)
    need(type(word_cap) is int and word_cap >= 1, 'positive word cap')
    import virtual_contraction_ladders_20261004 as virtual
    found = virtual.recognize(source)
    if found is None or found[0].j+1 > word_cap:
        return None
    chart, child = found
    record = macro('virtual:'+str(chart.s), source, (1,)*chart.j+(2*chart.s+1,))
    return finish(source, 'REDUCED', 'common_future', record.endpoint, (child,), (record,), ())


def lookup_ray_certificate(source):
    """Optional exact lookup in nineteen certified rays, not a residue verdict.

    This helper does not participate in the frozen policy comparison tables.
    It returns None if this supplied integer fails membership, even though the
    ray compiler can construct some other certified integer in the same class.
    """
    source_guard(source)
    import mod19_route_lifts_20261004 as rays
    ray = rays.BANK[source % 19]
    coefficient, offset, denominator = ray.affine_data()
    numerator = denominator*source+offset
    if numerator % coefficient:
        return None
    quotient = numerator//coefficient
    if quotient <= 0 or quotient & (quotient-1):
        return None
    exponent = quotient.bit_length()-1
    if exponent % 18:
        return None
    return rays.certificate_for(ray, exponent//18)


def complete_at_join(source, word, suffix):
    """A certified actual join need not be smaller than the original source."""
    codec = observer.ray_codec()
    codec.audit_certificate(suffix)
    path = replay(source, word)
    need(codec.expand(suffix) == path[-1], 'supplied actual join certificate')
    cert = suffix
    for i in reversed(range(len(word))):
        least = codec.kappa(path[i+1] % 9, path[i] % 3)
        need(word[i] >= least and (word[i]-least) % 6 == 0, 'certified-join row/block')
        cert = codec.extend(cert, path[i] % 3, (word[i]-least)//6)
    need(codec.expand(cert) == source, 'certified-join original source')
    return cert


def remember_suffixes(source, cert, certificates, provenance):
    """Store only actual suffixes of one already verified home certificate."""
    codec = observer.ray_codec()
    need(codec.expand(cert) == source, 'suffix-memory parent source')
    value, cut = source, 0
    while True:
        if value in certificates:
            need(certificates[value] == cert, 'canonical suffix identity')
        else:
            certificates[value] = cert
            provenance[value] = (source, cut)
        if cert == codec.ROOT:
            need(value == 1, 'suffix memory ends at root')
            break
        value, a = step(value)
        need(a == codec.exponent(cert), 'certified suffix cut, not an orbit-search guess')
        cert, cut = cert.parent, cut+1


def seed_closure(limit, budget, policy, learned=(), suffix_memory=False, statistics=None):
    """One increasing pass, seed1 only. No orbit-derived child certificates."""
    need(type(limit) is int and limit >= 1, 'positive finite universe limit')
    need(type(suffix_memory) is bool, 'explicit suffix-memory option')
    need(statistics is None or type(statistics) is dict, 'optional statistics dictionary')
    select(1, budget, policy, learned=learned)
    codec = observer.ray_codec()
    certificates = {1: codec.ROOT}
    provenance = {1: (1, 0)}
    pending, unresolved = [], []
    for n in range(3, limit+1, 2):
        if n in certificates:
            continue
        result = select(n, budget, policy, learned=learned)
        if suffix_memory:
            path = replay(n, result.word)
            cut = next((i for i in range(1, len(path)) if path[i] in certificates), None)
            if cut is not None:
                cert = complete_at_join(n, result.word[:cut], certificates[path[cut]])
                remember_suffixes(n, cert, certificates, provenance)
                continue
        if result.status == 'PENDING':
            pending.append(n)
            continue
        available = [a for a in result.dependencies if a in certificates]
        if not available:
            unresolved.append(n)
            continue
        child = min(available, key=lambda a: (codec.ranks(certificates[a])[0], a))
        cert = discharge(result, certificates[child])
        if suffix_memory:
            remember_suffixes(n, cert, certificates, provenance)
        else:
            certificates[n] = cert
    # Later certified suffixes can settle inputs that were pending on their turn.
    pending = [n for n in pending if n not in certificates]
    unresolved = [n for n in unresolved if n not in certificates]
    restricted = {n: cert for n, cert in certificates.items() if n <= limit}
    for value, (parent, cut) in provenance.items():
        node = certificates[parent]
        for _ in range(cut):
            need(node.parent is not None, 'recorded suffix cut lies on parent proof')
            node = node.parent
        need(node == certificates[value], 'retained parent/cut provenance audit')
    need(len(restricted)+len(pending)+len(unresolved) == (limit+1)//2,
         'finite closure partition including later certified suffixes')
    if statistics is not None:
        statistics.update(distinct_certified_sources=len(certificates),
                          outside_universe=len(certificates)-len(restricted),
                          separate_odd_edge_count=sum(codec.ranks(cert)[0] for cert in restricted.values()))
    return restricted, pending, unresolved


def main():
    print('ADAPTIVE BOUNDARY SELECTOR: sound bounded reductions, not global coverage')
    print('Inherited bank rows:', len(bank()), 'maximum literal macro length:',
          max(row['J']+1 for row in bank().values()))
    universe = tuple(range(3, 8192, 2))
    print('Universe: odd integers 3..8191; count', len(universe),
          '; per-macro word cap128; fixed policy orders; no warm cache')
    for budget in (1, 2, 4, 8):
        output, results = {}, {}
        for policy in ('literal', 'inherited', 'adaptive'):
            rows = {n: select(n, budget, policy) for n in universe}
            results[policy] = rows
            closed = [r for r in rows.values() if r.status == 'REDUCED']
            output[policy] = dict(reduced=len(closed), pending=len(rows)-len(closed),
                                  common_future=sum(r.kind == 'common_future' for r in closed),
                                  max_actual_odd_horizon=max((len(r.word) for r in closed), default=0))
        gains = [n for n in universe if results['inherited'][n].status == 'PENDING'
                 and results['adaptive'][n].status == 'REDUCED']
        losses = [n for n in universe if results['inherited'][n].status == 'REDUCED'
                  and results['adaptive'][n].status == 'PENDING']
        print('budget', budget, json.dumps(output, sort_keys=True))
        print('  adaptive-minus-inherited:', len(gains), 'first', gains[:12],
              '; reverse difference:', len(losses), 'first', losses[:12])

    repeats = 0
    for n in range(1, 32768, 2):
        predicted = (v2(11*n+19)-1)//4
        actual, x = 0, n
        for _ in range(8):
            word = []
            for _ in range(3):
                if x == 1:
                    break
                x, a = step(x)
                word.append(a)
            if word != [1, 1, 2]:
                break
            actual += 1
        need(predicted == actual, 'independent maximal112 repetition decoder')
        if predicted:
            repeated112(n, predicted)
            repeats += 1
    print('Independent 112 guard controls:16384 odd sources;', repeats, 'positive matches')

    codec = observer.ray_codec()
    compiled, common = 0, 0
    for n in range(3, 1024, 2):
        result = select(n, 8)
        if result.status == 'REDUCED':
            child = min(result.dependencies)
            certificate = discharge(result, codec.encode_source(child))
            need(certificate == codec.encode_source(n), 'independent canonical-route comparison')
            compiled += 1
            common += result.kind == 'common_future'
    print('Supplied-child route AST compilations:', compiled, '; common-future:', common)

    print('Seed-only ascending closure: odd1..1023; sole initial certificate1; no orbit-derived children')
    for budget in (1, 2, 4, 8):
        totals, sets = {}, {}
        for policy in ('inherited', 'adaptive'):
            certs, pending, unresolved = seed_closure(1023, budget, policy)
            sets[policy] = set(certs)
            totals[policy] = dict(certified_including_root=len(certs), pending=len(pending),
                                  reduced_but_child_unavailable=len(unresolved))
        gains = sorted(sets['adaptive']-sets['inherited'])
        losses = sorted(sets['inherited']-sets['adaptive'])
        print('  budget', budget, json.dumps(totals, sort_keys=True),
              '; gains', len(gains), gains[:12], '; losses', len(losses), losses[:12],
              '; adaptive pending first', pending[:15])

    learned = discover_descent(27)
    need(learned is not None, 'bounded least-pending discovery completed')
    p, q, b = compose(learned.word)
    print('Learned from least pending27: odd length/valuation', len(learned.word), sum(learned.word),
          '; exact cylinder', learned.residue, 'mod', learned.modulus,
          '; least endpoint', learned.endpoint_at_residue)
    for k in (0, 1, 2, 19, 100):
        n = learned.residue+learned.modulus*k
        endpoint = macro('learned_cylinder', n, learned.word).endpoint
        need(endpoint == learned.endpoint_at_residue+2*p*k < n, 'independent learned-cylinder lifts')
    for budget in (1, 2, 4, 8):
        base, _, _ = seed_closure(1023, budget, 'adaptive')
        added, pending, unresolved = seed_closure(1023, budget, 'adaptive', (learned,))
        gains, losses = sorted(set(added)-set(base)), sorted(set(base)-set(added))
        print('  learned/root-only budget', budget, 'certified', len(added), 'pending', len(pending),
              'unavailable child', len(unresolved), '; gains', len(gains), gains[:16],
              '; losses', len(losses), losses[:12])
    print('Certified suffix memory: same increasing pass and macro budget; exact parent/cut provenance')
    for budget in (1, 2, 4, 8):
        before, p0, u0 = seed_closure(1023, budget, 'adaptive', suffix_memory=True)
        after, p1, u1 = seed_closure(1023, budget, 'adaptive', (learned,), suffix_memory=True)
        print('  budget', budget, 'without learned', len(before), len(p0), len(u0),
              '; with learned', len(after), len(p1), len(u1),
              '; remaining uncertified first', sorted(p1+u1)[:16])
    rules = []
    print('Four bounded refinement rounds:budget8, root1 only, suffix memory, discovery cap128')
    for _ in range(4):
        memory_stats = {}
        current, pending, unresolved = seed_closure(1023, 8, 'adaptive', tuple(rules), True,
                                                     statistics=memory_stats)
        missing = sorted(set(range(1, 1024, 2))-set(current))
        need(bool(missing), 'this declared refinement round has a residual')
        point = missing[0]
        row = discover_descent(point, cap=128)
        need(row is not None, 'least residual admits checked whole-cylinder rule')
        rules.append(row)
        current, pending, unresolved = seed_closure(1023, 8, 'adaptive', tuple(rules), True,
                                                     statistics=memory_stats)
        print('  learned', point, 'length', len(row.word), 'valuation', sum(row.word),
              'least endpoint', row.endpoint_at_residue, '; certified', len(current),
              '; pending', pending, '; unavailable child', unresolved)
    need(len(current) == 512 and not pending and not unresolved, 'complete declared finite universe')
    print('Final logical proof-memory counts:', json.dumps(memory_stats, sort_keys=True),
          '; node count includesroot; separate count is odd edges, not bytes')
    for n, cert in current.items():
        need(codec.encode_source(n) == cert, 'independent final canonical route audit, after closure')

    import virtual_contraction_ladders_20261004 as virtual
    small = virtual_selection(79)
    need(small.dependencies == (67,) and small.word == (1, 1, 1, 3), 'explicit79 virtual rule')
    need(all(x > 79 for x in replay(79, small.word)[1:]), 'all actual79 prefixes grow')
    # Completed families supply their child AST directly, with no orbit search.
    transferred = 0
    for s in range(1, 5):
        source_ast, child_ast, b = virtual.completed_certificate(s, 0)
        cap = codec.bit_bounds(source_ast)[1]
        if cap <= 100000:
            n = codec.expand(source_ast, bit_cap=cap)
            chosen = virtual_selection(n)
            need(discharge(chosen, child_ast, bit_cap=cap+100) == source_ast,
                 'completed virtual family compiles without child orbit search')
            transferred += 1
    row = virtual.chart(5)
    n = row.r+row.modulus
    large = virtual_selection(n)
    need(len(large.macros) == 1 and len(large.word) == row.j+1 > 16,
         'one macro is not one odd-time unit')
    need(virtual_selection(n, word_cap=16) is None, 'declared macro expansion cap')
    print('Explicit virtual79:', small.dependencies, 'join', small.endpoint,
          '; completed-child AST handoffs:', transferred,
          '; s5 macro count/actual odd horizon:', len(large.macros), len(large.word))

    import mod19_route_lifts_20261004 as rays
    recognized = aliases = 0
    for ray in rays.BANK:
        coefficient, offset, denominator = ray.affine_data()
        for parameter in range(10):
            n = (coefficient*2**(18*parameter)-offset)//denominator
            cert = lookup_ray_certificate(n)
            need(cert == rays.certificate_for(ray, parameter) and codec.expand(cert) == n,
                 'supplied-source membership returns its actual home certificate')
            recognized += 1
            need((n+38) % 19 == n % 19 and lookup_ray_certificate(n+38) is None,
                 'same residue is not exact ray membership')
            aliases += 1
    print('Optional exact mod19-ray lookup:', recognized, 'certified supplied members;',
          aliases, 'same-residue nonmembers rejected; excluded from policy tables')

    # Restricted supplied-child test: keeping only the smallest sibling loses181.
    need(siblings(725) == (181, 45, 11), 'all intermediate sibling heights')
    first = literal_macro('step', 483, 1)
    join = literal_macro('sibling_join', 725, 1)
    result = finish(483, 'REDUCED', 'common_future', join.endpoint, siblings(725),
                    (first, join), ())
    need(discharge(result, codec.encode_source(181)) == codec.encode_source(483),
         'intermediate181 discharges483')
    need(181 != min(result.dependencies), 'minimum-only grammar loses this supplied proof')
    retained = select(893, 1)
    need(retained.endpoint == 335 and retained.dependencies == (223, 335),
         'smaller sibling must not erase already smaller actual join')
    need(discharge(retained, codec.encode_source(335)) == codec.encode_source(893),
         'actual join remains a usable discharge alternative')

    failures = 0
    for bad in (True, 1.0, 0, -1, 2):
        try:
            select(bad)
        except ValueError:
            failures += 1
        else:
            raise ValueError('invalid input accepted')
    need(select(1, 0).status == 'ROOT' and select(27, 0).status == 'PENDING', 'root versus zero budget')
    need(discover_descent(9) is None and discover_descent(17) is None,
         'failed root boundary returns an undispatched obligation')
    try:
        discharge(select(27, 0), codec.ROOT)
    except ValueError:
        failures += 1
    else:
        raise ValueError('pending search discharged')
    # Identical mod19 phase cannot replace the dyadic guard.
    need(observer.chart112_observer(7, 1)['coordinate'] == observer.chart112_observer(45, 1)['coordinate'],
         'same19 phase')
    need(repeated112(7, 1) is not None and repeated112(45, 1) is None, 'dyadic guard distinguishes phase aliases')
    print('Hostiles:483 intermediate181,893 retains223/335, root/zero budget,7/45 mod19 alias;',
          failures, 'rejected inputs/discharges')
    print('No adaptive-policy completeness, new basin, runtime speedup, or unbounded coverage claim.')


if __name__ == '__main__':
    main()
