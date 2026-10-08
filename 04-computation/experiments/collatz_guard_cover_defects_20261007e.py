"""Exact dyadic guard covers, finite head cuts and alternative proof defects.

Rule names are provenance labels, never ROOT premises. Literal graph nodes
are independently authenticated by the inherited first-hit consumer.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations
from collections import deque, Counter

import collatz_marked_completion_20261007c as proof


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def natural(n, minimum=0):
    need(type(n) is int and n >= minimum, "exact natural in the declared domain")


@dataclass(frozen=True, order=True)
class Cell:
    bits: int
    residue: int


def cell(c):
    need(type(c) is Cell, "exact dyadic Cell")
    natural(c.bits); natural(c.residue)
    need(c.residue < 1 << c.bits, "canonical dyadic residue")
    return c


def contains(c, parameter):
    cell(c); natural(parameter)
    return parameter % (1 << c.bits) == c.residue


def inside(small, large):
    cell(small); cell(large)
    return small.bits >= large.bits and contains(large, small.residue)


def count(c, lower, upper):
    cell(c); natural(lower); natural(upper)
    need(lower <= upper, "ordered finite interval")
    m = 1 << c.bits
    return (upper-1-c.residue)//m - (lower-1-c.residue)//m


def first(c, lower):
    cell(c); natural(lower)
    return lower+(c.residue-lower) % (1 << c.bits)


@dataclass(frozen=True)
class Rule:
    name: str
    guard: Cell
    minimum: int
    obligation: str


def rule(r):
    need(type(r) is Rule, "exact guard Rule")
    need(type(r.name) is str and bool(r.name), "named rule")
    need(type(r.obligation) is str and bool(r.obligation), "marked terminal obligation")
    cell(r.guard); natural(r.minimum)
    return r


@dataclass(frozen=True)
class Bank:
    chart: str
    rules: tuple


def bank(chart, rules):
    need(type(chart) is str and bool(chart), "source-owned parameter chart")
    need(type(rules) is tuple, "tuple of rules")
    for r in rules: rule(r)
    need(len({r.name for r in rules}) == len(rules), "rule names are unique")
    return Bank(chart, tuple(sorted(rules, key=lambda r: r.name)))


def audit_bank(b):
    need(type(b) is Bank and bank(b.chart, b.rules) == b, "canonical Bank")
    return b


def extend(b, addition):
    audit_bank(b); audit_bank(addition)
    need(b.chart == addition.chart, "do not mix source parameter charts")
    old = {r.name: r for r in b.rules}
    for r in addition.rules:
        need(r.name not in old or old[r.name] == r, "do not silently redefine a rule")
        old[r.name] = r
    return bank(b.chart, tuple(old.values()))


def available(b, parameter):
    audit_bank(b); natural(parameter)
    return tuple(r.name for r in b.rules
                 if parameter >= r.minimum and contains(r.guard, parameter))


@dataclass(frozen=True)
class Piece:
    guard: Cell
    names: tuple


def partition(rules):
    """All-label partition, with at most 1+sum(bits) trie leaves."""
    b = bank("internal partition", rules)
    result = []
    stack = [(Cell(0, 0), (), b.rules)]
    while stack:
        current, inherited, pending = stack.pop()
        here = tuple(r.name for r in pending if r.guard == current)
        names = tuple(sorted(inherited+here))
        descendants = tuple(r for r in pending if r.guard != current)
        if not descendants:
            result.append(Piece(current, names))
            continue
        for bit in (1, 0):
            child = Cell(current.bits+1, current.residue+(bit << current.bits))
            below = tuple(r for r in descendants if inside(r.guard, child))
            stack.append((child, names, below))
    return tuple(sorted(result, key=lambda p: p.guard))


@dataclass(frozen=True)
class Stratum:
    lower: int
    upper: object
    pieces: tuple


def strata(b):
    """Exact on nonnegative integer parameters, including every finite cut."""
    audit_bank(b)
    cuts = sorted({0} | {r.minimum for r in b.rules})
    result = []
    for i, lower in enumerate(cuts):
        upper = cuts[i+1] if i+1 < len(cuts) else None
        active = tuple(r for r in b.rules if r.minimum <= lower)
        result.append(Stratum(lower, upper, partition(active)))
    return tuple(result)


def tail_mass(b):
    return sum((F(1, 1 << p.guard.bits) for p in strata(b)[-1].pieces if p.names), F())


def new_mass(b, addition):
    return tail_mass(extend(b, addition))-tail_mass(b)


def uncovered_count(b, lower, upper):
    natural(lower); natural(upper)
    need(lower <= upper, "ordered finite interval")
    total = 0
    for s in strata(b):
        lo = max(lower, s.lower)
        hi = min(upper, s.upper) if s.upper is not None else upper
        if lo < hi:
            total += sum(count(p.guard, lo, hi) for p in s.pieces if not p.names)
    return total


def least_uncovered(b):
    for s in strata(b):
        candidates = [first(p.guard, s.lower) for p in s.pieces if not p.names]
        if candidates:
            answer = min(candidates)
            if s.upper is None or answer < s.upper:
                return answer
    return None


def verify_partition(b, claimed):
    """Canonical exact witness; equality follows only after exact type checks."""
    audit_bank(b)
    need(type(claimed) is tuple, "tuple of strata")
    for s in claimed:
        need(type(s) is Stratum, "exact Stratum")
        natural(s.lower)
        if s.upper is not None: natural(s.upper)
        need(type(s.pieces) is tuple, "tuple of pieces")
        for p in s.pieces:
            need(type(p) is Piece and type(p.names) is tuple, "exact labelled piece")
            cell(p.guard)
            need(all(type(n) is str for n in p.names), "exact rule identifiers")
    need(claimed == strata(b), "complete all-label partition and finite cuts")
    return True


@dataclass(frozen=True)
class Network:
    nodes: tuple


def network(nodes):
    """Retain every authenticated alternative at each explicitly named port."""
    need(type(nodes) is tuple, "tuple of literal proof nodes")
    for n in nodes: proof.audit(n)
    ports = {n.source for n in nodes}
    roots = {}
    for n in nodes:
        if type(n) is proof.Root:
            need(n.source not in roots or roots[n.source] == n, "unique actual ROOT word")
            roots[n.source] = n
        if type(n) is proof.Transport:
            need(n.child in ports, "retain an explicit child port")
    need(roots.get(1) == proof.Root(1, ()), "literal ROOT boundary")
    return Network(tuple(sorted(set(nodes), key=repr)))


def audit_network(net):
    need(type(net) is Network and network(net.nodes) == net, "canonical Network")
    return net


def grounded(net):
    """Least reachability over both lawful directions of every common future."""
    audit_network(net)
    selected = {n.source: n for n in net.nodes if type(n) is proof.Root}
    reverse = {}
    for n in net.nodes:
        if type(n) is proof.Transport:
            reverse.setdefault(n.child, []).append(n)
            flipped = proof.Transport(n.child,n.source,n.right,n.left,n.endpoint,'join')
            reverse.setdefault(n.source, []).append(flipped)
    queue = deque(sorted(selected))
    while queue:
        child = queue.popleft()
        for link in reverse.get(child, ()):
            if link.source not in selected:
                selected[link.source] = link
                queue.append(link.source)
    return selected


def defect_cut(net):
    audit_network(net)
    return frozenset(n.source for n in net.nodes)-frozenset(grounded(net))


def verify_cut(net, cut):
    audit_network(net)
    need(type(cut) is frozenset and all(type(n) is int for n in cut), "exact defect ports")
    need(cut <= {n.source for n in net.nodes}, "cut within this literal network")
    for n in net.nodes:
        if type(n) is proof.Root:
            need(n.source not in cut, "no ROOT seed behind the cut")
        if type(n) is proof.Transport:
            need((n.child in cut)==(n.source in cut), "BOTH directions of EVERY common future stay behind the cut")
    return True


def export(net, source):
    proof.odd(source)
    selected = grounded(net)
    if source not in selected:
        return None
    current, path = source, []
    while type(selected[current]) is proof.Transport:
        link = selected[current]
        path.append(link)
        current = link.child
    root = selected[current]
    for link in reversed(path): root = proof.discharge(link, root)
    return root


@dataclass(frozen=True)
class Instance:
    rule_name: str
    parameter: int
    node: object


def consume_instances(b, parameter, source, instances, support=()):
    """Exact pointwise consumer; guard labels alone cannot ground a cylinder."""
    audit_bank(b); natural(parameter); proof.odd(source)
    need(type(instances) is tuple and type(support) is tuple, "tuple of supplied evidence")
    names = available(b, parameter)
    nodes = [proof.Root(1, ())]
    for item in instances:
        need(type(item) is Instance and type(item.rule_name) is str, "exact rule instance")
        natural(item.parameter)
        need(item.parameter == parameter and item.rule_name in names, "same active parameter guard")
        proof.audit(item.node)
        need(type(item.node) in (proof.Root, proof.Transport) and item.node.source == source,
             "same literal source, actual receipt or ROOT word")
        nodes.append(item.node)
    for node in support: proof.audit(node)
    nodes.extend(support)
    ports = {n.source for n in nodes}
    needed = {source} | {n.child for n in nodes if type(n) is proof.Transport}
    nodes.extend(proof.Hole(n, "unsupplied exact child") for n in sorted(needed-ports) if n != 1)
    return network(tuple(nodes))


def forward(source, word):
    endpoint = proof.replay(source, word)
    if endpoint == 1: return proof.audit(proof.Root(source, word))
    return proof.audit(proof.Transport(source, endpoint, word, (), endpoint, "forward"))


def parameter_bank():
    """Authenticate the retained arithmetic rows; leave every child OPEN."""
    import collatz_parameter_cover_20261007e as arithmetic
    rows = arithmetic.load_rules()
    rules = []
    for i, row in enumerate(rows):
        arithmetic.audit_rule(row)
        rules.append(Rule('bank'+str(i).zfill(4),
            Cell(row.parameter_bits, row.parameter_residue), row.minimum_parameter,
            'OPEN M_(924745897+2^32*t-'+str(row.drop)+')'))
    new = bank('M_(924745897+2^32*t), t>=0', tuple(rules))
    old = bank(new.chart, (Rule('old D8',Cell(47,0),0,'OPEN M_(924745889+2^32*t)'),))
    return extend(old,new), old, rows


def add_signed(b):
    """A separate extension of the frozen positive bank, with reserve intact."""
    audit_bank(b)
    need(b.chart == 'M_(924745897+2^32*t), t>=0', 'same signed source chart')
    import collatz_signed_parameter_fill_20261007e as arithmetic
    rows, targeted = arithmetic.load()
    rules = []
    for i, row in enumerate(rows+(targeted,)):
        arithmetic.audit(row)
        rules.append(Rule('signed'+str(i).zfill(4),Cell(row.bits,row.residue),row.minimum,
                         'OPEN M_(924745897+2^32*t-'+str(row.drop)+')'))
    return extend(b,bank(b.chart,tuple(rules)))


def main():
    checks = Counter()
    def check(ok, category):
        checks[category] += 1
        if not ok: raise ArithmeticError((category, checks[category]))
    def rejects(fn):
        try: fn()
        except (ValueError, TypeError): check(True, "hostile")
        else: check(False, "hostile")

    # Complete small bank universe, independent direct pointwise comparisons.
    cells = tuple(Cell(b, r) for b in range(4) for r in range(1 << b))
    banks = 0
    for size in range(4):
        for chosen in combinations(cells, size):
            rows = tuple(Rule(str(i), c, (i*3)%5, 'OPEN') for i, c in enumerate(chosen))
            b = bank('small', rows); witness = strata(b)
            check(verify_partition(b, witness), 'partition')
            for s in witness:
                check(sum((F(1, 1 << p.guard.bits) for p in s.pieces), F()) == 1, 'kraft')
                check(len(s.pieces) <= 1+sum(r.guard.bits for r in rows if r.minimum <= s.lower), 'size')
            for t in range(24):
                direct = tuple(r.name for r in b.rules if t >= r.minimum and t%(1 << r.guard.bits)==r.guard.residue)
                matches = [p.names for s in witness if t >= s.lower and (s.upper is None or t < s.upper)
                           for p in s.pieces if contains(p.guard, t)]
                check(matches == [direct] and available(b,t) == direct, 'pointwise')
            check(uncovered_count(b,0,24) == sum(not available(b,t) for t in range(24)), 'finite_count')
            actual = next((t for t in range(24) if not available(b,t)), None)
            check(least_uncovered(b) == actual, 'least')
            check(tail_mass(b) == F(sum(any(t%(1 << r.guard.bits)==r.guard.residue for r in rows)
                                           for t in range(8)),8), 'density')
            banks += 1

    full_tail = bank('finite cut', (Rule('tail',Cell(0,0),1,'OPEN'),))
    check(tail_mass(full_tail)==1 and least_uncovered(full_tail)==0, 'finite_cut')
    for depth in range(1,65):
        b = bank('shrinking residual', tuple(Rule(str(k),Cell(k+1,1 << k),0,'OPEN') for k in range(depth)))
        check(tail_mass(b)==1-F(1,1 << depth) and least_uncovered(b)==0, 'limit')
        check(tuple(p.guard for p in strata(b)[-1].pieces if not p.names)==(Cell(depth,0),), 'limit')

    # Actual no-new-mass proof gain on N(t)=8t+3; final oddness requires t0mod8.
    broad = bank('N(t)=8t+3', (Rule('one',Cell(0,0),0,'endpoint open'),))
    narrow = bank(broad.chart, (Rule('one-four',Cell(3,0),0,'endpoint open'),))
    both = extend(broad,narrow)
    check(new_mass(broad,narrow)==0, 'retained_alternative')
    open_net = consume_instances(broad,0,3,(Instance('one',0,forward(3,(1,))),))
    check(export(open_net,3) is None and defect_cut(open_net)==frozenset((3,5)), 'retained_alternative')
    closed_net = consume_instances(both,0,3,(Instance('one',0,forward(3,(1,))),
        Instance('one-four',0,forward(3,(1,4)))))
    check(export(closed_net,3)==proof.Root(3,(1,4)), 'retained_alternative')
    check(export(closed_net,5)==proof.Root(5,(4,)), 'retained_alternative')
    rejects(lambda: forward(35,(1,4)))
    check(proof.replay(67,(1,4))==19 and available(both,8)==('one','one-four'), 'retained_alternative')
    # A supplied reverse marks local payment; the logical reverse was already lawful.
    reverse = proof.audit(proof.Transport(5,3,(),(1,),5,'paid'))
    cycle = network((proof.Root(1,()),forward(3,(1,)),reverse))
    check(defect_cut(cycle)==frozenset((3,5)) and verify_cut(cycle,defect_cut(cycle)), 'network')
    repaired = network(cycle.nodes+(proof.Root(5,(4,)),))
    check(export(repaired,3)==proof.Root(3,(1,4)) and export(repaired,5)==proof.Root(5,(4,)), 'network')
    rejects(lambda: verify_cut(repaired,frozenset((3,5))))
    # An alternative direct word and a cyclic choice coexist without forcing a bad selection.
    alternative = network(cycle.nodes+(proof.Root(3,(1,4)),))
    check(export(alternative,5)==proof.Root(5,(4,)), 'network')

    rejects(lambda: cell(Cell(True,0)))
    rejects(lambda: cell(Cell(2,4)))
    rejects(lambda: contains(Cell(0,0),True))
    rejects(lambda: bank('x',(Rule('a',Cell(0,0),0,'x'),Rule('a',Cell(1,0),0,'y'))))
    rejects(lambda: extend(broad,bank('changed source',narrow.rules)))
    rejects(lambda: consume_instances(both,1,11,(Instance('one-four',1,forward(3,(1,4))),)))
    rejects(lambda: consume_instances(both,0,3,(Instance('one',0,forward(11,(1,))),)))
    rejects(lambda: consume_instances(both,0,3,(Instance('one',True,forward(3,(1,))),)))
    rejects(lambda: verify_cut(cycle,frozenset((True,3,5))))
    s = strata(full_tail)
    rejects(lambda: verify_partition(full_tail,(replace(s[0],lower=False),)+s[1:]))

    actual_bank, old, rows = parameter_bank()
    actual_strata = strata(actual_bank)
    check(len(rows)==429 and len(actual_strata)==1, 'arithmetic_bank')
    actual_pieces = actual_strata[0].pieces
    covered = tuple(p for p in actual_pieces if p.names)
    residual = tuple(p for p in actual_pieces if not p.names)
    mass = sum((F(1,1 << p.guard.bits) for p in covered),F())
    import collatz_parameter_cover_20261007e as arithmetic
    check(mass==arithmetic.mass(rows), 'arithmetic_bank')
    check(len(actual_pieces)==27254 and len(residual)==26825, 'arithmetic_bank')
    check(max(p.guard.bits for p in actual_pieces)==238, 'arithmetic_bank')
    check(min(p.guard.residue for p in residual)==least_uncovered(actual_bank)==4, 'arithmetic_bank')
    check(next(p.guard for p in residual if contains(p.guard,4))==Cell(9,4), 'arithmetic_bank')
    check(residual[:2]==(Piece(Cell(7,36),()),Piece(Cell(7,116),())), 'arithmetic_bank')
    check(sum(count(p.guard,0,1024) for p in residual)==431, 'arithmetic_bank')
    check(len(covered)==429 and all(p.names for p in covered), 'arithmetic_bank')
    check(next(p.names for p in covered if contains(p.guard,0))==('bank0000','old D8'), 'arithmetic_bank')
    for t in range(1024):
        direct = any(arithmetic._covered((r,),t) for r in rows)
        check(bool(available(actual_bank,t))==direct, 'arithmetic_points')
    # Three independently supplied exact ordinary-source receipts from the new rows.
    # These test the receipt consumer, not ROOT claims about the astronomical family.
    for row in rows[:3]:
        macro = arithmetic.macro(row)
        e = max(2*sum(row.source_head), row.drop+1)
        source = arithmetic.context.fixed_source(macro,e,0)
        receipt = arithmetic.receipt(row,source)
        node = proof.paid(receipt.source, receipt.child, receipt.source_word,
                          receipt.child_word, receipt.endpoint)
        check(type(node) in (proof.Root,proof.Transport), 'actual_receipts')

    combined = add_signed(actual_bank)
    combined_pieces = strata(combined)[-1].pieces
    combined_residual = tuple(p for p in combined_pieces if not p.names)
    combined_mass = sum((F(1,1 << p.guard.bits) for p in combined_pieces if p.names),F())
    import collatz_signed_parameter_fill_20261007e as signed
    signed_rows,targeted = signed.load()
    direct_cells = [(r.parameter_residue,r.parameter_bits) for r in rows]
    direct_cells += [(r.residue,r.bits) for r in signed_rows+(targeted,)]
    check(combined_mass==signed.mass(direct_cells), 'signed_bank')
    check(combined_mass>mass and available(combined,6), 'signed_bank')
    check(least_uncovered(combined)==4, 'signed_bank')
    for p in combined_pieces:
        if p.names:
            check(set(available(combined,p.guard.residue))==set(p.names), 'signed_labels')
    combined_least_cell = next(p.guard for p in combined_residual if contains(p.guard,4))
    check(len(combined_pieces)==38259 and len(combined_residual)==37693, 'signed_bank')
    check(combined_least_cell==Cell(9,4), 'signed_bank')

    print('PROVED all-label dyadic partition and exact finite-cut strata; provenance is not ROOT evidence.')
    print('Small complete universe:',banks,'banks from15cells, up to3rules; parameters0..23.')
    print('PROVED exact Kraft mass, finite interval counts, and least uncovered parameter without common-modulus expansion.')
    print('Hostiles: full asymptotic mass misses finite t0; increasing mass1 limit can still miss t0.')
    print('Actual N(t)=8t+3: nested word14 adds zero guard mass but grounds source3; contained rules cannot be proof-pruned.')
    print('All-alternatives cut checks EVERY dependency; supplied ROOT5 breaks actual3/5 proof cycle.')
    print('Authenticated new bank429rows plus old0/2^47: routed mass',mass)
    print('Additional mass beyond old cell',mass-F(1,1 << 47))
    print('Exact frontier:27254labelled cells,26825uncovered; maximum238bits; least missingt4 liesin4mod512.')
    print('Coarsest uncovered cells36mod128,116mod128; sample431/1024uncovered. Child obligations remain OPEN.')
    print('Separate signed extension137rows: combined mass',combined_mass)
    print('Combined labelled leaves',len(combined_pieces),'uncovered',len(combined_residual),
          'least missing4 cell',combined_least_cell,'sample missing',sum(count(p.guard,0,1024) for p in combined_residual))
    print('Checks',dict(sorted(checks.items())))
    print('Total checks',sum(checks.values()))


if __name__ == '__main__': main()
