#!/usr/bin/env python3
"""Independent order-eleven finite audit: raw forward-pair Gaussian equations.

The only producer data read are found.json's seven literal certificates.
Every source and every old word is reconstructed independently from the
published source definition and displayed base ledgers in LOW_ORDERS_REVIEW.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product, permutations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'found.json'

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def pair(a, b):
    return (min(a, b), max(a, b))

def path(q):
    return {pair(a, b) for a, b in zip(q, q[1:])}

class Linear:
    def __init__(self, n):
        self.n = n
        self.pairs = list(combinations(range(n), 2))
        self.index = {p: i for i, p in enumerate(self.pairs)}
        self.m = len(self.pairs)
        self.one = 1 << self.m
        self.mask = self.one - 1
        self.basis = {}
        self.equations = []
        self.ok = True

    def copy(self):
        c = Linear(self.n)
        c.basis = self.basis.copy()
        c.equations = self.equations[:]
        c.ok = self.ok
        return c

    def arc(self, a, b, flips=frozenset()):
        p = pair(a, b)
        return (1 << self.index[p]) ^ (self.one if (a > b) ^ (p in flips) else 0)

    def reduce(self, row):
        for pivot in sorted(self.basis, reverse=True):
            if row & (1 << pivot):
                row ^= self.basis[pivot]
        return row

    def add(self, row):
        self.equations.append(row)
        row = self.reduce(row)
        if row & self.mask:
            self.basis[(row & self.mask).bit_length() - 1] = row
        elif row:
            self.ok = False

    def force(self, a, b, value):
        self.add(self.arc(a, b) ^ (self.one if value else 0))

    def anti(self, q, flips=frozenset()):
        r = dict(zip(q, reversed(q)))
        for a, b in self.pairs:
            self.add(self.arc(a, b, flips) ^ self.arc(r[b], r[a], flips))

    def dimension(self):
        need(self.ok, 'dimension of inconsistent family')
        return self.m - len(self.basis)

    def models(self):
        need(self.ok, 'models of inconsistent family')
        free = [v for v in range(self.m) if v not in self.basis]
        for choice in range(1 << len(free)):
            value = sum(((choice >> j) & 1) << v for j, v in enumerate(free))
            for pivot, row in sorted(self.basis.items()):
                parity = ((row & self.mask & value).bit_count() & 1) ^ bool(row & self.one)
                if parity:
                    value |= 1 << pivot
            need(all(((row & self.mask & value).bit_count() & 1) == bool(row & self.one)
                     for row in self.equations), 'full raw equation solution')
            rows = [0] * self.n
            for i, (a, b) in enumerate(self.pairs):
                u, v = (a, b) if value >> i & 1 else (b, a)
                rows[u] |= 1 << v
            yield tuple(rows)

def source(n, incoming, opened, chord):
    s = Linear(n)
    s.anti(tuple(range(n)))
    for i in range(n - 1):
        s.force(i, i + 1, 1)
    s.force(n - 1, 0, 1)
    for i in range(1, n - 1):
        s.force(0, i, i not in incoming)
    s.force(1, n - 2, opened)
    s.force(2, n - 2, chord)
    need(s.ok, ('source inconsistent', n, incoming, opened, chord))
    return s

def old_words(n, I, opened, chord):
    a, b, x, y, c, z = 0, 1, 2, n-3, n-2, n-1
    M = list(range(2, c))
    if opened:
        if y not in I:
            return 'open_release', [[b]+M+[c,z,a], [z]+M+[c,b,a], [a,z]+M+[c,b]]
        if c not in I:
            return 'open_penultimate', [[a,b,c,z]+M, [b,c,a,z]+M, [b]+M+[c,a,z]]
        h = next(iter(I-{y,c}))
        cut = y-int(h == y-1)
        L, R = list(range(2,cut)), list(range(cut,c))
        return 'open_endrun', [[a,b,c,z]+L+R, [b,c,z]+L+[a]+R, [b]+L+[a]+R+[c,z]]
    if c in I:
        need(y not in I, 'closed-endpair family is separate')
        cut = min({z-v for v in I}-{1})
        L, R = list(range(1,cut)), list(range(cut,c))
        return 'closed_last_movable', [[a]+L+[z]+R+[c], L+[z]+R+[a,c], L+R+[a,c,z]]
    if not chord:
        return 'closed_zero', [[c]+M+[b,z,a], [z,c]+M+[b,a], [a,z,b]+M+[c]]
    if x not in I:
        return 'closed_one', [[c,b,z,a]+M, [z,c,b,a]+M, [a,z,b]+M+[c]]
    if y not in I:
        cut = min(z-v for v in I if z-v > 2)
        name = 'closed_second_movable'
    elif 3 not in I:
        cut, name = y, 'closed_second_fixed'
    else:
        cut, name = n-4, 'closed_second_tail4'
    L, R = list(range(2,cut)), list(range(cut,c))
    return name, [[a,b]+L+[z]+R+[c], [a]+L+[z]+R+[b,c], [a]+L+R+[b,c,z]]

def symbolic_ledger(s, words):
    P = list(range(s.n))
    flips = set()
    blue_families = []
    for q in words:
        need(sorted(q) == P, 'symbolic word is a permutation')
        for a, b in zip(q, q[1:]):
            need(s.reduce(s.arc(a,b,flips)) == s.one, ('unforced join',q,a,b))
        t = s.copy()
        t.anti(q, flips)
        blue_families.append(t)
        flips ^= set(s.pairs) ^ path(q)
    need(flips == set(s.pairs) ^ path(P), 'complete symbolic labelled K target')
    need(path(words[0]) ^ path(words[1]) ^ path(words[2]) == path(P), 'edge XOR')
    return blue_families

def endpair_subcells(s, h, chord):
    """Literal complete complementary splits from the inherited h-family proofs.

    Return (name, guards, words). Guards branch only on genuine source arcs;
    each leaf's cardinality is checked below against the unsplit source cell.
    """
    n = s.n
    a,b,x,y,c,z = 0,1,2,n-3,n-2,n-1
    D = list(range(3,n-3))
    U = list(range(4,n-3))
    def forward(cut):
        L,R = list(range(3,cut)),list(range(cut,n))
        return [[a,b]+L+R+[x],[b]+L+[a]+R+[x],[b,x]+L+[a]+R]
    def backward(start):
        V = list(range(start,n))
        tail = list(range(3,start))
        return [[a]+V+[b,x]+tail,[x,a]+V+[b]+tail,[x,a,b]+tail+V]
    if h == 2:
        if chord:
            return [('h2_xc',[],[[a,b,x]+D+[z,y,c],
                                  [a,x]+D+[z,y,b,c],[a,x]+D+[y,b,c,z]])]
        return [('h2_cx_3c',[(3,c,1)],[[z,a,b,x]+D+[y,c],
                    [z,a,x]+D+[b,y,c],[a,x]+D+[b,y,c,z]]),
                ('h2_cx_c3',[(3,c,0)],[[b,x,a,3]+U+[y,c,z],
                    [b,x,a,3,c,z]+U+[y],[a,b,x,3,c,z]+U+[y]])]
    if h == 3:
        W = list(range(4,n-3))
        R = [y,c,z]
        return [('h3_forward',[(1,3,1)],forward(y)),
                ('h3_backward_4b',[(1,3,0),(1,4,0)],backward(4)),
                ('h3_backward_b4',[(1,3,0),(1,4,1)],
                 [[a,b]+W+R+[x,3],[b]+W+[a]+R+[x,3],[b,x,3]+W+[a]+R])]
    return [(f'h{h}_forward',[(1,3,1)],forward(h)),
            (f'h{h}_backward',[(1,3,0)],backward(5 if h == 4 else 4))]

def arc(rows, a, b):
    return (rows[a] >> b) & 1

def valid(rows):
    n = len(rows)
    need(all(0 <= row < 1 << n for row in rows), 'row bit range')
    need(all(not arc(rows,a,a) for a in range(n)), 'no loops')
    need(all(arc(rows,a,b) ^ arc(rows,b,a) == 1 for a,b in combinations(range(n),2)), 'tournament')

def mismatch(rows, q):
    r = dict(zip(q, reversed(q)))
    for a, b in combinations(range(len(rows)),2):
        if arc(rows,a,b) != arc(rows,r[b],r[a]):
            return {'ordered_pair':[a,b], 'reflected_ordered_pair':[r[b],r[a]],
                    'arc_values':[arc(rows,a,b),arc(rows,r[b],r[a])]}
    return None

def move(rows, q):
    E = path(q)
    n = len(rows)
    return tuple(sum((arc(rows,a,b) ^ (pair(a,b) not in E)) << b
                     for b in range(n) if a != b) for a in range(n))

def twins(rows):
    n = len(rows)
    return [[a,b] for a,b in combinations(range(n),2)
            if all(arc(rows,a,v) == arc(rows,b,v) for v in range(n) if v not in (a,b))]

def literal_source(rows, I, opened, chord):
    n = len(rows)
    valid(rows)
    need(all(arc(rows,v,v+1) for v in range(n-1)), 'source path arcs')
    need(arc(rows,n-1,0) and mismatch(rows,list(range(n))) is None, 'source closed blue')
    need({v for v in range(1,n-1) if arc(rows,v,0)} == set(I), 'exact source inports')
    need({v for v in range(n-1) if arc(rows,n-1,v)} == {0}|{n-1-v for v in I}, 'exact reflected outports')
    need(rows[0].bit_count() == n-5 and rows[-1].bit_count() == 4, 'endpoint scores')
    need(arc(rows,1,n-2) == opened and arc(rows,2,n-2) == chord, 'literal branch bits')

def literal_ledger(rows, words, all_black):
    n = len(rows)
    target = move(rows,list(range(n)))
    trace, flags, xor = [], [], set()
    for q in words:
        need(sorted(q) == list(range(n)), 'literal word permutation')
        need(all(arc(rows,a,b) for a,b in zip(q,q[1:])), 'literal directed joins')
        witness = mismatch(rows,q)
        if all_black:
            need(witness is not None, 'actual word must be black')
        flags.append(witness is None)
        trace.append({'word':q,'source_rows':rows,'black_witness':witness,
                      'actual_first_inports':[v for v in range(n) if v != q[0] and arc(rows,v,q[0])],
                      'actual_last_outports':[v for v in range(n) if v != q[-1] and arc(rows,q[-1],v)],
                      'literal_twin_pairs':twins(rows)})
        after = move(rows,q)
        valid(after)
        need((mismatch(after,q) is None) == (witness is None), 'blue invariance across own K edge')
        rows = after
        xor ^= path(q)
    need(xor == path(list(range(n))) and rows == target, 'literal complete target and XOR')
    return flags, trace, target

def solver_controls():
    """Exhaustively cross-check this elimination/decoder on small raw systems."""
    state = 104729
    consistent = inconsistent = 0
    for case in range(96):
        s = Linear(4)
        rows = []
        for _ in range(case % 13):
            state = (1103515245 * state + 12345) & 0x7fffffff
            row = (state >> 8) & ((1 << 7)-1)
            rows.append(row)
            s.add(row)
        brute = {value for value in range(64)
                 if all(((row & 63 & value).bit_count() & 1) == bool(row & 64) for row in rows)}
        need(s.ok == bool(brute), 'brute Gaussian consistency control')
        if s.ok:
            decoded = {sum(arc(A,a,b) << i for i,(a,b) in enumerate(s.pairs)) for A in s.models()}
            need(decoded == brute and len(decoded) == 1 << s.dimension(), 'brute exhaustive model control')
            consistent += 1
        else:
            inconsistent += 1
    need(consistent and inconsistent, 'both solver control polarities')
    return {'raw_systems':96,'brute_assignments_per_system':64,
            'consistent_systems':consistent,'inconsistent_systems':inconsistent}

def main():
    controls = solver_controls()
    data_bytes = INPUT.read_bytes()
    data = json.loads(data_bytes)
    table = {tuple(item['source_rows']):item for item in data}
    need(len(table) == len(data) == 7, 'seven distinct supplied literal repair sources')
    counts, dimensions = Counter(), Counter()
    failures, endpairs, replay = [], [], []
    exhausted = set()
    n = 11

    def audit_leaf(s, ports, opened, chord, name, qs, kind, guards):
        blue = symbolic_ledger(s,qs)
        modelsets = [set(t.models()) if t.ok else set() for t in blue]
        union = set().union(*modelsets)
        if not union:
            counts[kind+'_uniform_success'] += 1
            return
        counts[kind+'_failed_leaves'] += 1
        failures.append({'ports':ports,'middle_open':opened,'x_to_c':chord,
                         'old_ledger':name,'kind':kind,'guards':guards,
                         'source_dimension':s.dimension(),
                         'blue_loci':[{'step':i+1,'dimension':t.dimension(),'sources':len(ms)}
                                      for i,(t,ms) in enumerate(zip(blue,modelsets)) if t.ok],
                         'distinct_blue_union_sources':len(union)})
        for rows in sorted(union):
            need(rows in table, ('missing exhaustive replacement',ports,name,rows))
            need(rows not in exhausted, 'each literal source belongs to exactly one source leaf')
            item = table[rows]
            need(item['n'] == n and tuple(item['ports']) == ports
                 and item['middle_open'] == opened and item['x_to_c'] == chord,
                 'all literal port/branch labels match independent family')
            literal_source(rows,set(ports),opened,chord)
            for a,b,value in guards:
                need(arc(rows,a,b) == value, 'literal leaf guard')
            old_flags,_,_ = literal_ledger(rows,qs,False)
            need(old_flags == [rows in ms for ms in modelsets], 'actual old blue flags match exact loci')
            flags,trace,target = literal_ledger(rows,item['words'],True)
            need(flags == [False]*3, 'three actual black repairs')
            replay.append({'ports':ports,'source_rows':rows,'old_ledger':name,
                           'old_blue_flags':old_flags,'old_words':qs,'new_steps':trace,
                           'literal_target_rows':target,'target_twin_pairs':twins(target)})
            exhausted.add(rows)
            counts[kind+'_distinct_sources_repaired'] += 1

    for ports in combinations(range(2,n-1),3):
        I = set(ports)
        for opened,chord in product((False,True),repeat=2):
            s = source(n,I,opened,chord)
            dimensions[s.dimension()] += 1
            counts['all_source_cells'] += 1
            if not opened and {n-3,n-2} <= I:
                counts['closed_endpair_cells'] += 1
                h = next(iter(I-{n-3,n-2}))
                leaves = []
                for name,guards,qs in endpair_subcells(s,h,chord):
                    t = s.copy()
                    for a,b,value in guards:
                        t.force(a,b,value)
                    need(t.ok, ('nonempty exact endpair leaf',h,chord,name))
                    counts['endpair_leaves'] += 1
                    audit_leaf(t,ports,opened,chord,name,qs,'endpair',guards)
                    leaves.append({'name':name,'guards':guards,'dimension':t.dimension(),'words':qs})
                need(sum(1 << leaf['dimension'] for leaf in leaves) == 1 << s.dimension(),
                     'full endpair cardinality coverage')
                for first,second in combinations(leaves,2):
                    t = s.copy()
                    for a,b,value in first['guards']+second['guards']:
                        t.force(a,b,value)
                    need(not t.ok, 'disjoint endpair subleaves')
                endpairs.append({'ports':ports,'x_to_c':chord,'source_dimension':s.dimension(),'leaves':leaves})
            else:
                counts['base_cells'] += 1
                name,qs = old_words(n,I,opened,chord)
                audit_leaf(s,ports,opened,chord,name,qs,'base',[])
    need(exhausted == set(table), 'the seven proposed sources equal the entire blue-conditioned union')
    need(dimensions == {14:224}, 'unfiltered raw source dimensions')
    expected = {'all_source_cells':224,'closed_endpair_cells':12,'endpair_leaves':25,
                'endpair_uniform_success':23,'endpair_failed_leaves':2,'endpair_distinct_sources_repaired':3,
                'base_cells':212,'base_uniform_success':208,'base_failed_leaves':4,'base_distinct_sources_repaired':4}
    need(dict(counts) == expected, ('full order-eleven partition',dict(counts)))
    need(all(locus['dimension'] == 0 for failure in failures for locus in failure['blue_loci']),
         'all old blue loci are singletons')
    actual_failures = {(tuple(f['ports']),f['x_to_c'],f['old_ledger'],tuple(q['step'] for q in f['blue_loci'])) for f in failures}
    expected_failures = {((3,4,7),True,'closed_one',(3,)),((3,5,7),True,'closed_one',(3,)),
                         ((4,6,7),False,'closed_zero',(3,)),((5,6,7),False,'closed_zero',(3,)),
                         ((6,8,9),True,'h6_backward',(2,)),((7,8,9),False,'h7_forward',(2,3))}
    need(actual_failures == expected_failures, 'exact base and residual failures')
    clone_source = (126,156,440,1008,1248,1474,1926,1793,1555,1079,15)
    clone_record = next(item for item in replay if tuple(item['source_rows']) == clone_source)
    need([step['literal_twin_pairs'] for step in clone_record['new_steps']] == [[],[],[[4,5]]],
         'exact formation of transient literal clone pair')
    clone_state = clone_record['new_steps'][2]['source_rows']
    need(arc(clone_state,5,4) and not clone_record['target_twin_pairs'],
         'ordered clone pair (5,4) splits at exact target')
    # A completely exhausted known hostile must remain outside the claimed
    # threshold; it prevents misreading a failed template as the only obstacle.
    hostile = (6,4,8,19,103,79,15)
    literal_source(hostile,{3,4,5},False,False)
    H = move(hostile,list(range(7)))
    hp = [q for q in permutations(range(7)) if all(arc(H,a,b) for a,b in zip(q,q[1:]))]
    need(len(hp) == 9 and all(mismatch(H,q) is None for q in hp), 'pure-blue hostile retained')
    print(json.dumps({'status':'PASS','method':'independent raw forward-pair Gaussian; exact source partitions and complete singleton union; actual arc replay',
                      'input_sha256':sha256(data_bytes).hexdigest(),
                      'independent_solver_controls':controls,'counts':dict(sorted(counts.items())),
                      'source_dimensions':dict(dimensions),'literal_marked_sources':224*(1<<14),
                      'failed_old_leaves':failures,'endpair_partition':endpairs,'literal_repairs':replay,
                      'transient_clone_control':{'source':clone_source,'state_after_two_moves':clone_state,
                         'ordered_clone_block':[5,4],'target':clone_record['literal_target_rows'],
                         'target_twin_pairs':clone_record['target_twin_pairs']},
                      'order7_hostile':{'source':hostile,'target':H,'permutations_tested':5040,'directed_paths':9,'black_paths':0},
                      'scope':'FINITE-EXACT full marked directed closed blue endpoint-four closure at order11; combined with audited n>=12 yields n>=11. Optimal threshold and whole black-core connectivity remain open.'},indent=2,sort_keys=True))

if __name__ == '__main__':
    main()

