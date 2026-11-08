#!/usr/bin/env python3
"""Portable exact audit of scalar switches and the graft D-column lemma.
No producer imports. Every failed condition raises under both Python and -O.
"""
import itertools
import json
from collections import defaultdict
from functools import lru_cache


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def mm(a, b):
    return tuple(tuple(sum(x*y for x, y in zip(row, col)) for col in zip(*b)) for row in a)


def mv(a, v):
    return tuple(sum(x*y for x, y in zip(row, v)) for row in a)


I = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
O = ((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1))
Z = ((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1))
J = ((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))
L = ((1,-1,3,-2),(0,0,1,-1),(2,-1,4,-1),(1,-1,1,0))
C = ((3,-3,9,-6),(1,-1,2,-2),(3,-2,11,-6),(1,-1,3,-2))
OPS = dict(O=O, Z=Z, J=J, L=L, C=C)
GADGETS = dict(J=(0,), L=(2,0), C=(2,4,1))


def wm(word):
    out = I
    for c in word:
        out = mm(OPS[c], out)
    return out


def tournaments(n):
    edges = list(itertools.combinations(range(n), 2))
    for orient in range(1 << len(edges)):
        a = [0]*n
        for i, (u,v) in enumerate(edges):
            if orient >> i & 1:
                a[u] |= 1 << v
            else:
                a[v] |= 1 << u
        yield tuple(a)


def path_tables(a):
    n, full = len(a), (1 << len(a))-1
    ends = [[0]*n for _ in range(full+1)]
    for mask in range(1, full+1):
        for j in range(n):
            if not mask >> j & 1:
                continue
            prev = mask ^ (1 << j)
            ends[mask][j] = 1 if prev == 0 else sum(ends[prev][i] for i in range(n) if prev >> i & 1 and a[i] >> j & 1)
    H = [1] + [sum(row) for row in ends[1:]]
    return H, ends


def start_counts(a, root):
    n, full = len(a), (1 << len(a))-1
    starts = [0]*(full+1)
    @lru_cache(None)
    def f(mask, first):
        rest = mask ^ (1 << first)
        return 1 if not rest else sum(f(rest, j) for j in range(n) if rest >> j & 1 and a[first] >> j & 1)
    for mask in range(1,full+1):
        if mask >> root & 1:
            starts[mask] = f(mask, root)
    return starts


def data(a, root):
    H, ends = path_tables(a)
    start = start_counts(a, root)
    full = (1 << len(a))-1
    U = sum(start[s]*H[full ^ s] for s in range(full+1))
    return (H[full], ends[full][root], U, H[full ^ (1 << root)]), H, ends, start


def profile(a, root):
    return data(a, root)[0]


def submasks(mask):
    s = mask
    while True:
        yield s
        if not s:
            break
        s = (s-1) & mask


def graft_matrix(a, root=0):
    state, paths, ends, start = data(a, root)
    full = (1 << len(a))-1
    two = [sum(paths[s]*paths[t ^ s] for s in submasks(t)) for t in range(full+1)]
    H,e,alpha,q = state
    Q,s,b = two[full], start[full], two[full ^ (1 << root)]
    delta = sum(ends[t][root]*paths[full ^ t] for t in range(full+1))
    T = sum(start[t]*two[full ^ t] for t in range(full+1))
    return ((H,-H,Q-H,2*H-Q),(e,-e,delta-e,e-delta),
            (alpha,s-alpha,2*T-alpha,3*alpha-s-2*T),
            (1,0,0,0) if len(a)==1 else (q,-q,b-q,2*q-b))


def literal(a, root, c):
    n = len(a)
    if c == 'O':
        return tuple(row | (1 << n) for row in a)+(0,), root
    if c == 'Z':
        return tuple(row | ((1 << n) if i != root else 0) for i,row in enumerate(a))+(1 << root,), root
    b = GADGETS[c]
    bm = ((1 << len(b))-1) << n
    out = tuple(row | (bm if i != root else 0) for i,row in enumerate(a))
    out += tuple((row << n) | (1 << root) for row in b)
    return out, n


def apply_graph(a, root, word):
    for c in word:
        a, root = literal(a, root, c)
    return a, root


def main():
    ed = (1,0,0,0)
    check(mv(O,ed)==ed and mv(Z,ed)==ed, 'O and Z fix e_D')
    A, B = wm('JOLO'), wm('JJLO')
    # The E-column is checked below from basis multiplication, not assumed
    # by the E=0 switch theorem.
    restricted = (0,2,3)
    for name, M, scalar in [('JOLO',A,7),('JJLO',B,5)]:
        want = ((5,scalar,0),(0,0,0),(18,20,0),(3,3,0))
        got = tuple(tuple(row[j] for j in restricted) for row in M)
        check(got==want, name+' E=0 operator')
    A2,B2 = wm('JOZJO'),wm('JJJJO')
    for name,M,scalar in [('JOZJO',A2,7),('JJJJO',B2,6)]:
        check(tuple(tuple(row[j] for j in restricted) for row in M)==((5,scalar,0),(0,0,0),(12,14,0),(3,3,0)),name+' E=0 operator')
    for M,alpha,beta in ((A,6,20),(A2,4,14)):
        check(tuple(M[2][j]-alpha*M[3][j] for j in restricted)==(0,2,0), 'exact inverse U row')
        check(tuple(beta*M[3][j]-3*M[2][j] for j in restricted)==(6,0,0), 'exact inverse D row')
    side_det = 18*3-20*3
    check(side_det == -6, 'sidecar determinant')
    for c,b in GADGETS.items():
        check(graft_matrix(b)==OPS[c], c+' literal gadget statistics')
    inputs = [(a,0) for n in range(1,5) for a in tournaments(n)]
    direct_switches = 0
    for a,r in inputs:
        a,r = literal(a,r,'O')
        x = profile(a,r)
        check(x[1]==0, 'O creates E=0 input')
        for word in ('JOLO','JJLO','JOZJO','JJJJO'):
            b,s = apply_graph(a,r,word)
            check(profile(b,s)==mv(wm(word),x), 'literal '+word)
            direct_switches += 1
    column_gadgets = column_suffixes = 0
    for a,r in inputs:
        K = graft_matrix(a,r)
        check(mv(K,ed)==profile(a,r), 'graft D-column is gadget profile')
        column_gadgets += 1
        for c in ('','O','Z','J','L','C'):
            b,s = apply_graph(a,r,c)
            y = profile(b,s)
            check(mv(wm(c),mv(K,ed))==y, 'physical suffix D-column')
            check(y[0]%2==1, 'odd physical Hamilton count')
            check(y[2]>=1 and y[3]>=1, 'positive sidecars after graft suffix')
            column_suffixes += 1
    all_words = nontrivial = 0
    for length in range(6):
        for letters in itertools.product(OPS, repeat=length):
            word = ''.join(letters)
            col = mv(wm(word),ed)
            all_words += 1
            first = next((i for i,c in enumerate(word) if c in GADGETS),None)
            if first is None:
                check(col==ed, 'graft-free suffix preserves pure-D difference')
            else:
                a = GADGETS[word[first]]
                target = mv(wm(word[first+1:]),profile(a,0))
                check(col[0]%2==1, 'odd D-column multiplier')
                check(col==target and col[2]>0 and col[3]>0, 'D-column suffix identity')
                nontrivial += 1
    fibers = []
    states = {(1,0,2,1)}
    for depth in range(1,13):
        states = {mv(M,x) for x in states for M in (A,B)}
        by_sidecar = defaultdict(set)
        for D,E,U,d in states:
            by_sidecar[(E,U,d)].add(D)
            check(E==0 and d%3==0, 'switch output face and congruence')
            check(D in (5*d//3, U-13*d//3), 'sidecars explicitly recover two scalars')
            check(U-6*d>0, 'strict positive switch gap')
        check(max(map(len,by_sidecar.values()))<=2, 'fiber width')
        fibers.append(dict(depth=depth,states=len(states),sidecars=len(by_sidecar),max_width=max(map(len,by_sidecar.values()))))
    mixed_fibers=[]
    mixed={(1,0,2,1)}
    for depth in range(1,8):
        mixed={mv(M,x) for x in mixed for M in (A,B,A2,B2)}
        groups=defaultdict(set)
        for D,E,U,d in mixed:
            groups[(E,U,d)].add(D)
            check(U!=6*d, 'disjoint pair sidecar ranges')
            options=(5*d//3,U-13*d//3) if U>6*d else (U//2-d//3,U-7*d//3)
            check(D in options, 'mixed fiber branch detected by sidecars')
        check(max(map(len,groups.values()))<=2, 'mixed all-four fiber bound')
        mixed_fibers.append(dict(depth=depth,states=len(mixed),sidecars=len(groups),max_width=max(map(len,groups.values()))))
    # Independent exhaustive short-transfer discovery: retain all exact-length
    # matrices, rather than the producer's global shortest-word pruning.
    E0=tuple(tuple(int(i==j) for j in restricted) for i in range(4))
    level={E0}
    transfers=set()
    distinct_prefixes=set()
    first_seen=[]
    for depth in range(8):
        first_seen.append(len(level-distinct_prefixes))
        distinct_prefixes.update(level)
        transfers.update(mm(O,M) for M in level)
        if depth<7:
            level={mm(op,M) for op in OPS.values() for M in level}
    grouped=defaultdict(list)
    for M in transfers:
        grouped[(M[2],M[3])].append(M[0])
    def det3(M):
        return sum((1 if p in ((0,1,2),(1,2,0),(2,0,1)) else -1)*M[0][p[0]]*M[1][p[1]]*M[2][p[2]] for p in itertools.permutations(range(3)))
    pairs=0
    for side,rows in grouped.items():
        if len(rows)<2:
            continue
        for a,b in itertools.combinations(rows,2):
            check(a[:2]!=b[:2], 'no d-only switch in finite word universe')
            check(det3((a,*side))==0 and det3((b,*side))==0, 'no rank-three collision branch in finite word universe')
            pairs+=1
    short_scan=dict(prefix_alphabet='O,Z,J,L,C',maximum_prefix_length=7,mandatory_suffix='O',restricted_input='(D,0,U,d)',raw_prefix_word_count=sum(5**k for k in range(8)),new_transfers_by_minimum_prefix_length=first_seen,terminal_transfer_count=len(transfers),multiple_scalar_sidecar_buckets=sum(len(rows)>1 for rows in grouped.values()),unordered_switch_pairs=pairs,d_only_switch_pairs=0,rank_three_switch_branches=0)
    check(mv(A,(1,0,2,1))==(19,0,58,9), 'switch positive control')
    check(mv(B,(1,0,2,1))==(15,0,58,9), 'switch paired control')
    diff = tuple(u-v for u,v in zip(mv(A,(1,0,2,1)),mv(B,(1,0,2,1))))
    check(diff==(4,0,0,0), 'pure D control')
    check(mv(J,diff)==(4,4,4,4), 'later graft destroys common sidecars')
    return dict(status='PASS',scope='Exact algebra plus finite literal controls; no scalar coverage or universal adaptive grammar claim',
        word_convention='left to right',switch_matrices={'JOLO':A,'JJLO':B,'JOZJO':A2,'JJJJO':B2},short_scan=short_scan,sidecar_determinant=side_det,
        literal_switch_universe='All 75 labelled tournaments of orders 1..4, root 0, each prefixed by O; all four switch words; max final order 10',
        literal_switch_checks=direct_switches,gadget_universe='All 75 labelled tournaments of orders 1..4, root 0',
        graft_column_checks=column_gadgets,literal_suffixes=['','O','Z','J','L','C'],literal_column_suffix_checks=column_suffixes,
        algebraic_word_universe='Every word of length 0..5 over O,Z,J,L,C',algebraic_words=all_words,graft_containing_words=nontrivial,
        iteration_control_seed=[1,0,2,1],iteration_fibers=fibers,mixed_all_four_iteration_fibers=mixed_fibers,
        controls={'pure_D_difference':diff,'same_suffix_J_difference':mv(J,diff)})


if __name__ == '__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
