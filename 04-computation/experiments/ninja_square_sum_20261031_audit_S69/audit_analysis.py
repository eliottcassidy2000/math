"""S69 audit (independent referee, 2026-10-31): recomputation of every number in THM-4850 / THM-4851 / HYP-9447.

Uses only the audit's own programs (ss_enum.exe, ss_terminal_big.exe, built from this folder's .c files) and the
repository's frozen data for the prior-art cross-checks (threshold42 SMALL_BASES.json, the 2026-10-30 square bridge).

usage: python audit_analysis.py PATH_TO_ss_enum.exe PATH_TO_ss_terminal_big.exe  > audit_analysis.out
"""
import collections
import gzip
import json
import subprocess
import sys
import time
from math import isqrt
from pathlib import Path

T0 = time.time()
ENUM, BIG = sys.argv[1], sys.argv[2]
HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
RES = REPO / '05-knowledge' / 'results'


def run(exe, *args):
    return subprocess.run([exe, *map(str, args)], capture_output=True, text=True, check=True).stdout


def kv(line):
    return dict(t.split('=', 1) for t in line.split() if '=' in t)


def say(*a):
    print(*a, flush=True)


def sq(x):
    return isqrt(x) ** 2 == x


def head(t):
    say('')
    say('=' * 100)
    say(t)
    say('=' * 100)


def check_weak_B(S, L, start=(1, 3), ends=None):
    """Independent weak-B (or other normalisation) checker. ends(order) -> required last vertex or None."""
    n = len(S)
    assert len(L) == n + 1, 'lengths'
    for P in (S, L):
        assert sorted(P) == list(range(1, len(P) + 1)), 'permutation'
        assert tuple(P[:len(start)]) == tuple(start), 'prefix'
        if ends is not None:
            assert P[-1] == ends(len(P)), 'end'
        assert all(sq(a + b) for a, b in zip(P, P[1:])), 'square edge'
    posS = {v: i for i, v in enumerate(S)}
    posL = {v: i for i, v in enumerate(L)}
    assert all(posS[v] % 2 == posL[v] % 2 for v in range(1, n + 1)), 'positional parity'
    return True


B_END = (lambda order: 8 if order % 2 == 0 else 10)

# ------------------------------------------------------------------------------------------------ A census vs OEIS
head('A. Census of Q_N (own enumerator; paths rooted at every start, cycles rooted at the largest vertex) vs OEIS')
# fetched 2026-10-31 from https://oeis.org/A071983/list and https://oeis.org/A071984/list
A071983 = dict(zip(range(15, 59), [1, 1, 1, 0, 0, 0, 0, 0, 3, 0, 10, 12, 35, 52, 19, 20, 349, 392, 669, 4041, 17175,
                                   12960, 14026, 11889, 29123, 39550, 219968, 553694, 2178103, 5301127, 12220138,
                                   38838893, 68361609, 140571720, 280217025, 204853870, 738704986, 2368147377,
                                   5511090791, 9802605881, 21164463050, 47746712739, 68092497615, 123092214818]))
A071984 = dict(zip(range(32, 66), [1, 1, 11, 57, 31, 20, 25, 50, 64, 464, 1062, 4337, 10091, 21931, 69623, 115913,
                                   227893, 457707, 297126, 1051583, 3377189, 7618873, 12476654, 25832098, 55792448,
                                   75126741, 129180538, 357114149, 823402071, 3902161448, 20043267339, 131420398568,
                                   347422743997, 811591067418]))
census = {}
for line in run(ENUM, 'census', 1, 42).splitlines():
    d = kv(line)
    census[int(d['N'])] = {k: int(v) for k, v in d.items()}
    say('  ' + line)
for N in range(15, 43):
    assert census[N]['paths'] == A071983[N], ('A071983', N)
for N in range(1, 43):
    assert census[N]['cycles'] == A071984.get(N, 0), ('A071984', N)
assert [N for N in range(2, 43) if census[N]['paths']] == [15, 16, 17, 23] + list(range(25, 43))
say('  OK: paths = A071983 for 15..42 (reversal identified), cycles = A071984 for 32..42 (rotation+reversal),')
say('      no cycle below 32; Hamiltonian paths for N in {15,16,17,23} u [25,42] (N = 1 trivial)')

# ------------------------------------------------------------------------------------------------ B weak-B census
head('B. Weak B ("ninja") pairs, exhaustive for 2 <= n <= 53 (own enumerator, kind B)')
pairsB = {}
for fn in ('pairs_B_2_47.out', 'pairs_B_48_53.out'):
    for line in (HERE / fn).read_text().splitlines():
        d = kv(line)
        pairsB[int(d['n'])] = d
for n in range(2, 54):
    d = pairsB[n]
    if n >= 30:
        say(f"  n={n}: S={d['S_paths']} L={d['L_paths']} L_masks={d['L_distinct_masks']} S_ok={d['S_compatible']} "
            f"L_ok={d['L_compatible']} pairs={d['pairs']} exists={d['exists']}")
assert [n for n in range(2, 54) if pairsB[n]['exists'] == 'YES'] == list(range(42, 54))
S_ok = [int(pairsB[n]['S_compatible']) for n in range(42, 54)]
assert S_ok == [10, 4, 45, 3, 5, 11, 859, 108, 153, 230, 10710, 32242], S_ok
assert [int(pairsB[n]['pairs']) for n in range(42, 50)] == [20, 13, 221, 3, 11, 17, 3384, 276]
sess_S = {40: 11, 41: 1, 42: 139, 43: 136, 44: 946, 45: 1129, 46: 3864, 47: 6083, 48: 8212, 49: 19615, 50: 8591,
          51: 24777, 52: 153594, 53: 302246}
sess_L = {40: 1, 41: 139, 42: 136, 43: 946, 44: 1129, 45: 3864, 46: 6083, 47: 8212, 48: 19615, 49: 8591, 50: 24777,
          51: 153594, 52: 302246, 53: 481425}
for n in range(40, 54):
    assert int(pairsB[n]['S_paths']) == sess_S[n] and int(pairsB[n]['L_paths']) == sess_L[n], n
say('  OK: no weak B pair for n <= 41; pairs for every 42 <= n <= 53; all Shorter/Longer/compatible/pair counts of the')
say('      session table reproduced; pair counts for 50..53 (blank in the note): '
    + ', '.join(f"{n}: {pairsB[n]['pairs']}" for n in range(50, 54)))

head('B2. Other normalisations of a "nice pair" (same parity condition), for comparison')
for fns, kind in ((('pairs_G_14_45.out', 'pairs_G_46_50.out'), 'G'), (('pairs_F_14_36.out', 'pairs_F_37_41.out'), 'F')):
    rows = {}
    for fn in fns:
        rows.update({int(kv(l)['n']): kv(l) for l in (HERE / fn).read_text().splitlines()})
    ex = [n for n in sorted(rows) if rows[n]['exists'] == 'YES']
    say(f'  kind {kind}: pairs exist for n in {ex}  (range {min(rows)}..{max(rows)})')
    if kind == 'G':
        assert ex == list(range(41, 51))
        say('    G = Gerbicz 2018 normalisation (both start 1; odd order ends 3, even order ends 8): first pair at n = 41,')
        say('    exactly the lower end of his base block 41..2032 (session124_gerbicz_reconstruction_20261031)')
        G_END = (lambda order: 8 if order % 2 == 0 else 3)
        out = run(ENUM, 'pairdump', 41, 1, 'G').splitlines()
        Sg, Lg = (list(map(int, out[i].split(': ')[1].split(','))) for i in (0, 1))
        assert check_weak_B(Sg, Lg, start=(1,), ends=G_END)
        say('    witness at n = 41 (checked):  S = ' + ','.join(map(str, Sg)))
        say('                                  L = ' + ','.join(map(str, Lg)))
    else:
        assert ex[:5] == [15, 16, 25, 26, 27] and 28 not in ex and ex[5:] == list(range(29, 42))
        say('    F = no prefix / end condition (the "nice pair" of square_sum_B_base_interval_20261019, line 5):')
        say('    pairs at n = 15, 16, 25, 26, 27 and every 29..41 (28 has none); every weak B pair is one, so with')
        say('    THM-4850 free nice pairs exist exactly for n in {15,16,25,26,27} u [29, oo)')
        out = run(ENUM, 'pairdump', 15, 1, 'F').splitlines()
        Sf, Lf = (list(map(int, out[i].split(': ')[1].split(','))) for i in (0, 1))
        assert check_weak_B(Sf, Lf, start=())
        say('    witness at n = 15 (checked):  S = ' + ','.join(map(str, Sf)) + '   L = ' + ','.join(map(str, Lf)))

# ------------------------------------------------------------------------------------------------ C gate
head('C. The gate: odd-order Hamiltonian paths 1,3 -> 10')
odd10 = {N: census[N]['p13_10'] for N in range(15, 43, 2)}
say('  ' + str(odd10))
assert [N for N, c in odd10.items() if c] == [41] and odd10[41] == 1
gate = run(ENUM, 'list', 41, 1, 3, 10).split()
gate_path = list(map(int, gate[0].split(',')))
assert gate[-1] == 'total=1'
assert gate_path == [1, 3, 6, 30, 19, 17, 32, 4, 21, 15, 34, 2, 7, 18, 31, 33, 16, 9, 40, 41, 8, 28, 36, 13, 23, 26, 38,
                     11, 5, 20, 29, 35, 14, 22, 27, 37, 12, 24, 25, 39, 10]
say('  OK: the unique Q_41 path equals the one printed in THM-4850 (2); none at odd N in 15..39')
say(f"  even-order 1,3 -> 8 first at N = {min(N for N in range(2, 43) if census[N]['p13_8'])}")

# ------------------------------------------------------------------------------------------------ D cycles
head('D. Hamiltonian cycles of Q_32..Q_42, the unique 32- and 33-cycles')


def canon(c):
    i = c.index(1)
    c = c[i:] + c[:i]
    if c[1] > c[-1]:
        c = [c[0]] + c[1:][::-1]
    return tuple(c)


def edges(c):
    return {frozenset((c[i], c[(i + 1) % len(c)])) for i in range(len(c))}


cyc = {}
for N in range(32, 43):
    out = run(ENUM, 'cycles', N).split()
    lst = [canon(list(map(int, t.split(',')))) for t in out[:-1]]
    assert out[-1] == f'total={len(lst)}' and len(set(lst)) == len(lst)
    assert all(sorted(c) == list(range(1, N + 1)) and all(sq(c[i] + c[(i + 1) % N]) for i in range(N)) for c in lst)
    cyc[N] = lst
    assert len(lst) == A071984[N]
c32, c33 = cyc[32][0], cyc[33][0]
say('  Q_32: ' + ','.join(map(str, c32)))
say('  Q_33: ' + ','.join(map(str, c33)))
assert list(c32) == [1, 8, 28, 21, 4, 32, 17, 19, 30, 6, 3, 13, 12, 24, 25, 11, 5, 31, 18, 7, 29, 20, 16, 9, 27, 22, 14,
                     2, 23, 26, 10, 15]
dele = sorted(tuple(sorted(e)) for e in edges(c32) - edges(c33))
add = sorted(tuple(sorted(e)) for e in edges(c33) - edges(c32))
say(f'  deleted {dele}, added {add}')
assert dele == [(5, 31), (16, 20)] and add == [(5, 20), (16, 33), (31, 33)]
nb33 = [x for x in range(1, 33) if sq(33 + x)]
assert nb33 == [3, 16, 31] and not any(frozenset((a, b)) in edges(c32) for a in nb33 for b in nb33 if a < b)
for c in (c32, c33):
    assert {c[1], c[-1]} == {8, 15}
say('  OK: reversal of 31,18,7,29,20 + insertion of 33; 33 has neighbours 3,16,31, pairwise non-adjacent on the 32-cycle;')
say('      both cycles pass 8-1-15')
# prior art: the 2026-10-30 square bridge (32 -> 6242) records this very step as one insertion + one reversal
br = json.loads(gzip.decompress((RES / 'square_sum_spectrum_overnight_20261030' / 'square_bridge.json.gz').read_bytes()))
p = list(br['initial'])
assert canon(p) == c32
st = br['steps'][0]
p.insert(st['insert'] + 1, 33)
for a, b in st['reverse']:
    p[a + 1:b + 1] = p[a + 1:b + 1][::-1]
assert canon(p) == c33
say(f"  PRIOR ART: square_sum_spectrum_overnight_20261030/square_bridge.json.gz step 32->33 = {st} reproduces the")
say('      unique 33-cycle (one insertion + one reversal).')
direct = [N for N, s in zip(range(33, 43), br['steps'][:10]) if not s['reverse']]
say(f'      bridge steps 33..42 that are direct insertions (no reversal): {direct}')

# ------------------------------------------------------------------------------------------------ E terminal spectrum
head('E. Terminal-form cycles (8-1-3 for even N, 15-1-3 for odd N)')
term = {}
for line in run(ENUM, 'terminal', 8, 46, -1).splitlines():
    d = kv(line)
    N = int(d['N'])
    key = [k for k in d if k.startswith('terminal(')][0]
    term[N] = int(d[key])
for line in run(ENUM, 'terminal', 47, 62, 1).splitlines():
    d = kv(line)
    N = int(d['N'])
    key = [k for k in d if k.startswith('terminal(')][0]
    term[N] = int(d[key].rstrip('+'))
say('  exact counts N = 34..46: ' + str({N: term[N] for N in range(34, 47)}))
exists = [N for N in range(8, 63) if term[N]]
assert exists == [34, 35, 36, 37, 38] + list(range(40, 63))
for N in range(8, 47):   # cross-check with the census counts
    if N <= 42:
        assert term[N] == (census[N]['p13_8'] if N % 2 == 0 else census[N]['p13_15'])
# second, independent program (randomised Warnsdorff + forcing, complete when it terminates)
big_small = run(BIG, 30, 62).splitlines()
for line in big_small:
    d = kv(line.split(' path=')[0])
    N = int(d['N'])
    if 'EXHAUSTED' in line:
        assert N in (30, 31, 32, 33, 39), N
    else:
        P = list(map(int, line.split(' path=')[1].split(',')))
        T = 8 if N % 2 == 0 else 15
        assert sorted(P) == list(range(1, N + 1)) and P[:2] == [1, 3] and P[-1] == T
        assert all(sq(a + b) for a, b in zip(P, P[1:])) and sq(P[-1] + 1)
say('  OK: exists exactly for N in {34..38} u [40,62] (two own programs; complete searches at 30..33 and 39)')
rec = json.loads((HERE / 'terminal_big_63_2620_receipt.json').read_text())
assert rec['status'] == 'passed' and rec['range'] == [63, 2620] and rec['verified_cycles'] == 2558
say(f"  63..2620: all {rec['verified_cycles']} N re-found by ss_terminal_big and re-verified in Python "
    f"(receipt terminal_big_63_2620_receipt.json; max nodes {rec['max_nodes_last_attempt']})")
say(f'  Q_39: {len(cyc[39])} Hamiltonian cycles; through 8-1-3 or 15-1-3: '
    f"{sum(1 for c in cyc[39] if {c[1], c[-1]} in ({3, 8}, {3, 15}))}")
assert len(cyc[39]) == 50 and not any({c[1], c[-1]} in ({3, 8}, {3, 15}) for c in cyc[39])

# ------------------------------------------------------------------------------------------------ F insertion cost
head('F. Insertion cost: min over cycles C of Q_N of |E(C) - E(C\' - (N+1))|, per cycle C\' of Q_(N+1)')
mins = []
for N in range(32, 42):
    hist = collections.Counter()
    E1 = [edges(c) for c in cyc[N]]
    for c2 in cyc[N + 1]:
        E2 = {e for e in edges(c2) if N + 1 not in e}
        hist[min(len(e1 - E2) for e1 in E1)] += 1
    mins.append(min(hist))
    say(f'  {N} -> {N + 1}: {len(cyc[N])} -> {len(cyc[N + 1])} cycles; min {min(hist)}; histogram {dict(sorted(hist.items()))}')
assert mins == [2, 3, 2, 3, 3, 4, 4, 3, 2, 2], mins
say('  OK: minima 2,3,2,3,3,4,4,3,2,2 (histograms identical to the session .out)')

# ------------------------------------------------------------------------------------------------ G statistics
head('G. Third vertex and the 24-25 edge: ninja-pair members vs ALL paths with the same prefix/end (baseline)')


def paths_list(order, end):
    out = run(ENUM, 'list', order, 1, 3, end).split()
    assert out[-1].startswith('total=')
    return [tuple(map(int, t.split(','))) for t in out[:-1]]


def has2425(P):
    return any({P[i], P[i + 1]} == {24, 25} for i in range(len(P) - 1))


for n in range(42, 50):
    out = run(ENUM, 'pairdump', n, 10 ** 9)
    Ss = {tuple(map(int, l.split(': ')[1].split(','))) for l in out.splitlines() if l.strip().startswith('S: ')}
    Ls = {tuple(map(int, l.split(': ')[1].split(','))) for l in out.splitlines() if l.strip().startswith('L: ')}
    lines = out.splitlines()
    for i in range(0, len(lines) - 1, 2):      # every printed pair re-checked by the independent checker
        if lines[i].strip().startswith('S: '):
            check_weak_B(list(map(int, lines[i].split(': ')[1].split(','))),
                         list(map(int, lines[i + 1].split(': ')[1].split(','))), ends=B_END)
    mem = Ss | Ls
    allS, allL = paths_list(n, B_END(n)), paths_list(n + 1, B_END(n + 1))
    base = set(allS) | set(allL)
    t_mem = collections.Counter(P[2] for P in mem)
    t_base = collections.Counter(P[2] for P in base)
    e_mem, e_base = sum(map(has2425, mem)), sum(map(has2425, base))
    frac = lambda c, tot: f'{c}/{tot} = {100 * c / tot:.0f}%'
    t2233_mem = t_mem[22] + t_mem[33]
    t2233_base = t_base[22] + t_base[33]
    say(f'  n={n}: members {len(mem)} (S {len(Ss)}, L {len(Ls)}); third vertex {dict(sorted(t_mem.items()))}; '
        f'22/33 {frac(t2233_mem, len(mem))}; 24-25 {frac(e_mem, len(mem))}')
    say(f'        baseline all {len(base)} paths: third vertex {dict(sorted(t_base.items()))}; 22/33 '
        f'{frac(t2233_base, len(base))}; 24-25 {frac(e_base, len(base))}')

# ------------------------------------------------------------------------------------------------ H numeration
head('H. Numeration n = 49^k b + sum d_i 49^i (b in [B0, 49 B0 + 23], d_i in [24, 72])')


def numerals_upto(B0, X):
    hits = collections.Counter()
    top = 49 * B0 + 23
    layer = list(range(B0, top + 1))
    while layer:
        nxt = []
        for v in layer:
            if v <= X:
                hits[v] += 1
                for d in range(24, 73):
                    w = 49 * v + d
                    if w <= X:
                        nxt.append(w)
        layer = nxt
    return hits


for B0 in (42, 41, 53):
    X = 400000
    h = numerals_upto(B0, X)
    ok = all(h[v] == 1 for v in range(B0, X + 1)) and all(v >= B0 for v in h)
    say(f'  base block [{B0}, {49 * B0 + 23}]: every n in [{B0}, {X}] has exactly one numeral: {ok}')
    assert ok
S = [(85 * 49 ** l - 1) // 2 for l in range(4)]
T = [(107 * 49 ** l - 1) // 2 for l in range(4)]
G = [(83 * 49 ** l - 1) // 2 for l in range(4)]
assert S == [42, 2082, 102042, 5000082] and T[1] == 2621 and G[1] == 2033


def numeral(n, B0):
    ds = []
    while n > 49 * B0 + 23:
        m = (n - 24) // 49
        ds.append(n - 49 * m)
        n = m
    return n, ds[::-1]


assert all(numeral(S[l], 42) == (42, [24] * l) for l in range(4))
assert all(numeral(T[l], 42) == (53, [24] * l) for l in range(4))
assert all(numeral(G[l], 41) == (41, [24] * l) for l in range(4))
for n in (2082, 2130, 2131, 102102, 128453, 6294221):
    say(f'  {n} = {numeral(n, 42)}')
say('  threshold42 packet: S_l = (85*49^l - 1)/2 = [42; 24,...,24] = 42, 2082, 102042, 5000082 (first l-digit numerals)')
say('  session124 core:    T_l = (107*49^l - 1)/2 = [53; 24,...,24];  Gerbicz base 41: (83*49^l - 1)/2 = [41; 24,...,24]')

# ------------------------------------------------------------------------------------------------ I prior-art data
head('I. Prior art in the repository: threshold42 SMALL_BASES.json (2026-10-31 00:44) checked with the own checker')
small = json.loads((RES / 'session124_square_prefix_threshold42_20261031' / 'SMALL_BASES.json').read_text())
for n in range(42, 53):
    Sx, Lx = small[str(n)]
    assert len(Sx) == n and check_weak_B(Sx, Lx, ends=B_END)
say('  OK: the eleven threshold42 base pairs at n = 42..52 are weak B pairs (independent checker)')

# ------------------------------------------------------------------------------------------------ J Zeckendorf table
head('J. Zeckendorf / Fibonacci arithmetic quoted in the note (numerology table)')
F = [1, 2]
while F[-1] < 10 ** 4:
    F.append(F[-1] + F[-2])


def zeck(n):
    out = []
    for f in reversed(F):
        if f <= n:
            out.append(f)
            n -= f
    return out


claims = {32: [21, 8, 3], 33: [21, 8, 3, 1], 34: [34], 35: [34, 1], 39: [34, 5], 41: [34, 5, 2], 42: [34, 8],
          15: [13, 2], 16: [13, 3], 17: [13, 3, 1], 23: [21, 2], 25: [21, 3, 1]}
for n, z in claims.items():
    assert zeck(n) == z, (n, zeck(n))
Fib = [0, 1]
while len(Fib) < 15:
    Fib.append(Fib[-1] + Fib[-2])
assert (Fib[2], Fib[4], Fib[6], Fib[9]) == (1, 3, 8, 34)
assert sq(1 * 3 + 1) and sq(1 * 8 + 1) and sq(3 * 8 + 1)
say('  OK: all quoted Zeckendorf forms; 1,3,8 = F_2,F_4,F_6; 34 = F_9; {1,3,8} is a Diophantine triple')
say(f'\nALL AUDIT ASSERTIONS PASSED ({time.time() - T0:.1f} s)')
