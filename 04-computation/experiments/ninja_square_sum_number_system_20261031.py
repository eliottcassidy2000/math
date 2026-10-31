"""opus-2026-10-31-S69: the square-sum numbers start at 42 (THM-4850, THM-4851, HYP-9447) - driver with assertions.

Runs the exhaustive enumerator `ninja_pairs_square_sum_20261031.c` (compiled to ENUM) and reads the large-N terminal
search outputs `terminal_cycles_large_20261031.out` (63..1000, pre-registered) and `_window.out` (1001..2620).

  A  census N = 14..42: Hamiltonian paths and cycles (= OEIS A071983, A071984), terminal paths
  B  ninja (weak B) pairs n = 15..53, exhaustive: none for n <= 41, all n in 42..53
  C  the gate: the unique Q_41 Hamiltonian path from 1,3 to 10
  D  the unique cycles of Q_32 and Q_33 and the reversal-insertion relating them
  E  terminal-form spectrum {34..38} u [40, oo): exhaustive to 62, search 63..2620, session124 above
  F  insertion costs 32 -> 42
  G  first ninja pairs: the vertex after 1,3 (22 / 33) and the 24-25 edge
  H  the base-49 ninja numeration with digits 24..72 and base block [42, 2081]
  I  Zeckendorf forms of the threshold numbers (numerology table)

usage: python3 ninja_square_sum_number_system_20261031.py PATH_TO_ENUM_EXE
"""
import collections
import subprocess
import sys
import time
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent
ENUM = sys.argv[1]


def run(*args):
    return subprocess.run([ENUM, *map(str, args)], capture_output=True, text=True, check=True).stdout


def say(*a):
    print(*a, flush=True)


def kv(line):
    out = {}
    for tok in line.split():
        if '=' in tok:
            k, v = tok.split('=', 1)
            out[k] = v
    return out


# ---------------------------------------------------------------- A
say('=' * 78)
say('A. Census of Q_N: Hamiltonian paths and cycles (reversals identified), terminal paths')
say('=' * 78)
A071983 = {14: 0, 15: 1, 16: 1, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 3, 24: 0, 25: 10, 26: 12, 27: 35, 28: 52,
           29: 19, 30: 20, 31: 349, 32: 392, 33: 669, 34: 4041}
A071984 = {32: 1, 33: 1, 34: 11, 35: 57, 36: 31, 37: 20, 38: 25, 39: 50, 40: 64, 41: 464, 42: 1062}
census = {}
for line in run('paths', 14, 42).splitlines():
    d = kv(line)
    N = int(d['N'])
    census[N] = d
    say('  ' + line)
for N, v in A071983.items():
    assert int(census[N]['ham_paths']) == v, N
for N in range(14, 43):
    assert int(census[N]['ham_cycles']) == A071984.get(N, 0), N
paths_exist = [N for N in range(14, 43) if int(census[N]['ham_paths']) > 0]
assert paths_exist[:5] == [15, 16, 17, 23, 25] and paths_exist[5:] == list(range(26, 43))
say('  paths exist for N in {15,16,17,23} and every N >= 25 (to 42); first cycle N = 32 (unique), N = 33 (unique)')

# ---------------------------------------------------------------- B
say('')
say('=' * 78)
say('B. Ninja (weak B) pairs, exhaustive for 15 <= n <= 53')
say('=' * 78)
ninja = {}
for line in run('pairs', 15, 53).splitlines():
    n = int(line.split()[0].split('=')[1])
    ninja[n] = 'NINJA=YES' in line
    m = line.split('compatible_pairs(shorter paths matching some longer mask)=')
    cnt = int(m[1].split()[0])
    say(f'  n={n}: ninja={"YES" if ninja[n] else "no"}  matching Shorter paths={cnt}  | ' + line.split(' compatible')[0].split(' ', 1)[1])
assert [n for n in range(15, 54) if ninja[n]] == list(range(42, 54))
say('  => no ninja pair for any n <= 41; ninja pairs for every 42 <= n <= 53 (and 53..2620 certified, recursion above)')

# ---------------------------------------------------------------- C
say('')
say('=' * 78)
say('C. The gate: Hamiltonian paths of odd order from 1,3 to 10')
say('=' * 78)
odd_to10 = {N: int(census[N]['paths_13_to_10']) for N in range(15, 43) if N % 2 == 1}
say('  odd N: number of Hamiltonian paths 1,3 -> 10: ' + str(odd_to10))
assert [N for N, c in odd_to10.items() if c] == [41] and odd_to10[41] == 1
gate = run('show', 41, 1, 3, 10, 5)
say('  Q_41 unique path: ' + gate.splitlines()[0].strip())
assert gate.strip().endswith('total=1')

# ---------------------------------------------------------------- D
say('')
say('=' * 78)
say('D. The unique cycles of Q_32 and Q_33')
say('=' * 78)


def cycles_of(N):
    out = []
    for l in run('cycles', N, 10 ** 7).splitlines():
        l = l.strip()
        if l and l[0].isdigit():
            out.append(list(map(int, l.split(','))))
    return out


def edges(P):
    return {frozenset((P[i], P[(i + 1) % len(P)])) for i in range(len(P))}


c32, c33 = cycles_of(32), cycles_of(33)
assert len(c32) == 1 and len(c33) == 1
say('  Q_32: ' + ','.join(map(str, c32[0])))
say('  Q_33: ' + ','.join(map(str, c33[0])))
E32, E33 = edges(c32[0]), edges(c33[0])
deleted = sorted(tuple(sorted(e)) for e in E32 - E33)
added = sorted(tuple(sorted(e)) for e in E33 - E32)
say(f'  Q_33 cycle = Q_32 cycle minus {deleted} plus {added}: the segment 31..20 reversed and 33 inserted')
assert deleted == [(5, 31), (16, 20)] and added == [(5, 20), (16, 33), (31, 33)]
for P in (c32[0], c33[0]):
    i = P.index(1)
    assert {P[i - 1], P[(i + 1) % len(P)]} == {8, 15}
say('  both cycles pass 8-1-15 (8+1 = 9, 15+1 = 16): the two terminal closings of session124 at once')

# ---------------------------------------------------------------- E
say('')
say('=' * 78)
say('E. Terminal-form cycles (8-1-3 for even N, 15-1-3 for odd N)')
say('=' * 78)
term = {}
for line in run('exists', 30, 62).splitlines():
    d = kv(line)
    N = int(d['N'])
    term[N] = 'parity_terminal' in line and '=yes' in line.split('other')[0]
say('  exhaustive 30..62: terminal form exists for ' + str([N for N in range(30, 63) if term[N]]))
assert [N for N in range(30, 63) if term[N]] == [34, 35, 36, 37, 38] + list(range(40, 63))
big = {}
for fname in ('terminal_cycles_large_20261031.out', 'terminal_cycles_large_20261031_window.out'):
    for line in (HERE / fname).read_text().splitlines():
        d = kv(line)
        big[int(d['N'])] = ('FOUND' in line and 'verified=yes' in line)
assert sorted(big) == list(range(63, 2621)) and all(big.values())
say('  search 63..1000 (pre-registered HYP-9447) and 1001..2620: every N FOUND with an independently verified cycle')
say('  => terminal cycles exist exactly for N in {34,...,38} u [40, oo) (N >= 2621: session124 terminal theorem)')

# ---------------------------------------------------------------- F
say('')
say('=' * 78)
say('F. Adding a node: minimal number of old-cycle edges deleted, 32 -> 42')
say('=' * 78)
C = {N: [edges(P) for P in cycles_of(N)] for N in range(32, 43)}
mins = {}
for N in range(32, 42):
    hist = collections.Counter()
    for E2 in C[N + 1]:
        E2r = {e for e in E2 if N + 1 not in e}
        hist[min(len(E1 - E2r) for E1 in C[N])] += 1
    mins[N] = min(hist)
    say(f'  {N} -> {N + 1}: cycles {len(C[N])} -> {len(C[N + 1])}; min deletions {mins[N]}; histogram {dict(sorted(hist.items()))}')
assert all(v >= 2 for v in mins.values()) and mins[32] == 2
say('  no step from 32 to 42 admits a direct insertion (one deleted edge); 32 -> 33 needs exactly two')

# ---------------------------------------------------------------- G
say('')
say('=' * 78)
say('G. First ninja pairs: the vertex after 1,3 and the consecutive edge 24-25 (sum 49)')
say('=' * 78)


def has_edge(P, a, b):
    return any({P[i], P[i + 1]} == {a, b} for i in range(len(P) - 1))


for n in range(42, 50):
    out = run('pairdump', n, 10 ** 8)
    Ss = {tuple(map(int, l.split(': ')[1].split(','))) for l in out.splitlines() if l.strip().startswith('S: ')}
    Ls = {tuple(map(int, l.split(': ')[1].split(','))) for l in out.splitlines() if l.strip().startswith('L: ')}
    npairs = int(out.strip().splitlines()[-1].split('=')[-1])
    third = collections.Counter(p[2] for p in Ss | Ls)
    e = sum(has_edge(p, 24, 25) for p in Ss | Ls)
    say(f'  n={n}: path pairs {npairs}; distinct members {len(Ss | Ls)}; vertex after 1,3: {dict(sorted(third.items()))}; '
        f'members with 24-25: {e}')
    if n in (42, 43):
        assert set(third) <= {13, 22, 33} and third[22] + third[33] > third[13]

# ---------------------------------------------------------------- H
say('')
say('=' * 78)
say('H. The base-49 ninja numeration: n = 49^k b + sum d_i 49^i, b in [42, 2081], d_i in [24, 72]')
say('=' * 78)


def ninja_numeral(n):
    assert n >= 42
    digits = []
    while n > 2081:
        m = (n - 24) // 49
        d = n - 49 * m
        assert 24 <= d <= 72 and m >= 42
        digits.append(d)
        n = m
    return n, digits[::-1]


def value(b, digits):
    v = b
    for d in digits:
        v = 49 * v + d
    return v


prev = None
for n in range(42, 300000):
    b, ds = ninja_numeral(n)
    assert value(b, ds) == n and 42 <= b <= 2081
    if prev is not None and n > 2082:           # successor = odometer on the last digit
        pb, pds = prev
        if pds and pds[-1] < 72:
            assert (b, ds) == (pb, pds[:-1] + [pds[-1] + 1])
    prev = (b, ds)
for n in (42, 2081, 2082, 2130, 2131, 102102, 128453, 6294221):
    b, ds = ninja_numeral(n)
    say(f'  {n} = [{b}; ' + ', '.join(map(str, ds)) + ']')
say('  digits 24..72 form a complete residue system mod 49, so every n >= 42 has exactly one numeral; adding a node is')
say('  the odometer d -> d + 1, carrying 72 -> 24 into the next place (since 49 m + 72 + 1 = 49 (m + 1) + 24)')
say('  constructive: ninja_numeration_replay_20261031.py applies session124\'s 49 digit operators to the new pairs')
say('  m = 42..52 (539 outputs, all verified by session124\'s own check_pair)')

# ---------------------------------------------------------------- I
say('')
say('=' * 78)
say('I. Zeckendorf forms of the owner\'s numbers and the thresholds (NUMEROLOGY table; no claim)')
say('=' * 78)
F = [1, 2]
while F[-1] < 10 ** 4:
    F.append(F[-1] + F[-2])


def zeck(n):
    parts = []
    for f in reversed(F):
        if f <= n:
            parts.append(f)
            n -= f
    return parts


for n in (1, 3, 8, 15, 16, 17, 22, 23, 24, 25, 32, 33, 34, 35, 39, 40, 41, 42, 49, 72):
    say(f'  {n:3d} = ' + ' + '.join(map(str, zeck(n))))
assert zeck(32) == [21, 8, 3] and zeck(33) == [21, 8, 3, 1] and zeck(42) == [34, 8]
say('')
say(f'ALL ASSERTIONS PASSED ({time.time() - T0:.1f}s)')
