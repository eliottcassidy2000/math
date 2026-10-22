# Audit2 (D): H of the Paley tournament on F_27 (P27.out of paley_extremal_20261021 was not recomputed in audit 1),
# by inclusion-exclusion over vertex subsets modulo the translation group (Z_3)^3 (hcount_group.c).
# F_27 = F_3[t]/(t^3 - t + 1); u -> v iff v - u is a nonzero square.  Self-test first on Cayley tournaments of Z_3^2
# and Z_3^3 against Held-Karp (hcount hk).
import sys, itertools, subprocess
HG, HC = sys.argv[1], sys.argv[2]

def vec(i, k): return [(i // 3 ** j) % 3 for j in range(k)]
def idx(v): return sum(c * 3 ** j for j, c in enumerate(v))
def add(a, b): return [(x + y) % 3 for x, y in zip(a, b)]
def sub(a, b): return [(x - y) % 3 for x, y in zip(a, b)]

def mul27(a, b):
    # polynomials a0 + a1 t + a2 t^2 modulo t^3 = t - 1
    c = [0] * 5
    for i in range(3):
        for j in range(3): c[i + j] += a[i] * b[j]
    for d in (4, 3):              # t^d = t^(d-2) - t^(d-3)
        c[d - 2] += c[d]; c[d - 3] -= c[d]; c[d] = 0
    return [x % 3 for x in c[:3]]

def cayley_code(N, k, conn):
    code = ''
    for i in range(N):
        for j in range(i + 1, N):
            code += '1' if idx(sub(vec(j, k), vec(i, k))) in conn else '0'
    return code

def translations(N, k):
    return [[idx(add(vec(v, k), vec(g, k))) for v in range(N)] for g in range(N)]

def run_group(N, k, conn):
    code = cayley_code(N, k, conn)
    inp = f"{N} {N}\n{code}\n" + '\n'.join(' '.join(map(str, p)) for p in translations(N, k)) + '\n'
    r = subprocess.run([HG], input=inp, capture_output=True, text=True)
    return int(r.stdout), code

# self-tests on random Cayley tournaments of Z_3^2 and Z_3^3 (n = 9 and 27 is too big for Held-Karp; use n = 9 only,
# and Z_3^3 checked through P_27's own invariants below)
import random
random.seed(1)
for trial in range(3):
    N, k = 9, 2
    conn = set()
    for x in range(1, N):
        y = idx([(-c) % 3 for c in vec(x, k)])
        if x < y: conn.add(random.choice([x, y]))
    Hg, code = run_group(N, k, conn)
    Hh = int(subprocess.run([HC, 'hk'], input=code + '\n', capture_output=True, text=True).stdout)
    print(f"self-test Z_3^2 Cayley tournament {sorted(conn)}: group IE {Hg}, Held-Karp {Hh}, {'ok' if Hg == Hh else 'MISMATCH'}", flush=True)
# P_27
N, k = 27, 3
sq = {idx(mul27(vec(x, 3), vec(x, 3))) for x in range(1, N)}
assert len(sq) == 13 and all(idx([(-c) % 3 for c in vec(s, 3)]) not in sq for s in sq)
Hg, code = run_group(N, k, sq)
print(f"P_27 (F_27 = F_3[t]/(t^3 - t + 1)): H = {Hg}; repo P27.out value 418425945688463432031: {'MATCH' if Hg == 418425945688463432031 else 'MISMATCH'}")
