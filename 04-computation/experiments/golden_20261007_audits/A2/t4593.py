"""A2 audit of THM-4593: exact pair chains (Fraction), witnesses, -1 shadow, brute-force all-lags identity."""
from fractions import Fraction as Fr
import random, re, math, sys

def Tmap(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2

def step(state, beta):
    """THM-4581 pair-chain table; e in Z[1/3]; parity of e as 2-adic number."""
    k, e = state
    num, den = e.numerator, e.denominator   # den is a power of 3 (odd)
    sigma = num % 2                          # e mod 2 (den odd)
    if sigma == 0 and beta == 0:
        return (k, e / 2)
    if sigma == 0 and beta == 1:
        return (k, (3 * e + 1 - Fr(3) ** k) / 2)
    if sigma == 1 and beta == 0:
        return (k + 1, (3 * e + 1) / 2)
    return (k - 1, (e - Fr(3) ** (k - 1)) / 2)

def bits_of(y, t):
    out = []
    for _ in range(t):
        out.append(y % 2); y = Tmap(y)
    return out

# ---------- sanity: chain vs direct orbits on random integers ----------
random.seed(1)
bad = 0
for trial in range(300):
    y = random.getrandbits(200); r = random.randint(1, 40)
    st = (0, Fr(r)); u, v = y + r, y
    for t in range(150):
        b = v % 2
        st = step(st, b); u, v = Tmap(u), Tmap(v)
        k, e = st
        if Fr(u) != Fr(3) ** k * v + e: bad += 1; break
print("chain vs direct orbit mismatches:", bad)

# ---------- R = 2 witness y = 21 mod 32 ----------
bits = bits_of(21, 5)
sts = {r: (0, Fr(r)) for r in (1, 2)}
for t in range(5):
    for r in sts: sts[r] = step(sts[r], bits[t])
print("R=2 witness states at t=5:", sts)
for trial in range(5):
    q = random.getrandbits(100); y = 32 * q + 21
    a = [y + i for i in range(3)]
    for t in range(5): a = [Tmap(x) for x in a]
    assert a[1] == a[2] == 27 * q + 20 and a[0] == 3 * q + 2
print("direct: y,y+1,y+2 at t=5 = 3q+2, 27q+20, 27q+20 for random lifts: OK")

# ---------- all witnesses R<=32 ----------
lines = open(sys.argv[1]).read().strip().splitlines()
for ln in lines:
    m = re.match(r"R=(\d+): .* class y = (\d+) mod 2\^(\d+) \(partners coalesced at Terras time (\d+)", ln)
    R, rho, t0, tc = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
    bits = bits_of(rho, t0)
    sts = [(0, Fr(r)) for r in range(1, R + 1)]
    absorbed = False; first_all = None
    for t in range(t0):
        sts = [step(s, bits[t]) for s in sts]
        if any(s == (0, Fr(0)) for s in sts): absorbed = True
        if first_all is None and len(set(sts)) == 1: first_all = t + 1
    common = len(set(sts)) == 1
    # direct check on 3 random lifts
    okd = True
    for trial in range(3):
        y = rho + (random.getrandbits(64) << t0)
        vals = [y + r for r in range(R + 1)]
        for t in range(t0): vals = [Tmap(x) for x in vals]
        k, e = sts[0]
        if not (len(set(vals[1:])) == 1 and Fr(vals[1]) == Fr(3) ** k * vals[0] + e and vals[0] != vals[1]): okd = False
    print(f"R={R:2d} t0={t0:3d} none_absorbed={not absorbed} all_equal_at_t0={common} first_all_equal={first_all} state=(k={sts[0][0]}, e={str(sts[0][1])[:30]}) direct_ok={okd}")

# ---------- -1 shadow: y = -1 mod 2^L, partners r: k = O_L(r-1) - L ----------
L = 30
okc = 0
for r in range(1, 41):
    y = (1 << 200) * random.getrandbits(50) - 1   # y = -1 mod 2^200 (positive representative? use 2-adic -1 lift)
    y = (y % (1 << 200))  # positive integer = -1 mod 2^200
    bits = bits_of(y, L)
    assert all(b == 1 for b in bits)
    st = (0, Fr(r)); absorbed = False
    for t in range(L):
        st = step(st, bits[t]); absorbed |= (st == (0, Fr(0)))
    O = sum(bits_of(r - 1, L)) if r > 1 else 0
    if st[0] == O - L and not absorbed: okc += 1
print(f"-1 shadow: k_L = O_L(r-1) - L holds for {okc}/40 partners (L={L})")
