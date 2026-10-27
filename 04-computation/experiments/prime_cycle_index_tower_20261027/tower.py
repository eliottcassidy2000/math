"""The prime cycle-index tower: types(n) = non-rigid census at n + rigid rows from the ledger.

At each order n the species identity (session56, THM-4803; klein's THM-4814 (B)) determines everything
except the primes of order n themselves from the types of orders < n:
    A = SEQ(S),  S = X + sum_Q Z_{Aut Q}[A]                      (unlabelled classes)
    SC(n) = sc_nonstrong + sc_strong_nonprime + SCP(n)            (self-converse classes)
    c = MSET(c_s),  c_s = X + sum_{Q up to duality} Z_{Aut+-Q}[c] (intransitivity classes)
    1 + T_lab = 1/(1 - S_lab),  S_lab = x + sum_Q (Z_lab)          (labelled)
so with A(n), SC(n) (Burnside, burnside_counts.py; = OEIS A000568, A002785) it gives the number P(n) of
prime classes and SCP(n) of self-converse prime classes.  The non-rigid primes (nonrigid_primes.c)
give every row of the order-n type table with |Aut| > 1; the rigid rows are then
    rigid self-converse     = SCP(n) - (non-rigid self-converse),  anti-automorphism type 1^1 2^((n-1)/2) or 2^(n/2)
    rigid not self-converse = P(n) - SCP(n) - (non-rigid not self-converse).
Independent check at every order: the labelled prime count from 2^C(n,2) (exact, labelled ledger)
must equal sum over types of count * n!/|Aut|  (the Burnside defect identity).
No gentourng census is used; at n <= 10 the result is compared with klein's census (THM-4814 (A)).

Usage: python3 tower.py NMAX [NRDIR]   (NRDIR holds nr{n}.txt from nonrigid_primes, default out/)
"""
import sys, os, re
from fractions import Fraction
from math import factorial, comb
from burnside_counts import counts

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
NRDIR = sys.argv[2] if len(sys.argv) > 2 else "out"
KLEIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/data")

def parse_ms(s):
    res = []
    for tok in s.split():
        ct, mult = tok.rsplit(':', 1); d = {}
        for part in ct.strip(',').split(','):
            if part: l, e = part.split('^'); d[int(l)] = int(e)
        res.append((d, int(mult)))
    return res

def read_types(fn):
    T = []
    for line in open(fn):
        m = re.match(r'TYPE count=(\d+) k=(\d+) aut=(\d+) AUT:(.*?) SC=(\d) ANTI:(.*)$', line.strip())
        if not m: continue
        T.append(dict(k=int(m.group(2)), cnt=int(m.group(1)), G=int(m.group(3)), aut=parse_ms(m.group(4)),
                      sc=int(m.group(5)), anti=parse_ms(m.group(6)), line=line.strip()))
    return T

# ---------------- truncated power series over Fraction ----------------
def pmul(p, q, n):
    r = [Fraction(0)] * (n + 1)
    for i, a in enumerate(p):
        if a == 0 or i > n: continue
        for j in range(n + 1 - i):
            if q[j]: r[i + j] += a * q[j]
    return r

def subst_power(ser, l, n):
    """ser(x^l) truncated at degree n (ser[0] ignored -> taken 0)"""
    r = [Fraction(0)] * (n + 1)
    for j in range(1, n + 1):
        if l * j <= n: r[l * j] = Fraction(ser[j])
    return r

def cycle_index_coeff(ms, series_for_l, n):
    """[x^n] sum_{(ct, mult)} mult * prod_l series_for_l(l)(x^l)^{c_l}"""
    tot = Fraction(0)
    for ct, mult in ms:
        poly = [Fraction(0)] * (n + 1); poly[0] = Fraction(1)
        for l, c in ct.items():
            base = subst_power(series_for_l(l), l, n)
            for _ in range(c): poly = pmul(poly, base, n)
        tot += mult * poly[n]
    return tot

types = {}   # order -> list of type dicts
a = [Fraction(0)] * (NMAX + 1); s = [Fraction(0)] * (NMAX + 1)
sc = [Fraction(0)] * (NMAX + 1); scs = [Fraction(0)] * (NMAX + 1)
c = [Fraction(0)] * (NMAX + 1); cs = [Fraction(0)] * (NMAX + 1)
a[1] = s[1] = sc[1] = scs[1] = c[1] = cs[1] = Fraction(1)
# labelled
Tlab = [Fraction(2 ** comb(m, 2), factorial(m)) for m in range(NMAX + 1)]   # EGF of all labelled tournaments
inv = [Fraction(0)] * (NMAX + 1); inv[0] = Fraction(1)
for m in range(1, NMAX + 1): inv[m] = -sum(Tlab[j] * inv[m - j] for j in range(1, m + 1))
slab = [0] + [-inv[m] * factorial(m) for m in range(1, NMAX + 1)]             # labelled strong (A054946)
Plab = {}
Tplus = [Fraction(0)] + Tlab[1:]
L = [Fraction(0)] * (NMAX + 1); Ls = [Fraction(0)] * (NMAX + 1)
L[1] = Ls[1] = Fraction(1)

def seq_count(ser, n):
    """[x^n] 1/(1 - ser) for ser without constant term"""
    f = [Fraction(0)] * (n + 1); f[0] = Fraction(1)
    for m in range(1, n + 1): f[m] = sum(ser[j] * f[m - j] for j in range(1, m + 1))
    return f[n]

def seqp(ser, p, t):
    f = [Fraction(0)] * (t + 1); f[0] = Fraction(1)
    for _ in range(p):
        g = [Fraction(0)] * (t + 1)
        for i in range(t + 1):
            if f[i] == 0: continue
            for j in range(1, t + 1 - i): g[i + j] += f[i] * ser[j]
        f = g
    return f[t]

def mset_count(ser, n):
    """[x^n] prod_j (1 - x^j)^(-ser[j])  (Euler transform: multisets of ser[j] kinds of size j)"""
    f = [Fraction(0)] * (n + 1); f[0] = Fraction(1)
    for j in range(1, n + 1):
        k = int(ser[j])
        if k == 0: continue
        g = [Fraction(0)] * (n + 1)
        for m in range(n + 1):
            if f[m] == 0: continue
            t = 0
            while m + j * t <= n:
                g[m + j * t] += f[m] * comb(k + t - 1, t); t += 1
        f = g
    return f[n]

report = []
for n in range(2, NMAX + 1):
    An, SCn = counts(n)
    lower = [t for k in range(3, n) for t in types.get(k, [])]
    # unlabelled: strong non-prime and non-strong
    s_np = sum(t['cnt'] * cycle_index_coeff(t['aut'], lambda l: a, n) / t['G'] for t in lower)
    s_tmp = s[:]; s_tmp[n] = s_np
    a_np = seq_count(s_tmp, n)            # all classes except the order-n primes
    nonstrong = a_np - s_np
    P = An - a_np
    # self-converse
    ns = Fraction(0)
    for p in range(1, n // 2 + 1):
        if n % 2 == 0: ns += seqp(s, p, n // 2)
        for j in range(1, n - 2 * p + 1):
            if (n - j) % 2 == 0 and (n - j) // 2 >= p: ns += seqp(s, p, (n - j) // 2) * scs[j]
    snp = Fraction(0)
    for t in lower:
        if not t['sc']: continue
        orb = cycle_index_coeff(t['aut'], lambda l: a, n) / t['G']
        fixed_anti = cycle_index_coeff(t['anti'], lambda l: (sc if l % 2 == 1 else a), n)
        orbpm = (cycle_index_coeff(t['aut'], lambda l: a, n) + fixed_anti) / (2 * t['G'])
        snp += t['cnt'] * (2 * orbpm - orb)
    SCP = SCn - ns - snp
    # labelled primes
    tot = Fraction(0)
    for k in range(3, n):
        if k not in Plab: continue
        r = [Fraction(0)] * (n + 1); r[0] = Fraction(1)
        for _ in range(k): r = pmul(r, Tplus, n)
        tot += Plab[k] / factorial(k) * r[n]
    Plab_n = slab[n] - tot * factorial(n)
    # ---- types at order n: non-rigid census + rigid rows from P, SCP ----
    tn = []; defect_ok = None
    if n >= 3:
        fn = os.path.join(NRDIR, f"nr{n}.txt")
        if os.path.exists(fn):
            nr = read_types(fn)
            nr_sc = sum(t['cnt'] for t in nr if t['sc']); nr_nsc = sum(t['cnt'] for t in nr if not t['sc'])
            rig_sc = SCP - nr_sc; rig_nsc = P - SCP - nr_nsc
            assert rig_sc >= 0 and rig_nsc >= 0 and rig_nsc % 2 == 0, (n, rig_sc, rig_nsc)
            anti_t = f"1^1,2^{(n-1)//2}," if n % 2 else f"2^{n//2},"
            if rig_nsc: tn.append(dict(k=n, cnt=int(rig_nsc), G=1, aut=[({1: n}, 1)], sc=0, anti=[], line=f"TYPE count={int(rig_nsc)} k={n} aut=1 AUT: 1^{n},:1 SC=0 ANTI:"))
            if rig_sc: tn.append(dict(k=n, cnt=int(rig_sc), G=1, aut=[({1: n}, 1)], sc=1, anti=parse_ms(anti_t + ":1"), line=f"TYPE count={int(rig_sc)} k={n} aut=1 AUT: 1^{n},:1 SC=1 ANTI: {anti_t}:1"))
            tn += nr
            types[n] = tn
            plab_types = sum(Fraction(t['cnt'] * factorial(n), t['G']) for t in tn)
            defect_ok = (plab_types == Plab_n)
        else:
            defect_ok = None
    # update series with the order-n primes (if types known) or with the totals
    s[n] = s_np + P; a[n] = An
    sc[n] = SCn; scs[n] = SCn - ns
    Plab[n] = Plab_n
    # intransitivity classes
    tot_c = Fraction(0)
    for t in lower:
        if t['sc']: tot_c += t['cnt'] * cycle_index_coeff(t['aut'] + t['anti'], lambda l: c, n) / (2 * t['G'])
        else: tot_c += Fraction(t['cnt'], 2) * cycle_index_coeff(t['aut'], lambda l: c, n) / t['G']
    cs_np = tot_c
    cs[n] = cs_np + (P + SCP) / 2            # each prime class of order n up to duality, with singleton modules
    c[n] = mset_count(cs, n)
    # labelled intransitivity structures
    Ls_n = Fraction(0)
    for k in range(3, n + 1):
        if k in Plab:
            r = [Fraction(0)] * (n + 1); r[0] = Fraction(1)
            for _ in range(k): r = pmul(r, L, n)
            Ls_n += Fraction(1, 2) * Plab[k] / factorial(k) * r[n]
    Ls[n] = Ls_n
    # 1 + L = exp(Ls): L[n] via the standard recurrence for exp
    E = [Fraction(0)] * (n + 1); E[0] = Fraction(1)
    for m in range(1, n + 1): E[m] = sum(Fraction(j, m) * Ls[j] * E[m - j] for j in range(1, m + 1))
    L[n] = E[n]
    report.append((n, An, nonstrong, s_np, P, SCn, ns, snp, SCP, cs[n], c[n], slab[n], Plab_n, L[n] * factorial(n), defect_ok))

print("n | A(n) | non-strong | strong non-prime | P(n) | SC(n) | SC non-strong | SC strong non-prime | SCP(n) | c_s(n) | c(n) | s_lab(n) | P_lab(n) | l(n) | defect check")
for row in report:
    print(" | ".join(str(int(x)) if isinstance(x, Fraction) else str(x) for x in row))
for n in sorted(types):
    with open(os.path.join(NRDIR, f"t{n}.txt"), "w") as f:
        for t in types[n]: f.write(t['line'] + "\n")
# compare with klein's census for n <= 10
for n in sorted(types):
    fn = os.path.join(KLEIN, f"prime_types_n{n}.txt")
    if not os.path.exists(fn): continue
    mine = sorted(t['line'] for t in types[n]); his = sorted(l.strip() for l in open(fn) if l.startswith("TYPE"))
    print(f"n={n}: tower types {'==' if mine == his else '!='} klein census ({len(mine)} rows)")
