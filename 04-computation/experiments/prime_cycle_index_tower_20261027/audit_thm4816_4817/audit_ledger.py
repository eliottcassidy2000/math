"""Auditor's independent ledger (exact rationals) for THM-4816 / THM-4817.

Inputs: type tables (klein census files n<=10, the auditor's own non-rigid census n=11,12[,13], optionally the
producer's tables), A(n) and SC(n) from audit_burnside.py.  Rigid rows at orders >= 11 are closed by the ledger
(THM-4816 A3/A4), re-derived here.  Everything is Fraction / int; no floats anywhere.

Identities (re-derived by the auditor):
  strong:  S(x) = x + sum_Q Z_{Aut Q}[A](x);   1 + A = 1/(1 - S)
  SC:      SC(n) = SC_ns(n) + SC_snp(n) + SCP(n),
           SC_ns  = [x^n] (S(x^2)/(1 - S(x^2))) (1 + Sc(x)),       Sc = strong self-converse classes
           SC_snp = [x^n] sum_{Q sc} (1/|Aut|) sum_{psi anti} prod_{cycles l} (SC(x^l) if l odd else A(x^l))
  intrans: c_s = x + sum_{Q non-sc} (1/2) Z_Aut[c] + sum_{Q sc} Z_{Aut+-}[c];  1 + c = exp(sum_k c_s(x^k)/k)
  labelled: T = sum 2^C(n,2) x^n/n!, S_lab = 1 - 1/T, P_lab(n) = s_lab(n) - n! sum_{3<=k<n} P_lab(k)/k! [x^n](T-1)^k
           l_s = x + (1/2) sum_k P_lab(k) l^k/k!,  1 + l = exp(l_s)
"""
from fractions import Fraction as Fr
from math import factorial, comb, gcd
from collections import Counter
import re, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../..'))
KLEIN = os.path.join(ROOT, '05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/data')
PROD = os.path.join(HERE, '..', 'out')


def parse_ct(s):
    d = Counter()
    for part in re.split(r'[,.]', s):
        if part:
            l, e = part.split('^')
            d[int(l)] += int(e)
    return tuple(sorted(d.items()))


def parse_ms(s):
    res = Counter()
    for tok in s.split():
        ct, mult = tok.rsplit(':', 1)
        res[parse_ct(ct)] += int(mult)
    return tuple(sorted(res.items()))


def parse_types(fn):
    rows = []
    for line in open(fn):
        line = line.strip()
        if not line.startswith('TYPE'):
            continue
        m = re.match(r'TYPE count=(\d+) (?:k|n)=(\d+) aut=(\d+) AUT:(.*?) SC=(\d) ANTI:(.*)$', line)
        assert m, line
        cnt, k, G = int(m.group(1)), int(m.group(2)), int(m.group(3))
        aut, sc, anti = parse_ms(m.group(4)), int(m.group(5)), parse_ms(m.group(6))
        assert sum(x for _, x in aut) == G, line
        if sc:
            assert sum(x for _, x in anti) == G, line
        else:
            assert not anti, line
        rows.append((k, cnt, G, aut, sc, anti))
    return rows


def rowkey(r):
    return (r[0], r[2], r[3], r[4], r[5])


def table_dict(rows):
    d = Counter()
    for r in rows:
        d[rowkey(r)] += r[1]
    return d


# ---------------- power series helpers ----------------
def mul(a, b, N):
    r = [Fr(0)] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                if b[j]:
                    r[i + j] += x * b[j]
    return r


def subst(F, l, N):
    """F(x^l) truncated"""
    r = [Fr(0)] * (N + 1)
    for j in range(1, N // l + 1):
        r[l * j] = Fr(F[j])
    return r


def cycle_prod(part, series_for_len, N):
    prod = [Fr(0)] * (N + 1)
    prod[0] = Fr(1)
    for l, m in part:
        F = series_for_len(l)
        for _ in range(m):
            prod = mul(prod, F, N)
    return prod


def coef_ms(ms, series_for_len, n):
    """[x^n] sum_{(part, mult)} mult * prod_cycles series(l)(x^l)"""
    tot = Fr(0)
    for part, mult in ms:
        tot += mult * cycle_prod(part, series_for_len, n)[n]
    return tot


def burnside_table(nmax):
    sys.path.insert(0, HERE)
    import audit_burnside
    A = [0] * (nmax + 1)
    SC = [0] * (nmax + 1)
    for n in range(1, nmax + 1):
        A[n], SC[n] = audit_burnside.counts(n)
    return A, SC


class Ledger:
    def __init__(self, A, SC, N):
        self.N = N
        self.A = [Fr(x) for x in A[:N + 1]]
        self.SC = [Fr(x) for x in SC[:N + 1]]
        self.types = {}  # order -> list of rows
        self.S = [Fr(0)] * (N + 1)
        self.NS = [Fr(0)] * (N + 1)
        self.SNP = [Fr(0)] * (N + 1)
        self.P = [Fr(0)] * (N + 1)
        self.Sc = [Fr(0)] * (N + 1)
        self.SCns = [Fr(0)] * (N + 1)
        self.SCsnp = [Fr(0)] * (N + 1)
        self.SCP = [Fr(0)] * (N + 1)
        self.S[1] = self.A[1] = Fr(1)
        self.Sc[1] = Fr(1)

    def rows_below(self, n):
        for k in sorted(self.types):
            if k < n:
                for r in self.types[k]:
                    yield r

    def close_order(self, n):
        """compute NS, S, SNP, P, SC parts at order n from types of orders < n"""
        N = n
        # non-strong: [x^n] (S/(1-S) - S) with S known below n
        Sser = [Fr(0)] + [self.S[j] if j < n else Fr(0) for j in range(1, N + 1)]
        # sequences of length >= 2
        seq = [Fr(0)] * (N + 1)
        power = Sser[:]
        for k in range(2, n + 1):
            power = mul(power, Sser, N)
            seq = [seq[i] + power[i] for i in range(N + 1)]
        self.NS[n] = seq[n]
        self.S[n] = self.A[n] - self.NS[n]
        Aser = lambda l: subst(self.A, l, N)
        snp = Fr(0)
        for (k, cnt, G, aut, sc, anti) in self.rows_below(n):
            snp += cnt * coef_ms(aut, Aser, n) / G
        self.SNP[n] = snp
        self.P[n] = self.S[n] - snp
        # self-converse
        S2 = subst([Fr(0)] + [self.S[j] if j <= n // 2 else Fr(0) for j in range(1, N + 1)], 2, N)
        seq2 = [Fr(0)] * (N + 1)
        power = S2[:]
        for k in range(1, n // 2 + 1):
            seq2 = [seq2[i] + power[i] for i in range(N + 1)]
            power = mul(power, S2, N)
        onepSc = [Fr(1)] + [self.Sc[j] if j < n else Fr(0) for j in range(1, N + 1)]
        self.SCns[n] = mul(seq2, onepSc, N)[n]
        self.Sc[n] = self.SC[n] - self.SCns[n]
        mixed = lambda l: subst(self.SC if l % 2 == 1 else self.A, l, N)
        scsnp = Fr(0)
        for (k, cnt, G, aut, sc, anti) in self.rows_below(n):
            if sc:
                scsnp += cnt * coef_ms(anti, mixed, n) / G
        self.SCsnp[n] = scsnp
        self.SCP[n] = self.Sc[n] - scsnp

    def rigid_rows(self, n, nonrigid_rows):
        nr_sc = sum(r[1] for r in nonrigid_rows if r[4])
        nr_nsc = sum(r[1] for r in nonrigid_rows if not r[4])
        P, SCP = self.P[n], self.SCP[n]
        assert P.denominator == 1 and SCP.denominator == 1
        rig_sc = int(SCP) - nr_sc
        rig_nsc = int(P) - int(SCP) - nr_nsc
        idt = (((1, n),), 1),
        anti = ((((1, 1), (2, (n - 1) // 2)) if n % 2 else ((2, n // 2),)), 1),
        rows = []
        if rig_nsc:
            rows.append((n, rig_nsc, 1, idt, 0, ()))
        if rig_sc:
            rows.append((n, rig_sc, 1, idt, 1, anti))
        return rows, rig_nsc, rig_sc


def intransitivity(ledger_types, P, SCP, N):
    """c_s(n), c(n) for n <= N; P[n], SCP[n] needed for the order-n primes themselves"""
    c = [Fr(0)] * (N + 1)
    cs = [Fr(0)] * (N + 1)
    for n in range(1, N + 1):
        if n == 1:
            cs[1] = Fr(1)
        else:
            cser = lambda l: subst(c, l, n)
            tot = Fr(0)
            for k in sorted(ledger_types):
                if k >= n:
                    continue
                for (kk, cnt, G, aut, sc, anti) in ledger_types[k]:
                    if sc:
                        tot += cnt * (coef_ms(aut, cser, n) + coef_ms(anti, cser, n)) / (2 * G)
                    else:
                        tot += Fr(cnt, 2) * coef_ms(aut, cser, n) / G
            if n >= 3:
                tot += (P[n] + SCP[n]) / 2
            cs[n] = tot
        # c(n) = [x^n] exp(sum_k cs(x^k)/k) - 1 ; use the standard Euler-transform recursion
        # n c(n) = sum_{j=1}^n b(j) c(n-j) with c(0)=1, b(j) = sum_{d|j} d cs(d)
        cc = [Fr(1)] + c[1:n]
        b = [Fr(0)] * (n + 1)
        for j in range(1, n + 1):
            b[j] = sum(d * cs[d] for d in range(1, j + 1) if j % d == 0)
        c[n] = sum(b[j] * cc[n - j] for j in range(1, n + 1)) / n
    return cs, c


def labelled(N):
    T = [Fr(2 ** comb(m, 2), factorial(m)) for m in range(N + 1)]
    inv = [Fr(0)] * (N + 1)
    inv[0] = Fr(1)
    for m in range(1, N + 1):
        inv[m] = -sum(T[j] * inv[m - j] for j in range(1, m + 1))
    slab = [0] + [-inv[m] * factorial(m) for m in range(1, N + 1)]
    Tm1 = [Fr(0)] + T[1:]
    pw = {1: Tm1[:]}
    for k in range(2, N + 1):
        pw[k] = mul(pw[k - 1], Tm1, N)
    Plab = [Fr(0)] * (N + 1)
    npl = [Fr(0)] * (N + 1)
    for n in range(3, N + 1):
        npl[n] = factorial(n) * sum(Plab[k] / factorial(k) * pw[k][n] for k in range(3, n))
        Plab[n] = slab[n] - npl[n]
    # l: l_s = x + 1/2 sum_k P_lab(k) l^k/k!, 1 + l = exp(l_s)
    L = [Fr(0)] * (N + 1)
    Ls = [Fr(0)] * (N + 1)
    for n in range(1, N + 1):
        if n == 1:
            Ls[1] = Fr(1)
        else:
            tot = Fr(0)
            for k in range(3, n + 1):
                lk = [Fr(1)] + [Fr(0)] * n
                for _ in range(k):
                    lk = mul(lk, L[:n] + [Fr(0)], n)
                tot += Fr(1, 2) * Plab[k] / factorial(k) * lk[n]
            Ls[n] = tot
        # L[n] = [x^n] exp(Ls): use f' = Ls' f, f = 1 + L : n f_n = sum_{j=1}^n j Ls_j f_{n-j}
        f = [Fr(1)] + L[1:n]
        L[n] = sum(j * Ls[j] * f[n - j] for j in range(1, n + 1)) / n
    lab = [L[n] * factorial(n) for n in range(N + 1)]
    return slab, npl, Plab, lab


def fmt(x):
    if isinstance(x, Fr):
        if x.denominator != 1:
            return f"{x.numerator}/{x.denominator}"
        x = x.numerator
    return f"{x:,}"
