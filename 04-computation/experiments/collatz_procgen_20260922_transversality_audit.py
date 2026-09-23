#!/usr/bin/env python3
"""Independent audit of Theorem R / Theorem S (transversality lane, collatz-procgen-20260922).

Independent code path: the Bernstein number is computed from Bernstein's closed formula
    Phi(w) = - sum_{l>=1} 2^(d_l) 3^(-l)    (d_1<d_2<... the positions of the 1s),
not from the affine recursion R_(ze)=m_e R_z + r_e 2^|z| used by the lane; the periodic
approximants Phi(u v^inf) are summed as exact geometric series with Fractions.
Checks, for the shortcut map T(x)=x/2, (3x+1)/2:
  (1) the residue X = Phi(s) mod 2^N really has parity prefix s[0,N) (T iterated mod 2^(N-t));
  (2) Lemma J: s(i+q)=s(i) whenever i, i+1 are not in J; J-gaps >= q_(n+1); min J < q_n+q_(n+1);
  (3) the isometry: v_2(Phi(s) - P/Q) equals the common-prefix length lambda(s, u v^inf);
  (4) the proof's lambda lower bounds for Options A and B, and the gain
      g = lambda - log2(|P|+|Q|) (reduced fraction) against q_(n+1)/mu - (mu-1) q_n.
"""
import sys, math
from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 90
LOG32 = Decimal(2).ln() / Decimal(3).ln()          # log_3 2 = 0.6309...

def mech(alpha, rho, n):
    """lower mechanical word s(i)=floor((i+1)a+rho)-floor(i a+rho), i<n."""
    out = []
    prev = int((Decimal(0) * alpha + rho).to_integral_value(rounding='ROUND_FLOOR'))
    for i in range(n):
        cur = int((Decimal(i + 1) * alpha + rho).to_integral_value(rounding='ROUND_FLOOR'))
        out.append(cur - prev); prev = cur
    return out

def margin_ok(alpha, rho, n, eps):
    # certify every letter: frac(i a + rho) stays > eps from 0 and 1 (eps >> precision error);
    # the exact integer point i=0, rho=0 is exact and skipped
    for i in range(n + 1):
        if i == 0 and rho == 0: continue
        f = (Decimal(i) * alpha + rho) % 1
        if f < eps or 1 - f < eps:
            return False, i
    return True, None

def bernstein_mod(w, N, qm=3, r=1):
    """Phi_T(w) mod 2^N for T(x)=x/2 (even), (qm x + r)/2 (odd): Phi = -r sum_l 2^(d_l) qm^(-l)."""
    M = 1 << N; invq = pow(qm, -1, M); x = 0; l = 0; pq = 1
    for d, b in enumerate(w[:N]):
        if b:
            l += 1; pq = pq * invq % M
            x = (x - r * (1 << d) * pq) % M
    return x

def parity_prefix(x, N, qm=3, r=1):
    bits = []
    for t in range(N):
        b = x & 1; bits.append(b)
        x = (qm * x + r) >> 1 if b else x >> 1
        x &= (1 << (N - t - 1)) - 1 if N - t - 1 > 0 else 0
    return bits

def phi_ep(u, v, qm=3, r=1):
    """Phi_T(u v^inf) exactly, from the closed formula summed as a geometric series."""
    s = Fraction(0); l = 0
    for d, b in enumerate(u):
        if b:
            l += 1; s -= Fraction(r * 2 ** d, qm ** l)
    ku = l; p = len(v); kv = sum(v)
    blk = Fraction(0); i = 0
    for d, b in enumerate(v):
        if b:
            i += 1; blk -= Fraction(r * 2 ** (len(u) + d), qm ** (ku + i))
    ratio = Fraction(2 ** p, qm ** kv)          # 2-adically |ratio| < 1 (p >= 1)
    return s + blk / (1 - ratio)

def v2_diff(X, fr, N):
    P, Q = fr.numerator, fr.denominator
    assert Q % 2 == 1
    M = 1 << N
    d = (X - P * pow(Q, -1, M)) % M
    if d == 0: return N
    return (d & -d).bit_length() - 1

def common_prefix(w, u, v, L):
    n = 0
    while n < L:
        t = n
        c = u[t] if t < len(u) else v[(t - len(u)) % len(v)]
        if c != w[t]: break
        n += 1
    return n

def cf_convergent_dens(alpha, K):
    a = []; x = alpha
    for _ in range(K):
        ai = int(x); a.append(ai); fr = x - ai
        if fr == 0: break
        x = 1 / fr
    q = [1, a[1]] if len(a) > 1 else [1]
    for k in range(2, len(a)):
        q.append(a[k] * q[-1] + q[-2])
    return a, q

ROWS = [
    # (label, qm, r, alpha(Decimal), rho)
    ("3x+1 critical log_3 2, rho=0",      3, 1, LOG32, Decimal(0)),
    ("3x+1 critical log_3 2, rho=1/2",    3, 1, LOG32, Decimal(1) / 2),
    ("3x+1 critical, rho=(3-sqrt5)/2",    3, 1, LOG32, (3 - Decimal(5).sqrt()) / 2),
    ("3x+1 supercritical sqrt2/2, rho=0", 3, 1, Decimal(2).sqrt() / 2, Decimal(0)),
    ("3x-1 supercritical sqrt2/2, rho=1/3", 3, -1, Decimal(2).sqrt() / 2, Decimal(1) / 3),
    ("3x+1 subcritical 1/phi, rho=0",     3, 1, (Decimal(5).sqrt() - 1) / 2, Decimal(0)),
    ("5x+1 alpha=1/phi (mu=1.435), rho=0", 5, 1, (Decimal(5).sqrt() - 1) / 2, Decimal(0)),
    ("5x+1 alpha=sqrt2/2 (mu=1.642), rho=0", 5, 1, Decimal(2).sqrt() / 2, Decimal(0)),
    ("CONTROL 5x+1 alpha=0.9+ (mu~2.1 > thresholds), rho=0", 5, 1, Decimal(9) / 10 + Decimal(2).sqrt() / 1000, Decimal(0)),
]

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    L = 3 * N
    print(f"Independent audit of Theorem R/S; N={N} bits of Phi, words of length {L}.")
    print("gain = lambda - log2(|P|+|Q|) (reduced P/Q); bound = q_(n+1)/mu - (mu-1) q_n (proof, up to -O(log q_(n+1)));")
    print("'cap' = lambda reached N (true lambda larger). A row passes if every check is True.\n")
    for (lab, qm, r, alpha, rho) in ROWS:
        a, qs = cf_convergent_dens(alpha, 30)
        mu = max(1.0, float(alpha) * math.log2(qm))
        ok, bad = margin_ok(alpha, rho, L, Decimal('1e-60'))
        s = mech(alpha, rho, L)
        X = bernstein_mod(s, N, qm, r)
        pp = parity_prefix(X, N, qm, r)
        print(f"== {lab}: mu={mu:.4f} mu(mu-1)={mu*(mu-1):.4f} (<phi? {mu*(mu-1) < 1.6180339887}); "
              f"pq={a[1:12]}; letters certified={ok}; parity prefix of Phi residue == word: {pp == s[:N]}")
        allok = ok and pp == s[:N]
        best = []
        for qi, q in enumerate(qs):
            if q < 2 or q > N // 2 or qi + 1 >= len(qs): continue
            qn1 = qs[qi + 1]
            pn = int((Decimal(q) * alpha).to_integral_value())
            delta = Decimal(q) * alpha - pn
            J = [j for j in range(L) if
                 int((Decimal(j) * alpha + rho + delta).to_integral_value(rounding='ROUND_FLOOR')) !=
                 int((Decimal(j) * alpha + rho).to_integral_value(rounding='ROUND_FLOOR'))]
            Js = set(J)
            per_ok = all(s[i + q] == s[i] for i in range(L - q - 1) if i not in Js and (i + 1) not in Js)
            gap_ok = all(J[k + 1] - J[k] >= qn1 for k in range(len(J) - 1))
            if not J:
                print(f"   q_n={q}: J empty in window (q_n+1={qn1} > window)"); continue
            j0 = J[0]; j1 = J[1] if len(J) > 1 else None
            first_ok = j0 < q + qn1
            uA, vA = [], s[0:q]
            frA = phi_ep(uA, vA, qm, r)
            lamA = common_prefix(s, uA, vA, L); lamA2 = v2_diff(X, frA, N)
            gA = min(lamA, N) - math.log2(abs(frA.numerator) + frA.denominator)
            uB, vB = s[0:j0 + 1], s[j0 + 1:j0 + 1 + q]
            frB = phi_ep(uB, vB, qm, r)
            lamB = common_prefix(s, uB, vB, L); lamB2 = v2_diff(X, frB, N)
            gB = min(lamB, N) - math.log2(abs(frB.numerator) + frB.denominator)
            iso = (lamA2 == min(lamA, N)) and (lamB2 == min(lamB, N))
            lamok = lamA >= j0 - 1 + q and (j1 is None or lamB >= j1 - 1 + q)
            bound = qn1 / mu - (mu - 1) * q
            cap = " cap" if max(lamA, lamB) >= N else ""
            rowok = per_ok and gap_ok and first_ok and iso and lamok
            allok = allok and rowok
            print(f"   q_n={q:6d} q_n+1={qn1:7d} j0={j0:6d}: LemmaJ(period,gaps,first)={per_ok},{gap_ok},{first_ok} "
                  f"lambda A/B={lamA}/{lamB} bounds ok={lamok} isometry={iso} | gain={max(gA, gB):9.1f} bound={bound:9.1f}{cap}")
        print(f"   ROW {'PASS' if allok else 'FAIL'}\n")

if __name__ == '__main__':
    main()
