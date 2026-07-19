#!/usr/bin/env python3
"""
death-star-2026-07-19-S59 -- HYP-7885 part 2: the SINGLE-FAR ABSORPTION ATLAS.

THE LEMMA (far-element absorption, single-far form; elementary):
  Let B be a finite set of positive integers, theta in (0, 1/2), and let
  l(B,theta) = length of the longest closed interval inside
  S(B,theta) = {t in [0,1] : ||u t|| >= theta for all u in B}.
  If x >= 2*theta / l(B,theta)  then M(B u {x}) >= theta.
  Proof: {t : ||x t|| < theta} is a disjoint union of OPEN arcs of width
  2theta/x centered at k/x, separated by gaps of width (1-2theta)/x > 0.
  A closed interval I of length >= 2theta/x cannot lie inside one open arc;
  if it meets two arcs it contains the whole safe gap between them; if it
  meets at most one arc it contains points outside every arc. Either way I
  contains t with ||x t|| >= theta.  Apply with I = the longest interval of
  S(B,theta).  QED

APPLICATION: N speeds, first-gap window (1/(N+1), 2/(2N+1)), theta = 2/(2N+1)
(the window's TOP -- M >= theta excludes the open window).  Base
B = B_{N,i} = {1..N} \ {i} has N-1 <= 12 speeds, so by SETTLED LRC(<=13)
M(B) >= 1/N > 2/(2N+1) = theta, hence S(B,theta) has nonempty interior and
X0 = ceil(2 theta / l) is finite.  Members of the single-defect single-far
stratum {1..N}\{i} u {x} therefore all have x <= X0(N,i): a FINITE CHECK
decides the stratum COMPLETELY.  This upgrades the stratum from 'searched'
(mac-mini-S26 ~9000 families empirical) to 'classified'.

Also emitted: the resonance atlas (exact M(B_i u {x}) for all x in the finite
range, with maximizer witness (q,a) and active pair) -- the data for the
general binder law extending THM-633 (i=N at N=12) and HYP-4516 (canonical
i=N-1, w=3(N-1)).

At N=13 we additionally flag members inside (1/14, 3/41) -- the interior of
opus THM-1240's interval -- to decide that question's single-far stratum.
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
from itertools import combinations
import sys, time

# ---------- exact M with witness ----------
def cand_denoms(S):
    Q = set()
    for v in S: Q.add(2*v)
    for x, y in combinations(S, 2):
        Q.add(x+y); Q.add(abs(x-y))
    Q.discard(0)
    return sorted(Q)

def M_exact_wit(S):
    """returns (M, q, a) exact; integer hot loop."""
    S = sorted(set(S))
    bn, bd, ba = 0, 1, 0
    for q in cand_denoms(S):
        for a in range(1, q//2 + 1):
            mn = q
            for v in S:
                r = (v*a) % q
                if r > q - r: r = q - r
                if r < mn:
                    mn = r
                    if mn * bd <= bn * q: break
            if mn * bd > bn * q:
                bn, bd, ba = mn, q, a
    return F(bn, bd), bd, ba

def M_exact(S, stop_above=None):
    S = sorted(set(S))
    bn, bd = 0, 1
    sa = stop_above
    for q in cand_denoms(S):
        for a in range(1, q//2 + 1):
            mn = q
            for v in S:
                r = (v*a) % q
                if r > q - r: r = q - r
                if r < mn:
                    mn = r
                    if mn * bd <= bn * q: break
            if mn * bd > bn * q:
                bn, bd = mn, q
                if sa is not None and bn * sa.denominator >= sa.numerator * bd:
                    return F(bn, bd)
    return F(bn, bd)

# ---------- exact safe-interval computation ----------
def safe_intervals(B, theta):
    """S(B,theta) = [0,1] minus union over u in B, k of open arcs
    ((k-theta)/u, (k+theta)/u).  Returns list of closed intervals [a,b] (a<=b,
    Fractions) and the max length.  Exact rational sweep."""
    events = []   # (pos, +1 open-danger / -1 close-danger)
    for u in B:
        for k in range(0, u+1):
            lo = (F(k) - theta) / u
            hi = (F(k) + theta) / u
            events.append((max(lo, F(0)), 1))
            events.append((min(hi, F(1)), -1))
    events.sort()
    out = []
    depth = 0
    prev = F(0)
    for pos, d in events:
        if depth == 0 and pos > prev:
            out.append((prev, pos))
        depth += d
        if depth == 0:
            prev = pos
        elif depth < 0:
            raise RuntimeError("sweep imbalance")
    if prev < F(1):
        out.append((prev, F(1)))
    lens = [b - a for a, b in out]
    lmax = max(lens) if lens else F(0)
    return out, lmax

def primitive(S):
    g = reduce(gcd, S)
    return tuple(sorted(v//g for v in S))

def main():
    log = lambda s="": print(s, flush=True)
    log("== HYP-7885 part 2: single-far absorption atlas (death-star-S59) ==")
    log("Lemma check + per-(N,i) X0 + complete finite classification.\n")

    # empirical lemma spot-check first (belt and braces)
    B = [v for v in range(1, 13) if v != 5]
    th = F(2, 25)
    iv, l = safe_intervals(B, th)
    x0 = ceil(F(2)*th / l)
    log(f"  sanity: N=12 i=5 base safe-set: {len(iv)} intervals, l = {l} "
        f"= {float(l):.5f}, X0 = {x0}")
    viol = 0
    for x in range(x0, x0 + 40):
        if x in B: continue
        m = M_exact(B + [x], stop_above=th)
        if m < th: viol += 1
    log(f"  lemma spot-check above X0 (40 x's): violations = {viol} (want 0)\n")

    grand = {}
    for N in range(6, 14):
        theta = F(2, 2*N+1)
        lo = F(1, N+1)
        members = []
        log(f"-- N = {N}: window ({lo}, {theta}) --")
        for i in range(1, N+1):
            Bset = [v for v in range(1, N+1) if v != i]
            iv, l = safe_intervals(Bset, theta)
            if l == 0:
                log(f"   i={i:>2}: SAFE SET HAS NO INTERVAL (unexpected; M(B)<=theta?)")
                continue
            X0 = ceil(F(2)*theta / l)
            # finite complete check + margin band for empirical cross-validation
            found_i = []
            hi_check = max(X0, 12*N) + 25
            for x in range(N+1, hi_check+1):
                if x in Bset: continue
                Sfam = Bset + [x]
                m = M_exact(Sfam, stop_above=theta)
                if lo < m < theta:
                    Mw, q, a = M_exact_wit(Sfam)
                    inside = (x <= X0)
                    found_i.append((x, Mw, q, a, inside))
            tag = ""
            if found_i:
                for x, Mw, q, a, inside in found_i:
                    members.append((i, x, Mw, q))
                    star = "" if inside else "  ** ABOVE X0 (lemma violation!)"
                    log(f"   i={i:>2}: x={x:>4}  M={Mw} (q={q}, a={a}){star}")
            # compact line per i
            log(f"   i={i:>2}: l={float(l):.5f} X0={X0:>4} checked x<= {hi_check}"
                f"  members: {len(found_i)}")
        grand[N] = members
        if N == 13:
            innermost = [(i,x,Mw) for (i,x,Mw,q) in members if F(1,14) < Mw < F(3,41)]
            log(f"   N=13 members inside (1/14, 3/41): {innermost if innermost else 'NONE'}")
        log(f"   => N={N} single-far stratum members (complete): "
            f"{[(i, x, str(Mw)) for i, x, Mw, q in members] if members else 'EMPTY'}\n")

    log("== SUMMARY: complete single-far classification, N = 6..13 ==")
    for N in sorted(grand):
        ms = grand[N]
        s = ", ".join(f"(i={i},x={x},M={Mw})" for i, x, Mw, q in ms) if ms else "EMPTY"
        log(f"  N={N:>2}: {s}")

if __name__ == "__main__":
    main()
