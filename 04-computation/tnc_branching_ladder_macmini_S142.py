#!/usr/bin/env python3
"""
The TNC BRANCHING LADDER: factor at each break, recurse, prune by nondegeneracy
                                                        (mac-mini-S142)
================================================================================
THM-1610 found that the linear peel handles exactly the r_j with M | j and stalls at
M !| j.  HYP-8440(a) asked whether the first broken equation, being a single quadratic
in one unknown, could be solved explicitly to re-linearise downstream.

IT CAN.  At (M,N) = (2,2) the m=3 equation is  [u^6]R^3 = -3 r3 (r1^3 - r3),  which
FACTORS into two linear branches, and BOTH die:
    r3 = 0      => r4 = -r1 r3 = 0, contradicting deg R = d (r_d != 0);
    r3 = r1^3   => r4 = -r1^4, and then [u^8]R^4 = -12 r1^8, forcing r1 = 0,
                   which sends r3, r4 to 0 -- same contradiction.
That is a complete elementary proof of TNC(2,2).

This file automates the method:
  peel linearly while possible; at a break, FACTOR the equation and branch on its
  irreducible factors; prune any branch forcing r_d = 0; declare the bidegree CLOSED
  when every branch dies.

Also recorded: the DUALITY.  Lambda(1/u) = R(1/u) u^M = u^{M-d} R*(u) = R*(u) u^{-N}
with R* the reversed polynomial, and CT is invariant under u -> 1/u.  Hence
    TNC(M,N)  <=>  TNC(N,M)   with R <-> R*,
so WLOG M <= N -- which is exactly why min(M,N) is the right index (cf. THM-1595).
"""
import sympy as sp
from fractions import Fraction as F

def polymul(a, b, cap):
    out = [0]*(min(len(a)+len(b)-1, cap+1))
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i+j > cap: break
            out[i+j] += x*y
    return out

def cMM(rs, M, m):
    """[u^{Mm}] R^m for R = 1 + sum r_i u^i."""
    cap = M*m; R = [1] + list(rs); P = [1]
    for _ in range(m): P = polymul(P, R, cap)
    return P[cap] if cap < len(P) else 0

def close_bidegree(M, N, MMAX=None, maxdepth=40, verbose=True):
    """Try to close TNC(M,N) by the branching ladder.  Returns (closed?, log)."""
    d = M + N
    r = sp.symbols(f'r1:{d+1}')
    if MMAX is None: MMAX = 2*d + 4
    eqs = [sp.expand(cMM(list(r), M, m)) for m in range(1, MMAX+1)]
    log = []
    # a branch is a substitution dict; process a stack
    stack = [({}, 0)]          # (subs, next equation index)
    alive = []
    steps = 0
    while stack:
        steps += 1
        if steps > 4000: return None, log + ["ABORTED: too many branches"]
        subs, idx = stack.pop()
        # non-degeneracy: r_d must not be forced to 0
        rd = sp.simplify(sp.sympify(r[d-1]).subs(subs))
        if rd == 0:
            log.append(f"    branch pruned (r_{d} = 0, degenerate)")
            continue
        if idx >= len(eqs):
            alive.append(subs); continue
        e = sp.expand(eqs[idx].subs(subs))
        if e == 0:
            stack.append((subs, idx+1)); continue
        free = sorted(e.free_symbols, key=lambda s: int(str(s)[1:]))
        if not free:
            log.append("    branch pruned (nonzero constant equation)")
            continue
        top = free[-1]
        deg = sp.Poly(e, top).degree()
        if deg == 1:
            sol = sp.solve(e, top)[0]
            ns = dict(subs); ns[top] = sp.simplify(sol)
            ns = {k: sp.simplify(v.subs(subs)) for k, v in ns.items()}
            stack.append((ns, idx+1))
        else:
            facs = sp.factor_list(e)[1]
            got = False
            for (fac, _mult) in facs:
                if top not in fac.free_symbols: continue
                sols = sp.solve(fac, top)
                for s in sols:
                    ns = dict(subs); ns[top] = sp.simplify(s)
                    ns = {k: sp.simplify(v.subs(subs)) for k, v in ns.items()}
                    stack.append((ns, idx+1)); got = True
            if not got:
                # equation involves top only nonlinearly with no rational roots
                log.append(f"    STALL at m={idx+1}: deg {deg} in {top}, no rational branch")
                alive.append(subs)
    return (len(alive) == 0), log

print("=" * 78)
print("THE DUALITY (free, and it explains the min(M,N) indexing)")
print("=" * 78)
print("  Lambda(1/u) = R(1/u) u^M = u^{M-d} R*(u) = R*(u) u^{-N},  R* the reversed poly.")
print("  CT is invariant under u -> 1/u, so CT(Lambda^m) = CT((Lambda*)^m) and")
print("      TNC(M,N)  <=>  TNC(N,M)  with R <-> R*.   WLOG M <= N.")
print()
# verify the duality numerically
def CTpow_lam(rs, M, m):
    return cMM(rs, M, m)
ok = True
for (M, N, rs) in ((2, 3, [F(1), F(2), F(-1), F(3), F(1)]),
                   (2, 4, [F(1), F(0), F(2), F(-1), F(1), F(2)]),
                   (3, 4, [F(2), F(1), F(-1), F(1), F(3), F(1), F(1)])):
    d = M+N
    star = list(reversed([F(1)] + rs))          # R*(u) coefficients, constant first
    c0 = star[0]
    starn = [c/c0 for c in star[1:]]            # normalise R*(0) = 1
    for m in range(1, 7):
        a = cMM(rs, M, m)
        b = cMM(starn, N, m) * c0**m
        if sp.simplify(a - b) != 0: ok = False
print(f"  duality CT(Lambda^m) = CT(Lambda*^m) verified numerically: {ok}")

print()
print("=" * 78)
print("THE BRANCHING LADDER, applied bidegree by bidegree (WLOG M <= N)")
print("=" * 78)
print(f"{'(M,N)':>8} {'d':>3} {'closed?':>9}   notes")
for (M, N) in ((1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(2,5),(3,3),(3,4)):
    try:
        closed, log = close_bidegree(M, N)
    except Exception as ex:
        closed, log = None, [f"ERROR {type(ex).__name__}"]
    stalls = [l for l in log if "STALL" in l]
    note = ("every branch dies" if closed
            else ("ABORTED" if closed is None else f"{len(stalls)} stall(s): "
                  + (stalls[0].strip() if stalls else "branch survived")))
    print(f"{str((M,N)):>8} {M+N:>3} {str(closed):>9}   {note}")

print()
print("=" * 78)
print("THE (2,2) PROOF, written out")
print("=" * 78)
r1, r2, r3, r4 = sp.symbols('r1 r2 r3 r4')
E = [sp.expand(cMM([r1, r2, r3, r4], 2, m)) for m in range(1, 5)]
print(f"  m=1: {E[0]} = 0            =>  r2 = 0")
E2 = sp.expand(E[1].subs({r2: 0}))
print(f"  m=2: {E2} = 0   =>  r4 = -r1*r3")
E3 = sp.factor(sp.expand(E[2].subs({r2: 0, r4: -r1*r3})))
print(f"  m=3: {E3} = 0   -- FACTORS")
print(f"       branch A: r3 = 0     =>  r4 = 0, contradicting r4 != 0 (deg R = 4)")
E4 = sp.factor(sp.expand(E[3].subs({r2: 0, r3: r1**3, r4: -r1**4})))
print(f"       branch B: r3 = r1^3  =>  r4 = -r1^4, and then m=4 gives {E4} = 0")
print(f"                 => r1 = 0 => r3 = r4 = 0, same contradiction.")
print("  Both branches die.  TNC(2,2) PROVED, elementarily, from m = 1,2,3,4 only.")
