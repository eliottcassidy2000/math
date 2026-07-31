#!/usr/bin/env python3
"""Exact (big-integer) solver for the AMM 12592 within-block deadline profile.

THE SYSTEM (THM-3008).  On the dyadic shell [m,2m), a balanced rule deciding
n = m+k by flip m+k+1+a_k exists iff

    sum_{k,i} p_{k,i} C(L_k, j-1-i) = C(m,j)/2,   1 <= j <= m-1,
    0 <= p_{k,i} <= C(a_k,i),      L_k = m-1-k-a_k.                     (*)

TWO STRUCTURAL FACTS make this solvable exactly and fast.

(1) EVERY EQUATION IS EXACTLY CENTRED.  Summing the boxes at level
    t := m-1-j gives  sum_k C(a_k+L_k, t-k+1) = sum_k C(m-1-k, t-k+1)
    = C(m, t+1)  by the hockey-stick identity -- exactly TWICE the target.
    So p = C/2 is always the (generally non-integral) centre, and (*) is a
    pure integrality/granularity question, never a capacity question.

(2) VARIABLES ENTER WITH COEFFICIENT 1.  Reading the equations from the top,
    t = 0,1,2,..., the variable p_{k,i} first appears at level
    t_0 = a_k + k - 1 - i with coefficient C(L_k, L_k) = 1.  The sole
    exception is p_{0,a_0}, which appears at t = 0 with coefficient L_0.
    Hence the system is TRIANGULAR: at each level one linear equation is
    closed and a fresh batch of unit-coefficient variables is introduced.

So a level-by-level solve is exact and complete up to the choice of how each
residual is split among the fresh variables; we try several splitting
policies and back off when a level goes out of range.

NECESSARY CONDITION (certified lower bounds).  Relaxing all equations but
one, level t is a bounded subset-sum: coins C(L_k, t-k+1-s) with
multiplicities C(a_k,s), target C(m,t+1)/2.  If the target misses the
achievable set, the profile is infeasible.  The classic case is t = 0, where
the coins are {L_0 (mult 1), 1 (mult a_0), 1 (mult 1)}: the achievable set is
[0,a_0+1] u [L_0, L_0+a_0+1], and the target m/2 lands in the gap unless
a_0 >= m/2 - 1.  That is exactly THM-2160 S6.2, in one line.
"""
import sys
from fractions import Fraction
from math import comb

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1000000)


# --------------------------------------------------------------- the system

def levels(m, a, tmax=None):
    """For each level t = 0..min(tmax,m-2) return the list (k, i, coef, box)."""
    L = [m - 1 - k - a[k] for k in range(m)]
    out = []
    hi = m - 1 if tmax is None else min(tmax, m - 1)
    for t in range(hi):
        row = []
        for k in range(m):
            for i in range(a[k] + 1):
                e = m - 2 - t - i
                if 0 <= e <= L[k]:
                    row.append((k, i, comb(L[k], e), comb(a[k], i)))
        out.append(row)
    return out, L


def target(m, t):
    return comb(m, t + 1) // 2


def check_centred(m, a):
    """Fact (1): every level's box total is exactly twice its target."""
    rows, L = levels(m, a)
    for t, row in enumerate(rows):
        if sum(c * b for _, _, c, b in row) != 2 * target(m, t):
            return False
    return True


# ------------------------------------------------- exact triangular solve

def _proportional(R, boxes):
    """x_j as close as possible to R * b_j / cap, exactly summing to R."""
    cap = sum(b for _, b in boxes)
    if cap == 0:
        return [0] * len(boxes)
    base = [(R * b) // cap for _, b in boxes]
    rem = R - sum(base)
    # hand the remainder to the largest fractional parts, respecting boxes
    order = sorted(range(len(boxes)),
                   key=lambda j: (R * boxes[j][1]) % cap, reverse=True)
    j = 0
    while rem > 0:
        idx = order[j % len(order)]
        if base[idx] < boxes[idx][1]:
            base[idx] += 1
            rem -= 1
        j += 1
    return base



def solve(m, a, policy="proportional", p0_choices=(1, 0)):
    """Exact level-by-level solve.  Returns (True, p) or (False, reason).

    Maintains the accumulated polynomial sum_k P_k(u)(1+u)^(L_k) incrementally:
    fixing p_{k,i} = v adds v * u^i * (1+u)^(L_k), an O(L_k) update.  Level t
    reads off the coefficient of u^(m-2-t)."""
    L = [m - 1 - k - a[k] for k in range(m)]
    rowpow = {}
    for k in range(m):
        if L[k] not in rowpow:
            rowpow[L[k]] = [comb(L[k], e) for e in range(L[k] + 1)]
    binom_a = [[comb(a[k], i) for i in range(a[k] + 1)] for k in range(m)]

    fresh = {t: [] for t in range(m - 1)}
    special = None
    for k in range(m):
        for i in range(a[k] + 1):
            t0 = a[k] + k - 1 - i
            if t0 < 0:
                special = (k, i)            # only ever (0, a_0)
            elif t0 <= m - 2:
                fresh[t0].append((k, i))

    for p0 in p0_choices:
        p = {}
        acc_poly = [0] * (m + 1)

        def add(k, i, v):
            if not v:
                return
            row = rowpow[L[k]]
            for e, c in enumerate(row):
                acc_poly[i + e] += v * c

        if special is not None:
            p[special] = p0
            add(special[0], special[1], p0)
        ok = True
        for t in range(m - 1):
            acc = acc_poly[m - 2 - t]
            R = target(m, t) - acc
            boxes = [(kv, binom_a[kv[0]][kv[1]]) for kv in fresh[t]]
            cap = sum(b for _, b in boxes)
            if R < 0 or R > cap:
                ok = False
                break
            # split R across the fresh unit-coefficient variables.
            # The exact centre p = C/2 solves every equation, so the policy
            # that preserves feasibility is to keep each variable at the same
            # fraction R/cap of its own box (proportional = stay centred).
            if policy == "proportional":
                assign = _proportional(R, boxes)
            elif policy == "greedy-first":
                assign, left = [], R
                for _, b in boxes:
                    take = min(b, left); assign.append(take); left -= take
            elif policy.startswith("alt"):
                # The extracted optima sit at BOX EXTREMES with alternating
                # signs (P_k ~ +-(-1)^y, the alternating atom), not near the
                # box centre.  Aim each fresh variable at +-its box with the
                # sign alternating in the level index, then correct to hit R.
                sgn = 1 if policy.endswith("+") else -1
                pref = []
                for idx, (kv, b) in enumerate(boxes):
                    s_i = sgn * (1 if (t + kv[0]) % 2 == 0 else -1)
                    pref.append(b if s_i > 0 else 0)
                need = R - sum(pref)
                assign = pref[:]
                j = 0
                while need != 0 and j < 4 * len(boxes) + 4:
                    for idx in range(len(boxes)):
                        b = boxes[idx][1]
                        if need > 0 and assign[idx] < b:
                            step = min(b - assign[idx], need)
                            assign[idx] += step; need -= step
                        elif need < 0 and assign[idx] > 0:
                            step = min(assign[idx], -need)
                            assign[idx] -= step; need += step
                        if need == 0:
                            break
                    j += 1
                if need != 0:
                    ok = False
                    break
            elif policy == "greedy-last":
                assign, left = [0] * len(boxes), R
                for idx in range(len(boxes) - 1, -1, -1):
                    take = min(boxes[idx][1], left); assign[idx] = take; left -= take
            else:
                raise ValueError(policy)
            assert sum(assign) == R
            for (kv, b), take in zip(boxes, assign):
                assert 0 <= take <= b
                p[kv] = take
                add(kv[0], kv[1], take)
        if ok:
            for kv in [(k, i) for k in range(m) for i in range(a[k] + 1)]:
                p.setdefault(kv, 0)
            return True, p
    return False, "range"


def verify(m, a, p):
    """Independent exact re-check of (*)."""
    L = [m - 1 - k - a[k] for k in range(m)]
    for k in range(m):
        for i in range(a[k] + 1):
            if not (0 <= p[(k, i)] <= comb(a[k], i)):
                return False
    for j in range(1, m):
        s = 0
        for k in range(m):
            for i in range(a[k] + 1):
                e = j - 1 - i
                if 0 <= e <= L[k]:
                    s += p[(k, i)] * comb(L[k], e)
        if s != comb(m, j) // 2:
            return False
    return True


# ----------------------------------------- necessary condition (lower bounds)

def sumset_hits(coins, tgt):
    """Can sum_{c,mult} x_c * c  equal tgt with 0 <= x_c <= mult?
       Exact when the small coins make the set gap-free; otherwise a sound
       (necessary-condition) test using the completeness criterion."""
    items = sorted(coins)                      # (value, mult)
    reach = 0                                  # [0, reach] fully achievable
    for val, mult in items:
        if val > reach + 1:
            # gap: everything in (reach, val) is unreachable from small coins
            if reach < tgt < val:
                return False
        reach += val * mult
    return 0 <= tgt <= reach


def necessary(m, a, tmax=12):
    """Level-by-level subset-sum necessary condition (top tmax levels)."""
    rows, L = levels(m, a, tmax)
    if tmax is None:
        tmax = len(rows)
    for t in range(min(tmax, len(rows))):
        coins = {}
        for _, _, c, b in rows[t]:
            coins[c] = coins.get(c, 0) + b
        if not sumset_hits(list(coins.items()), target(m, t)):
            return False, t
    return True, None


# ------------------------------------------------------------------- rho

def profile(m, C):
    a = []
    for k in range(m):
        n = m + k
        cap = Fraction(C) * n - n - 1
        if cap < 0:
            return None
        a.append(min(m - 1 - k, int(cap)))
    return a


def rho_bounds(m, policies=("alt+", "alt-", "proportional", "greedy-last", "greedy-first")):
    cands = sorted({Fraction(m + k + 1 + aa, m + k)
                    for k in range(m) for aa in range(m - k)})

    def feas(C):
        a = profile(m, C)
        if a is None:
            return False, None          # no legal depth profile at this ratio
        ok_nec, _ = necessary(m, a)
        if not ok_nec:
            return False, a
        for pol in policies:
            ok, p = solve(m, a, policy=pol)
            if ok:
                assert verify(m, a, p), ("witness failed exact check", m, a, pol)
                return True, a
        return None, a                    # undecided: no witness, no refutation

    lo, hi = 0, len(cands) - 1
    upper = None
    while lo <= hi:                       # least C with a verified witness
        mid = (lo + hi) // 2
        got, a = feas(cands[mid])
        if got is True:
            upper = (cands[mid], a)
            hi = mid - 1
        else:
            lo = mid + 1
    lower = None
    for C in cands:                       # greatest C certified infeasible
        got, a = feas(C)
        if got is False:
            lower = C
        else:
            break
    return lower, upper


if __name__ == "__main__":
    ms = [int(x) for x in sys.argv[1:]] or [4, 8, 16, 32, 64, 128, 256]
    print("m     rho lower bound      rho upper bound (verified witness)")
    for m in ms:
        lo, up = rho_bounds(m)
        lo_s = "-" if lo is None else f"{lo} = {float(lo):.6f}"
        if up is None:
            up_s = "none found"
        else:
            up_s = f"{up[0]} = {float(up[0]):.6f}"
        print(f"{m:5d}  > {lo_s:22s}  <= {up_s}")
        if up is not None:
            T = [m + k + 1 + up[1][k] for k in range(m)]
            print(f"       T(n) head: {T[:12]}{' ...' if m > 12 else ''}")
        sys.stdout.flush()
