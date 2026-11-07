# audit_H1 follow-up: independent checks of THM-4900(7): weights from the exponential walk (Session159 c_(r,q), not the recursion (R)), (b)1-(b)2 for e <= 200, r <= 100, the proof's Pascal split and slack, the (d) sandwiches, and the (c) error bound on random tournament pairs (orders 4-8) and on signed synthetic coefficient vectors.
import math
import random
import sys
from fractions import Fraction
from audit_tournament_lib import rand_tour, transitive, c3, cover_poly, ham_paths

RMAX, EMAX = 100, 202
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        if len(fails) <= 10:
            print("FAIL:", msg)


# ---------- weights U_(e,r) = sum_q c_(r,q) q^e from M e^(qy) = q e^((q+1)y) + q^2 e^(qy) - q e^((q-1)y) ----------
def walk(rmax):
    W = {1: 1}
    out = [dict(W)]
    for r in range(rmax):
        new = {}
        for q, c in W.items():
            if q == 0:
                continue                      # M 1 = 0
            for qq, wt in ((q + 1, q), (q, q * q), (q - 1, -q)):
                new[qq] = new.get(qq, 0) + wt * c
        W = new
        out.append(W)
    return out


WL = walk(400)
U = []
for r in range(RMAX + 1):
    row = [0] * (EMAX + 1)
    for q, c in WL[r].items():
        p = 1                                 # q^0 = 1, including 0^0 = 1 (needed for d = 0)
        for e in range(EMAX + 1):
            row[e] += c * p
            p *= q
    U.append(row)
check(U[0] == [1] * (EMAX + 1), "U_(e,0) = 1")

# consistency with the first-step recursion (R) (a different derivation)
badR = 0
for r in range(1, RMAX + 1):
    for e in range(0, EMAX - 1):
        v = U[r - 1][e + 2] + sum(2 * math.comb(e, 2 * k + 1) * U[r - 1][e - 2 * k] for k in range(e // 2 + 1) if 2 * k + 1 <= e)
        badR += v != U[r][e]
check(badR == 0, "(R) mismatch")
print(f"U_(e,r) from the exponential walk agrees with (R) for 1 <= r <= {RMAX}, 0 <= e <= {EMAX - 2}: mismatches {badR}")

# ---------- (b)1, (b)2 ----------
v1 = v2 = 0
minrat = {}
for r in range(0, RMAX + 1):
    best = None
    for e in range(2, 201):
        if not U[r][e] >= (r + 1) * U[r][e - 2]:
            v2 += 1
        rat = Fraction(U[r][e], (r + 1) * U[r][e - 2])
        if best is None or rat < best[0]:
            best = (rat, e)
        if r >= 1 and e + 2 <= EMAX and not U[r][e] >= U[r - 1][e + 2]:
            v1 += 1
    minrat[r] = best
check(v1 == 0 and v2 == 0, "(b)")
print(f"(b)1 U_(e,r) >= U_(e+2,r-1) and (b)2 U_(e,r) >= (r+1) U_(e-2,r), 2 <= e <= 200, 0 <= r <= {RMAX}: "
      f"violations {v1}, {v2}")
print("     min_e U_(e,r)/((r+1)U_(e-2,r)) [attained at e]:",
      "; ".join(f"r={r}: {float(minrat[r][0]):.3f} [e={minrat[r][1]}]" for r in (0, 1, 2, 3, 5, 10, 20, 50, 100)))
# equality cases of (b)2
eq = [(e, r) for r in range(0, RMAX + 1) for e in range(2, 201) if U[r][e] == (r + 1) * U[r][e - 2]]
print(f"     equality U_(e,r) = (r+1)U_(e-2,r) only at: {eq[:6]}{' ...' if len(eq) > 6 else ''} ({len(eq)} cases)")

# ---------- the proof's split, step by step (r >= 1, e >= 2) ----------
bad_split = bad_hyp = bad_slack = 0
for r in range(1, RMAX + 1):
    lam = r
    for e in range(2, 199):
        A1 = U[r - 1][e + 2] + sum(2 * math.comb(e - 2, 2 * k + 1) * U[r - 1][e - 2 * k] for k in range(e // 2 + 1))
        S = sum((4 * math.comb(e - 2, 2 * k) + (2 * math.comb(e - 2, 2 * k - 1) if k >= 1 else 0)) * U[r - 1][e - 2 * k]
                for k in range(e // 2 + 1))
        bad_split += (A1 + S != U[r][e])                    # Pascal split of (R)
        bad_hyp += not (A1 >= lam * U[r][e - 2])             # hypothesis applied at heights e+2 and e-2k >= 3
        bad_slack += not (S >= U[r][e - 2])                  # slack dominates U_(e-2,r)
        # heights where the hypothesis is used: e + 2 and e - 2k with C(e-2,2k+1) != 0
        used = [e + 2] + [e - 2 * k for k in range(e // 2 + 1) if math.comb(e - 2, 2 * k + 1) != 0]
        check(min(used) >= 2, "hypothesis used below height 2")
check(bad_split == bad_hyp == bad_slack == 0, "split")
print(f"proof steps for 1 <= r <= {RMAX}, 2 <= e <= 198: Pascal split exact (bad {bad_split}); "
      f"A1 >= r U_(e-2,r) (bad {bad_hyp}); S >= U_(e-2,r) (bad {bad_slack}); hypothesis only used at heights >= 3 or e+2")

# ---------- (d) ----------
bd1 = bd2 = 0
for r in range(1, RMAX + 1):
    for d in range(0, 61):
        x = Fraction(U[r][d], U[r - 1][d + 2])
        if not (1 <= x <= 1 + Fraction(2 ** d, r)):
            bd1 += 1
        if d >= 2:
            y = Fraction(U[r][d - 2] * U[r - 1][d + 2], U[r][d] * U[r - 1][d])     # rho_(d,r)/rho_(d+2,r-1)
            if not (1 / (1 + Fraction(2 ** d, r)) <= y <= 1 + Fraction(2 ** (d - 2), r)):
                bd2 += 1
check(bd1 == bd2 == 0, "(d)")
print(f"(d) 1 <= U_(d,r)/U_(d+2,r-1) <= 1 + 2^d/r (0 <= d <= 60) and 1/(1+2^d/r) <= rho_(d,r)/rho_(d+2,r-1) <= 1 + 2^(d-2)/r "
      f"(2 <= d <= 60), 1 <= r <= {RMAX}: violations {bd1}, {bd2}")

# ---------- (c) ----------
def Ud(d, r):
    return sum(c * (q ** d if d else 1) for q, c in WL[r].items())


Ucache = {}
def Uc(d, r):
    if (d, r) not in Ucache:
        Ucache[(d, r)] = Ud(d, r)
    return Ucache[(d, r)]


def check_pair(swT, swS, cT, cS, n, r):
    """returns (|E|, bound) for cores with signed cover coefficients sw (sw[j] = coefficient of q^(n-2j))."""
    sT = sum(abs(x) for x in swT[1:])
    sS = sum(abs(x) for x in swS[1:])
    HT = sum(swT[j] * Uc(n - 2 * j, r) for j in range(len(swT)))
    HS = sum(swS[j] * Uc(n - 2 * j, r) for j in range(len(swS)))
    rho = Fraction(Uc(n - 2, r), Uc(n, r))
    E = Fraction(HT, HS) - 1 - 2 * (cT - cS) * rho
    B = rho * (sT + sS + 2 * abs(cT - cS) * sS) / (r + 1 - sS)
    return abs(E), B


rng = random.Random(4907)
pool = {}
for n in range(4, 9):
    lst = [transitive(n)] + [rand_tour(n, rng) for _ in range(14 if n < 8 else 8)]
    if n % 2 == 1:
        lst.append([[1 if (j - i) % n in range(1, (n + 1) // 2) else 0 for j in range(n)] for i in range(n)])   # rotational
    pool[n] = [(A, cover_poly(A), c3(A)) for A in lst]
cnt = bad = 0
worst = Fraction(0)
sig_ok = True
for n, lst in pool.items():
    for A, co, cc in lst:
        sw = [co[n - 2 * j] for j in range(n // 2 + 1)]
        sig_ok &= (sum(sw[1:]) == ham_paths(A) - 1) and min(sw) >= 0
    for (A, coT, cT), (B_, coS, cS) in ((x, y) for x in lst for y in lst):
        swT = [coT[n - 2 * j] for j in range(n // 2 + 1)]
        swS = [coS[n - 2 * j] for j in range(n // 2 + 1)]
        sS = sum(swS[1:])
        rmin = sS                     # need r + 1 > sigma(T')
        for r in sorted({rmin, rmin + 1, rmin + 3, 2 * rmin + 1, min(400, 3 * rmin + 7), 400}):
            if r > 400 or r + 1 <= sS:
                continue
            e, b = check_pair(swT, swS, cT, cS, n, r)
            cnt += 1
            if e > b:
                bad += 1
            elif b:
                worst = max(worst, e / b)
check(bad == 0, "(c) tournaments")
print(f"(c) |E| <= rho[sigma(T)+sigma(T')+2|dc3|sigma(T')]/(r+1-sigma(T')) on {cnt} (T, T', r) triples "
      f"(random/transitive/rotational cores of orders 4-8, r from the threshold r+1 > sigma(T') to 400): violations {bad}, "
      f"max |E|/bound = {float(worst):.4f}")
print(f"    sigma(X) = H(X) - 1 and all sw_j >= 0 on the pool: {sig_ok}")

# signed synthetic coefficient vectors: the derivation only uses |sw_j|
cnt2 = bad2 = 0
worst2 = Fraction(0)
for trial in range(400):
    n = rng.randint(4, 9)
    m = n // 2
    def vec():
        v = [1] + [rng.randint(-12, 12) for _ in range(m)]
        return v
    swT, swS = vec(), vec()
    cT, cS = Fraction(swT[1], 2), Fraction(swS[1], 2)
    sS = sum(abs(x) for x in swS[1:])
    for r in sorted({sS, sS + 2, 2 * sS + 3, 300}):
        if r + 1 <= sS or r > 400:
            continue
        e, b = check_pair(swT, swS, cT, cS, n, r)
        cnt2 += 1
        if e > b:
            bad2 += 1
        elif b:
            worst2 = max(worst2, e / b)
check(bad2 == 0, "(c) signed")
print(f"(c) same bound with signed random sw_j (|sw_j| <= 12, c3 := sw_1/2) on {cnt2} triples: violations {bad2}, "
      f"max |E|/bound = {float(worst2):.4f}")

# ---------- fitted slope eps of c3_first_order_law vs 2 rho ----------
for n, N, eps_fit in ((4, 402, 3.9452e-4), (5, 401, 3.9612e-4), (6, 402, 3.9443e-4), (7, 401, 3.9595e-4)):
    r = (N - n) // 2
    rho = Fraction(Uc(n - 2, r), Uc(n, r))
    print(f"n = {n}, N = {N}: 2 rho_(n,r) = {float(2 * rho):.4e} vs the fitted eps = {eps_fit:.4e} "
          f"(relative difference {abs(eps_fit / float(2 * rho) - 1):.1e})")

print("FAILURES:", len(fails))
sys.exit(1 if fails else 0)
