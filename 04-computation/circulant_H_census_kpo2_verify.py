# circulant_H_census_kpo2_verify.py
# ADVERSARIAL VERIFIER for session kind-pasteur-2026-06-10-S2, thread A.
# Independent re-verification of claims A7, A8, A9, A11 (H counts, R ratios,
# circulant census n = 7..19).
#
# FRESH method:
#  * Hamiltonian-path counts via a numpy popcount-LAYERED Held-Karp
#    (dp[(S,v)] = #paths starting at vertex 0 covering S ending at v),
#    fixed start justified INDEPENDENTLY by vertex-transitivity of circulants
#    (rotation x->x+1 is an automorphism, so #paths starting at v is the same
#    for every v; total H = n * H_from_0). NOT the worker's code.
#  * validated against a from-scratch pure-Python ALL-STARTS subset DP at
#    n = 7, 9, 11 (no fixed-start trick, no symmetry), and against canon
#    values H_paley(7)=189, H_paley(11)=95095, A038375(9)=3357 (THM-329/064).
#  * census symmetry reduction by my own orbit enumeration under the unit
#    action S -> u*S; orbit-invariance of H itself adversarially re-checked by
#    computing H on EVERY member of every orbit for n <= 13 and on one full
#    orbit at n = 17.
# Exact integer arithmetic (int64 safe: max H ~ 1.2e12 << 2^63).
# verifier: kind-pasteur-o2-verify

import sys, math
from math import gcd, factorial
import numpy as np

CHECKS = []
def check(name, ok, detail=""):
    CHECKS.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))

def is_tournament_set(S, n):
    return all(((d in S) + ((n - d) in S)) == 1 for d in range(1, n))

def all_circ_tournaments(n):
    out = []
    for mask in range(1 << (n - 1)):
        S = frozenset(d for d in range(1, n) if (mask >> (d - 1)) & 1)
        if is_tournament_set(S, n):
            out.append(S)
    return out

def orbits_under_units(sets_list, n):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    seen = set()
    orbs = []
    for S in sets_list:
        if S in seen: continue
        orb = {frozenset((u * d) % n for d in S) for u in units}
        seen |= orb
        orbs.append(sorted(orb, key=lambda s: sorted(s)))
    return orbs

# ---------------- pure-Python ALL-STARTS Held-Karp (validator; no tricks)
def H_full_python(S, n):
    adj = [[1 if ((b - a) % n) in S else 0 for b in range(n)] for a in range(n)]
    out_nb = [[b for b in range(n) if adj[a][b]] for a in range(n)]
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(1, n):
        ndp = {}
        for (Sm, v), c in dp.items():
            if bin(Sm).count("1") != size: continue
            for w in out_nb[v]:
                if not (Sm >> w) & 1:
                    k = (Sm | (1 << w), w)
                    ndp[k] = ndp.get(k, 0) + c
        # merge layers (keep all for popcount filter trick)
        for k, c in ndp.items():
            dp[k] = dp.get(k, 0) + c
    full = (1 << n) - 1
    return sum(c for (Sm, v), c in dp.items() if Sm == full)

# ---------------- numpy layered fixed-start Held-Karp
_layer_cache = {}
def layers_with_bit0(n):
    if n in _layer_cache: return _layer_cache[n]
    size = 1 << n
    a = np.arange(size, dtype=np.int64)
    pc = np.zeros(size, dtype=np.int8)
    for b in range(n):
        pc += ((a >> b) & 1).astype(np.int8)
    has0 = (a & 1).astype(bool)
    layers = [a[(pc == k) & has0] for k in range(n + 1)]
    _layer_cache[n] = layers
    return layers

def H_circulant(S, n):
    """total #directed Hamiltonian paths = n * (#paths starting at vertex 0)"""
    layers = layers_with_bit0(n)
    out_nb = [[w for w in range(n) if ((w - v) % n) in S] for v in range(n)]
    dp = np.zeros((1 << n, n), dtype=np.int64)
    dp[1, 0] = 1
    for k in range(1, n):
        sub = layers[k]
        for v in range(n):
            colv = dp[sub, v]
            nzm = colv != 0
            if not nzm.any(): continue
            src = sub[nzm]; vals = colv[nzm]
            for w in out_nb[v]:
                m2 = ((src >> w) & 1) == 0
                if not m2.any(): continue
                tgt = src[m2] | (1 << w)
                dp[tgt, w] += vals[m2]
    full = (1 << n) - 1
    h0 = int(dp[full, :].sum())
    return n * h0

def paley_set(p):
    assert p % 4 == 3
    return frozenset({(x * x) % p for x in range(1, p)})

def rot_set(n):
    return frozenset(range(1, (n - 1) // 2 + 1))

# =================================================================
print("=" * 72)
print("VALIDATION: numpy fixed-start*n vs pure-Python all-starts DP")
print("=" * 72)
for n, S in [(7, paley_set(7)), (7, rot_set(7)), (9, rot_set(9)),
             (9, frozenset({1, 2, 3, 5})), (11, paley_set(11)), (11, rot_set(11))]:
    hf = H_full_python(S, n)
    hc = H_circulant(S, n)
    check(f"validator n={n} S={sorted(S)}", hf == hc, f"full={hf} fixed-start*n={hc}")

check("canon H_paley(7)=189", H_circulant(paley_set(7), 7) == 189)
check("canon H_paley(11)=95095", H_circulant(paley_set(11), 11) == 95095)
check("canon/THM-064 H({1,2,3,5}@9)=3357", H_circulant(frozenset({1, 2, 3, 5}), 9) == 3357)

# =================================================================
print("=" * 72)
print("CENSUS n = 7..19: orbits under unit action, exact H per class")
print("=" * 72)
claimed_nclasses = {7: 2, 9: 4, 11: 4, 13: 6, 15: 16, 17: 16, 19: 30}
results = {}
for n in [7, 9, 11, 13, 15, 17, 19]:
    sets_list = all_circ_tournaments(n)
    orbs = orbits_under_units(sets_list, n)
    check(f"census n={n}: 2^((n-1)/2) tournaments = {2**((n-1)//2)}, #classes",
          len(sets_list) == 2 ** ((n - 1) // 2) and len(orbs) == claimed_nclasses[n],
          f"#classes={len(orbs)} (claimed {claimed_nclasses[n]})")
    rows = []
    for orb in orbs:
        if n <= 13:
            hs = {H_circulant(S, n) for S in orb}   # H on EVERY member
            ok_inv = (len(hs) == 1)
            h = hs.pop()
            if not ok_inv:
                check(f"census n={n}: H constant on orbit {sorted(orb[0])}", False)
        else:
            h = H_circulant(orb[0], n)
        rows.append((h, sorted(orb[0]), len(orb)))
    if n <= 13:
        check(f"census n={n}: H constant on every orbit (all members computed)", True)
    rows.sort(reverse=True)
    results[n] = rows
    for h, rep, osz in rows[:4]:
        tag = ""
        if frozenset(rep) == rot_set(n) or any(frozenset(rep) == frozenset((u * d) % n for d in rot_set(n)) for u in range(1, n) if gcd(u, n) == 1):
            tag += " [rot-class]"
        if n in (7, 11, 19) and any(frozenset(rep) == frozenset((u * d) % n for d in paley_set(n)) for u in range(1, n) if gcd(u, n) == 1):
            tag += " [paley-class]"
        print(f"  n={n}: H={h:>15,}  rep={rep} (orbit size {osz}){tag}")
# one full orbit at n=17 as spot-check of unit invariance at larger n
orb17 = orbits_under_units([rot_set(17)], 17)[0]
hs17 = {H_circulant(S, 17) for S in orb17}
check("census n=17: H constant on the full rotation orbit", len(hs17) == 1, f"orbit size {len(orb17)}")

# =================================================================
print("=" * 72)
print("A9: maximizers and the n=19 top-3")
print("=" * 72)
def class_of(S, n):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    return {frozenset((u * d) % n for d in S) for u in units}

for n, expect_max_is in [(7, "paley"), (11, "paley"), (13, "rot"), (15, "rot"), (17, "rot"), (19, "rot")]:
    hmax, rep, _ = results[n][0]
    cls = class_of(frozenset(rep), n)
    if expect_max_is == "paley":
        ok = paley_set(n) in cls
    else:
        ok = rot_set(n) in cls
    check(f"A9 n={n}: circulant H-maximizer is {expect_max_is}", ok, f"max H={hmax:,} rep={rep}")

check("A9 n=9 maximizer class contains {1,2,3,5}, H=3357",
      frozenset({1, 2, 3, 5}) in class_of(frozenset(results[9][0][1]), 9) and results[9][0][0] == 3357,
      f"max={results[9][0]}")
check("A9 H_paley(7)=189 is census max at n=7", results[7][0][0] == 189)
check("A9 H_paley(11)=95095 is census max at n=11", results[11][0][0] == 95095)
check("A9 H_rot(13)=3711175 (MISTAKE-011b value) is census max at n=13", results[13][0][0] == 3711175,
      f"max={results[13][0][0]:,}")

h19 = results[19]
print("  n=19 top 5 classes:")
for h, rep, osz in h19[:5]:
    print(f"    H={h:>15,}  rep={rep}")
hrot19 = H_circulant(rot_set(19), 19)
hpaley19 = H_circulant(paley_set(19), 19)
hmid19 = H_circulant(frozenset({1, 2, 3, 4, 5, 6, 7, 8, 10}), 19)
check("A9 H_rot(19) = 1,184,212,824,763", hrot19 == 1184212824763, f"{hrot19:,}")
check("A9 H({1..8,10}@19) = 1,178,609,421,219", hmid19 == 1178609421219, f"{hmid19:,}")
check("A9 H_paley(19) = 1,172,695,746,915", hpaley19 == 1172695746915, f"{hpaley19:,}")
check("A9 ordering rot > {1..8,10} > paley at n=19", hrot19 > hmid19 > hpaley19)
check("A9 gap rot-paley = 11,517,077,848", hrot19 - hpaley19 == 11517077848, f"{hrot19-hpaley19:,}")
top3 = [h for h, _, _ in h19[:3]]
check("A9 census top-3 at n=19 = (rot, {1..8,10}-class, paley) i.e. Paley ranks THIRD",
      top3 == [hrot19, hmid19, hpaley19],
      f"top3={[f'{h:,}' for h in top3]}")

# =================================================================
print("=" * 72)
print("A7/A8: exact R_rot(n) = H * 2^(n-1)/n! vs naive exp(tanh 1)")
print("=" * 72)
hrot = {}
for n in range(5, 20, 2):
    hrot[n] = H_circulant(rot_set(n), n)
check("A7 H_rot(17) = 13,689,269,499", hrot[17] == 13689269499, f"{hrot[17]:,}")
check("A7 H_rot(19) = 1,184,212,824,763", hrot[19] == 1184212824763, f"{hrot[19]:,}")
claimed_R = {5: 2.00000, 7: 2.22222, 9: 2.30476, 11: 2.38646, 13: 2.44113,
             15: 2.48658, 17: 2.52227, 19: 2.55197}
naive = math.exp(math.tanh(1.0))
print(f"  exp(tanh 1) = {naive:.6f}  (worker claims 2.14169)")
check("A7 exp(tanh 1) = 2.14169 (5 dp)", abs(naive - 2.14169) < 6e-6, f"{naive:.6f}")
prev = 0.0
ok_all = True
for n in range(5, 20, 2):
    R = hrot[n] * 2 ** (n - 1) / factorial(n)
    ok = abs(R - claimed_R[n]) < 6e-6
    ok_all &= ok
    print(f"  n={n:>2}: H_rot={hrot[n]:>15,}  R={R:.5f}  claimed {claimed_R[n]:.5f}  {'ok' if ok else 'MISMATCH'}")
    assert R > prev; prev = R
check("A7 all exact R_rot(n) n=5..19 match claimed values (5 dp) and increase", ok_all)
check("A7 R_rot exceeds naive exp(tanh 1) from n=9 on",
      hrot[7] * 2 ** 6 / factorial(7) < naive < hrot[9] * 2 ** 8 / factorial(9))
check("A8 status check: R_rot(19)=2.55197 < e=2.71828, gap still " +
      f"{math.e - hrot[19]*2**18/factorial(19):.4f}", hrot[19] * 2 ** 18 / factorial(19) < math.e)

# =================================================================
print("=" * 72)
print("A11: H equalities for the relabeled sets")
print("=" * 72)
h_block9 = H_circulant(frozenset({1, 2, 5, 6}), 9)
h_rot9 = hrot[9] if 9 in hrot else H_circulant(rot_set(9), 9)
check("A11 H({1,2,5,6}@9) = H(rot_9) = 3267", h_block9 == h_rot9 == 3267, f"{h_block9}, {h_rot9}")
h_block13 = H_circulant(frozenset({1, 2, 3, 7, 8, 9}), 13)
check("A11 H({1,2,3,7,8,9}@13) = H(rot_13) = 3711175",
      h_block13 == hrot[13] == 3711175, f"{h_block13}, {hrot[13]}")

print("=" * 72)
nf = sum(1 for _, ok in CHECKS if not ok)
print(f"TOTAL {len(CHECKS)} checks, {nf} failures")
sys.exit(1 if nf else 0)
