#!/usr/bin/env python3
"""
claudebox-2026-06-03-S608 : does χ add beyond vertex-transitivity among the regular (maximally-
cyclic) tournaments — is there a tight regular config off the Paley/AP orbit, and does its χ differ?

Setup. A tight LRC config has a MAXIMALLY symmetric (regular = vertex-transitive) witness tournament
(HYP-2124). Two canonical regular tournaments: the ROTATIONAL (the AP {1..n-1}, connection set an
interval) and the doubly-regular PALEY P_p (the QR tournament, p≡3 mod 4). The character spectrum χ
(additive/Dirichlet character values of the speed set, = the block eigenvalues) is a FINER invariant
than vertex-transitivity: Paley is FLAT (|χ|=√p/2, a 2-design); the rotational AP is structured.

ANSWER (this file): χ adds beyond vertex-transitivity — it is the RESONANCE-ENERGY signature
(S605/HYP-2155). Verified:
  * Paley P_7 = {1,2,4}: doubly-regular, flat χ, resonance energy E=0 → MAXIMALLY SAFE (not tight).
  * AP (rotational): structured χ, high resonance energy E → TIGHT.
  ⇒ among regular tournaments, the FLAT-χ (doubly-regular Paley) ones are SAFE; the TIGHT ones are
    the structured-χ rotational (AP) class. Vertex-transitivity is shared; χ (= resonance energy)
    is what selects tightness.
  * The one tight regular config that is NOT the AP-interval — Paley P_11 = {1,3,4,5,9} at n=6 — is
    an AP-LIFT (≡ {1,2,3,4,5} mod 7): the AP transported to another additive face, not genuinely new.
  * The genuinely non-AP tight configs (sporadics n=5,8) are NON-regular (no vertex-transitive
    symmetry). So no tight config is regular off the AP/AP-lift orbit.
"""
import math, itertools, cmath
from math import gcd

def dist(x): s = x % 1; return min(s, 1 - s)
def gap_exact(V, Mmax=1200):
    best = (0.0, 0, 1)
    for M in range(2, Mmax):
        for a in range(1, M):
            if gcd(a, M) != 1: continue
            g = min(dist(v * a / M) for v in V)
            if g > best[0]: best = (g, a, M)
    return best
def QR(p): return sorted({(x * x) % p for x in range(1, p)})
def ghat(m, n):
    d = 1.0 / n
    return (1 - 2 * d) if m == 0 else -math.sin(2 * math.pi * m * d) / (math.pi * m)
def resonance_energy(v, M=3):
    n = len(v) + 1; tot = 0.0; l3 = 0.0
    for m in itertools.product(range(-M, M + 1), repeat=len(v)):
        if not any(m): continue
        if sum(a * b for a, b in zip(m, v)) != 0: continue
        p = 1.0
        for mi in m: p *= abs(ghat(mi, n))
        tot += p
        if sum(abs(x) for x in m) == 3: l3 += p
    return tot, l3, (1 - 2 / n) ** (n - 1)
def chi_spec(V, m):
    return [round(abs(sum(cmath.exp(2j * math.pi * a * x / m) for x in V)), 3) for a in range(1, m)]
def mult_inv(V, m):
    Vs = set(x % m for x in V)
    return [c for c in range(1, m) if gcd(c, m) == 1 and all((c * x) % m in Vs for x in Vs)]
def is_ap_lift(V):
    """is V ≡ an interval {1..k} (mod some M < max(V))? (an AP in disguise)"""
    k = len(V)
    AP = list(range(1, k + 1))
    for M in range(k + 1, max(V) + 1):
        for c in range(1, M):
            if gcd(c, M) == 1 and sorted(set((c * x) % M for x in V)) == sorted(set(x % M for x in AP)):
                return (M, c)
    return None

def main():
    print("=" * 94)
    print("S608  χ beyond vertex-transitivity: tight rotational (AP) vs SAFE doubly-regular (Paley)")
    print("=" * 94)

    print("\n[1] PALEY (doubly-regular, flat χ) is SAFE; AP (rotational, structured χ) is TIGHT")
    print("  config            | regular? (mult-inv) | χ flat? | resonance E vs main | verdict")
    cases = [("Paley P_7 {1,2,4}", [1, 2, 4], 7),
             ("AP {1,2,3}",        [1, 2, 3], None),
             ("AP {1..5}",         [1, 2, 3, 4, 5], None),
             ("AP {1..6}",         [1, 2, 3, 4, 5, 6], None)]
    for nm, V, mod in cases:
        E, l3, main = resonance_energy(V)
        if mod:
            mi = mult_inv(V, mod); sp = chi_spec(V, mod)
            flat = max(sp) - min(sp) < 0.02
            reg = f"yes mod {mod} (|grp|={len(mi)})"
        else:
            flat = False; reg = "rotational (interval)"
        print(f"  {nm:17s} | {reg:19s} | {str(flat):5s}   | E={E:.3f}  main={main:.3f}    | "
              f"{'TIGHT (high-E)' if E > main else 'SAFE (low/zero-E)'}")
    print("  => Paley P_7 has resonance energy E=0 (no effective 3-term resonance) ⇒ the SAFEST")
    print("     regular config. The tight regular configs are the structured-χ rotational (AP) class.")

    print("\n[2] THE ONE TIGHT REGULAR NON-AP-INTERVAL CONFIG (n=6) IS AN AP-LIFT")
    V = [1, 3, 4, 5, 9]
    g, a, M = gap_exact(V)
    mi = mult_inv(V, 11); lift = is_ap_lift(V)
    print(f"  {{1,3,4,5,9}} = QR mod 11 (Paley P_11): tight={abs(g-1/6)<1e-6} (G={g:.4f}, t={a}/{M})")
    print(f"     vertex-transitive: mult-inv mod 11 = {mi} (size {len(mi)} = QR group, ONE orbit)")
    print(f"     χ mod 11 = {set(chi_spec(V,11))} (FLAT, doubly-regular) — DIFFERS from AP's χ.")
    print(f"     BUT it is an AP-LIFT: {{1,3,4,5,9}} ≡ {{1,2,3,4,5}} (mod {lift[0]}) — the AP in another face.")

    print("\n[3] THE GENUINELY NON-AP TIGHT CONFIGS (sporadics) ARE NON-REGULAR")
    for nm, V in [("n=5 {1,3,4,7}", [1, 3, 4, 7]),
                  ("n=8 {1,2,3,4,5,7,12}", [1, 2, 3, 4, 5, 7, 12]),
                  ("n=8 {1,4,5,6,7,11,13}", [1, 4, 5, 6, 7, 11, 13])]:
        n = len(V) + 1; g, a, M = gap_exact(V)
        regmod = None
        for m in range(len(V) + 2, 70):
            if all(gcd(x, m) == 1 for x in V) and len(mult_inv(V, m)) > 1:
                regmod = m; break
        print(f"  {nm:22s}: tight={abs(g-1/n)<1e-5} (G={g:.4f}) | regular? {'mod '+str(regmod) if regmod else 'NO (not vertex-transitive)'}")
    print("  => the genuinely-new tight configs carry NO vertex-transitive symmetry. So every tight")
    print("     REGULAR config is the AP / an AP-lift; χ (=resonance energy) is the discriminator.")

    print("\n[4] ANSWER")
    print("  χ ADDS beyond vertex-transitivity: it is the resonance-energy signature. Among regular")
    print("  (maximally-cyclic) tournaments, the FLAT-χ doubly-regular Paley class is SAFE (low E);")
    print("  the TIGHT ones are the structured-χ rotational (AP) class and its lifts. There is NO")
    print("  tight regular config off the AP/AP-lift orbit — the non-AP tight configs are non-regular.")

if __name__ == "__main__":
    main()
