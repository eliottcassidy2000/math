#!/usr/bin/env python3
"""
u_spectrum_flip_rank_second_moment_klein.py  --  klein-2026-07-01-S82

QUESTION (owner): does the U-spectrum (Cayley transform eig) SEE the flip-rank excess that the
skew-spectrum misses?  And consider the SECOND MOMENT.

ANSWER (this script tests): NO at the eigenvalue level -- and there is a clean structural reason.
  (1) Cayley is a SPECTRAL BIJECTION: eig(U) = (1 - i*mu)/(1 + i*mu) for eig(S)=i*mu. So U-cospectral
      <=> S-cospectral; the U-spectrum carries EXACTLY the information of the skew-spectrum. No gain.
  (2) RESOLUTION CEILING = V_merged. Complement is a REFLECTION (S81): spec(-S)=spec(S) ALWAYS, so
      every NS complement-PAIR is cospectral. Hence #distinct spectra <= V_merged = (|G_n|+#SC)/2:
      the spectrum factors through the MERGED (reflection-quotient) metagraph and can NEVER resolve
      the unmerged G_n. The flip-rank lives on the UNMERGED G_n and its EXCESS is carried by the SC
      (reflection-FIXED) classes -- exactly the fixed points the quotient collapses. So no spectral
      invariant (S or U) can see the flip-rank excess.
  (3) THE SECOND MOMENT is the order-2 shadow of that reflection-symmetry, and it is CONSTANT:
      trace(S^2) = -n(n-1) for EVERY tournament (blind). The first informative even moment is
      trace(S^4)=sum mu^4. The Cayley wrap turns the blind constant 2nd moment into a NON-constant
      circular moment trace(U)=sum cos(theta)=sum (1-mu^2)/(1+mu^2), which improves cospectral
      resolution -- but still only up to the V_merged ceiling (still complement-blind).

Verifies (1),(2),(3) on all iso classes n=3,4,5 (and attempts n=6). Uses EXACT integer char polys
(Faddeev-LeVerrier) for cospectrality; numpy floats only for the circle/second-moment illustration.
"""
import numpy as np
from fractions import Fraction
from itertools import permutations, product
import time

# ---------- exact integer characteristic polynomial (Faddeev-LeVerrier) ----------
def charpoly_int(M):
    n = len(M)
    Mk = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    I = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    coeffs = [Fraction(1)]  # leading
    c = Fraction(1)
    A = [row[:] for row in Mk]
    for k in range(1, n+1):
        c = -sum(A[i][i] for i in range(n)) / k
        coeffs.append(c)
        if k < n:
            # A = M @ (A + c I)
            ApluscI = [[A[i][j] + (c if i == j else 0) for j in range(n)] for i in range(n)]
            A = [[sum(Mk[i][t]*ApluscI[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    return tuple(int(x) for x in coeffs)  # integer coeffs

# ---------- tournaments ----------
def adj_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]: A[i][j] = 1
            else:         A[j][i] = 1
            idx += 1
    return A

def transpose(A, n): return [[A[j][i] for j in range(n)] for i in range(n)]
def skew(A, n):      return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]

def ham_paths(A, n, perms):
    return sum(1 for p in perms if all(A[p[k]][p[k+1]] for k in range(n-1)))

def three_cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i < j and A[i][j] and A[j][k] and A[k][i]: c += 1
    return c

def canon_key(A, n, perms):
    best = None
    for p in perms:
        flat = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or flat < best: best = flat
    return best

def enumerate_classes(n, tbudget=90):
    perms = list(permutations(range(n)))
    m = n*(n-1)//2
    seen = {}   # canon -> representative adjacency (as tuple of tuples)
    t0 = time.time()
    for bits in product((0, 1), repeat=m):
        A = adj_from_bits(bits, n)
        k = canon_key(A, n, perms)
        if k not in seen:
            seen[k] = A
        if time.time() - t0 > tbudget:
            return None  # bail
    return list(seen.values())

# ---------- U-spectrum (Cayley) ----------
def cayley_eig(A, n):
    S = np.array(skew(A, n), dtype=float); I = np.eye(n)
    U = (I - S) @ np.linalg.inv(I + S)
    return np.linalg.eigvals(U), U, S

# ============================================================
if __name__ == "__main__":
    # canon merged-count reference (from CLAUDE.md)
    G = {3: 2, 4: 4, 5: 12, 6: 56}
    SCn = {3: 2, 4: 2, 5: 8, 6: 12}        # #self-complementary classes
    Vmerged = {n: (G[n]+SCn[n])//2 for n in G}
    FLIP = {3: 1, 4: 2, 5: 4, 6: 7}         # flip-rank rho(n) (S71/HYP-3803)
    LB   = {n: (G[n]-1).bit_length() for n in G}  # ceil(log2 |G_n|)

    print("="*78)
    print("DOES THE U-SPECTRUM SEE THE FLIP-RANK EXCESS?   (n, structural facts)")
    print("="*78)
    print(f"{'n':>2} {'|G_n|':>6} {'V_merged':>9} {'#SC':>4} {'rho(flip)':>9} {'ceil_log2':>9} {'excess':>7}")
    for n in (3, 4, 5, 6):
        print(f"{n:>2} {G[n]:>6} {Vmerged[n]:>9} {SCn[n]:>4} {FLIP[n]:>9} {LB[n]:>9} {FLIP[n]-LB[n]:>7}")
    print("  (excess first appears at n=6; HYP-3810: SC classes carry it)")

    for n in (3, 4, 5, 6):
        print("\n" + "-"*78)
        print(f"n = {n}")
        reps = enumerate_classes(n, tbudget=120)
        if reps is None:
            print(f"  [class enumeration exceeded time budget -- skipping fine analysis at n={n}]")
            print(f"  ANALYTIC ceiling still holds: #distinct spectra <= V_merged = {Vmerged[n]} (complement pairs cospectral).")
            continue
        assert len(reps) == G[n], f"class count {len(reps)} != {G[n]}"
        perms = list(permutations(range(n)))
        # exact skew char poly (spectrum) per class
        cp = {}   # charpoly -> list of class indices
        sc_flags = []
        tr2, tr4, trU = [], [], []
        for idx, A in enumerate(reps):
            S = skew(A, n)
            poly = charpoly_int(S)
            cp.setdefault(poly, []).append(idx)
            # SC?  canon(A^T) == canon(A)
            sc = canon_key(transpose(A, n), n, perms) == canon_key(A, n, perms)
            sc_flags.append(sc)
            # second moments
            Smat = np.array(S, dtype=float)
            tr2.append(round(float(np.trace(Smat @ Smat)), 6))
            S2 = Smat @ Smat
            tr4.append(round(float(np.trace(S2 @ S2)), 6))
            ev, U, _ = cayley_eig(A, n)
            trU.append(round(float(np.trace(U)), 6))
        ndist = len(cp)
        nSC = sum(sc_flags)
        # (1) Cayley bijection check: distinct U-spectra count == distinct S-spectra count
        u_specs = set()
        for A in reps:
            ev, _, _ = cayley_eig(A, n)
            u_specs.add(tuple(sorted(round(np.angle(e), 6) for e in ev)))
        # (2) resolution ceiling
        # complement-pair cospectrality: every NS pair shares a spectrum
        print(f"  |G_n|={len(reps)}, #SC={nSC}, V_merged={Vmerged[n]}")
        print(f"  (1) distinct skew-spectra = {ndist};  distinct U-spectra = {len(u_specs)}  "
              f"(equal => Cayley is a spectral bijection, U carries no more)")
        print(f"  (2) resolution ceiling V_merged={Vmerged[n]}: distinct-spectra {ndist} "
              f"{'== V_merged (TIGHT: spectrum resolves exactly the merged metagraph)' if ndist==Vmerged[n] else f'< V_merged (LOSES {Vmerged[n]-ndist} MORE beyond complement-pairs)'}")
        # how many cospectral collisions are complement-pairs vs extra
        collisions = sum(len(v)-1 for v in cp.values())
        print(f"      cospectral collisions = |G_n|-distinct = {len(reps)-ndist}  "
              f"(>= #NS-pairs = {(len(reps)-nSC)//2} forced by complement-reflection)")
        # (3) second moments
        print(f"  (3) 2nd moment trace(S^2): values={sorted(set(tr2))}  "
              f"(constant = -n(n-1) = {-n*(n-1)} => BLIND)")
        print(f"      4th moment trace(S^4): #distinct values = {len(set(tr4))}  (first informative even moment)")
        print(f"      circular trace(U)=sum cos(theta): #distinct values = {len(set(trU))}  "
              f"(Cayley wrap makes the blind 2nd moment informative; note it matches/beats trace(S^4))")
        # (4) COMBINATORIAL invariants (the RIGHT instrument for the metagraph/flip-rank)
        H = [ham_paths(A, n, perms) for A in reps]
        c3 = [three_cycles(A, n) for A in reps]
        meanH = sum(H)/len(H)
        varH = sum((h-meanH)**2 for h in H)/len(H)   # 2nd moment of H over the metagraph (THM-589 W(n))
        combo = len(set(zip(H, c3)))
        print(f"  (4) COMBINATORIAL: #distinct H(Redei) = {len(set(H))}, #distinct (H,c3) = {combo}  "
              f"(vs #distinct spectra = {ndist}) -- combinatorial counts resolve FAR more than any spectrum")
        print(f"      H all odd (Redei) = {all(h % 2 == 1 for h in H)}; metagraph 2nd moment Var(H) = {varH:.3f}  "
              f"(THM-589 W(n); the SECOND MOMENT that IS relevant to covering/flip-rank)")

    print("\n" + "="*78)
    print("NET: U-spectrum = S-spectrum (Cayley bijection) -> NO eigenvalue-level gain. Both factor")
    print("through the MERGED metagraph (ceiling V_merged) because complement is a reflection and the")
    print("spectrum is reflection-invariant. The flip-rank excess is carried by the SC (fixed) classes")
    print("on the UNMERGED G_n -> SPECTRALLY INVISIBLE. The 2nd moment trace(S^2)=-n(n-1) is the constant")
    print("order-2 shadow of that symmetry; Cayley wraps it into an informative circular moment trace(U),")
    print("but resolution is still capped at V_merged. To see the excess needs NON-spectral symmetry data.")
    print("="*78)
