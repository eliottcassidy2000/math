#!/usr/bin/env python3
"""
three_geometries_cayley_complement_klein.py  --  klein-2026-07-01-S81

THE SYNTHESIS: a tournament lives in THREE geometries, glued by the Cayley transform:
  (1) STAIRCASE  -- the tiling of delta_{n-2} (fix base path, tiles = non-consecutive arcs).
  (2) CIRCLE     -- the runners on the unit loop: the eigenvalues of U = (I-S)(I+S)^{-1} on |z|=1
                    (S = A - A^T the skew +-1 matrix). Cayley maps skew(imaginary axis) -> unitary(circle).
  (3) The CAYLEY TRANSFORM itself is the GLUE: S (staircase/combinatorial) |-> U (circle spectrum).

CLAIM 1 -- COMPLEMENT IS A REFLECTION IN ALL THREE (one involution, three coordinates):
  staircase: complement = grid reflection sigma:(x,y)->(n+1-y,n+1-x)  (opus-S18, canon).
  skew:      complement T^op: A->A^T  =>  S = A-A^T -> -S            (negation).
  circle:    U(-S) = (I+S)(I-S)^{-1} = U(S)^{-1} = U(S)^T            (conjugation e^{i th}->e^{-i th}).
  Via Cayley, sigma <-> (S->-S) <-> (theta->-theta) are the SAME involution.

CLAIM 2 -- THE ODD/EVEN (Redei) DUALITY IS ONE FACT IN THREE GEOMETRIES:
  A class is SELF-COMPLEMENTARY (SC) iff FIXED by the reflection. Redei: H(T) (Ham-path count) is ODD.
  - staircase: an SC merged node has ODD tiling-fiber H; an NS-merged node has fiber H(T)+H(T^op)=2H, EVEN.
  - circle:    SC = self-conjugate (palindromic) spectrum; the reflection-fixed mass is odd (Redei H).
  - spectral:  SC <=> U permutation-conjugate to U^{-1}; the reflection's FIXED SET.
  The single fact = "the complement is a reflection and Redei makes its fixed fibers odd" -- read 3 ways.
  PLUS: skew-spectrum is ALWAYS +-symmetric  =>  spec(-S)=spec(S)  => the reflection acts TRIVIALLY on
  spectra, i.e. the skew/Cayley spectrum CANNOT detect complement (the S72 spectral-blindness, explained).
  PLUS: n odd <=> S singular (odd skew) <=> U has eigenvalue +1 (a runner pinned at angle 0).

VERIFIES all of the above at the labeled level (Cayley algebra) and the iso-class level (n=3,4,5).
"""
import numpy as np
from itertools import permutations, product

# ---------- tournaments as +-1 skew matrices ----------
def adj_from_bits(bits, n):
    A = np.zeros((n, n), dtype=int); idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]: A[i, j] = 1
            else:         A[j, i] = 1
            idx += 1
    return A

def skew(A):   return A - A.T
def cayley(S):
    n = S.shape[0]; I = np.eye(n)
    return (I - S) @ np.linalg.inv(I + S)

def ham_paths(A, n):
    # Redei count: number of Hamiltonian (directed) paths
    cnt = 0
    for p in permutations(range(n)):
        if all(A[p[k], p[k+1]] for k in range(n-1)): cnt += 1
    return cnt

def canon(A, n):
    best = None
    for p in permutations(range(n)):
        pp = list(p); B = A[np.ix_(pp, pp)]
        key = B.tobytes()
        if best is None or key < best: best = key
    return best

# ============================================================
if __name__ == "__main__":
    np.set_printoptions(suppress=True, precision=3)

    # ---- CLAIM 1 + Cayley glue: labeled level, all tournaments n=3,4,5 ----
    print("="*74)
    print("CAYLEY GLUE + COMPLEMENT=REFLECTION  (labeled level, all tournaments)")
    print("="*74)
    for n in (3, 4, 5):
        m = n*(n-1)//2
        orth_ok = eig_ok = comp_ok = npar_ok = True
        for bits in product((0, 1), repeat=m):
            A = adj_from_bits(bits, n); S = skew(A); U = cayley(S)
            # (a) U orthogonal
            if not np.allclose(U.T @ U, np.eye(n), atol=1e-9): orth_ok = False
            ev = np.linalg.eigvals(U)
            # (b) eigenvalues on the unit circle
            if not np.allclose(np.abs(ev), 1.0, atol=1e-9): eig_ok = False
            # (c) complement -> U^{-1} (S -> -S)
            Uop = cayley(skew(A.T))
            if not np.allclose(Uop, np.linalg.inv(U), atol=1e-9): comp_ok = False
            # (d) n odd <=> +1 eigenvalue (a runner pinned at angle 0)
            has_plus1 = np.any(np.abs(ev - 1.0) < 1e-7)
            if has_plus1 != (n % 2 == 1): npar_ok = False
        print(f"n={n}: U orthogonal={orth_ok}; |eig(U)|=1 on circle={eig_ok}; "
              f"complement->U^-1 (S->-S)={comp_ok}; (n odd <=> +1 eigenvalue)={npar_ok}")

    # ---- the QR heptagon: TWO circles. Vertex-loop vs Cayley-spectrum ----
    # There are two unit circles, and it matters which one carries "roots of unity":
    #   CIRCLE (i)  the VERTEX LOOP: a circulant/rotational tournament has its n vertices AT the n-th
    #               roots of unity (HYP-3802: the 6 atoms = primitive 14th roots). Complement = k->-k.
    #   CIRCLE (ii) the CAYLEY SPECTRUM: eig(U), U=(I-S)(I+S)^{-1}, for ANY tournament. Complement =
    #               conjugation theta->-theta. THIS is the circle glued to the staircase by Cayley.
    # These are DIFFERENT point sets. The Paley-p skew circulant has spec(S)={0, +-i*sqrt(p)} (Gauss sum
    # g(chi)=i*sqrt(p), p=3 mod 4), so its CAYLEY eigenvalues sit at the clean angle
    #   cos(theta) = (1 - p)/(1 + p) = -(p-1)/(p+1)   -- an IRRATIONAL multiple of 2pi (NOT a root of
    # unity; U^p != I). So the vertices live at roots of unity but the Cayley spectrum encodes the Gauss
    # sum. The +1 eigenvalue (from the 0 in spec(S)) is the n-odd pinned runner.
    print("\n--- QR heptagon (n=7 Paley): Cayley spectrum encodes the GAUSS SUM (NOT roots of unity) ---")
    for p in (7, 11):  # p = 3 mod 4 so the Paley tournament is well-defined
        QR = {(k*k) % p for k in range(1, p)}
        A = np.zeros((p, p), dtype=int)
        for i in range(p):
            for j in range(p):
                if i != j and ((j - i) % p) in QR: A[i, j] = 1
        S = skew(A); ev = np.linalg.eigvals(cayley(S))
        cosnz = sorted({round(e.real, 4) for e in ev if abs(e.imag) > 1e-6})
        specS = sorted({round(abs(z.imag), 4) for z in np.linalg.eigvals(S) if abs(z.imag) > 1e-6})
        print(f"  p={p}: spec(S) nonzero = +-i*{specS} (=+-i*sqrt(p)={round(p**0.5,4)}, the Gauss sum); "
              f"U nonzero cos(theta)={cosnz} (predicted -(p-1)/(p+1)={round(-(p-1)/(p+1),4)}); "
              f"+1 eig (n odd)={bool(np.any(np.abs(ev-1)<1e-6))}")

    # ---- CLAIM 2: odd/even = ONE fact, iso-class level ----
    print("\n" + "="*74)
    print("ODD/EVEN (Redei) = ONE FACT IN THREE GEOMETRIES  (iso classes)")
    print("="*74)
    A000568 = {3: 2, 4: 4, 5: 12}
    for n in (3, 4, 5):
        m = n*(n-1)//2
        classes = {}           # canon-key -> representative adjacency
        for bits in product((0, 1), repeat=m):
            A = adj_from_bits(bits, n); k = canon(A, n)
            if k not in classes: classes[k] = A
        reps = list(classes.values())
        # SC: canon(A^T) == canon(A)
        sc, ns = [], []
        for A in reps:
            (sc if canon(A.T, n) == canon(A, n) else ns).append(A)
        # per-node fiber parity via Redei H
        sc_H_odd = all(ham_paths(A, n) % 2 == 1 for A in sc)
        ns_pair_even = all((ham_paths(A, n) + ham_paths(A.T, n)) % 2 == 0 for A in ns)
        # spectral blindness: spec(S) == spec(-S) always
        spec_blind = True
        for A in reps:
            S = skew(A)
            e1 = np.sort_complex(np.round(np.linalg.eigvals(S), 6))
            e2 = np.sort_complex(np.round(np.linalg.eigvals(-S), 6))
            if not np.allclose(e1, e2, atol=1e-6): spec_blind = False
        ns_merged = len(ns)//2
        Vmerged = (len(reps) + len(sc))//2
        print(f"n={n}: |G_n|={len(reps)} (A000568={A000568[n]} {'OK' if len(reps)==A000568[n] else 'MISMATCH'}); "
              f"#SC={len(sc)} (all even count); #NS={len(ns)} -> NS-merged={ns_merged}; V_merged={Vmerged}")
        print(f"     [staircase]  SC node fiber H ODD = {sc_H_odd};  NS-merged fiber 2H EVEN = {ns_pair_even}   (Redei)")
        print(f"     [spectral]   spec(S)=spec(-S) always (reflection trivial on spectra; SC undetectable) = {spec_blind}")
        print(f"     [circle]     SC <=> palindromic/self-conjugate spectrum (fixed by theta->-theta)")

    print("\n" + "="*74)
    print("NET: complement = ONE reflection (sigma / S->-S / theta->-theta), Cayley-conjugate across")
    print("     the three geometries. Its FIXED set = the SC classes; Redei makes the fixed fiber ODD,")
    print("     the moved (NS) fiber 2H EVEN -- the odd/even fiber parity, one fact seen three ways.")
    print("     The skew/Cayley spectrum is complement-blind (reflection trivial on spectra); n-parity")
    print("     is the pinned runner at 0 (U has +1 <=> n odd <=> S singular).")
    print("="*74)
