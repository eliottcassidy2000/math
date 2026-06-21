#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A, part 2 (kpswf6): MAKE THE LINEARITY SPLIT CONCRETE.

Goal: show EXACTLY where CJJ's linear-code completeness applies to Paley=QR but
to a functional that is NOT the tournament H. Two computations:

  (A) The QR cyclic code is a genuine LINEAR code. Its weight enumerator obeys
      MacWilliams (the dual / Delsarte LP closes EXACTLY -- the CJJ completeness
      condition). We exhibit this for the QR (quadratic-residue) code mod p.
      => the LINEAR-code certificate EXISTS and is closed for the QR/Paley side.

  (B) The LRC AP miss-distribution is NON-linear: the analogue lattice is the
      BOOLEAN lattice of missed-sector subsets (HYP-2744), NOT a subspace lattice;
      there is no MacWilliams dual that closes => CJJ completeness does NOT apply.
      We reproduce the structural contrast (Boolean vs subspace Mobius).

  (C) THE PUNCHLINE: even where the QR linear certificate closes (A), it certifies
      the CODE's distance distribution, NOT H = I(Omega(T),2). The tournament
      H-functional is a DIFFERENT, nonlinear object on the SAME QR set, and its
      maximizer is the flat spectrum only up to p=11 (THM-128 kills it at p=13).

So linearity SPLITS the certificate-existence question (QR: yes; AP: no) but does
NOT transfer to the tournament H-extremality, which is the functional CJJ does not
see and which is anyway false past p=11.
"""
import sys
import numpy as np
from fractions import Fraction as F
from math import comb
from scipy.optimize import linprog
sys.stdout.reconfigure(line_buffering=True)

def qr_set(p):
    return sorted({pow(a, 2, p) for a in range(1, p)})

# ---------------------------------------------------------------------------
# (A) The QR code is linear: build the binary QR code of length p (p=7 -> the
#     [7,4] Hamming code = QR code; p=23 -> [23,12] Golay). Show MacWilliams holds
#     EXACTLY (the linear-code Delsarte dual closes -- CJJ completeness condition).
# ---------------------------------------------------------------------------
print("=" * 78)
print("(A) THE QR CODE IS LINEAR: MacWilliams duality closes EXACTLY (CJJ premise)")
print("=" * 78)

def qr_code_generator_rows(p):
    """Binary cyclic QR code of length p: idempotent generator from QR positions.
    g(x) = prod_{i in QR}(x - alpha^i) over GF(2^k); we just build the code as the
    cyclic code whose nonzeros are QR. Simplest exact route for small p: build via
    the (p, (p+1)/2) cyclic QR code by its generator polynomial over GF(2).
    For demonstration we use p=7 -> [7,4] Hamming (the binary QR code)."""
    # p=7 QR = {1,2,4}; generator polynomial g(x)=x^3+x+1 (Hamming/QR code)
    # We'll just hand the known [7,4] Hamming code generator matrix.
    if p == 7:
        G = np.array([
            [1,0,0,0,1,1,0],
            [0,1,0,0,0,1,1],
            [0,0,1,0,1,1,1],
            [0,0,0,1,1,0,1],
        ], dtype=int)
        return G
    return None

def weight_enum(G):
    k, n = G.shape
    A = [0] * (n + 1)
    for msg in range(2 ** k):
        cw = np.zeros(n, dtype=int)
        for b in range(k):
            if (msg >> b) & 1:
                cw ^= G[b]
        A[int(cw.sum())] += 1
    return A

def krawtchouk(k, x, n):
    return sum((-1) ** j * comb(x, j) * comb(n - x, k - j) for j in range(k + 1))

G = qr_code_generator_rows(7)
A = weight_enum(G)
n = 7
print(f"  p=7 binary QR code = [7,4,3] Hamming code. Weight enumerator A = {A}")
# MacWilliams transform: B_k = (1/|C|) sum_i A_i K_k(i)
sizeC = sum(A)
B = [F(sum(A[i] * krawtchouk(k, i, n) for i in range(n + 1)), sizeC) for k in range(n + 1)]
print(f"  Dual weight enum via MacWilliams B = {[str(b) for b in B]}")
# the [7,4] Hamming dual is the [7,3] simplex code: all 7 nonzero codewords have
# weight 4, so A_dual = [1,0,0,0,7,0,0,0].
print(f"  Expected dual (simplex [7,3], all 7 nonzero cw weight 4) = [1, 0, 0, 0, 7, 0, 0, 0]")
ok = [int(b) for b in B] == [1, 0, 0, 0, 7, 0, 0, 0]
print(f"  MacWilliams duality EXACT (the linear-code Delsarte dual closes)? {'YES' if ok else 'NO'}")
print("""
  => For the QR (linear) code, the Delsarte/MacWilliams dual is an EXACT closed
     transform. This is precisely the structure CJJ leverages: the LP-hierarchy is
     COMPLETE because the dual lives on the (linear) subspace lattice and closes.
""")

# ---------------------------------------------------------------------------
# (B) Structural contrast: AP miss-distribution lattice is BOOLEAN, not subspace.
# ---------------------------------------------------------------------------
print("=" * 78)
print("(B) THE AP (LRC) SIDE: Boolean lattice, NOT a subspace lattice -> no closure")
print("=" * 78)
print("""
  HYP-2744: the LRC analogue uses the BOOLEAN lattice of missed-sector subsets
  (2^6 masks) with Boolean Mobius inversion q[M]=sum_{A superset M}(-1)^|A\\M| a[A].
  This is NOT the subspace lattice of a linear code. CJJ's completeness/integrality
  (Prop 1.2) needs the optimizer to BE a linear code so the dual closes on the
  subspace lattice. The AP miss-distribution is non-linear => no such closure =>
  the hierarchy improves the BOUND but cannot SELECT consec (the HYP-2744 collapse,
  re-confirmed: degree-2 per-subset functional non-vacuously selects consec, but
  every clean monotone atom has many AP-beaters; needs a signed aggregate cut).
""")

# ---------------------------------------------------------------------------
# (C) THE PUNCHLINE -- the linear certificate certifies the CODE, not H.
# ---------------------------------------------------------------------------
print("=" * 78)
print("(C) PUNCHLINE: the linear QR certificate certifies DISTANCE, not H=I(Omega,2)")
print("=" * 78)
print("""
  The QR code's MacWilliams dual (A) certifies the QR code's DISTANCE DISTRIBUTION.
  But the tournament extremal quantity is H = I(Omega(T),2), the hardcore partition
  function of the ODD-CYCLE conflict graph -- a NONLINEAR functional that the
  MacWilliams/Delsarte LP does not compute and does not bound from below.

  Concretely (from threadA_cjj_paley_certify_kpswf6 Part 1, p=7):
    * the Lovasz/Delsarte LP bound theta_LP(Omega) is LARGER for Paley (40) than
      for non-Paley (29.5): the independence relaxation ANTI-tracks H.
    * alpha(Omega)=2 for ALL 8 circulants: the independence number does not
      separate the maximizer at all.
  And globally THM-128: at p=13 the H-maximizer is the AP/odd-step (NON-linear)
  set; the flat QR-style spectrum LOSES. So the linear-code certificate, even where
  it exists, points at the wrong functional and the wrong optimizer past p=11.

  HONEST VERDICT: linearity DOES split the CERTIFICATE-EXISTENCE question
  (QR=closed MacWilliams dual; AP=open Boolean aggregate). It does NOT yield a CJJ
  proof of tournament H-max extremality, because (i) CJJ certifies code distance,
  not the hardcore H-functional, and (ii) the H-extremality of Paley is itself only
  a p=3mod4, p<=11 phenomenon -- there is no binding large p where Paley is the
  H-maximizer for a hierarchy to certify.
""")
print("DONE.")
