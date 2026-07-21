#!/usr/bin/env python3
"""
positivity_direction_markovkrein_deathstar_S77.py

The DECISIVE part of the "does GMC(2) => LRC(14)" positivity-transfer test, using
Agent-2's exact objects. Two facts break the naive transfer, one survives.

Agent-2 verbatim structure:
 * The object that IMPLIES LRC(14) is the moment FLOOR B_d (THM-661): a
   1-D Markov-Krein moment problem on a SCALAR W = sum_i (g_i - 1/7)_+ in [0,6/7],
   B_d = max{ sum c_i E[W^i] : sum c_i w^i <= 1_{w>0} on [0,6/7] }.  B_2 = Paley-Zygmund.
 * The object SHARED with GMC is disc_v (THM-731/732) = Parseval energy of the
   arc-endpoint jump-sum S_m = sum_j sign_j e^{-2pi i m x_j}  (= GMC reconstruction sum).
   LRC needs disc_v <= UPPER bound; GMC positivity gives LOWER bounds.

Demonstrates:
 (E) DIRECTION MISMATCH + the energy-conservation rescue. Parseval:
     sum_{k!=0}|ghat(k)|^2 = |G'|/r - (|G'|/r)^2 is FIXED. disc_v is the v-dilated
     sub-sum. So  "disc_v <= B"  <=>  "complement energy >= fixed - B": an UPPER
     bound on disc_v IS a LOWER bound on complementary energy -- which a positive
     (sum-of-squares / Fejer) certificate CAN provide. This is the one surviving
     positivity idea, and it does NOT need GMC(2) as a statement.
 (F) The moment FLOOR is 1-D Markov-Krein: its Hankel is PSD for FREE (moments of
     a real RV); the floor depends only on the moment VALUES E[W^i]. GMC(2) fixes
     no such values -> cannot discharge the LRC crux. Show B_2 = E[W]^2/E[W^2].
"""
from fractions import Fraction as Fr

try:
    import numpy as np
    HAVE_NP = True
except Exception:
    HAVE_NP = False


def sep(t): print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


def min_eig(mat):
    if HAVE_NP:
        return float(np.min(np.linalg.eigvalsh(np.array(mat, float))))
    n = len(mat); A = [[float(mat[i][j]) for j in range(n)] for i in range(n)]
    for k in range(n):
        p = A[k][k] - sum(A[k][t] ** 2 for t in range(k))
        if p <= 0: return -1.0
        p = p ** .5
        for i in range(k + 1, n):
            A[i][k] = (A[i][k] - sum(A[i][t] * A[k][t] for t in range(k))) / p
        A[k][k] = p
    return 1.0


# ----------------------------------------------------------------------
# (E) DIRECTION MISMATCH and the energy-conservation rescue.
#     Concrete good set G' subset Z/r. ghat(k) = (1/r) sum_{x in G'} e(-2pi i k x/r).
#     Total nonzero energy is FIXED by Parseval; disc_v is a sub-sum of it.
# ----------------------------------------------------------------------
sep("(E) DIRECTION: LRC needs disc_v <= UPPER bound; positivity gives LOWER bounds")
import cmath, math


def fourier_energy(Gp, r):
    """returns (ghat0_sq, total_nonzero_energy, per-k |ghat(k)|^2 list)."""
    e = []
    for k in range(r):
        s = sum(cmath.exp(-2j * math.pi * k * x / r) for x in Gp) / r
        e.append(abs(s) ** 2)
    total_nz = sum(e[1:])
    return e[0], total_nz, e


def disc_v(e, r, v):
    # disc_v = sum_{m != 0} |ghat(m v mod r)|^2  (m=1..r-1), indices mod r
    return sum(e[(m * v) % r] for m in range(1, r) if (m * v) % r != 0)


r = 14                     # LRC14 modulus flavor
Gp = [1, 2, 3, 5, 8]       # a toy good set of residues mod 14
g0, tot, e = fourier_energy(Gp, r)
print(f"good set G'={Gp} mod {r}:  |G'|/r={Fr(len(Gp),r)},  Parseval nonzero energy={tot:.5f}")
print(f"   check: |G'|/r - (|G'|/r)^2 = {len(Gp)/r - (len(Gp)/r)**2:.5f}  (= total nonzero energy)")
for v in (3, 5, 7):
    dv = disc_v(e, r, v)
    print(f"   v={v}: disc_v={dv:.5f},  complement energy = total - disc_v = {tot-dv:.5f}")
print("""
  READING:
   * LRC covering route B needs  disc_v <= (upper bound)   [|eps_v|^2 <= (6/49) disc_v].
   * Positive-definiteness / Hankel-PD (GMC's escape) yields LOWER bounds, never upper.
     -> the shared object disc_v is bounded the WRONG WAY by GMC positivity: MISMATCH.
   * RESCUE (energy conservation): total nonzero energy is FIXED (Parseval), so
       disc_v <= B   <=>   complement energy >= total - B.
     An UPPER bound on disc_v is a LOWER bound on the complementary (non-v-resonant)
     energy -- a LOWER bound, which a sum-of-squares / Fejer positive certificate CAN
     supply. THIS is the one surviving positivity idea; it is a METHOD, provable
     directly, and needs NO assumption of GMC(2).""")

# ----------------------------------------------------------------------
# (F) The LRC(14) moment FLOOR is a 1-D Markov-Krein problem: Hankel PSD is FREE,
#     the floor depends only on the moment VALUES. GMC(2) fixes no values.
# ----------------------------------------------------------------------
sep("(F) The moment FLOOR is 1-D Markov-Krein: positivity FREE, crux = the VALUES")
# toy scalar W supported in [0, 6/7]: take W in {0, 1/2, 6/7} with probs.
support = [Fr(0), Fr(1, 2), Fr(6, 7)]
probs = [Fr(1, 2), Fr(1, 3), Fr(1, 6)]      # a genuine prob measure => moments PSD free
def mom(i): return sum(p * (w ** i) for p, w in zip(probs, support))
M = [mom(i) for i in range(5)]              # E[W^0..4]
print(f"toy W on [0,6/7], moments E[W^0..4] = {[str(x) for x in M]}")
Hk = [[float(mom(i + j)) for j in range(3)] for i in range(3)]   # 3x3 Hankel
print(f"3x3 Hankel of E[W^i]: min eig = {min_eig(Hk):.4e}  -> PSD (FREE: moments of a real RV)")
mu = sum(p for p, w in zip(probs, support) if w > 0)            # P(W>0)
B2 = float(M[1] ** 2 / M[2])                                    # Paley-Zygmund floor
print(f"P(W>0) = {mu} ;  Paley-Zygmund floor B_2 = E[W]^2/E[W^2] = {B2:.4f}  <= P(W>0)={float(mu):.4f}: {B2<=float(mu)}")
print("""
  READING:
   * The Hankel of E[W^i] is PSD automatically (W is a real random variable) -- the
     LRC floor does NOT need any positivity theorem; positivity is free here.
   * The floor B_d is an LP over the moment cone; its value is fixed by the moment
     VALUES E[W^i]. The open LRC crux is computing/bounding those values (the
     decorrelation tail) -- pure harmonic analysis of the SPEED SET (arithmetic).
   * GMC(2) is a statement about a DIFFERENT measure (Gaussian) and says NOTHING
     about the values E[W^i] of the LRC scalar W. It cannot discharge this crux.""")

sep("VERDICT")
print("""GMC(2) => LRC(14) does NOT go through, for three now-sharp reasons:

 1. TWO OBJECTS, neither transfers.
    - The object that IMPLIES LRC(14) (moment floor B_d) is 1-D Markov-Krein:
      positivity FREE, crux = moment VALUES (decorrelation) -> GMC(2)-blind.
    - The object SHARED with GMC (disc_v = jump-sum energy) needs an UPPER bound;
      GMC positivity gives LOWER bounds. DIRECTION MISMATCH.
 2. LAYER MISMATCH. LRC covering is TORAL (constant-term) = DvdK = the PROVED half
    of GMC(2)=L_s o CT_u; GMC(2)'s OPEN half is the perpendicular RADIAL/Gaussian
    layer. Assuming the open half is assuming the LRC-irrelevant half.
 3. NOT A POSITIVITY. GMC(2) is a nullcone STRUCTURE (variety) conjecture, not a
    Hankel positivity; the positive reformulation is methodological, so even the
    first lemma of the transfer (GMC(2) <=> H_GMC >= 0) is not in hand.

WHAT SURVIVES (the honest, creative residue):
 * The ONE machine-verified link (disc_v = Parseval of the GMC jump-sum) + energy
   conservation gives a concrete NEW positivity target for LRC alone: bound disc_v
   from ABOVE by lower-bounding the complementary (non-v-resonant) energy with a
   sum-of-squares / Fejer certificate. Method transfer, not statement implication.
 * LRC-AP's extremal config IS a genuine power-sum moment-nullcone (roots of unity,
   depth 14) -- so LRC sits ADJACENT to the THM-1775 template (shares the jump-sum
   object + the positivity manoeuvre) but is NOT an instance (its proof uses a
   Markov-Krein floor + a certificate discrepancy, not a nullcone vanishing).
 * The strongest TRUE conditional: assuming GMC(2) closes the repo's own GMC(2)
   program (radial bridge); it does NOT give JC (needs GMC for ALL n) and does NOT
   give LRC(14).""")
