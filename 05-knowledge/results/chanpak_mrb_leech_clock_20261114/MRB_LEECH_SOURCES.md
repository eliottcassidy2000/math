# MRB, the Leech lattice, and which transport laws really lift

**Scope.** The primary papers and earlier repository proofs are inputs,
not new MRB or Leech theorems. The new statements below are exact
rank-two consequences for corridor transfers. They neither prove
arithmetic properties of the MRB constant nor identify a tournament
lattice with the Leech lattice.

## 1. What the supplied MRB papers actually provide

The archived Session 31 input includes and was reread in the original
supplied texts:

- Marvin Ray Burns, *The Many Worlds View of the MRB Constant:
  Multiplicative Growth, Logarithmic Gain, and Alternating Aggregation*,
  §§2–3,9–12 (15-page supplied PDF; [author's description](https://www.linkedin.com/posts/marvin-ray-burns-27a39329_the-many-worlds-view-of-the-mrb-constant-activity-7504251711335763969-HHGo)).
  The gain law is \(k^{1/k}=e^{\log k/k}\); the signed adjacent
  pairing is distinct from generating those gains. Subtract
  \(\log k/k\) before summing to remove the leading slow layer;
  the nonlinear remainder is absolutely summable. The signed pairs
  are absolutely summable, whereas the individual nonzero-sign
  subseries diverge in magnitude. Pairing without reordering recovers
  the fixed value; rearranging unpaired summands is not a construction
  for that same fixed value.
- Burns, *What Lebsir's Prime Progressions Gave the MRB Problem*
  ([DOI](https://doi.org/10.29322/IJSRP.16.10.2026.p17806)),
  pp. 20–23 (supplied PDF). The square block
  \([n^2,(n+1)^2)\) has exact even center \(n(n+1)\).
  Centered second differences supply a curvature correction, while
  odd moments cancel by reflection. The first displayed curvature
  identity needs the [Session 31 sign repair](https://github.com/eliottcassidy2000/math/blob/095ef14e1e88351526125b56f2679ed5d7d40cbe/05-knowledge/results/session31_mirror_mrb_20261022/square-blocks/AUDIT.md):
  the correction has \((-1)^d\), not \((-1)^{n+d}\).
  That repair also gives a positive integral kernel and rigorously
  bounds the finite-block remainder. No prime-progression conjecture
  is needed for the symmetry.
- Burns, [*The Analytic Structure of the MRB Constant*](https://www.academia.edu/171694177/The_Analytic_Structure_of_the_MRB_Constant_Abel_Plana_Summation_Reciprocal_Residues_and_a_Universal_Greedy_Expansion),
  and [*A Convergent Greedy Residue with Geometry*](https://www.academia.edu/175280204/A_Convergent_Greedy_Residue_with_Geometry_Ramanujan_and_Michele_Nardellis_Version),
  are primary background for the contour/residue lane. See the
  [Session 31 contour audit](https://github.com/eliottcassidy2000/math/blob/095ef14e1e88351526125b56f2679ed5d7d40cbe/05-knowledge/results/session31_mirror_mrb_20261022/mrb/MRB_CONTOUR_AUDIT.md)
  and the [Session 144 strengthening](https://github.com/eliottcassidy2000/math/blob/095ef14e1e88351526125b56f2679ed5d7d40cbe/05-knowledge/results/session144_jacobian_eth_mrb_20261103/mrb/RESEARCH.md).
  The branch jump gives the actual residue \(A(8)\). The exact
  CM identities concern that residue, not the difference
  \(A(8)-C_{\rm MRB}\); in fact \(0<A(8)-C_{\rm MRB}<q^4\),
  \(q=e^{-\pi\sqrt2}\). Explicit entire perturbations retain
  the branch jumps, any finite prefix, and all uniform
  right-half-plane inverse-power asymptotics while changing the
  alternating functional. They do not retain the complete defining
  nth-root function. The operative lesson is to retain the data on
  which the desired **readout** depends.

A prior [MRB/square-probe synthesis](https://github.com/eliottcassidy2000/math/blob/095ef14e1e88351526125b56f2679ed5d7d40cbe/05-knowledge/results/mrb_square_probe_spectrum_20261103/README.md)
implemented this lesson physically: a period-four MRB filter detects
the difference of two named singleton-birth counts on a tournament
source. Reversing a triangle can agree on *all* one-birth probes
but differ on a two-birth word. Hence pairing, resonance, and phase
filters can expose missing coordinates; a scalar obtained from them
does not create the next legal source. The all-modulus terrace
clock in the [companion proof](ODD_MONODROMY.md) is a new, exact
version of that warning.

## 2. The CM value occurs in a genuine conjugacy, but outside the physical cone

The prior corridor analysis bounded contraction in log odds by
\(3-2\sqrt2\). Independently, the MRB contour audit's exact AGM
calculation gives the elliptic modulus \(\sqrt2-1\) at nome
\(e^{-\pi\sqrt2}\), hence
\(\lambda(i\sqrt2)=(\sqrt2-1)^2=3-2\sqrt2\).
This is an exact shared algebraic value, and it has a direct
operator interpretation.

At a one-zero seam the odds transfer with next width-power \(y\) is
\[
 t\longmapsto y(2+2t)/(1+2t).
\]
Set \(y=1/2\) formally. The map \(F(t)=(1+t)/(1+2t)\)
fixes \(t=\pm1/\sqrt2\), and conjugating by
\(z=(t-1/\sqrt2)/(t+1/\sqrt2)\) gives
\[
 z(F(t))=-(3-2\sqrt2)\,z(t).
\]
Indeed the two fixed lines imply multiplication in \(z\), and
the multiplier is \(F'(1/\sqrt2)\). The analytic continuation
to \(y=1/2\) is **not** a legal positive gap width; no physical
tournament transfer or MRB sum is transported by this conjugacy.

On the physical cone \(t\ge2\), the sharp uniform log-odds
contraction is instead **\(2/15\)**. Every input and outgoing
physical state satisfies \(t\ge2\), since \(f_q\ge e_q\) and
\(g_q\ge f_q\) and \(y\ge2\). For \(q=1\) the absolute
derivative is
\[
 2t/((2+2t)(1+2t)).
\]
It decreases for \(t>1/\sqrt2\), so its maximum occurs at \(t=2\),
a physical input (seed width one). For \(q\ge2\) we have
\(e_q/f_q\ge3/4\) and \(g_q/f_q\le4/3\); both hold for
\(q=2,3,4\) and then propagate by the positive recurrence.
For any \(t>0\), the general log-odds derivative maximum is
\[
 \frac{|\sqrt{e_qg_q}-f_q|}{\sqrt{e_qg_q}+f_q}
 \le 7-4\sqrt3 < 2/15.
\]
Its elementary proof minimizes the denominator divided by \(t\).
The two reset lengths have zero derivative. This sharper physical
bound can improve projective-state interval estimates; it does not
imply that approximated integer counts have any missing residues.

## 3. Two honest centered corrections

**An exact integer analogue of dominant-term subtraction.**
Repeat one-zero seams, with every physical positive gap width
equal to one. Starting at \((1,2)^t\) the state matrix is
\(\left(\begin{smallmatrix}1&2\\4&4\end{smallmatrix}\right)\).
The whole counts \(H_k\), with \(k\) seams, satisfy
\[
 H_0=3,\quad H_1=17,\quad H_{k+2}=5H_{k+1}+4H_k.
\]
For \(\mu_\pm=(5\pm\sqrt{41})/2\) and
\(\beta=(3\sqrt{41}-19)/(2\sqrt{41})\),
\[
 H_k=(3-\beta)\mu_+^k+\beta\mu_-^k.
\]
Here \(0<\beta<1/50\) and \(-1<\mu_-<0\).
Thus exact nearest-integer rounding of the dominant term recovers
the ordinary count at every depth. The residual alternates and
has a uniform strict rounding certificate, unlike a fit inferred
from main-term agreement. This is a genuine positive-width count
family, but the sequence is exponentially sparse.

**A square-center-style block check that exposes gaps, not fills them.**
For a symmetric fixed two-wing corridor with exponent sum \(s=a+b\),
write \(a=s/2+h\), \(b=s/2-h\). The exact center difference is
\[
 H(h)-H(0)=2 f_q 2^{s/2}(\cosh(h\log2)-1).
\]
It is even, positive for nonzero \(h\), and has no odd Taylor layers.
For integer widths, the displacement \(h\) is integral when \(s\)
is even and half-integral when \(s\) is odd. At even \(s\ge4\),
the first positive step from the balanced center increases
\(H\) by \(f_q2^{s/2-1}\). So this genuine centered curvature makes
a **growing** integer gap at fixed center. An MRB-style small
normalized analytic error can coexist with a large unnormalized
integer target deficit. A different frame, not a more refined
real Taylor expansion, is needed to fill that deficit.

## 4. What the six integers mean on the Leech side

The requested \(4,2,5,3,7,13\) labels were recovered in the
[historical valuation chain](https://github.com/eliottcassidy2000/math/blob/095ef14e1e88351526125b56f2679ed5d7d40cbe/05-knowledge/results/session32_square_monster_20261023/VALUATION_CYCLE_AND_QUOTIENT.md):
\(13\to2\to5\to3\to7\to4\); a five-state cycle requires explicitly
identifying 4 with 13. It is not six iterations of a globally
defined number-theoretic map. The active physical corridor uses
4 and deliberately does not identify it with singular 13.

There is a precise Leech-side role for all six numbers.
For a pure \(A\)-type Niemeier root system \(A_{h-1}^{\,k}\),
\(k(h-1)=24\); hence the Coxeter number \(h\) runs through
\(2,3,4,5,7,9,13,25\). All six requested labels are present;
the omitted labels \(9,25\) come from divisors \(8,24\) of 24.
At \(h=13\), the type is \(A_{12}^2\), whose associated deep
hole has 26 vertices in the repository's description.
See [Borcherds's Leech paper](https://math.berkeley.edu/~reb/papers/leech/leech.pdf)
for the deep-hole correspondence and
[THM-4758](https://github.com/eliottcassidy2000/math/blob/1b91742ce/01-canon/theorems/THM-4758-leech-deep-hole-swap-theorem-and-the-deep-hole-theta-identity.md)
for the full pure-\(A\) formulas. There is no asserted map from
these Coxeter numbers to the Tribonacci determinant zeros.
Indeed the earlier [24-port audit](https://github.com/eliottcassidy2000/math/blob/095ef14e1e88351526125b56f2679ed5d7d40cbe/05-knowledge/results/session291_cluster_ports_20261113/README.md)
found its natural paired-port code has distance 2, not Golay's 8.
Counting 24 named objects does not supply Leech glue.

## 5. A real two-dimensional gluing test at the reset

Put \(w=(47,56)\), \(w^\perp=(-56,47)\), and
\(Q=w^tw=5345=5\cdot1069\); \(1069\) is prime by trial division.
The reset is \(C_{13}=ww^t\).
The integral lattice \(L=\mathbb Zw+\mathbb Zw^\perp\) has
index \(Q\) in the standard integer lattice and Gram matrix \(QI\).
Write the two scalar boundary probes
\[
 U=w^tz,\qquad V=(w^\perp)^tz .
\]
Then the image of integral vectors is **exactly**
\[
 V\equiv 2387 U\pmod{5345},\qquad 2387^2\equiv-1\pmod{5345}.
\]
In fact the inverse coordinates are
\((47U-56V,\;56U+47V)/Q\), so the congruence follows in both
directions by invertibility of 47 and 56 modulo \(Q\).
In dual coordinates \(L^*/L\cong(\mathbb Z/Q)^2\), the graph
of this congruence is its actual gluing subgroup into
\(\mathbb Z^2\). It is isotropic for the discriminant *bilinear*
form; the lattice is odd, so no even quadratic-form or Leech
identification is being made. Orthogonality over the reals
does **not** mean that the two probes are independently realizable
over the integers.

The reset retains \(U\) and loses the exact value of \(V\).
A fixed \(U\) determines the residue of \(V\) but leaves the
integer lift free, corresponding to shifts \(z\mapsto z+jw^\perp\).
If the second physical entry must be divisible by \(2^r\),
differences of allowed shifts must satisfy \(2^r\mid j\)
because 47 is odd. Positivity restricts the remaining interval;
the actual prefix grammar must still be checked. This is a
concrete and testable "retain the transverse coordinate" rule.

A hostile for a norm-only decoder is
\(5345=4^2+73^2=47^2+56^2\). The alternative vector \((4,73)\)
gives the same Gram matrix and index (from the two
Gaussian factorizations of \(5\cdot1069\)), but a different
glue graph and a different retained linear form. It is not
a corridor reset direction. This makes the Leech analogy
operational at rank two: the gluing subgroup, not just rank,
determinant, or a named norm shell, must enter the transition
rule.
