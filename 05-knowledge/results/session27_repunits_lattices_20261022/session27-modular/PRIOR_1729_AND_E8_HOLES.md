# Bounded addendum: recovered1729 mechanisms and the E8 hole analogue

Status: prior-work recovery plus a proved null-quotient consequence. No change to the frozen main report. No repository edits or new theorem IDs. Source ref for recovered scripts:54bf96642934c6cdf01f02e364f3bba00892ae22.

## A.1729 was already more than one historical coincidence

The remote search found three genuinely different prior lanes. Their status must not be merged.

1. **Exact recurrence/near-miss lane.** `04-computation/ramanujan_near_miss_kpc1.py` (June10,2026; header routes HYP-2370) and its verify/verify2 companions recover Ramanujan's rational generating functions

   A(x)=(1+53x+9x²)/D(x), B(x)=(2−26x−12x²)/D(x), C(x)=(2+8x−10x²)/D(x),
   D(x)=1−82x−82x²+x³.

   Every sequence follows u_(n+3)=82u_(n+2)+82u_(n+1)−u_n and has a canonical integral extension to all n∈Z. At n=−1 the triple is(−9,12,10), so a_n³+b_n³−c_n³=−1 becomes9³+10³=12³+1=1729. Negative indices are the reversible-recurrence extension, not negative coefficients of the Taylor series at0.

   The inherited script proves the all-n identity by a27-dimensional tensor-cube companion recurrence: the cubic expressions and(−1)^n satisfy the same monic reversible recurrence, and27 consecutive exact zero residuals suffice by Cayley–Hamilton. The full fetched script, copied unchanged to `ramanujan_near_miss_inherited.py`, reran successfully; `ramanujan_near_miss_replay.txt` records ALL CHECKS PASS. This is a reusable algebraic mechanism, not a single matching integer. It is inherited work, not a Session27 discovery.

   **Replay does not validate all printed commentary.** The inherited script calls X³+Y³=Z³+W³ a “quartic surface”; the homogeneous equation is cubic. Its claim that the finite Eisenstein unit group makes±1 the only possible near-miss is also unjustified: unit-valued residuals are a chosen subproblem, and arbitrary integral residuals exist, e.g.1³+1³−2³=−6. The exact cubic recurrence identity does not rely on either sentence.

2. **Arithmetic-geometry lane, CITED and freshly inspected.** Ono–Trebat-Leder's corrected arXiv v5, https://arxiv.org/html/1510.00735v5 , Theorems1.1/1.3 give

   k(T)=63(3T²−3T+1)(T²+T+1)(T²−3T+3),
   X³+Y³=k(T).

   This elliptic curve has rank2 over Q(T), rank4 over Q(sqrt(−3))(T), and its smooth minimal elliptic surface is K3 with geometric Picard number18. The geometric proof uses six type-IV bad fibers and the Shioda–Tate relationρ=rank+14. The two rational sections are written explicitly in the paper. This is an authentic1729/arithmetic-geometry connection; it does not identify a tournament invariant with a K3 period. The2017 erratum https://link.springer.com/article/10.1007/s40993-017-0076-8 supplies omitted Ramanujan/Hirschhorn historical attribution, rather than retracting these theorems. Distinguish the raw cubic surface X³+Y³=Z³+W³ from the associated elliptic K3 construction. The family has j=0, so it is not the j(i)=1728 numerical lane.

3. **Tournament ratio and asymptotic lane.** `04-computation/ham_excess_spectra_kps_S134.py` (July26,2026; HYP-9028) separately records H(Paley11)/|Aut|=1729 and checks Korselt's criterion, while the p=19 ratio is not Carmichael. It also samples H(T)/(n!/2^(n−1)) for rotational and Paley tournaments. Its e-limit is explicitly consistency with finite data, not a theorem. The source correctly distinguishes H(R_n)/n from full automorphism orbit count when only rotations have been accounted for. No large run was repeated in this addendum.

4. **Metaphorical modular lane.** `04-computation/ramanujan_1729_bridge_s133.py`, `bridge_1728_1729_s122.py`, and `moonshine_plus_one_s18g.py` combine j(i)+1, factorization, OCR denominators and Monster factors. The fetched S133 header asserts a “tournament discriminant as modular form” but supplies no valid modular source-to-target map. Preserve numerical identities and provenance; do not promote this rhetoric to a theorem. In particular the same additive+1 in many formulas does not define an intertwining operation.

Read-only recovered URLs follow the pattern:
https://github.com/eliottcassidy2000/math/blob/54bf96642934c6cdf01f02e364f3bba00892ae22/04-computation/ramanujan_near_miss_kpc1.py
and the corresponding filenames above. The remote search catalogue remains in the repunit worker's deliverable; this bounded pass read the near-miss file in full and the other two specified files through line120.

## B. The lower-rank Lorentzian analogue is exact, but not another Niemeier census

**PROVED, importing the classical uniqueness of positive even unimodular rank8 E8.** For L=II_(9,1)≅E8⊕U and any primitive isotropic v∈L, the quotient v-perp/Zv is E8.

Proof: unimodularity and primitivity give w∈L with(v,w)=1. Since L is even, replace w by w−(w,w)v/2; this remains integral, pairs1 with v, and has norm0. Thus Zv+Zw is an integral unimodular hyperbolic plane U, splits orthogonally from L, and its complement is positive even unimodular rank8. It is E8 by the rank8 classification. The complement identifies with v-perp/Zv. Two such decompositions are isometric by mapping their hyperbolic bases and their E8 complements; therefore there is one primitive isotropic-vector orbit under the full orthogonal group (equivalently one ray cusp after the usual cone convention). ∎

The same proof in II_(25,1) produces rank24 complements; there the24 Niemeier isometry classes give24 cusp types. The difference is the classification of positive-definite complements, not subtraction26−2 versus10−2 alone. The rank8 uniqueness/classification input is recorded in the historical introduction of Borcherds, https://math.berkeley.edu/~reb/papers/leech/leech.pdf . The null-quotient argument here is written out to avoid relying on an unstated analogy.

**CITED E8 geometry.** At minimum squared norm2, E8 has covering radius1. A deep hole c=e1 has nearest lattice points e1±e_i (i=1,…,8), forming a16-vertex cross-polytope. Its affine lattice-isometry stabilizer is W(D8). This is not the Leech's23-type deep-hole/Niemeier correspondence. The primary treatment is *Polarization of lattices: Stable cold spots and spherical designs*, Lemma3.6 and AppendixA.3.1, https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/polarization-of-lattices-stable-cold-spots-and-spherical-designs/16B869EA8216480737BB2BD00B115F9D . The reference also distinguishes deep and shallow E8 holes; do not infer every hole is a deep hole. This addendum does not need a new numerical orbit census.

**Elementary hostile/normalization test.** In E8=D8∪((1/2,…,1/2)+D8), every listed e1±e_i belongs to D8 and is at distance1 from e1. No integer lattice point can lie closer since e1 itself has odd coordinate sum and is not in D8; every half-integer point has squared distance≥8/4=2. Hence these16 points are precisely the nearest set at c. In the Leech normalization the deep-hole radius is sqrt(2), not1, and its rootless lattice has minimum squared norm4, not2. One cannot transport the hole diagram merely by copying the label8→24.

Preserved/lost ledger: null quotient preserves evenness, unimodularity and induced positive metric but discards the chosen hyperbolic-plane embedding and Lorentzian coordinates. A hole→nearest-shell construction preserves distances about the marked center but is not the same operation as the null quotient. A lawful comparison must retain the center, normalization, lattice embedding and group action. None of these constructions currently sends forbidden H-values7,21 to E8/Leech classes.
