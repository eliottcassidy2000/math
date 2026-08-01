---
id: THM-3000
title: "Fixed-edge cumulant-curvature universality and bounded-jet transfer"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
source: codex-gmc-fixed-edge-curvature-2026-07-31 (reservation); klein-S428 (proof, mechanism, third-edge closed form, hostile)
depends_on: []
related:
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity
  - THM-2988-first-gap-self-curvature-negativity
  - THM-2994-first-gap-hurwitz-hermite-biehler-prefix
  - MISTAKE-339
script: 04-computation/gmc_fixed_edge_curvature_universality_and_bounded_jet_transfer_thm3000.py
output: 05-knowledge/results/gmc_fixed_edge_curvature_universality_and_bounded_jet_transfer_thm3000.out
script_sha256: 6cd895b12fa53a2b772f2d03a3b014a3d0eea0f8e6ce9f4559ce1df086a578cc
output_sha256: 4c5cb347fe3fc5c12cb4b0b3bd90e002e6311aff72945f30614ffd979871d92b
hash_basis: LF-normalized bytes
---

# THM-3000 -- fixed-edge cumulant-curvature universality and bounded-jet transfer

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is an abstract coefficient theorem about normalized Newton ratios of an
arbitrary real polynomial.  It has **no** proved dependency on any first-gap,
Macaulay, resultant, wall, or GMC object, and it proves **no** GMC(2), ULC,
arbitrary-radial, or all-width family conclusion by itself.  Its use for the
first-gap family is isolated in section 7 and is explicitly CONDITIONAL.

The reservation stub was opened by `codex-gmc-fixed-edge-curvature-2026-07-31`
with exactly this intended statement; this file discharges it.

## 1. Universe and coordinates

Let

    N(n)=sum_(i=0)^d a_i n^i,  a_d>0,

and, for a fixed edge index `k>=2`, assume `a_(d-1),...,a_(d-k-1)>0`.  As in
THM-2997 (1) put

    h_j=a_(d-j)/(a_d binom(d,j)),  h_0=1,
    R_j=h_j^2/(h_(j-1)h_(j+1)).                        (1)

Define the **formal log jets** by

    log(N(n)/(a_d n^d))=sum_(j>=1)(ell_j/j)n^(-j),
    ell_1=u=a_(d-1)/a_d,                                (2)

extending THM-2997 (4) to all orders; equivalently `ell_j=(-1)^(j-1)p_j`,
where `p_j` is the `j`-th power sum of the roots of `N` (`N=a_d prod(n+r_m)`).
Define the **normalized log jets** and the **root-measure moments**

    J_j=ell_j/u^j,                                      (3)
    m_j=p_j/d,  x=m_2/m_1^2=-d ell2/u^2,  z=m_3/m_1^3=d^2 ell3/u^3,
    w=m_4/m_1^4=-d^3 ell4/u^4.                          (4)

`x,z` are exactly THM-2997 (7); `w` is the new third-edge coordinate.  All of
`x,z,w,J_j` are invariant under `n->cn`, `N->cN`.

By THM-2997 (2), `log(R_k/R_(k-1))=-Delta^3(log h)_(k-2)`, so every statement
below is a statement about the third difference of the normalized
log-coefficient sequence.

## 2. Statement (universality)

Grade the ring `Q[J_2,J_3,...][d^(-1)]` by

    wt(J_j)=j-1,    wt(d^(-1))=1.                       (5)

**Theorem.**  For every fixed `k>=2` there is a weight-graded expansion

    log(R_k/R_(k-1))=sum_(r>=2) W_r^(k),                (6)

with `W_r^(k)` homogeneous of weight `r`, and

    W_2^(k)=3J_2^2-2J_3-d^(-2)                          (7)

**independent of `k`**.  For `r>=3`, `W_r^(k)` is a polynomial in `k` of degree
exactly `r-2`.  Moreover the jet `J_j` first occurs in weight `r=j-1`, with the
exact coefficient

    [J_j] W_(j-1)^(k)=-(j-1)! binom(k-2,j-3).           (8)

After scaling to `m_1=1`, if the relative moments
`q_j=m_j/m_1^j` are bounded (equivalently `J_j=O(d^(1-j))`), weight `r` is
order `d^(-r)` and (7) reads

    log(R_k/R_(k-1))=(3x^2-2z-1)d^(-2)+O_k(d^(-3)).     (9)

The leading coefficient is **literally** `3x^2-2z-1`: not merely proportional
to it, and with no `k`-dependent constant.  Weight two is the unique order at
which the edge index drops out, by the degree statement.

Verified symbolically for `k=2..9`, `r<=5`, `j<=6`.

## 3. Cumulant form

Let `kappa_1,kappa_2,kappa_3` be the first three cumulants of the empirical
root measure (`kappa_2` variance, `kappa_3` third central moment).  Then

    3x^2-2z-1=(3kappa_2^2-2kappa_1 kappa_3)/kappa_1^4.  (10)

Hence the universal curvature is **positive iff `3kappa_2^2>2kappa_1kappa_3`**.
Consequences, all sharp:

- equal roots give `kappa_2=kappa_3=0`, hence one exact zero of the curvature
  (the boundary that THM-2988's self-curvature negativity probes from the
  inside); globally the zero set is the full surface in `(10)`, not just this
  point;
- any symmetric root measure has `kappa_3=0` and curvature
  `3kappa_2^2/kappa_1^4>=0`, with equality only at equal roots.  This is the
  precise sense in which the certificate is "Gaussian": a Gaussian root profile
  can never fail it;
- failure needs strictly positive skew past the exact threshold
  `kappa_3>3kappa_2^2/(2kappa_1)`.

## 4. Mechanism (why `k` cancels)

Write `h_j=j!e_j/(d)_j` with `(d)_j` the falling factorial, so
`log h_j=log(j!e_j)-log((d)_j)` and

    log(R_k/R_(k-1))=-Delta^3 log(k! e_k)+Delta^3 log((d)_k)   at k-2.

**The `-1` is combinatorial.**  `log((d)_k)=k log d-binom(k,2)/d-S_2(k)/(2d^2)-...`
with `S_2(k)=sum_(i<k)i^2`.  `Delta^3` kills the linear and quadratic parts;
`Delta^3 S_2=2`, so the falling-factorial normalization alone contributes
exactly `-d^(-2)`.  This term has nothing to do with the coefficients of `N`.

**The `3m_2^2-2m_3` is a Moebius/set-partition residue.**  By inclusion-exclusion
over set partitions of `[k]` (`k!e_k=sum_pi mu(pi) prod_B p_(|B|)`, `m_1=1`),

    k!e_k=d^k[1-binom(k,2)m_2/d+(3binom(k,4)m_2^2+2binom(k,3)m_3)/d^2+O(d^(-3))],

so the `d^(-2)` bracket of `log(k!e_k)` is
`3binom(k,4)m_2^2+2binom(k,3)m_3-binom(k,2)^2m_2^2/2`.  The two `m_2^2` pieces
collapse:

    3binom(k,4)-binom(k,2)^2/2=k(k-1)(3-2k)/4,          (11)

a **cubic** in `k` with `Delta^3=-3`; and `Delta^3[2binom(k,3)]=2`.  Therefore

    -Delta^3 log(k!e_k)|_(d^(-2))=3m_2^2-2m_3,

with every quartic and quadratic `k`-dependence cancelling inside (11).  This
is the weight-two content of universality: a quartic-in-`k` term and a
squared-quadratic-in-`k` term differ by a cubic, and `Delta^3` of a cubic is a
constant.

Here is the all-order argument. Scale to `m_1=1` and put
`A_k=k!e_k/d^k`. View its set-partition formula as a labelled hard-core
polymer gas: a non-singleton block `B` is a polymer of excess `|B|-1`, while
unused labels are singleton blocks. The coefficient of `d^(-r)` has total
polymer excess `sum_B(|B|-1)=r`. The standard Mayer expansion of `log A_k`
retains only clusters whose polymer-intersection graph is connected. For such
a cluster, adjoining each polymer along an intersection shows that its union
has at most `1+sum_B(|B|-1)=r+1` labels. Hence its contribution is a polynomial
in `k` of degree at most `r+1`. The single polymer of size `r+1` contributes

    (-1)^r r! binom(k,r+1) q_(r+1),

and no other term contains `q_(r+1)`. Hence the degree is exactly `r+1` in the
formal moment ring. The falling-factorial term has degree at most `r+1` and
contains no `q_(r+1)`, so it cannot cancel this top-moment coefficient.
Applying `-Delta^3` at `k-2` gives exact degree `r-2`. Since
`Delta^3 binom(k,r+1)=binom(k,r-2)`, and
`q_j/d^(j-1)=(-1)^(j-1)J_j`, the same calculation gives `(8)` for every `j`.
This proves the all-order degree and first-occurrence statements; the finite
symbolic checks are corroboration, not their proof.

## 5. Exact closed forms (no asymptotics)

With `e_1=u`, `e_2=(u^2+ell2)/2`, `e_3=(u^3+3u ell2+2ell3)/6`,
`e_4=(u^4+6u^2 ell2+3ell2^2+8u ell3+6ell4)/24` and `h_j=j!e_j/(d)_j`:

**Second edge** (recovers THM-2997 (8) exactly):

    h_2^3-h_1^3h_3=G_2 u^6/(d^6(d-1)^3(d-2)),
    G_2=d^2(3x^2-2z-1)+d(-x^3-6x^2+3x+4z)+2x^3-2z.      (12)

**Third edge** (new):

    h_1h_3^3-h_2^3h_4=G_3 u^10/(d^10(d-1)^4(d-2)^3(d-3)),   (13)

    G_3=d^6(3x^2-2z-1)
      +d^5(6w+x^3-12x^2-12xz+9x+8z)
      +d^4(-18wx-24w-15x^4-4x^3+30x^2z-15x^2+48xz+12z^2-14z)
      +d^3(18wx^2+72wx+24w+3x^5+60x^4+8x^3z+31x^3-120x^2z-36xz^2-12xz-48z^2)
      -2d^2(3wx^3+36wx^2+36wx+6x^5+30x^4+16x^3z-33x^2z-72xz^2-4z^3-18z^2)
      +4d(6wx^3+18wx^2+3x^5+8x^3z-27xz^2-8z^3)
      -24(wx^3-z^3).                                    (14)

Since `d>0` and `(d)_4>0` for `d>3`, and `h_1^3h_3,h_2^3h_4>0` under the
positivity hypothesis, `sign(log(R_2/R_1))=sign G_2` and
`sign(log(R_3/R_2))=sign G_3`.  Both leading coefficients are the **same**
`3x^2-2z-1`, at `d^2` and `d^6` respectively, matching (9) after dividing by
`d^4` and `d^8`.

**Top-jet coefficients factor.**  In (12) and (14),

    [z]G_2=-2(d-1)^2,       [w]G_3=+6(d-x)^3(d-2)^2.    (15)

These are the `j=k+1` instances of (8): `-2!binom(0,0)` and `+3!binom(1,1)`
after the weight normalization, with the alternating sign `(-1)^j`.

## 6. Hostile: bounded jets are indispensable, and the sharp condition is graded

`w` enters `G_3` at `d^5` while the curvature sits at `d^6`.  So the third edge
is controlled by the curvature **only while `w=o(d)`**.  This is attained.

**Hostile family.**  Fix `x,z` with `c:=3x^2-2z-1>0`, fix `alpha>c/6`, and for
each `d` set `m_1=1` (`u=d`) and

    ell2=-x u^2/d,  ell3=z u^3/d^2,  ell4=alpha u^4/d^2   (i.e. w=-alpha d).

Then `a_d,...,a_(d-4)>0` for all large `d` (the `ell4` sign makes `e_4` larger,
not smaller), so the family is a legitimate positive-coefficient universe, and
the remaining coefficients are free.  By (12) and (14),

    log(R_2/R_1) -> +c d^(-2),
    log(R_3/R_2) -> (c-6alpha) d^(-2)<0.

**Exact instance in the output:** `x=1`, `z=1/4`, `c=3/2`, `alpha=1/2`, so
`6alpha=3>3/2`.  For `d=200,400,800,1600,3200` the exact rational ratios give
`R_1<R_2` and `R_2>R_3` at every `d`, with `d^2 log(R_2/R_1)->3/2` and
`d^2 log(R_3/R_2)->-3/2`.

**Consequence.**  Uniform boundedness of the normalized log jets is *not* a
technical convenience: without it the second-edge certificate does not transfer
to the third edge, and the two edges can have opposite signs at every width.
The sharp higher-jet hypothesis is **graded, not uniform**:

    m_j/m_1^j=o(d^(j-3))   for 4<=j<=k+1.               (16)

Together with bounded `x,z` and curvature separated from zero, `(16)` suffices
for `(9)` to determine the sign of edge `k`, by `(8)`. Jet `j=3` belongs to the
leading curvature and is not a remainder hypothesis. Uniform boundedness of
all relative moments is a stronger sufficient special case, because
`O(1)=o(d^(j-3))` for every `j>=4`; the `j=4` hostile with `w asymp d` shows
the first threshold is attained. The alternating signs in `(8)` mean the
failure direction alternates with `j`.

## 7. CONDITIONAL transfer to the first-gap family

**This section is CONDITIONAL and is not used by sections 2--6.**

THM-2997 sections 5--6 prove, for the encoded continuation wall invoice, the
rational box

    A>0,  129/100<=x<=2,  0<=z<=39/20,  d>=701,         (17)

for every `M>=43`, and on that box `3x^2-2z-1>=923/10000`.  By section 2, the
**same** box gives the same positive leading curvature at **every** fixed edge
`k>=2`.  What is missing per edge is only the jet invoice (16).

For `k=3` the invoice is a single one-sided bound. By `(15)`, `[w]G_3>0` on
`(17)`, and uniformly on that compact box

    G_3/d^6=(3x^2-2z-1)+6w/d+O((1+|w|/d)/d).            (18)

Consequently `w=o(d)` suffices. More generally, the sharp uniform asymptotic
condition is

    liminf w/d > -923/60000,                            (19)

or, using `w=-d^3 ell4/u^4`,

    limsup d^2 ell4/u^4 < 923/60000.

A convenient strict sufficient form is
`ell4 <= (923/60000-epsilon)u^4/d^2` eventually for some `epsilon>0`.
The strict margin is essential; the former finite threshold table and its
non-strict boundary were false and are recorded in MISTAKE-339.

For the first-gap family THM-2997 (24) gives `u/M^2->131/12` and
`d/M->62/3`. Thus the scale on the right is asymptotic to `0.511...M^6`,
while any core whose roots have modulus `O(M)` has
`|p_4|=O(dM^4)=O(M^5)`. Equivalently
`d^2|p_4|/u^4=O(M^(-1))`, which clears `(19)` by a factor of `M` regardless
of the sign of `p_4`.

**So the third-edge invoice for the first-gap family is not a delicate exact
sign computation: a crude `|p_4|=O(M^5)` bound suffices asymptotically.**
The precise obligation handed to the resultant lane is therefore:

> supply the fourth reduced resultant log-jet numerator `P_4` and the wall's
> fourth power sum `w_4` (extending THM-2997 (21)--(22) by one index), or any
> `O(M^5)` a-priori bound on the core's fourth power sum.

This does not prove full ratio no-return, ULC, or GMC(2), and it does not
remove THM-2989/2997's continuation wall hypothesis.  It replaces a per-edge
Macaulay/resultant recomputation by a per-edge scalar jet bound.

## 8. Boundaries, quotient losses, and what is NOT claimed

- (9) is an asymptotic in `d` at **fixed** `k`.  It says nothing about `k`
  growing with `d`; the `k`-degree `r-2` statement shows the expansion is not
  uniform in `k`, so a `k asymp d` statement is a different theorem and is
  **not** proved here.
- The curvature is degree-`0` homogeneous, so it forgets the overall scale of
  the roots.  Restoration sidecar: `u` and `d` separately.
- The curvature vanishes on a codimension-one surface
  (`3kappa_2^2=2kappa_1kappa_3`).  On that surface every edge is decided at
  weight three, where `k`-dependence is already present (degree `1`).  So there
  is **no** universal tie-break: on the zero set the edges genuinely separate.
- Sign transfer requires `|3x^2-2z-1|` bounded away from `0` relative to the
  weight-three remainder; (17) supplies `923/10000` for the first-gap family.
- Nothing here supplies THM-2842's missing radial-variance multipliers, and
  nothing here is a moment decoder; MISTAKE-211 still applies.

## 9. Independent hostile audit

The independent audit rederived the connected-partition proof above, checked
the conversion between raw moments and normalized log jets, and replayed both
exact closed forms. It found and repaired two pre-promotion boundary defects.
First, `(16)` formerly started at `j=3`, although `J_3` is part of `W_2`; the
remainder condition starts at `j=4`. Second, the former finite box-threshold
table was not a safe uniform bound. At the exact corner

    d=701, x=129/100, z=39/20, w=-149/20,

the old claim applied (`-7.45>=-7.4613`) but the exact numerator is

    G_3=-114191274399994230172453/10000000000<0.

The repaired strict asymptotic statement `(18)--(19)` survives and is exactly
what the first-gap `O(M^5)` argument needs. The audit also checked the hostile
family's top coefficients directly and confirmed that it inhabits the stated
arbitrary-real-polynomial universe; it makes no real-rootedness claim.

## 10. Reproduction

    python3 04-computation/gmc_fixed_edge_curvature_universality_and_bounded_jet_transfer_thm3000.py

The script performs seven independent checks (A universality, B cumulant
identity, C set-partition mechanism, D/G weight grading and the first-occurrence
binomial law, E exact closed forms with a match against THM-2997 (8), a positive
control on four real-rooted families in exact rational arithmetic, and F the
hostile).  All report `True`.  Sections 2 and 5 are two logically independent
derivations (series in `1/d` from Newton's identities; exact rational-function
algebra from the exponential formula) that agree on the leading coefficient and
on (8) for `j=4`.
