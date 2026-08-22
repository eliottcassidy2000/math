---
id: THM-3501
title: "Rule 30 universal-cover Green potential and slack-holonomy seam"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A literal Green-path
  slack potential on the physical
  arrival-phase cover differentiates to the terminal current.  On every hard
  unwrapped depth k<P_k, its cyclic arc gives a coefficientwise physical-image
  slack--phase lift.  Universal endpoint arcs differ from cyclic arcs by an
  exact slack-holonomy seam whose low phase-Hasse moments are the complete
  descent obstruction.  An explicit restart contraction and either a
  raw-selected top carrier or a support-minimal Frobenius carrier give a
  degree-(k-4) lift of the THM-3492 inward axis.  No Rule 30 prize consequence
  is claimed.
source: root/rule30-sharp-unlocks/universal-cover-slack-phase-coupling/2026-08-16
audit: >
  PASS (2026-08-16), adversarial independent proof and replay audit.  The
  auditor rederived the Duhamel primitive, cut/seam identities, complete
  low-Hasse descent obstruction, restart contraction, both exact repairs, and
  the typed THM-3500 scalar bridge.  Every displayed depth-five/depth-six
  control and current preimage is explicitly gated.  Ordinary, optimized, and
  stored transcripts are byte-identical.  Universal scope rests on the
  displayed algebraic proofs, not on the finite universe k=5,...,16.
depends_on:
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas
  - THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
  - THM-3489-rule30-packed-restart-and-pointed-pascal-face
  - THM-3492-rule30-fiber
related:
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
script: 04-computation/rule30_universal_cover_slack_holonomy_thm3501.py
output: 05-knowledge/results/rule30_universal_cover_slack_holonomy_thm3501.out
script_sha256: 58f4f3a48f5eb7afecb88ae394540f7d82262116de7ac8ea1a4c3ba3ff05667e
output_sha256: 804e6f376e2c81b0ced3f1e2b9bb343d82b4bf06d3e23180beb62ce33a1961a5
hash_basis: raw bytes
---

# THM-3501 -- Rule 30 universal-cover Green potential and slack-holonomy seam

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3492 proves that the complete inward slack polynomial and the complete
terminal phase profile have a large family of algebraic common refinements.
It deliberately does not call any of those sections physical.  This theorem
starts one level earlier, at literal Rule 30 collision events and
their Rule-150 Green paths.  The resulting universal-cover potential has an
exact terminal-current derivative.  A physical phase cut then exposes two
different objects:

1. universal endpoint arcs, which retain the literal endpoint potential; and
2. cyclic current arcs, which descend coefficientwise through the physical
   arc operator.

Their difference is exactly one slack-holonomy seam.  After an explicit
slack-chain contraction, the cyclic object gives a degree-minimal common
refinement of the two THM-3492 axes.  The construction selects a mixed-kernel
class from raw spacetime data; it does not prove that the two marginal axes
alone select that class.

## 1. Inheritance, scope, and event conventions

All algebra below is over `F_2`.  Use centered Rule 30 rows

```text
a_0(0)=1,       a_0(j)=0 for j!=0,
f(l,c,r)=l+c+r+cr,
e_s(j)=a_s(j)a_s(j+1),                               (1)
```

and the Rule-150 Laurent operator and Green shell

```text
Q(y)=y^(-1)+1+y,
H_d(n)=[y^d]Q(y)^n
      =[x^(n+d)](1+x+x^2)^n.                        (2)
```

Put `H_d(n)=0` when `d<0`, `n<0`, or `d>n`.  THM-3468's Duhamel identity is

```text
CalA_n(y)=sum_j a_n(j)y^j,       CalE_s(y)=sum_j e_s(j)y^j,
CalA_n=Q^n CalA_0+sum_(0<=s<n) Q^(n-1-s)CalE_s.      (3)
```

The calligraphic row polynomials in (3) are distinct from the cyclic arc
operator `A_k` below.

For the terminal depth `k`, retain

```text
T_k(h)=a_(k+h)(h),
T_k(h+1)+T_k(h)=Q_k(h),
F_k=A_k Q_k,       A_k=1+S+...+S^(k-1).             (4)
```

The inheritance pass is:

1. closest proved mechanism: THM-3471's Green-slack current and THM-3489's
   pointed terminal arc;
2. canonical hostile: THM-3492's depth-six two-section ambiguity;
3. corrected near miss: form products before scalarization, but distinguish
   the phase **cover** from its cyclic descent; and
4. least-used sidecars: the physical cut `h=-k`, the slack holonomy across
   that cut, and the raw Green-path carrier of the monic inward face.

The main chain-refinement theorem is stated in the hard unwrapped regime

```text
k>=5,       p=P_k,       k<p.                        (5)
```

This is exactly the unwrapped line left open in THM-3489.  The raw potential
itself is valid for every physical `h>=-k`.  Strict wraps require a second
or further holonomy crossing and are not included in the cyclic-refinement
conclusion below.

## 2. The raw universal-cover Green potential

For `h>=-k`, put `n=k+h` and define

```text
E_k(h;q)
 =H_|h|(n) q^(n-|h|)
  +sum_(0<=s<n) sum_(j in Z)
    e_s(j) H_|h-j|(n-1-s)
      q^(n-1-s-|h-j|).                              (6)
```

Terms with distance greater than horizon are zero and are omitted.  Formula
(6) has a literal event interpretation.  The first term is the initial seed
atom.  A later atom consists of a physical occupied bond `(s,j)` and a
ternary Rule-150 path from that bond to `(n,h)`; `H` retains the parity of all
such paths, and the exponent is their excess time over absolute displacement.

### Theorem 2.1 (Duhamel primitive, degree, and raw monicity)

For every `k>=1` and `h>=-k`,

```text
boxed:
 E_k(h;1)=T_k(h),
 deg_q E_k(h;q)<=k,
 E_k(-k;q)=0,
 [q^k]E_k(0;q)=1.                                  (7)
```

In particular the raw center polynomial

```text
R_k^raw(q):=E_k(0;q)                                (8)
```

is monic of degree exactly `k`.

### Proof

Set `q=1` in (6) and take the coefficient at `y^h` in (3).  Palindromy of
`Q^m` gives exactly the two Green coefficients in (6), proving the first
identity.  At `h=-k`, the target is the exterior seed cell `(0,-k)`, so every
term is zero.

For the degree bound, the initial exponent is `k+h-|h|<=k`.  A live collision
has `|j|<=s`.  If `h>=0` and `j<=h`, its slack is

```text
k+h-1-s-(h-j)=k-1-s+j<=k-1.                        (9)
```

If `j>h`, then `h<s` and the same bound is strict.  For `h<0`, the source
exponent is at most `k+h-1<=k-1`.  At `h=0`, sign reversal pairs every
nonconstant central ternary path, leaving `H_0(k)=1`; hence the initial atom
is `q^k`, while every collision has smaller degree.  This proves (7).
`square`

The scalar current orientation follows directly from the local rule.  With
`n=k+h`,

```text
T_k(h+1)+T_k(h)
 =a_(n+1)(h+1)+a_n(h)
 =a_n(h+1) OR a_n(h+2)
 =Q_k(h).                                            (10)
```

Thus `E_k(h+1;q)+E_k(h;q)` is an exact slack-graded primitive of the physical
terminal current.

## 3. Cut-cyclic arcs and universal endpoint arcs

Assume (5) and choose the physical cut

```text
h_r=-k+r,       0<=r<=p.                             (11)
```

Modulo `p`, the cut origin is `h_0=-k=p-k=:h_*`, exactly the marked phase of
THM-3492 in the hard regime.

Define the `p` cut-current entries

```text
d_r(q)=E_k(h_(r+1);q)+E_k(h_r;q),       0<=r<p,     (12)
```

and its cyclic length-`k` arc

```text
Z_k^cyc(q,r)=sum_(i=0)^(k-1)d_(r+i mod p)(q).        (13)
```

Also retain the literal universal-cover endpoint difference

```text
Z_k^univ(q,r)=E_k(h_r+k;q)+E_k(h_r;q).               (14)
```

### Theorem 3.1 (physical current lift)

The cyclic table has

```text
boxed:
 Z_k^cyc(1,r)=F_k(h_r),
 Z_k^cyc(q,0)=R_k^raw(q),
 [q^m]Z_k^cyc in im(A_k) for every m.                (15)
```

### Proof

Equations (7) and (10) give `d_r(1)=Q_k(h_r)`.  Applying the cyclic arc gives
the first line of (15).  Since `k<p`, the marked arc `r=0,...,k-1` does not
cross the cut, so (12) telescopes to

```text
E_k(h_k;q)+E_k(h_0;q)=E_k(0;q)+E_k(-k;q)=R_k^raw(q). (16)
```

Finally (13) is literally `A_k d` coefficient by coefficient.  `square`

This is already more than the marginal fiber product: the phase and raw
slack axes now share one declared spacetime/path base.  It is not yet the
THM-3492 inward normalization.

## 4. Slack holonomy is the exact seam, and only its low tomography obstructs descent

For `0<=a<k`, define the slack holonomy word

```text
Omega_a(q)=E_k(h_0+p+a;q)+E_k(h_0+a;q).              (17)
```

THM-3489's packed restart gives only its scalar value:

```text
Omega_a(1)=epsilon_k                 for every a.    (18)
```

Put `W=Z^univ+Z^cyc`.

### Theorem 4.1 (exact one-seam formula)

For `0<=r<p`,

```text
boxed:
 W_r(q)=0,                              r<p-k;
 W_(p-k+a)(q)=Omega_a(q)+Omega_0(q),    0<=a<k.      (19)
```

In particular `W_r(1)=0` for every phase, but `W` need not vanish as a
slack-graded table.

### Proof

Before the arc crosses the cut, both (13) and (14) are the same ordinary
telescope.  If `r=p-k+a`, the cyclic telescope breaks into the two intervals
`r,...,p-1` and `0,...,a-1`.  Comparing their four endpoints with (14) leaves

```text
E(h_0+p+a)+E(h_0+p)+E(h_0+a)+E(h_0)
 =Omega_a+Omega_0.                                  (20)
```

Equation (18) proves scalar disappearance.  `square`

The seam is not itself the descent obstruction.  Put

```text
a_0=nu_2(k),       ell=2^a_0-1.                     (21)
```

THM-3481 identifies the cyclic arc image as

```text
I_k=im(A_k)=(X+1)^ell V_p.                           (22)
```

Let `M_j^r` denote the phase Hasse moment in the cut coordinate `r`.

### Theorem 4.2 (complete low-Hasse seam tomography)

For every `0<=j<ell`, as an identity of `q`-polynomials,

```text
boxed:
 M_j^r(W)=sum_(a=0)^(k-1) binom(a,j)Omega_a(q).      (23)
```

Consequently

```text
Z_k^univ is coefficientwise in I_k
 iff sum_a binom(a,j)Omega_a(q)=0
     for every 0<=j<ell.                             (24)
```

When `k` is odd, `ell=0` and the condition is vacuous.  Thus a nonconstant
slack-holonomy word alone is **not** an obstruction.

### Proof

Apply `M_j` to (19):

```text
M_j(W)=sum_a binom(p-k+a,j)(Omega_a+Omega_0).        (25)
```

Both `p` and `k` are divisible by `2^a_0`, so for `j<2^a_0`, Lucas gives
`binom(p-k+a,j)=binom(a,j)`.  The coefficient of `Omega_0` is

```text
sum_(a=0)^(k-1)binom(a,j)=binom(k,j+1)=0,            (26)
```

because `j+1<2^a_0` in the required range.  This proves (23).  Since
`Z^cyc` is already in `I_k`, and (22) says that the first `ell` Hasse moments
are the complete quotient, (24) follows.  `square`

Equation (23), not the raw nonconstancy of `Omega`, is the exact
Radon/Pascal-style seam detector.

## 5. Raw time zero, the time-two restart, and the inward gauge

Retain THM-3471/3488's notation

```text
B_k=H_1(k-2),
D_k(q)=B_k+sum_(u>=3)[z^k]R_u(z,q),
N=k-4.                                               (27)
```

Define the time-two full polynomial and the three-strip boundary current

```text
Dtilde_k(q)=B_k+sum_(u>=0)[z^k]R_u(z,q),
G_k(q)=sum_(u=0)^2[z^k]R_u(z,q).                     (28)
```

Then

```text
Dtilde_k=D_k+G_k,
G_k(1)=0,
deg G_k<=k-5=N-1.                                   (29)
```

Indeed THM-3471's unmarked three-strip circuit gives `G_k(1)=0`, while the
minimal live `(u,d)` pairs `(0,2),(1,2),(2,1)` in
`v=k-u-1-2d` give the displayed degree bound.

The generating series of `G_k` is exactly THM-3471 equation (14):

```text
sum_k G_k(q)z^k
 =z(1+q)W_q^2((1+q)W_q+z)
   /((1+qz)(1+W_q^2)).                              (30)
```

The truly raw time-zero potential has one additional, completely explicit
restart chain.

### Theorem 5.1 (exact restart and inward normalization)

For every `k>=5`,

```text
boxed:
 R_k^raw(q)+Dtilde_k(q)
  =Gamma_k(q)
  :=q^k+q^(k-2)+B_k(q^(k-3)+1),                     (31)

 R_k^raw(q)+D_k(q)=Gamma_k(q)+G_k(q)=:K_k(q),       (32)

 K_k(1)=0.                                          (33)
```

Write `G_k=(q+1)Ghat_k`.  Then the raw-to-inward homotopy is explicit:

```text
K_k=(q+1)L_k,
L_k=q^(k-2)(q+1)
    +B_k(1+q+...+q^(k-4))+Ghat_k.                   (33a)
```

Thus raw-to-inward normalization is an exact slack chain homotopy, not only
an equality after evaluating at one.

### Proof

In the time-zero Duhamel formula, the initial seed contributes `q^k`.
There is no time-zero collision.  At time one,

```text
e_1(-1)=e_1(0)=1.                                   (34)
```

The first bond contributes `B_k q^(k-3)` and the second contributes
`q^(k-2)`.  Restarting at time two replaces those three terms by the single
unmarked linear response `B_k`: indeed the time-two row is `11001`, and its
two distance-two Green terms cancel, leaving `H_1(k-2)`.

For completeness, every later center source has the exact change of variables

```text
d=|j|,       u=s-d,       v=k-s-1-d,
k=u+1+2d+v,
H_d(k-1-s)=K(d+v,v).                                (34a)
```

The last equality is ternary palindromy.  Folding the two signs of `j` gives
THM-3471's `alpha_u(d)`, and conversely every live `(u,d,v)` at target `k`
comes from one of these `s>=2` sources (the declared folded source is zero
when `s=u+d<2`).  Hence the entire later-source sum is exactly
`sum_(u>=0)[z^k]R_u(z,q)` in both expansions.  This proves (31).
Equations (28)--(29) give (32), and specialization gives (33).  `square`

The first-three-strip current is therefore only part of the raw/inward
mismatch.  The additional term `Gamma_k` is the exact time-zero/time-two
restart debt.  Calling the whole mismatch “the three strips” would be false.

## 6. A degree-minimal chain refinement of the THM-3492 axes

Define the pointed contraction

```text
C_N(q^m)=q^m,       m<=N;
C_N(q^m)=1,         m>N,                              (35)
```

and extend linearly.  It preserves evaluation at one, commutes with every
phase-linear operator, and

```text
f+C_N(f) is divisible by q+1.                        (36)
```

Write `P_N={f in F_2[q]:deg f<=N}`.
More precisely, on monomials define

```text
H_N(q^m)=0,                       m<=N;
H_N(q^m)=1+q+...+q^(m-1),        m>N.
```

Then the operator identity

```text
id+C_N=(q+1)H_N                                     (36a)
```

is the promised explicit contraction homotopy.

Every nonzero exponent in `Gamma_k` other than its optional constant is
strictly above `N=k-4`.  Therefore
`C_N(q^k+q^(k-2))=1+1=0` and
`C_N(B_k(q^(k-3)+1))=B_k(1+1)=0`, so `C_N(Gamma_k)=0`.
Also `deg Dtilde_k<=N` by (29) and THM-3488, so (31) gives

```text
C_N(R_k^raw)=Dtilde_k.                               (37)
```

Contract the physical cut-current and arc:

```text
dbar=C_N(d),
Zbar=A_k dbar=C_N(Z_k^cyc).                          (38)
```

Then

```text
deg_q Zbar<=N,
Zbar(1,r)=F_k(h_r),
Zbar(q,0)=Dtilde_k(q),
[q^m]Zbar in I_k.                                   (39)
```

There are two exact inward repairs.

### 6.1 Raw-selected top carrier

Put

```text
eta_k(r)=[q^N]dbar(q,r),
gamma_k(r)=[q^N]Zbar(q,r),       gamma_k=A_k eta_k.  (40)
```

THM-3488's monicity and (29) give `gamma_k(0)=1`; coefficientwise physical
descent gives `gamma_k in I_k`.  Therefore

```text
boxed:
 Z_k^can(q,r)=Zbar(q,r)+G_k(q)gamma_k(r)
             =A_k(dbar+G_k eta_k)(q,r).              (41)
```

satisfies

```text
Z_k^can(1,r)=F_k(h_r),
Z_k^can(q,0)=D_k(q),
Z_k^can in P_N tensor I_k.                           (42)
```

This carrier is selected by the raw chain itself.  It need not be one phase-
Hasse monomial; away from depth six its Hasse support can be broad.

### 6.2 Support-minimal Frobenius carrier

In the cut coordinate, let

```text
phi_k(X)=1,                         a_0=0;
phi_k(X)=1+X^(2^a_0)=(X+1)^(2^a_0), a_0>=1.          (43)
```

Identify this polynomial with its cut-coordinate coefficient profile,
`phi_k(r)=[X^r]phi_k(X)`.  Then `phi_k(r=0)=1` and `phi_k in I_k`.  For even
`k`, no one-point profile is
in `I_k` because its zeroth Hasse moment is one, so the two-point carrier in
(43) has minimum possible support.  For odd `k`, `A_k` is invertible and the
one-point carrier is minimum.  Hence

```text
boxed:
 Z_k^Fr(q,r)=Zbar(q,r)+G_k(q)phi_k(r)                (44)
```

also satisfies all three conclusions in (42), with support cost one or two.

For an explicit current preimage, write

```text
A_k(X):=X^(-(k-1))(1+X+...+X^(k-1))
       =(X+1)^ell U_k(X),       U_k(1)!=0             (45)
```

The unit `U_k` is invertible in `F_2[X]/(X^p-1)`.  A preimage of `phi_k` is
`U_k^(-1)` when `a_0=0` and `(X+1)U_k^(-1)` when `a_0>=1`.  Thus (44) is not
merely an image-membership assertion: it is an explicit cyclic-current chain
repair.

Equations (41) and (44) are exact common refinements beyond THM-3492's
marginal fiber product.  They do not make the refinement unique.  The raw
top carrier and the minimum-support Frobenius carrier can give different
mixed-kernel classes with the same two axes.

Since the marked inward polynomial has degree exactly `N`, no common
refinement with that marked axis can have smaller `q`-degree.  Thus both
repairs are degree-minimal in the literal sense used here.

## 7. Minimal positive and hostile controls

### 7.1 Depth five: holonomy nonconstancy is not an obstruction

At `k=5`,

```text
p=8,       ell=0,       N=1,
R_5^raw=q+q^3+q^5,
Dtilde_5=D_5=q,
Gamma_5=q^3+q^5,
G_5=0.                                               (46)
```

The exact seam word is nonzero at two phases, but `A_5` is invertible.  There
are no low Hasse conditions.  This is the cheapest positive control against
the false statement “nonconstant slack holonomy prevents cyclic descent.”

### 7.2 Depth six: the first genuine coefficientwise seam hostile

At `k=6`, in the cut order `h=-6,-5,...,1`,

```text
p=8,       ell=1,       N=2,
R_6^raw=q+q^2+q^4+q^6,
Dtilde_6=q+q^2,
D_6=1+q^2,
Gamma_6=q^4+q^6,
G_6=1+q,
K_6=1+q+q^4+q^6.                                    (47)
```

The universal endpoint seam has

```text
boxed: M_0^r(W)=q^2+q^4 !=0.                         (48)
```

Thus the literal universal endpoint table is not coefficientwise in
`I_6=(X+1)V_8`, even though it scalarizes to the correct phase profile.
The cyclic current and the seam tomography identify the exact first failed
implication.

The raw-selected carrier has cut support `{0,6}` and gives

```text
Z_6^can=
(1+q^2, 0, q, 1+q, 1, q, q^2, 1+q).                (49)
```

The minimum-support Frobenius carrier has cut support `{0,2}` and gives

```text
Z_6^Fr=
(1+q^2, 0, 1, 1+q, 1, q, 1+q+q^2, 1+q).           (50)
```

Both tables have degree two, every coefficient profile has zeroth phase
moment zero, their marked column is `D_6`, and their scalar values in cut
order are

```text
(0,0,1,0,1,1,1,0).                                  (51)
```

Rotating back to phases `0,...,7` gives exactly THM-3492's

```text
F_6=(1,0,0,0,1,0,1,1).                              (52)
```

The two repaired tables are a physical-arc-image **parity current-chain**
hostile: even the contracted raw chain does not remove section dependence
unless a carrier policy is declared.  They are not asserted to be literal
nonnegative or raw event-distribution tables.

## 8. Connection and loss ledger

| object / map | preserves | destroys or changes | needed sidecar |
|---|---|---|---|
| collision bond + Green paths -> `E_k` | physical source, endpoint, exact slack parity | individual path order after `H` parity | source time and target phase |
| `E(h+1)+E(h)` -> `d` | exact terminal current at `q=1` | scalar owner potential | physical cut `h=-k` |
| universal endpoint arc -> cyclic arc | phase profile and marked no-wrap raw polynomial | literal off-cut endpoint slack | full holonomy word `Omega` |
| seam -> low Hasse tomography | complete obstruction modulo `I_k` | high image-internal holonomy modes | orders `0,...,ell-1` |
| raw time zero -> time two | scalar center and all later collisions | grades `k,k-2,k-3` | explicit `Gamma_k` chain |
| time-two full -> inward | monic inward face | three boundary strips | `G_k` |
| `C_N` | evaluation at `q=1`, phase image, marked time-two full axis | exact raw slack above `N` | pointed grade zero and quotient chain |
| top or Frobenius carrier repair | full phase axis, inward marked axis, coefficientwise arc image | uniqueness of mixed correlation | carrier policy |

This is the Rule 30 analogue of THM-3466's pointed boundary primitive: a
basepoint makes the current exact, while the clutch across the cut retains
the missing information.  The analogy supplies no factorial-conjecture or
Keller transfer.  `B_k=H_1(k-2)` enters the restart chain explicitly; it is
not a free scalar backbone after regrading.

THM-2538's common-ancestry lesson is realized here by actual spacetime bonds
and Green paths before contraction.  The construction is nevertheless over
`F_2`; it is not a nonnegative transportation coupling.  The LRC owner/current
analogy ends at the cut and holonomy ledger: there is no loneliness predicate,
chronological arrival theorem, or positive-Haar current.  Equation (23) is a
finite Radon/Pascal tomography statement only for the seam quotient.

There is also an exact scalar bridge to the dyadic section cut studied by the
provisional THM-3500 candidate.  If `p_i^(n)` is the length-`n` Mealy section
mask, then for `k>=2` and `0<=t<n`, direct reversal gives

```text
p_(k-1)^(n)(n-1-t)
 =b_(k-1)(t) OR b_(k-2)(t)=:Q_k^edge(t).             (53)
```

For `n=2s`, define the upper/lower section defect and its reversal by

```text
Delta_(k-1)^(2s)(i)
 =p_(k-1)^(2s)(i)+p_(k-1)^(2s)(s+i),
Delta_k^rev(t)=Delta_(k-1)^(2s)(s-1-t),
                                      0<=i,t<s.
```

Then

```text
Delta_k^rev(t)
 =Q_k^edge(t+s)+Q_k^edge(t)
 =H_k^(s)(t+1)+H_k^(s)(t),
H_k^(s)(t)=b_k(t+s)+b_k(t).                          (54)
```

In particular, telescoping (54) gives
`M_0(Delta_k^rev)=b_k(2s)` because `b_k(0)=0`.
Equation (53) is the factor-by-factor reversal of THM-3463's section compiler;
(54) then uses only the packed bit recurrence
`b_k(t+1)+b_k(t)=b_(k-1)(t) OR b_(k-2)(t)`.

Our phase current is the translate `Q_k(h)=Q_k^edge(k+h)`, while

```text
Omega_a(1)=b_k(p+a)+b_k(a)=H_k^(p)(a).               (55)
```

Thus the two constructions share the same translation-coboundary complex at
the unmarked level.  At `s=p=P_k`, (55) is the constant `epsilon_k`, so its
ordinary coboundary vanishes.  The polynomial word `Omega_a(q)` is a strictly
finer slack grading and is not recoverable from the Boolean section defect.
Equations (53)--(55) are a typed comparison only; THM-3500 is not a dependency
of any conclusion above while its own status remains provisional.

## 9. Exact companion and scope boundary

Run

```bash
python3 04-computation/rule30_universal_cover_slack_holonomy_thm3501.py
python3 -O 04-computation/rule30_universal_cover_slack_holonomy_thm3501.py
```

The transcripts are byte-identical.  The companion uses explicit `check`
gates and verifies:

1. independent centered-row and packed-`Phi` Rule 30 engines through the
   declared row universe;
2. independent bit-polynomial and centered-Green recurrence engines;
3. the Mealy-section/current reversal and Boolean bit-holonomy coboundary
   bridge (53)--(55) on the declared finite section universe;
4. every universal-cover potential used by all cut, endpoint, current, and
   holonomy tables for each hard unwrapped depth `5<=k<=16`;
5. the Duhamel scalarization, degree bound, monicity, current orientation,
   cyclic arc, universal endpoint, seam, and every low Hasse identity;
6. the raw/time-two restart, first-three-strip split, inward monicity, and
   exact normalization `K`;
7. contraction/current commutation, raw-selected and support-minimal
   Frobenius repairs, marked axes, phase axes, degree bounds, and physical
   image conditions; and
8. the complete depth-five positive and depth-six hostile controls
   (46)--(52), including both cut and standard phase orientations.

The proofs above are finite algebraic identities for arbitrary hard
unwrapped `k`; the finite universe audits the implementation and controls.
The independent audit rederived the universal identities and replayed the
ordinary and optimized companions byte-for-byte against the pinned output.

No Rule 30 center nonperiodicity, balance, density, unpredictability, or
random-access lower bound follows.  The refinement is parity-valued rather
than positive, and the carrier policy remains an explicit sidecar.  No Rule
30 prize and no literature novelty is claimed.
