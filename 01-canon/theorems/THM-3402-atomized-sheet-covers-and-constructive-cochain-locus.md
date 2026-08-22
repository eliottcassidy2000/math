---
id: THM-3402
title: "Atomized sheet covers and constructive cochain locus"
status: >
  PROVED analytic all-q atomized-cover/complete-cochain equivalence + PROVED
  constructive exact-locus reconstruction + PROVED base-star and finite-mode
  compressions + FINITE-EXACT controls + INDEPENDENTLY PROOF- AND
  EVENT-AUDITED.  Repeated atoms of one owner are retained until the cochain
  glues them, then compress exactly to one consecutive phase mode.  The q=8
  singleton/domino mode criterion is analytic; corrected full-clutter counts
  remain separately FINITE-EXACT.  Speeds divisible by q are core speeds in
  body-relative uses.  The theorem determines B_q(U), while
  B_q(U) subset A_C and B_q(U) minus A_C subset Gamma_D remain distinct
  sidecars.  No refined-ledger decrement, physical transport, or LRC(14)
  conclusion follows.
source: root-2608-crouzeix-puzzle-2026-08-15
audit: independent direct event geometry, repeated-owner tuples, q-divisible core typing, generalized CRT reconstruction, strict real-line Helly, base-star compression, constructive full-locus equality, q=6 hostiles, current MISTAKE-383 clutter scope, and ordinary/optimized replay
depends_on:
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
  - THM-3395-small-sheet-typed-cover-star-cochain
related:
  - THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks
  - THM-3388-three-sheet-phase-triangle-cover-clutter
  - THM-3389-four-sheet-typed-cover-clutter
  - THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight
script: 04-computation/lrc_atomized_sheet_cover_constructive_locus_thm3402.py
output: 05-knowledge/results/lrc_atomized_sheet_cover_constructive_locus_thm3402.out
script_sha256: b969212a5bfc918098a8e9366c14d6bbaa24729f90a5e8ea6d18e62f3e255bea
output_sha256: 20a1ca6aec5a74f1cdc10227025429a6924d54d7d042d9b108ed573faa31a017
semantic_sha256: 39aade82ac505a78e51d779e500acfc1d1b916c5350632311e9274b727e55aec
audit_script: 04-computation/lrc_atomized_sheet_cover_constructive_locus_independent_audit_thm3402.py
audit_output: 05-knowledge/results/lrc_atomized_sheet_cover_constructive_locus_independent_audit_thm3402.out
audit_script_sha256: e346798c975d6980d5488e9089956ffd3e55de8919d11042008c60737b8f7131
audit_output_sha256: 5435b9463795b6f39163876e266c8d9bb7c0e98cdc06d2c15cc316290ac54396
audit_semantic_sha256: bb394b88aa71ebe9a892c370a071227a3e016adb8d550163efc73dba70995ddc
hash_basis: LF-normalized bytes
---

# THM-3402 -- atomized sheet covers and constructive cochain locus

**PROVED analytic all-`q` atom/cochain and constructive-locus theorem +
FINITE-EXACT controls + INDEPENDENTLY PROOF- AND EVENT-AUDITED.**

## 1. Inheritance and connection contract

THM-3387 identifies the exact transverse obstruction as the full cyclic-sheet
cover locus `B_q(U)`.  THM-3398 proves the exact all-sheet consecutive-mode
compiler, while THM-3395 proves its one-coset face through seven sheets.  Those
theorems are load-bearing, not superseded: they supply the owner-compressed
target, the exact `q<=7` face, the strict `q=7` endpoint convention, and the
first pairwise-compatible but globally incoherent hostiles.

The new coordinate beyond seven sheets is **owner multiplicity**.  One owner
may fire several multiplication-kernel cosets at the same phase.  Quotienting
those cosets to one unlabelled owner before gluing loses both their individual
sheet footprints and the requirement that they come from the same tooth
window.  Retaining them as repeated-owner atoms restores an exact theorem for
every `q>=2`.

The inheritance pass is therefore:

- closest proved mechanism: THM-3395's typed complete cochain and compatible
  star;
- canonical hostile: at `q=6`, `(2,8,14)` has a covering coset partition and
  nonempty pair fibres but no closing cochain;
- sharp new boundary: at `q=7` two adjacent unit atoms only touch, while at
  `q=8` they overlap strictly;
- least-used sidecar: the grouping of several atoms by one owner, retained
  through the complete integral cochain;
- corrected near miss: independent blocker phases, scalar capacities, and
  pairwise mode overlaps do not reconstruct one common physical phase.

| field | all-sheet connection |
|---|---|
| source | the same located sheet set `X_q=Z/qZ` for every owner |
| target | owner-labelled kernel-coset atoms with one complete affine integer cochain |
| map | `(u,k)` retains the coset `k+ker(u:X_q->X_q)` and its tooth interval |
| preserved | source identity, owner identity, sheet footprint, strict overlap, common phase |
| lost by scalar capacity | coset placement and every overlap/correlation relation |
| lost by independent phases | the common parameter `t` |
| lost by owner quotient | which several atoms must occur in one owner window |
| required sidecars | complete cochain, core typing, and the body/deleted-grid inclusions |
| cheapest decisive tests | repeated-owner q=8 tuples and the q=6 `(2,8,14)` hostile |

The pair observable is symmetric integer compatibility, not an intrinsic
orientation.  This is not a tournament encoding.

## 2. Exact atomization and capacity

Fix `q>=2` and a positive speed `u`.  Put

```text
X_q=Z/qZ,                g_u=gcd(u,q),
m_u=q/g_u,               a_u=u/g_u.                    (1)
```

Multiplication by `u` on `X_q` has kernel

```text
K_u={0,m_u,2m_u,...,(g_u-1)m_u}.                       (2)
```

For `k in X_(m_u)`, define the owner-labelled atom

```text
Q_(u,k)={k+j m_u mod q:0<=j<g_u}.                      (3)
```

At first-sheet coordinate `t in T=R/Z`, the sheets blocked by `u` are

```text
E_u^*(t)={ell in X_q:||u(t+ell/q)||<1/14}.             (4)
```

Expression `(4)` is constant on every atom `(3)`, because increasing `ell`
by `m_u` changes the phase by the integer `u/g_u`.  Hence `E_u^*(t)` is the
disjoint union of its firing atoms.

The `m_u` image phases are equally spaced.  With

```text
c_u=ceil(m_u/7),                                       (5)
```

the exact number `s_u(t)` of firing atoms obeys

```text
s_u(t) in {c_u-1,c_u},          max_t s_u(t)=c_u,
|E_u^*(t)|=g_u s_u(t).                                 (6)
```

Both values occur.  When `7|m_u`, the lower value occurs at the strict
endpoint alignment; no closed-arc replacement is valid.  Thus THM-3391's
uniform capacity is recovered atom by atom:

```text
lambda(u)=g_u ceil(m_u/7).                             (7)
```

For `q<=7` and `q` not dividing `u`, one has `c_u=1`, so a firing owner
contributes exactly one atom.  This is precisely the load-bearing THM-3395
regime.  At `q=8`, a unit owner can contribute two atoms for the first time.

For a nonnegative weight on a located support `S subset X_q`, replace the
size of `(3)` by the weight of `S intersect Q_(u,k)`.  Full positive-weight
cover means that the chosen atom footprints contain the positive support;
the analytic gluing theorem below is unchanged.

## 3. Atomic families and admissible cochains

Let `U` be a finite set of positive owner speeds.  An **atomic family**
`A` is a finite set of distinct pairs

```text
alpha=(u_alpha,k_alpha),
u_alpha in U,          k_alpha in X_(m_(u_alpha)).      (8)
```

Different atoms may have the same owner.  A positive witness never needs
more than `c_u` atoms of owner `u`; the theorem itself forces this, while
`(5)` is the finite-search cutoff.

For oriented atoms `alpha,beta`, with speeds `u,v` and labels `k,l`, assign
an integer `p_(alpha,beta)` satisfying

```text
p_(beta,alpha)=-p_(alpha,beta),                        (9)

p_(alpha,beta)
 ==(l-k)u v                  (mod q gcd(u,v)),          (10)

14|p_(alpha,beta)|<q(u+v).                             (11)
```

For every three distinct atoms `alpha,beta,gamma`, of speeds `u,v,h`, require

```text
h p_(alpha,beta)
+u p_(beta,gamma)
+v p_(gamma,alpha)=0.                                  (12)
```

Call `(9)`--`(12)` an **admissible complete affine cochain**.

## 4. The all-sheet atomic cover theorem

For every `q>=2` and every finite positive owner set `U`,

```text
B_q(U) is nonempty
iff there is an atomic family A such that
    union_(alpha in A) Q_alpha=X_q
and A admits an admissible complete affine cochain.     (13)
```

For a retained support `S`, define its own cover locus

```text
B_(q,S)(U)={t in T:S subset union_(u in U)E_u^*(t)}.   (14a)
```

Then `(13)` holds with `B_(q,S)(U)` on the left and the footprint condition

```text
S subset union_(alpha in A)Q_alpha.                    (14)
```

on the right. Thus `B_q(U)=B_(q,X_q)(U)`. For a positive weighted support,
take `S` to be its positive-support set; weights affect capacity bounds, not
the pointwise cover predicate.

Minimal atomic covers and minimal owner covers are different objects when one
owner supplies several atoms.  An owner-level edge is obtained only after
grouping every selected atom by owner.

### 4.1 Necessity

Suppose the owners cover at a common first-sheet time `t`.  Select enough of
the firing atoms to cover `X_q`.  For each selected atom choose an integer
tooth `b_alpha` and lifted centre

```text
x_alpha=b_alpha/u_alpha-k_alpha/q                     (15)
```

whose open interval of radius `1/(14u_alpha)` contains `t`.  Put

```text
p_(alpha,beta)
 =q u_alpha u_beta(x_alpha-x_beta).                    (16)
```

Expanding gives

```text
p_(alpha,beta)
=q(u_beta b_alpha-u_alpha b_beta)
 +(k_beta-k_alpha)u_alpha u_beta.                      (17)
```

The ideal generated by `u_alpha,u_beta` proves `(10)`, common strict overlap
proves `(11)`, and telescoping the actual centre differences proves `(12)`.

### 4.2 Triangle closure gives rational potentials

Normalize

```text
delta_(alpha,beta)
=p_(alpha,beta)/(q u_alpha u_beta).                    (18)
```

Then `(12)` is exactly

```text
delta_(alpha,beta)+delta_(beta,gamma)
+delta_(gamma,alpha)=0.                                (19)
```

Fix a base atom `alpha_0`, set `z_(alpha_0)=0`, and put

```text
z_alpha=-delta_(alpha_0,alpha).                        (20)
```

Equation `(19)` gives

```text
z_alpha-z_beta=delta_(alpha,beta).                     (21)
```

Thus the complete rational cochain is a coboundary.  Pairwise gap existence
without `(19)` does not produce these vertex potentials.

### 4.3 The congruences give one generalized-CRT shift

We seek one rational `s` such that every translated potential is an actual
tooth centre:

```text
z_alpha+s in (1/u_alpha)Z-k_alpha/q.                   (22)
```

The two progressions for atoms `alpha,beta` are compatible exactly when

```text
(-k_alpha/q-z_alpha)-(-k_beta/q-z_beta)
 in (1/u_alpha)Z+(1/u_beta)Z
 =gcd(u_alpha,u_beta)/(u_alpha u_beta) Z.              (23)
```

Substitution of `(21)` turns `(23)` into congruence `(10)`.  Clearing every
speed and offset denominator changes `(22)` into ordinary integer
congruences; their pairwise compatibility is `(23)`.  The generalized CRT
therefore supplies one common `s`.  No distinct-owner or coprimality
assumption is used, so repeated-owner atoms are included literally.

### 4.4 Coherent real lifts and strict Helly

Put

```text
x_alpha=z_alpha+s.                                    (24)
```

These are coherent real lifts, not independently chosen circular
representatives.  By `(11)` and `(18)`,

```text
|x_alpha-x_beta|
 <1/(14u_alpha)+1/(14u_beta).                         (25)
```

Hence the open real intervals centred at `x_alpha` with the displayed radii
intersect pairwise.  If `L` is the largest left endpoint and `R` the smallest
right endpoint, the two intervals attaining `L,R` intersect, so `L<R`.
Their full intersection is nonempty.  Reducing a common point modulo one
fires every selected atom and proves sufficiency in `(13)`.

The argument occurs on the universal cover after CRT alignment.  There is no
circular-wrap shortcut and no total-arc-length hypothesis.

## 5. Exact base-star compression

The complete graph need not be searched edge by edge.  Fix the base atom
`alpha_0` of speed `u_0` and choose only

```text
p_(alpha_0,alpha) in P_(alpha_0,alpha)                 (26)
```

from the pair fibres `(10)`--`(11)`.  For two nonbase atoms `alpha,beta` of
speeds `u_alpha,u_beta`, triangle closure forces

```text
p_(alpha,beta)
= [u_alpha p_(alpha_0,beta)
   -u_beta p_(alpha_0,alpha)]/u_0.                    (27)
```

The star is compatible iff every numerator in `(27)` is divisible by `u_0`
and every resulting integer lies in its required pair fibre.  Those checks
imply all triangle laws, including triangles with no base vertex.  Thus

```text
admissible complete cochain  iff  compatible base star. (28)
```

Formula `(27)` remains valid for equal owner speeds attached to distinct
atoms.  This is an exact reduction in search dimension, not a relaxation.

## 6. Constructive reconstruction of the full cover locus

The theorem determines the entire open set `B_q(U)`, not only whether it is
empty.  Define the first-sheet cover locus

```text
C_q(U)={t in T:union_(u in U)E_u^*(t)=X_q}.            (29)
```

With THM-3387's base coordinate `y`,

```text
C_q(U)=(t |-> qt)^(-1)(B_q(U)),
B_q(U)={q t mod 1:t in C_q(U)}.                        (30)
```

The `q` preimages of one base point differ by `1/q`; they merely permute the
sheet labels, so the full-cover predicate is unchanged.

Fix an atomic cover `A` and an admissible cochain `p`.  Choose the potentials
`z_alpha` from `(20)` and put

```text
S(A,p)=intersection_alpha
  [(1/u_alpha)Z-k_alpha/q-z_alpha].                    (31)
```

The CRT proof shows that `(31)` is nonempty.  If

```text
g_A=gcd_(alpha in A)u_alpha,                           (32)
```

then the intersection of the underlying lattices is `(1/g_A)Z`; hence, for
one solution `s_0`,

```text
S(A,p)=s_0+(1/g_A)Z.                                  (33)
```

For `s in S(A,p)`, define the strict common interval

```text
J(A,p;s)=intersection_(alpha in A)
 (z_alpha+s-1/(14u_alpha),
  z_alpha+s+1/(14u_alpha)).                            (34)
```

It is nonempty by `(25)`, and common translation gives

```text
J(A,p;s_0+n/g_A)=J(A,p;s_0)+n/g_A.                    (35)
```

Therefore the exact first-sheet and THM-3387 base loci are

```text
C_q(U)
= union_(A atom-covers X_q)
  union_(p admissible on A)
  union_(0<=n<g_A)
  [J(A,p;s_0)+n/g_A] mod Z,                            (36)

B_q(U)={q t mod 1:t belongs to the union in (36)}.     (37)
```

Every interval on the right fires its selected cover, so it lies in the left.
Conversely every `t in C_q(U)` selects firing atoms and its actual centre
cochain, so it occurs in `(36)`.  This proves equality, including repeated
owners and strict endpoints.  The independent audit reconstructs `(36)` on
`7,741` atom families at `204,684` event and midpoint samples, with `1,059`
cochains and `4,200` repeated-owner families.

## 7. Equivalent consecutive-mode compression for all q

The atomic theorem refines, and compresses exactly back to, THM-3398's current
finite-mode theorem.  In the image phase grid of owner `u`, atom `k` has residue

```text
r=a_u k mod m_u.                                      (38)
```

At one phase the firing residues form a cyclic consecutive block of size
`c_u-1` or `c_u`. A selected mode certificate is any nonempty cyclic
consecutive subblock of those residues. Starting at `r` and containing `s`
residues, such a certificate has

```text
1<=s<=ceil(m_u/7).                                    (39)
```

Starting with an admissible atom cochain, Section 4 first supplies an actual
common `t` in the strict interval `J`. For each owner, replace its selected
atoms by the **entire atom set actually firing at that same `t`**. In the image
coordinate `(38)` this is one cyclic consecutive block of size `c_u-1` or
`c_u`; adding genuinely firing atoms preserves the common phase and can only
enlarge the footprint. Thus necessity may enlarge to the full actual firing
block. Conversely, sufficiency needs only a selected consecutive subblock
satisfying `(39)`; its strict source interval may force extra unselected atoms
to fire, which again only enlarges the footprint. Every selected mode expands
back into owner-labelled atoms. No shortest-cyclic-hull convention is used.

For a block `r,r+1,...,r+s-1`, set

```text
h=-g_u(2r+s-1) mod 2q,
w=g_u[m_u-7(s-1)].                                    (40)
```

Its exact common source-time interval has centre lattice and radius

```text
H_(u,r,s)=(1/u)(Z+h/(2q)),
rho_(u,r,s)=w/(14 q u)
            =[1/14-(s-1)/(2m_u)]/u.                   (41)
```

The integer `w` is positive precisely in the lawful range `(39)`.  For two
modes `i,j`, with speeds `u_i,u_j`, define

```text
P_ij=2q u_i u_j(x_i-x_j).                             (42)
```

Their exact pair conditions are

```text
P_ij == h_i u_j-h_j u_i       (mod 2q gcd(u_i,u_j)),
7|P_ij|<w_i u_j+w_j u_i.                              (43)
```

The triangle law is

```text
u_k P_ij+u_i P_jk+u_j P_ki=0.                         (44)
```

The same coboundary, generalized-CRT, and real-line Helly proof shows:

```text
B_q(U) is nonempty
iff one mode per selected owner has footprints covering X_q
and the modes admit (43)--(44), equivalently a compatible star. (45)
```

Thus atomization is the lossless fine representation, while modes are the
exact owner-compressed representation.  The hard finite step is mode-cover
enumeration; no new analytic gluing principle appears at larger `q`.

## 8. The analytic q=8 singleton/domino corollary

At `q=8`, `(40)`--`(41)` reduce exactly to four mode species:

```text
gcd(u,8)=1:  8 singleton modes (w=8)
             + 8 adjacent-image domino modes (w=1);
gcd(u,8)=2:  4 antipodal-pair modes (w=8);
gcd(u,8)=4:  2 parity-quadruple modes (w=8);
8|u:         one all-sheet core mode.                  (46)
```

For an odd owner, “adjacent” means adjacent after multiplication by `u` in
the phase grid; the two sheet labels need not be numerically consecutive.
Singleton source radius is `1/(14u)` and domino source radius is `1/(112u)`.
Writing the mode centre as `(Z+h/16)/u`, equations `(42)`--`(44)` become

```text
P_ij=16u_i u_j(x_i-x_j),
P_ij == h_i u_j-h_j u_i       (mod 16 gcd(u_i,u_j)),
7|P_ij|<w_i u_j+w_j u_i,                              (47)
```

together with zero triangle circulation.  Covering all eight sheets plus
`(47)` is an exact analytic iff.  It proves the criterion, not a literal
finite-clutter enumeration.

This distinction is load-bearing after MISTAKE-383.  The current independent
full-subset audit gives the true `q=8` minimal-edge rank profile
`(4:15,5:17,6:6)`; the older number `32` was only the rank-at-most-five part.
Likewise the current `q=15` full profile is `(6:157,7:16,8:6)`.  Nothing in
the atom/cochain proof depends on either finite profile.  The theorem explains
why a correctly enumerated mode bank is exact, but never promotes a bounded
census into an analytic clutter classification.

## 9. Small-sheet and q=6 faces

When `q<=7` and `q` divides no transverse owner, every selected owner has one
atom.  Equations `(10)`--`(12)` and the base-star formula `(27)` are exactly
THM-3395.  Its q=2--7 literal classifications and body counts remain sourced
there; THM-3402 does not rebrand them as new computations.

At `q=6`, the CRT chart

```text
X_6 -> X_2 x X_3                                    (49)
```

types atoms as cells, columns, or rows according as `gcd(u,6)` is `1,2,3`.
There are `23` irredundant abstract covers: `11` capacity-equality partitions,
six weight-seven overlap covers, and six weight-eight overlap covers.  The
hostile `(2,8,14)` has a covering column partition and pair gaps but no
cochain; `(2,8,10)` is the matched positive.  This separates scalar capacity,
Boolean cover, and common-phase correlation.

## 10. q-divisible core typing

The analytic theorem permits `q|u`.  Then

```text
m_u=1,             Q_(u,0)=X_q,                       (50)
```

so the owner has sheet-independent status.  This is mathematically correct
inside `(13)` and is independently tested, including repeated q-divisible
owners.

In THM-3387's body decomposition, however, every such speed is a **core**
speed and must be routed to `A_C`; it must not remain in the transverse set
`U`.  Every body-relative invocation of THM-3402 therefore assumes

```text
q does not divide u for all u in U.                    (51)
```

Otherwise a full atom in the transverse clutter merely disguises core danger
and destroys the body typing.

## 11. Body-relative and deleted-grid sidecars

THM-3402 reconstructs `B_q(U)` by `(36)`--`(37)`.  It does not silently turn
transverse cover into a body-relative obstruction.  For

```text
F=qC disjoint union U,
A_C=union_(c in C)D_c,
Gamma_D=D^(-1)Z/Z,                                    (52)
```

THM-3387 remains the exact composition law:

```text
pointwise safe-image equality
 iff B_q(U) subset A_C,                                (53)

aligned open-cell completion
 iff B_q(U) minus A_C subset Gamma_D.                  (54)
```

These predicates can differ at strict endpoint handoffs.  A valid atomic or
mode cover can be harmless because its entire locus lies in core danger, and
an endpoint-only residual can violate `(53)` while preserving `(54)`.  The
constructive formula supplies the exact `B_q(U)` needed to test the sidecars;
it does not remove them.

## 12. Verification, provenance, and scope

The primary standard-library verifier declares finite exact universes through
`q=10`.  It checks `1,116` multiplier/phase cells, `407,566` atom tuples,
`225,883` repeated-owner tuples, and `135` owner-cover problems.  There are
`53,248` positive atom tuples and `6,469` positive repeated-owner tuples.  It
also freezes the strict `q=7/q=8/q=14` walls and repeated-owner full-cover
controls at `q=8,9,10`.

The independently written audit uses a separate event sweep and cochain
enumerator.  It checks `3,318` capacity phases, `11,704` pullbacks, `50,063`
atom tuples, `28,118` repeated-owner tuples, `957` core-typed tuples plus
`1,914` extra q-divisible controls, `105` owner covers, base-star compression,
the constructive locus, and the q=6 cover/hostile packet.  Both implementations
use exact integers and `Fraction`, with no floating literals or
optimization-dependent assertions.

Reproduce after installing the canonical paths with

```text
python3 04-computation/lrc_atomized_sheet_cover_constructive_locus_thm3402.py
python3 -O 04-computation/lrc_atomized_sheet_cover_constructive_locus_thm3402.py
python3 04-computation/lrc_atomized_sheet_cover_constructive_locus_independent_audit_thm3402.py
python3 -O 04-computation/lrc_atomized_sheet_cover_constructive_locus_independent_audit_thm3402.py
```

Ordinary and optimized outputs must LF-normalized-byte-match their stored
outputs.  The unbounded statements rest on Sections 4--8, not on finite
extrapolation.

This theorem proves an exact structural compiler for cyclic-sheet cover loci.
It does not decrement the refined ledger, turn a projected screen into a
physical cover, transport a reflected phase, close an arbitrary projected
sector, or prove LRC(14).
