---
id: THM-3395
title: "Small-sheet typed covers are exact star cochains"
status: >
  PROVED analytic q<=7 typed-coset cover criterion + FINITE-EXACT literal
  body-relevant q=2,...,7 audit + INDEPENDENTLY EVENT- AND PROOF-AUDITED.
  Below eight sheets, every
  firing transverse speed blocks one kernel coset.  A selected family covers
  iff those typed cosets cover Z/qZ and its affine integer gap fibres contain
  one complete closed cochain, equivalently one compatible star.  In the
  literal q=5 pool the clutter has 231 rank-five edges and independence
  profile (1,12,66,220,495,561,268,45,1,0,...); at q=6 it has 39 edges of
  ranks (3:3,4:29,5:7) and profile
  (1,12,66,217,441,515,304,76,5,0,...).  These give 1,617 and 1,471 globally
  transverse-safe six-body rows; the two and seven already-known core rescues
  recover THM-3387's exact totals 1,619 and 1,478.  No refined-ledger
  decrement or LRC(14) follows.
source: root-2608-compatibility-obstruction-2026-08-14
audit: exact integer construction of every relevant q=2,...,7 literal subset, independent THM-3387 event comparison on 6540 subsets, independent q=6 convergence artifact, strict q=7/q=8 boundary, q=5/q=6 closure hostiles and positives
depends_on:
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
related:
  - THM-3383-terminal-monomial-cone-polynomiality-fork
  - THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks
  - THM-3388-three-sheet-phase-triangle-cover-clutter
  - THM-3389-four-sheet-typed-cover-clutter
external:
  - Artem Chernikov, "SOP2 = SOP3", arXiv:2608.13291v1 (CITED PREPRINT MOTIVATION ONLY)
script: 04-computation/lrc14_small_sheet_typed_cover_star_cochain_thm3395.py
output: 05-knowledge/results/lrc14_small_sheet_typed_cover_star_cochain_thm3395.out
script_sha256: 0953401d98d62fd3118bd4a7bbeb50bd43a459a04dc120578ca6af355cada067
output_sha256: f7ed05e16fdd3660741aa8a79600cf9920bbebd8087c8d25a252ecca0dbc1ce5
semantic_sha256: a47bf86457583259aa3712bfd5ca849328a7fa7b2db8e194345a910abb6610c1
hash_basis: LF-normalized bytes
---

# THM-3395 -- small-sheet typed covers are exact star cochains

**PROVED analytic `q<=7` typed-coset cover criterion + FINITE-EXACT literal
`q=2,...,7` audit + INDEPENDENTLY EVENT- AND PROOF-AUDITED.**

## 1. Inheritance and connection contract

THM-3387 proves that the exact transverse obstruction is a full cover of the
cyclic sheet fibre.  THM-3388 and THM-3389 solve its three- and four-sheet
slices by retaining phase gaps.  The present theorem supplies one uniform
criterion for every small sheet degree needed by the literal six-body atlas,
including the previously unproved `q=5,6,7` structural slices.

| field | connection |
|---|---|
| source | distinct transverse speeds on the `q` sheets over one base phase |
| target | a typed set-cover clutter on `Z/qZ` decorated by an affine integral complete `1`-cochain |
| map | retain each owner speed, its blocked kernel coset, and its cleared signed centre gaps |
| preserved | simultaneous sheet ownership, strict tooth overlap, common source phase, cyclic labels |
| lost by the pair shadow | which pair witnesses use the same teeth, gap sign and magnitude, triangle closure |
| required sidecar | one blocked-coset label per owner and one compatible star of integer gaps |
| cheapest hostile | `q=6`, speeds `(2,8,14)`, labels `(0,1,2)`: every pair fibre is nonempty but the star fibre is empty |

This is a theorem about the THM-3387 sheet-cover obstruction.  It is not a
tournament encoding: the pair observable is symmetric, blockers can own more
than one sheet, and compatibility is carried by integer values rather than
edge orientations.

## 2. One firing speed is one kernel coset

Write `T=R/Z`, put

```text
D_u={x in T: ||ux||<1/14},
E_u(y)={ell in Z/qZ:(y+ell)/q in D_u}.                 (1)
```

Let `2<=q<=7`, let `q` not divide `u`, and set

```text
g_u=gcd(u,q),                       m_u=q/g_u.          (2)
```

As `ell` runs through `Z/qZ`, the phases

```text
u(y+ell)/q mod 1                                      (3)
```

take `m_u` equally spaced values, each with multiplicity `g_u`.  Their
spacing is `1/m_u>=1/7`, while the strict danger arc has length `1/7`.
Consequently `(1)` is either empty or exactly one coset

```text
K_u(k)={ell in Z/qZ:u(ell-k)=0 mod q}                  (4)
```

of the kernel of multiplication by `u`; it has size `g_u`.  Equality at
`m_u=7` causes no double hit because the danger arc is open.

This is the first sharp small-sheet boundary.  At `q=7,u=1,y=1/2`, sheets
zero and six lie exactly at `+/-1/14`; the strict hit count is zero although
the closed hit count would be two.  At `q=8,u=1`, the two sheets centred at
`-1/16` and `1/16` are both strictly dangerous.  Thus the one-coset typing
fails immediately at eight sheets.

## 3. Affine gap fibres

Take distinct transverse owners `u_0,...,u_(r-1)` and choose a representative
`k_i` for one kernel coset `(4)` of each owner.  A tooth blocking that coset
has a lifted source-time centre

```text
x_i=a_i/u_i-k_i/q,                         a_i in Z,    (5)
```

and radius `1/(14u_i)`.  For an oriented pair define

```text
p_ij=q u_i u_j(x_i-x_j)
    =q(u_j a_i-u_i a_j)+(k_j-k_i)u_i u_j.              (6)
```

Accordingly the exact finite pair-gap fibre is

```text
P_ij(k_i,k_j)={p in Z:
  p == (k_j-k_i)u_i u_j (mod q gcd(u_i,u_j)),
  14|p|<q(u_i+u_j)}.                                  (7)
```

The congruence in `(7)` is precisely the solvability condition for the two
tooth centres in `(5)`.  The inequality is precisely strict overlap of their
two source-time intervals.  Neither condition may be replaced by a floating
distance or a non-strict endpoint test.

Call antisymmetric integers `p_ij=-p_ji` an **admissible complete cochain**
when `p_ij` lies in `(7)` for every `i<j` and

```text
u_h p_ij+u_i p_jh+u_j p_hi=0                           (8)
```

for every triple of distinct owners.  Equivalently, the normalized gaps

```text
delta_ij=p_ij/(q u_i u_j)                              (9)
```

form an exact rational `1`-cochain: `delta_ij=z_i-z_j` for rational vertex
potentials `z_i`.

Only a star must be searched.  With owner zero as anchor, choose
`p_0i in P_0i`.  Then every remaining value is forced by

```text
p_ij=(u_i p_0j-u_j p_0i)/u_0.                         (10)
```

The star is compatible exactly when every value in `(10)` is integral and
lies in its required fibre `(7)`.  Equations `(8)` then hold automatically.
Thus “complete cochain” and “compatible star” are equivalent descriptions,
not different sufficient tests.

## 4. Exact typed-cover theorem

Let `U` be a finite set of at most six distinct positive integers, with `q`
dividing no member.  Then

```text
B_q(U) is nonempty
iff there are owners I subset U, labels (k_i), and an admissible cochain
    such that union_(i in I) K_(u_i)(k_i)=Z/qZ.         (11)
```

Here `B_q(U)` is THM-3387's full transverse-cover locus.  The statement is
literal: the same owner supplies every sheet in its kernel coset, and the
same selected tooth participates in all of its pair gaps.

### Necessity

Suppose `y` lies in `B_q(U)` and put `t=y/q`.  Keep the speeds that fire at
this base phase.  Section 2 says each has one blocked coset, and those cosets
cover `Z/qZ`.  Choose its representative `k_i` and its dangerous tooth
`a_i`; then the intervals centred at `(5)` all contain `t`.  Formula `(6)`
gives the congruence and strict bound `(7)`.  The actual centre differences
telescope, and multiplying

```text
(x_i-x_j)+(x_j-x_h)+(x_h-x_i)=0
```

by `q u_i u_j u_h` gives `(8)`.

### Sufficiency and the CRT gluing step

Conversely suppose the right side of `(11)` holds.  Triangle closure gives
rational potentials `z_i` satisfying `(9)`; explicitly take `z_0=0` and
`z_i=-delta_0i`, then use `(8)`.  Put

```text
L_i=(1/u_i)Z-k_i/q.                                   (12)
```

The congruence in `(7)` is equivalent to

```text
delta_ij in L_i-L_j.                                  (13)
```

Hence the rational lattice cosets `L_i-z_i`, which are the allowed common
shifts, meet pairwise.  Choose an integer `M` divisible by all `u_i` and by
the denominators of their offsets.  Multiplication by `M` turns these into
ordinary integer congruences, with moduli `M/u_i`.  Condition `(13)` is their
pairwise generalized-CRT compatibility.  The generalized CRT therefore
gives one common shift `s`, so

```text
x_i=z_i+s in L_i                                      (14)
```

simultaneously for every owner.  This is the missing global-witness step;
nonempty pair fibres alone do not give `(14)`.

The strict bounds in `(7)` say that the selected danger intervals around
the centres `(14)` meet pairwise.  Ordinary one-dimensional Helly now gives
a common lifted source time.  Intrinsically on the circle, the same step has
an important six-owner guardrail: a pairwise-intersecting circular-arc family
with empty total intersection covers the circle, while these arcs have total
length

```text
sum_i 1/(7u_i)<=r/7<=6/7<1.                           (15)
```

Thus no wrapped circular false positive is possible.  Their common source
time blocks every coset in the cover, proving `B_q(U)` nonempty.

The two bounds are deliberately visible: `q<=7` makes every firing blocker a
single kernel coset, and the six-owner LRC scope makes `(15)` strict.  The
`q=8` two-sheet hostile above forbids extrapolating the typed statement.

Finally, every inclusion-minimal coset cover has at most `q` owners: give
each indispensable owner a private sheet.  Thus the theorem reduces the
literal `q<=6` clutters to a finite search through ranks at most `q`.

## 5. Exact `q=6` quotient faces

The criterion contains the proved `q=2` and `q=3` laws as exact faces of the
six-sheet object, not as analogies.

If `u_i=2a_i` with `3` not dividing `a_i`, its block is a two-sheet fibre of
`Z/6Z -> Z/3Z`.  Substitution

```text
u_i=2a_i,                    p_ij=4r_ij                 (16)
```

turns `(7)`--`(8)` exactly into the `q=3` gap congruence, strict bound, and
closure equation of THM-3388.  In particular

```text
(2,8,14)=2(1,4,7)                                      (17)
```

with cosets `{0,3},{1,4},{2,5}` is the pullback of the ternary hostile:

```text
P_01={4},               P_02={-4},              P_12={-8,4}, (18)
```

yet `(10)` has no compatible star.  Replacing `14` by `10` gives the positive
cochain `(p_01,p_02,p_12)=(4,4,-4)`.

Similarly, if `u_i=3a_i` with `a_i` odd, then `p_ij=9r_ij` turns the size-three
block face into THM-3387's `q=2` gcd graph.  The mixed `q=6` clutter is the
genuine compatibility fibre product of these binary and ternary faces; it is
not determined by either quotient separately.

## 6. Literal `q=5,6,7` classification

For `q=5`, every transverse speed in the literal pool

```text
V_5={1,2,3,4,6,7,8,9,11,12,13,14}                     (19)
```

owns a singleton sheet.  The exact clutter has `231` inclusion-minimal
rank-five edges and independence profile

```text
(I_0,...,I_12)=(1,12,66,220,495,561,268,45,1,0,0,0,0). (20)
```

Pair feasibility is far from sufficient.  Assign labels `(0,2,3,1,4)` to
speeds `(1,2,3,4,9)`: all ten pair fibres are nonempty, but no compatible
star exists.  With the same labels, `(1,2,3,4,7)` is positive; one complete
cochain in lexicographic pair order is

```text
(-1,-1,-1,-2, 1,2,3, 1,1,-1).                        (21)
```

For `q=6`, the literal pool is

```text
V_6={1,2,3,4,5,7,8,9,10,11,13,14}.                    (22)
```

Its exact clutter has `39` inclusion-minimal edges, with rank profile and
independence profile

```text
(3:3,4:29,5:7),
(I_0,...,I_12)=(1,12,66,217,441,515,304,76,5,0,0,0,0). (23)
```

For `q=7` all literal transverse blockers are singletons.  A body row with a
nonempty core has at most five transverse speeds, so it cannot cover seven
sheets.  The relevant independent counts are therefore simply

```text
(I_0,...,I_5)=(1,12,66,220,495,792).                  (24)
```

The two core clocks for `q=5,6,7` give the globally transverse-safe counts

```text
q=5: 2I_5+I_4=2*561+495=1617,
q=6: 2I_5+I_4=2*515+441=1471,
q=7: 2I_5+I_4=2*792+495=2079.                         (25)
```

THM-3387's independent event atlas has respectively two, seven, and zero
core rescues.  Adding only those already-proved rescues gives exact row totals

```text
q=5:1619,                    q=6:1478,                 q=7:2079. (26)
```

No number in `(25)` or `(26)` is an additive refined-ledger decrement.

## 7. Typed-witness provenance and the JC boundary

Artem Chernikov's very recent preprint *SOP2 = SOP3*,
[arXiv:2608.13291v1](https://arxiv.org/abs/2608.13291) (submitted
2026-08-13), is **CITED motivation only**.  Its Lemma 2.3 keeps vertices as
witness-parameter pairs `(v_i,p_i)` and uses emptiness of shared-witness
fibres when constructing a directed three-cycle obstruction.  That typing
suggested the faithful LRC vertex `(u_i,k_i,p_0i)`: a speed, its coset, and its
star gap must remain attached.  The proof above is elementary arithmetic and
does not interpret a first-order theory, establish SOP for an LRC structure,
or use the preprint as a dependency.  “SOP2/SOP3” is not a synonym for the
cochain criterion.

THM-3383 supplies the adjacent guardrail.  There, exact rational torsor
decoding at `|g-ae|=1` still loses one polynomial target because the source-
ring intersection and boundary divisibility have been forgotten.  Here,
pair-gap existence likewise becomes a genuine source-time witness only after
the coset labels, closed star cochain, generalized CRT, and strict interval
geometry are restored.  The shared meta-pattern is **existential fibres must
be glued before quotienting**; no LRC-to-JC map or JC consequence is claimed.

## 8. Exact verification and scope

The standard-library companion constructs the typed clutter without calling
the event sweep and then compares it to THM-3387's independent rational-event
implementation on every literal relevant subset:

```text
q=2,...,7,                 0<=|U|<=5,
6540 subsets,              zero discrepancies.                         (27)
```

For `q<=6` it searches through the full possible minimal-cover rank `q` and
freezes the complete profiles; for `q=7` it audits exactly the ranks relevant
to a six-body row with nonempty core.  The `q=5` and `q=6` edge and profile
digests, both closure hostiles and positives, the independent `q=6`
convergence artifact, and the strict `q=7/q=8` boundary are all pinned in the
stored output.  There are no floating literals and no optimization-sensitive
assertions.

Reproduce with

```text
python 04-computation/lrc14_small_sheet_typed_cover_star_cochain_thm3395.py
python -O 04-computation/lrc14_small_sheet_typed_cover_star_cochain_thm3395.py
```

Ordinary and optimized outputs LF-normalized-byte-match the stored output.
This theorem closes the structural `q<=7` typed-cover classification needed
by the literal THM-3387 body atlas.  It does not close any new projected row,
transport a reflected phase to a physical tail, decrement the refined ledger,
or prove LRC(14).
