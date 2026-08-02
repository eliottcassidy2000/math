---
id: THM-3149
title: "Depth-three selector persistence and cross-support wall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  At support (1,3), bank I2, no probability law on all 134 physically
  multiplicity-valid pole prefixes of length at most three is Hasse-positive
  in every degree 5 through 11.  The four THM-3144 upset facets still suffice,
  with a new primitive positive integer combination strictly negative on all
  134 states.  The old depth-two combination has exactly one depth-three
  escape, (1,1,2), so the new certificate is a genuine repair.  At support
  (1,2), bank I2, two upset facets strictly separate all 42 depth-at-most-three
  states.  This is a finite averaged virtual-prefix theorem, not a sequential
  stopping process or an original-response decomposition.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
audit: >
  An independent immutable audit reconstructed both physical pole multisets
  and state censuses, rechecked the four-row 134-state separator and the
  two-row 42-state cross-support separator, and recovered every declared
  range and equality state.  It separately verified all six length-singleton
  descriptions, the Farkas implication, and the distinction between the full
  probability simplex and realizable sequential stopping laws.  Fresh normal
  and optimized runs byte-match the stored 1,309-byte transcript and both
  declared LF-normalized hashes.
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3144-mixed-depth-selector-persistence-death-barcode
related:
  - THM-3136-one-sided-fixed-reference-elementary-tail-hasse-no-go
  - THM-3147-length-singleton-endpoint-jet-facet-observer
script: 04-computation/gmc_depth_three_selector_persistence_wall_scout.py
output: 05-knowledge/results/gmc_depth_three_selector_persistence_wall_scout.out
script_sha256: 51d7b86393903b3498fe0178b95351fbd8bd5fe826ba7228c951669603b79c57
output_sha256: 02761eb027c8340b29c46f86a5964d39f741427699b6da93cbef6d7bee02efb5
hash_basis: LF-normalized bytes
---

# THM-3149 -- depth-three selector persistence and cross-support wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3144 proves a genuine selector death barcode on the 41 physical pole
prefixes of length one or two: a common positive selector survives through
degree 10, and the degree-11 fibre is separately nonempty, but no one selector
persists through all degrees 5 through 11.  Allowing every physical prefix of
length three does not restore a common section.  The obstruction is more
rigid than a generic larger-state Farkas wall: the same four upset facets
suffice after their coefficients are changed.

The same experiment at the earlier support `(1,2)` portability wall is even
smaller.  Two upset facets exclude all physical prefixes through depth three.

## 1. Physical states and response rows

Fix invariant `I=1` and bank `I2`.  For support `(a,b)`, let `P(a,b)` be the
reduced common-denominator pole multiset from THM-3120.  A physical state of
depth `d` is an unordered size-`d` submultiset of `P(a,b)`, respecting its
multiplicities.  If the state is `sigma`, subtract the corresponding virtual
alphabet from every signed bank atom and from the distinguished residual
alphabet `Q`.  As in THM-3137 and THM-3144, write

```text
G_N^sigma(mu)
 =Phi^sigma(h_N)m_mu[Q^sigma]
  -Phi^sigma(m_mu)h_N[Q^sigma].                            (1)
```

Every current `(1)` has total mass zero.  For a probability law `lambda` on
a finite state bank `S`, put

```text
G_N(lambda)=sum_(sigma in S) lambda_sigma G_N^sigma.        (2)
```

For a coarsening upset `U` in the partition lattice `P_N`, define the response
row

```text
R_(N,U)(sigma)=sum_(mu in U)G_N^sigma(mu).                  (3)
```

THM-3127 says that `(2)` is a nonnegative fine-to-coarse Hasse boundary only
if

```text
sum_sigma lambda_sigma R_(N,U)(sigma)>=0                   (4)
```

for every `U`.  Consequently, if finitely many such rows have positive
coefficients `c_i` and

```text
H(sigma)=sum_i c_i R_i(sigma)<0             for every sigma in S, (5)
```

then no probability law on `S` can satisfy all the corresponding Hasse
conditions: averaging `(5)` against `lambda` contradicts the positive
combination of `(4)`.  Both results below are explicit instances of `(5)`.

## 2. The 134-state support `(1,3)` wall

At support `(1,3)`, the reduced pole multiset is

```text
P(1,3)=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                (6)
```

The multiplicity-valid state counts are

```text
depth 1: 8,            depth 2: 33,            depth 3: 93,
total: 134.                                                   (7)
```

Use exactly the four upsets from the THM-3144 persistence certificate:

```text
U_8       =P_8 \ {(1^8)},
U_10      =P_10 \ {(1^10)},
U_11^(3)  =P_11 \ {(1^11),(2,1^9),(2,2,1^7)},
U_11^(1)  =P_11 \ {(1^11)}.                               (8)
```

Write their response rows in the displayed order as

```text
R_8, R_10, R_11^(3), R_11^(1).                            (9)
```

The old THM-3144 coefficient vector

```text
(461295,3948,1,22)                                         (10)
```

is strictly negative on all 41 depth-at-most-two states and on 92 of the 93
depth-three states, but it has the unique positive coordinate

```text
H_old(1,1,2)=744310486097342352.                           (11)
```

Thus simply citing the old certificate after enlarging the state bank would
be false.  The escape is real and uniquely localized.

Now take the primitive positive integer vector

```text
c_8  =3627302205077683560821210804482107686117129010,
c_10 =  29864988664447894013528399014705164256913403,
c_3  =      7514363586001882980386000831487329663540,
c_1  =    314610406969924544691100890629234040301424.      (12)
```

The exact companion evaluates

```text
H_3=c_8 R_8+c_10 R_10+c_3 R_11^(3)+c_1 R_11^(1)           (13)
```

on all 134 states.  Every coordinate is strictly negative, with exact range

```text
-674800650195958551767268010800617689469061950339284516178880
 <= H_3(sigma) <=
-1166613115072139157133415139738043609298077113840435285120. (14)
```

The upper equality set is exactly

```text
(1), (7), (1,1,2), (1,2,2),                               (15)
```

and the lower equality occurs only at `(5,5,6)`.  Applying `(5)` proves:

```text
there is no law on all physical depth-<=3 states at support (1,3), I2,
whose averaged current is Hasse-positive for every N=5,...,11.             (16)
```

The large coefficients in `(12)` are not a numerical approximation.  A
floating LP was used only to discover the active four-row support.  Equating
the four coordinates in `(15)` leaves a one-dimensional rational nullspace;
clearing denominators and dividing by the common gcd gives `(12)`.  The proof
is the subsequent exact 134-coordinate verification of `(13)--(14)`.

## 3. The 42-state support `(1,2)` wall

At support `(1,2)`, the reduced pole multiset is

```text
P(1,2)=(5,4,3,3,2,2,2,1,1,1,1).                          (17)
```

The depth-one, depth-two, and depth-three state counts are respectively

```text
5, 13, 24,                         total 42.                (18)
```

Only two upset rows are needed:

```text
V_6=P_6 \ {(1^6)},
V_9=P_9 \ {(1^9),(2,1^7),(3,1^6)}.                        (19)
```

For their response rows `S_6,S_9`, the primitive positive combination

```text
14903617223 S_6+641337 S_9                                (20)
```

is strictly negative on all 42 states.  Its exact range is

```text
-108248637454886456 <= (20) <= -2009645062627112.          (21)
```

The upper equality set is exactly `(1),(3,4,5)`, and the lower equality is
unique at `(2,3,4)`.  Therefore no probability law on all physical
depth-at-most-three states at support `(1,2)`, bank `I2`, is Hasse-positive
simultaneously in degrees 5 through 11.

This strengthens the one-pole portability wall in THM-3137: neither mixing
stopping depths nor adding every physical triple repairs that support.

## 4. Length-singleton compression of every active facet

The six rows used above look like partition-lattice data, but none requires
the full partition coordinate system.  In the notation of THM-3147, with
`ell=ell(mu)` and `m_1=m_1(mu)`, their upsets are exactly

```text
U_8       ={ell<=7},
U_10      ={ell<=9},
U_11^(1)  ={ell<=10},
U_11^(3)  ={ell<=8} union {ell=9 and m_1=8},

V_6       ={ell<=5},
V_9       ={ell<=6} union {ell=7 and m_1<=5}.              (22)
```

Indeed, the two degree-11 partitions of length nine are `(3,1^8)` and
`(2,2,1^7)`, with singleton counts eight and seven; the two degree-nine
partitions of length seven are `(3,1^6)` and `(2,2,1^5)`, with singleton
counts six and five.  All longer shapes are precisely the explicitly
excluded chains in `(8)` and `(19)`.

Therefore every coordinate in both Farkas certificates factors through the
THM-3147 length-singleton profile `Z_X(t,u,w)`.  Pure length detects four of
the six rows; one singleton coordinate distinguishes the remaining two.
This is a proved compression of the *dual obstruction*, not a claim that the
length-singleton profile determines the whole Hasse cone.  It identifies a
much smaller possible route to a depth-uniform recurrence for the dual
coefficients.

## 5. Exact verification

Run

```text
python 04-computation/gmc_depth_three_selector_persistence_wall_scout.py
python -O 04-computation/gmc_depth_three_selector_persistence_wall_scout.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_depth_three_selector_persistence_wall_scout.out.
```

The companion uses only integer and rational arithmetic.  It independently
reconstructs the two pole multisets; enumerates all `134` and `42` legal
states; reconstructs the six displayed upset rows from monomial symmetric
functions; verifies primitivity and positivity of both coefficient vectors;
checks every combined coordinate; and identifies both exact ranges and all
equality states.  The immutable LF-normalized hashes are recorded in the
frontmatter.

## 6. Scope and stopping boundary

The theorem concerns probability averages of the derived fixed-`Q` virtual
prefix currents `(1)`.  A state is physically multiplicity-valid, but an
arbitrary probability law on the state bank is not thereby the stopping law
of a sequential pole-removal process.  Proving nonexistence on this larger
simplex is nevertheless a valid obstruction to every such more structured
law.

The result does not prove:

- nonexistence after prefixes of length four or larger;
- a common selector obstruction at every support or bank;
- positivity or a decomposition of the original product-Gamma response;
- a Markov, martingale, or stopping-time realization of any selector;
- arbitrary-radial NC2 or the Gaussian Moment Conjecture.

The preserved structure is unusually small: from depth two to depth three,
the active upset support `(8)` is unchanged.  What changes is the positive
dual section over those facets.  The missing next sidecar is either a
depth-uniform dual certificate or a lawful compactness/recurrence mechanism
controlling how these dual coefficients evolve with prefix depth.

QED.
