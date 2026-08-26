---
id: THM-4169
title: "Prime-parent one-vertex augmentation and quartic Johnson transfer"
status: >
  PROVED ELEMENTARY PRIME-AUGMENTATION CLASSIFICATION + ROOT-PRESERVING
  BURNSIDE FORMULA + QUADRATIC-CAPACITY/QUARTIC-JOHNSON TRANSFER + CITED
  SCHMERL--TROTTER PRIME-DELETION INPUT + FINITE-EXACT
  ORDER-TEN PRIME-PARENT CENSUS AND CRITICAL CONTROLS + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. The 7,921,973,322 rows form a class-marked
  presentation cover, not a rooted-orbit or unrooted-isomorphism-class count.
  An exact root/card incidence corollary now gives the weighted quotient from
  rooted augmentation orbits to unrooted noncritical classes, including the
  asymmetric 1/p(T) card weight; its numerical total is not evaluated.
  Rational centrality uses the strict gate 2|C|<D; exact Johnson cosets retain
  the complete (H,c) sidecar. This theorem does not prove order-eleven
  Johnson centrality.
source: codex-frontier-synthesis-creative-20260826at
depends_on:
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4137-strong-tournament-centrality-complete-order-ten
  - THM-4144-order-eleven-large-homogeneous-module-johnson-centrality
  - THM-4163-order-eleven-homogeneous-pair-johnson-centrality
  - THM-4168-prime-order-eleven-nontrivial-automorphism-johnson-centrality
external_input: >
  CITED: J. H. Schmerl and W. T. Trotter, "Critically indecomposable
  partially ordered sets, graphs, tournaments and other binary relational
  structures," Discrete Mathematics 113 (1993), 191--205,
  DOI 10.1016/0012-365X(93)90516-V. Only the classification of critical
  prime tournaments, and hence the prime-deletion consequence, is imported.
audit_script: 04-computation/tournament_prime_parent_quartic_transfer_thm4169_audit.cpp
independent_audit_script: 04-computation/tournament_prime_parent_quartic_transfer_thm4169_independent_audit.py
rebuild_script: 04-computation/tournament_prime_parent_quartic_transfer_thm4169_rebuild.sh
selftest_output: 05-knowledge/results/tournament_prime_parent_quartic_transfer_thm4169_selftest.out
full_q10_output: 05-knowledge/results/tournament_prime_parent_quartic_transfer_thm4169_full_q10_poly.out
independent_audit_output: 05-knowledge/results/tournament_prime_parent_quartic_transfer_thm4169_independent_audit.out
order10_census_output: 05-knowledge/results/tournament_prime_parent_quartic_transfer_thm4169_order10_prime_census.out
unrooted_output: 05-knowledge/results/tournament_prime_parent_quartic_transfer_thm4169_unrooted.out
audit_script_sha256: 4c09fff12656bcaea2c6f8ee913262e99c3a75273f278dac938ce2601ccdf178
independent_audit_script_sha256: a91151b78fc38e9b1b2c27a584aa78eba10881f1b2f85138598b76cb2e4c7ef0
rebuild_script_sha256: 12d64cea5fcf68c56feb19479aaee3708ac7c6333f8a70e38b469e7fd294a77b
selftest_output_sha256: 2f6c6d9a204dc218e00124bef41befbbfdc7bc6259525b8f68516d64c75fb6bd
full_q10_output_sha256: bc4ae33043394739f3adbee66b7b6be0c71949f19dab0fb1173f577c3499b432
independent_audit_output_sha256: 476832e572e4df42122394c8af92abf80f23248ff079f94040e67c6d75b06184
order10_census_output_sha256: f336723f870e48c965450086f1cf6a5d7a2f87db2bea46bec9ef971df690b4b5
unrooted_output_sha256: 807b5fc55ad8af2b293fc5fcdc08184c6b316fc93f9bc0b8c4a9ce74fb87a2a1
gentourng_sha256: 89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110
shortg_sha256: 4eff129bdd76ac58d9f0856856c6bffc3e1d45614ad620a5b74bb5916bbc2bcb
hash_basis: >
  Raw file bytes; tracked text artifacts use LF line endings. The gentourng
  and shortg hashes are hashes of raw executable bytes.
primary_audit: >
  PASS. The warning-clean C++ selftest byte-matches under Apple Clang 17.0.0
  at -O0/-O3 and under ASan/UBSan. The -O3 binary additionally verifies every
  q10 transfer value, scans the full order-ten census, and regenerates the
  rooted-to-unrooted hostile.
independent_audit: >
  ACCEPT. A self-contained Python endpoint-convolution implementation shares
  no transfer code with C++, exhausts small parent universes, and matches all
  1,024 H values, 56,320 capacities, and 4,096 quartic packet coordinates in
  normal, optimized, and fixed-hash modes.
---

# THM-4169 -- prime-parent augmentation and quartic Johnson transfer

**PROVED ELEMENTARY + CITED SCHMERL--TROTTER INPUT + FINITE-EXACT +
VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Prime one-vertex augmentations

Let `Q` be a prime tournament on `q>=3` vertices. For
`z in {0,1}^{V(Q)}`, let `T_z=Q+x_z`, where

```text
x_z -> v  iff  z_v=1.
```

Then the `2^q` attachment patterns split disjointly as follows.

1. Exactly two are nonstrong: `z=0` makes `x` a sink and `z=1` makes `x`
   a source.
2. Exactly `2q` are strong and decomposable. They are indexed by
   `(v,epsilon) in V(Q) x {0,1}`: for `u!=v`, set
   `z_u=1` iff `v->u`, while the mutual bit `z_v=epsilon` is free. Then
   `{x,v}` is a homogeneous pair.
3. The remaining

   ```text
   2^q-2-2q
   ```

   patterns are prime, hence strong.

The boundary `q=3` is exact: the prime remainder is empty.

### Proof

A prime tournament of order at least three is strong. If `T_z` is not
strong, the strong subtournament `Q` lies in one strong component. Since the
only other vertex is `x`, that vertex uniformly dominates or is uniformly
dominated by `Q`, giving the two uniform patterns.

Let `M` be a proper nontrivial module of `T_z`. If `x` is not in `M`, then
`M` is a module of `Q`; primeness leaves only `M=V(Q)`, which occurs precisely
for a uniform extension. If `x` is in `M`, then `M-{x}` is a module of `Q`.
It is neither empty nor all of `Q`, so it is a singleton `{v}`. Thus every
remaining decomposable extension is one of the two clone patterns at `v`.
A pattern cannot clone two different old vertices, because those vertices
would form a module of `Q`. A clone pattern cannot be uniform, because a
source or sink in `Q` would make its complement a nontrivial module. This
proves the partition and its counts. QED.

## 2. Root-preserving Burnside count

The action of `Aut(Q)` on attachment patterns is the correct quotient for
rooted objects `(T_z,x)`. For `g in Aut(Q)`, let `c(g)` be its number of
cycles on `V(Q)` and `f(g)` its number of fixed vertices. The number of
root-preserving isomorphism orbits of prime attachments is

```text
N_Q=(1/|Aut(Q)|) sum_(g in Aut(Q)) (2^{c(g)}-2-2f(g)).
```

A pattern fixed by `g` is constant on the cycles of `g`, giving `2^{c(g)}`
fixed patterns. The source and sink are fixed. The clone patterns are
naturally `(v,epsilon)`, and exactly `2f(g)` of them are fixed. Subtracting
these disjoint nonprime types and applying Burnside proves the formula. QED.

This is a **root-preserving** count. An unrooted child isomorphism may move
`x` to an old vertex and is not controlled by `Aut(Q)`.

## 3. Quadratic capacities and the quartic Johnson packet

Write `H(z)=H(T_z)`. In the integer-capacity convention of THM-4114 and
THM-4128, let

```text
c_ij(z)=2w_ij(T_z)
```

for each unordered child edge, oriented by `T_z`. Then:

- `H(z)` is a multilinear Boolean polynomial of degree at most two;
- every old--old `c_ij(z)` has degree at most two; and
- every new edge `c_ix(z)` is affine and independent of the mutual bit
  `z_i`.

For every child vertex `i`, including `x`, put

```text
d_i(z)=sum_(j!=i)c_ij(z),
h_i(z)=sum_(i->j)c_ij(z)-sum_(j->i)c_ij(z),
C(z)=sum_i d_i(z)h_i(z),
D(z)=sum_(e<f, e disjoint f)c_e(z)c_f(z),
G_+(z)=D(z)+2C(z),
G_-(z)=D(z)-2C(z).
```

The four polynomials `C,D,G_+,G_-` have Boolean degree at most four. At
`q=10` each has at most

```text
sum_(k=0)^4 binom(10,k)=386
```

monomial coefficients, so one subset-zeta transform evaluates it at all
1,024 attachment patterns.

The complete cut-field carrier is the 56-coordinate family

```text
(H(z), (c_e(z))_(e in E(K_11))).
```

Indeed, for every `S subseteq V(T_z)`,

```text
F_(T_z)(S)=H(z)+sum_(i in S,j notin S,i->j)c_ij(z).
```

Consequently THM-4123 and the generic capacity-to-coset evaluator in
THM-4162 recover every exact Johnson-layer lattice, anchor, and coset floor.
The quartic rational packet is not a replacement for this `(H,c)` sidecar.

At order eleven the strict THM-4128 rational-centrality gate is

```text
G_+(z)>0 and G_-(z)>0,
```

equivalently

```text
2|C(z)|<D(z).
```

Equality is a tie with an outer layer and is not part of the strict claim.

### Degree proof

THM-4114 gives

```text
H(z)=H(Q)+sum_(i->j in Q)c^Q_ij z_i(1-z_j),
```

which is Boolean-quadratic. For any fixed permutation of `Q+x+y`, the
Hamilton-path indicator uses at most two literals involving the attachment
bits of `x`, because `x` has at most one predecessor and one successor.
Thus every one-ear value of `T_z`, and hence every four-value finite
difference defining an old--old capacity, has degree at most two.

For `c_ix`, contract the two possible adjacent blocks `ix` and `xi` in the
exposed-gap formula. At most one further old arc is incident with `x`, and
summing the two block orders removes the orientation of the mutual edge.
Thus `c_ix` is affine and independent of `z_i`.

Unsigned degrees have degree at most two. For an old vertex `i`, the new-edge
term in the signed degree is `(1-2z_i)c_ix(z)` and still has degree at most
two; at `x`, the signed degree has degree at most two while the unsigned
degree is affine. Hence every term of `C` has degree at most four. Every term
of `D` is a product of two quadratic capacities, so it also has degree at
most four. QED.

## 4. Cited deletion reduction and finite presentation cover

**CITED INPUT (Schmerl--Trotter, not proved here).** At each odd order at
least five, the only critical prime tournaments are
`T_(2k+1),U_(2k+1),W_(2k+1)`; there are no even-order critical tournaments.
Therefore every prime order-eleven tournament outside
`T_11,U_11,W_11` has a prime order-ten deletion.

**FINITE-EXACT.** Nauty's one-representative-per-isomorphism-class order-ten
stream contains

```text
9,733,056 tournament classes, of which 7,906,161 are prime.
```

Each prime order-ten representative has exactly

```text
2^10-2-20=1,002
```

prime attachment patterns. Hence

```text
7,906,161 * 1,002 = 7,921,973,322
```

class-marked labelled prime-attachment presentations form a surjective cover
of all noncritical prime order-eleven isomorphism classes.

This product is **not** a root-preserving orbit count and **not** an unrooted
class count. Parent automorphisms identify rooted presentations, and one
unrooted child can have multiple prime deletion cards. The exact total of
root-preserving orbits would be `sum_[Q] N_Q`; it is not computed here.

## 5. Critical controls

**FINITE-EXACT.** In the integer-capacity scale, with central layers
`t=+-1`, the three cited exceptional families have:

| row | `H` | `C` | `D` | prime cards | strong cards | exact-coset central margin | actual central margin |
|---|---:|---:|---:|---:|---:|---:|---:|
| `T_11` | `93,027` | `0` | `811,847,516,832` | `0` | `11` | `28,632` | `29,268` |
| `U_11` | `15,611` | `0` | `57,104,620,276` | `0` | `9` | `7,370` | `6,090` |
| `W_11` | `459` | `0` | `581,082,528` | `0` | `8` | `732` | `924` |

For each row, the rational, exact-coset, and actual optimizing imbalance sets
are exactly `{+1,-1}`. Each row is self-converse: use `i |-> -i` for `T_11`,
`i |-> 5-i (mod 11)` for the stated `U_11` convention, and reverse the
ten-vertex spine while fixing the apex for `W_11`. Converse negates `C` and
fixes `D`; self-converseness explains `C=0`.

## 6. Quotient losses and hostile controls

For the prime parent

```text
Q=111111101111111111111101111110110111110111111
```

the exact audit gives `|Aut(Q)|=1` and 1,002 prime attachment patterns, hence
1,002 root-preserving prime orbits. Nevertheless `shortg` finds only 1,000
unrooted child classes, with the two doubletons

```text
336 ~ 432,      368 ~ 400.
```

Thus even an asymmetric parent does not make all 1,002 unrooted children
pairwise nonisomorphic; forgetting the distinguished root is the lost
coordinate.

### Exact root/card incidence quotient

Let `Q` range over one representative of each prime order-ten isomorphism
class, and let `A_Q` be its prime attachment patterns. There is a canonical
bijection

```text
disjoint_union_[Q] Aut(Q)\A_Q
  <--> rooted prime pairs [(T,x)] with T-x prime.        (24)
```

Indeed, choosing an isomorphism `T-x -> Q` records the attachment pattern of
`x`; changing that identification acts by `Aut(Q)`. Conversely an attachment
recovers the rooted pair, and rooted isomorphisms give exactly the same orbit.

For a prime order-eleven child `T`, put

```text
P(T)={v in V(T):T-v is prime},
a(T)=|Aut(T)\P(T)|.                                    (25)
```

The fibre of the forget-root map over `[T]` has exactly `a(T)` elements.
Therefore weighting each rooted class by `1/a(T)` gives the exact unrooted
noncritical class count:

```text
boxed:
N_noncritical
 =sum_[Q] sum_([z] in Aut(Q)\A_Q) 1/a(T_z)
 =sum_[Q] 1/|Aut(Q)| sum_(z in A_Q)
      |Stab_(Aut(Q))(z)|/a(T_z).                        (26)
```

The second equality is orbit--stabilizer, orbit by orbit. The stabilizer is
canonically `Aut(T_z,x)`. On the asymmetric-child stratum it is trivial and
`a(T)=p(T):=|P(T)|`, so

```text
boxed:
N_(asymmetric,noncritical)
 =sum_[Q] 1/|Aut(Q)|
    sum_(z in A_Q, Aut(T_z)=1) 1/p(T_z).                (27)
```

This is an exact quotient formula, not an evaluation. For the hostile
doubletons `336~432` and `368~400`, all four children are asymmetric and
have exactly six prime deletion cards. The missing root/card weight is thus
literally `1/6` on each rooted presentation in these fibres.

For the same parent, the source pattern `1023` is nonstrong and has

```text
2|C|=3,818,407,904,
D=2,735,733,720,
D-2|C|=-1,082,674,184.
```

It is a positive hostile for the strict rational gate but is outside the
strong/prime target. Pattern `0` has positive margin `209,423,320`, and all
1,022 strong extensions of this particular parent pass. No all-parent
centrality conclusion follows.

Finally, prime child pattern `308` has exact layer lattices, for
`m=1,...,10`,

```text
(4,2,2,2,2,2,2,2,2,2).
```

The rational quartic gate therefore cannot replace the exact Johnson-coset
sidecar.

## 7. Reproduction

```text
bash 04-computation/tournament_prime_parent_quartic_transfer_thm4169_rebuild.sh
```

The recorded rebuild used `/usr/bin/clang++`, Apple clang version 17.0.0
(`clang-1700.0.13.5`), target `arm64-apple-darwin24.3.0`. The script resolves
`clang++` through `PATH` and does not pin a compiler binary or version. Its
selftest runs at `-O0/-O3` and under ASan/UBSan; the `-O3` binary runs the full
q10 cube, order-ten census, and `shortg` hostile. A separate Python
endpoint-convolution audit runs in normal, optimized, and fixed-hash modes.
Frozen outputs and tool hashes are bound above.

## 8. Scope

The source is a chosen prime order-ten representative together with an
attachment bit pattern. The target is the child carrier `(H,c)` and its
quartic packet `(C,D,G_+,G_-)`. Relabeling preserves the Johnson predicates,
but quotienting by parent automorphisms loses presentation labels, and
forgetting the root can merge further children. Root/card incidence is the
sidecar for class questions; `(H,c)` and the Johnson anchors/lattices are the
sidecar for exact-coset questions.

This theorem proves a complete finite reduction. It does **not** evaluate the
7.92-billion-presentation gate/coset bank, prove order-eleven Johnson
centrality, count order-eleven rooted or unrooted classes, or locate actual
maximizers outside the three critical controls. Combined with THM-4168, the
remaining open prime stratum is asymmetric; the number of still-open
presentations or classes has not been computed. QED.
