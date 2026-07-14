---
id: THM-770
title: Full-residue endpoint-owner rigidity through lift height twelve
status: PROVED (finite-exact owner-CSP at H=12; dilation and gcd descent are elementary)
source: codex-2026-07-14-S3 (endpoint/splice continuation)
depends_on:
  - THM-769   # only for the bounded shallow-tight corollary C
related:
  - THM-762   # signed pair decks at small rational periods
  - THM-763   # global finite height for primitive tight instances
  - THM-765   # hereditary primitivity / gcd-deck context
  - THM-766   # component-tooth cones and endpoint alignment
  - HYP-6775 # open tight twelve-speed rigidity
  - HYP-6800 # open n=12 sporadic branch
  - HYP-6820 # q<=25 and n=12 uniformity audit
verification: 04-computation/lrc13_full_residue_endpoint_owner_h12_codex_S3.cpp
  (+ 05-knowledge/results/lrc13_full_residue_endpoint_owner_h12_codex_S3.out)
---

# THM-770 — Full-residue endpoint-owner rigidity through lift height twelve

## Statement

For a finite set `W` of positive integer speeds and an integer `q>=3`, put

```text
G_q(W) = {t in R/Z : min_(w in W) ||w t|| > 1/q},
chi_q(W) = number of connected components of G_q(W).
```

Thus `chi_q(W)=0` exactly when there is no strict `1/q` witness.  In the
endpoint language of the companion audit, whenever the open danger union is
not the whole circle,

```text
chi_q(W)=kappa_q(W)-P_q(W),
```

where `kappa_q` counts open danger components and `P_q` counts protected
end/start splices.

### A. Dilation and gcd descent

For every positive integer `c`,

> **`chi_q(cW)=c chi_q(W)`.**                                      (1)

If the open danger union is not the whole circle, the two endpoint terms
separately scale:

> **`kappa_q(cW)=c kappa_q(W)` and `P_q(cW)=c P_q(W)`.**             (2)

Consequently, if `g=gcd(W)` and `U=W/g`, then

> **`chi_q(W)=g chi_q(U)`.**                                       (3)

If `q=13` and `W` is a complete nonzero residue system modulo `13`, then
`13` does not divide `g`, and `U` is again a complete nonzero residue system
after multiplication of its residue labels by `g^(-1) mod 13`.  Zero defect
therefore descends exactly to the primitive quotient.

### B. Finite-exact full-residue classification through height twelve

For a lift vector

```text
k=(k_1,...,k_12),       0<=k_r<=12,
```

define the labelled full-residue packet

```text
W(k)={r+13 k_r : 1<=r<=12}.                              (4)
```

Among all

```text
13^12 = 23,298,085,122,481
```

packets in (4),

> **`chi_13(W(k))=0` if and only if**
>
> **`W(k)=c*{1,2,...,12}`**
>
> for
>
> **`c in {1,2,...,12,14}`.**                             (5)

These are precisely the dilated arithmetic progressions contained in the
box whose scale is a unit modulo `13`.  In particular:

> **the unique primitive zero-defect packet in the height-twelve box is**
>
> **`{1,2,...,12}`.**                                     (6)

Equivalently, every other primitive full nonzero residue transversal in
`{1,...,168}` has a strict `1/13` witness.  Because every full-residue packet
has clearance exactly `1/13` at each `j/13`, `1<=j<=12`, zero defect here is
also exactly the equality `M(W)=1/13`.

### C. Bounded shallow tightness corollary

Let `A` be a primitive tight twelve-speed set with `max(A)<=168`.  If no
speed in `A` is divisible by `13`, then

> **`A={1,2,...,12}`.**                                   (6a)

Indeed, THM-769 proves that the no-multiple-of-13 (shallow) branch of a tight
twelve-set is exactly the complete nonzero residue branch.  The height bound
then places `A` in (4), and Part B applies.  Thus every other bounded tight
candidate must lie in THM-769's deep binding-scale branch.

## Proof

### 1. Scaling is a circle-cover identity

Let

```text
phi_c : R/Z -> R/Z,       phi_c(t)=ct.
```

Then, directly from the definitions,

```text
t in G_q(cW)
  iff ||cwt||>1/q for every w in W
  iff phi_c(t) in G_q(W).
```

Hence

```text
G_q(cW)=phi_c^(-1)(G_q(W)).                               (7)
```

Every component of `G_q(W)` is a proper open arc: `0` never belongs to the
strict safe set.  The inverse image of a proper open arc under the `c`-fold
circle covering `phi_c` is a disjoint union of exactly `c` open arcs.  This
proves (1).

The open danger union satisfies the same pullback identity.  When it is not
the whole circle, all its components are proper arcs, so their number
multiplies by `c`.  A protected splice is a local end/start incidence with no
continuing open tooth; a covering map is a local homeomorphism and gives each
such incidence exactly `c` inverse images.  This proves (2).  Taking `c=g`
and `W=gU` proves (3).

For a full nonzero residue packet modulo `13`, its gcd cannot be divisible by
`13`.  Division by `g` multiplies every residue by the unit `g^(-1)`, which
permutes the twelve nonzero classes.  This proves the full-residue descent
claim.

### 2. Universal protected points in the full-residue branch

At `t=j/13`, multiplication by the unit `j` permutes the nonzero residue
classes.  Thus

```text
{||w j/13|| : w in W}
```

has minimum `1/13`.  Exactly the two runners whose residues are
`+j^(-1)` and `-j^(-1)` bind; every other runner has clearance at least
`2/13`.  These are protected complementary-pair splices.  In particular the
danger union is not the whole circle, and

```text
chi_13(W)=0  iff  M(W)=1/13                              (8)
```

inside the full-residue branch.

### 3. Exact atomic-cell reduction

For a candidate speed `w`, its open danger teeth have endpoints

```text
L_(w,m)=(13m-1)/(13w),
R_(w,m)=(13m+1)/(13w),          m in Z/wZ.                (9)
```

Take the union `E` of the endpoints (9) for all 156 possible speeds
`r+13k`, `1<=r<=12`, `0<=k<=12`.  The components of `(R/Z)\E` are the
**atomic cells**.  The truth of

```text
||w t||<1/13
```

is constant on every atomic cell for every candidate speed `w`.

A packet has a strict witness if and only if it leaves an atomic cell
uncovered.  One direction is immediate: the rational midpoint of an
uncovered cell is strictly safe, since equality can occur only at an endpoint
in `E`.  Conversely, strict safety is an open condition.  Even if a strict
witness were itself in `E`, a small neighbourhood of it would contain an
atomic cell and remain strictly safe.  Thus the finite predicate is exactly
set cover on the atomic cells, with one choice required from each residue
class.

This is also an endpoint-Farey frontier formulation.  Comparing the left
endpoint of the `(v,n)` tooth with the right endpoint of the `(u,m)` tooth is
the integer inequality

```text
L_(v,n) <= R_(u,m)
  iff 13(un-vm) <= u+v.                                  (10)
```

Equality in (10) forces `13|(u+v)` and hence complementary nonzero residue
owners.  Protection is the additional assertion that no third open tooth
continues through that equality.  Formula (10) is the local edge; the atomic
cell/owner incidence hypergraph retains the simultaneous global predicate.

### 4. Lossless unique-owner recursion

At a node of the search, let `U` be the set of atomic cells not covered by the
choices already made, and let `R` be the remaining residue groups.  For each
`r in R`, let `C_r` be the union of the cells covered by all thirteen options
in residue group `r`.  Define

```text
F_r = U \ union_(s in R\{r}) C_s.                        (11)
```

Every cell in `F_r` is coverable by no remaining group except `r`.  Therefore
any successful completion must choose an option in group `r` which covers
all of `F_r`.  Deleting the other options is logically exact.  If no option
covers `F_r`, the node has no completion.  Otherwise one may branch on any
remaining group; choosing the group with the fewest surviving options changes
only the running time.

Repeated application of (11) enumerates every and only zero-defect packets.
The exact implementation uses integer rational midpoints, bit-vector cell
incidences, and `__int128` cross-products.  Its height-twelve certificate is

```text
distinct endpoints       24,008
atomic cells              24,008
64-bit words per cell set    376
search nodes          32,708,254
owner prunes          18,459,588
complete leaves               13
uncovered leaves                0
zero-defect rows                13.
```

Literal endpoint sweeps of the thirteen leaves give

```text
kappa_13=P_13=12c,       chi_13=0,
```

and sorting each leaf gives `c*{1,...,12}` for exactly the scales in (5).
There is one gcd-one leaf, at `c=1`.  This proves (5) and (6).

As an independent reduction audit, the program literally endpoint-sweeps all
`3^12=531,441` height-two packets and compares their zero/nonzero
classification with a separate atomic-cell owner search.  Both find exactly
the three rows with scales `1,2,3`.  Among the 531,375 primitive non-AP rows,
the minimum positive endpoint defect is `2`.  This checks the cell reduction
without relying on the height-twelve DFS itself.  ∎

## What this proves—and what remains open

Part A turns the global equality problem into a clean descent:

```text
full-residue zero defect
  -> divide by the common gcd
  -> primitive full-residue zero defect.                 (12)
```

The exact next theorem is therefore the **primitive descent trigger**

> `W` full nonzero modulo `13` and `chi_13(W)=0`
> implies `W={1,...,12}` or `gcd(W)>1`.

If proved, repeated use of Part A would force every zero-defect full-residue
packet to be a dilation of `{1,...,12}`.  THM-770 does **not** prove this
unbounded trigger.  It proves it exactly through lift height twelve, a much
larger scale-normal box than the earlier height-one audit.  THM-769 makes the
remaining split precise: Part C closes the bounded shallow branch, while the
deep binding-scale branch (`s>=2`) and shallow packets above height twelve
remain outside this finite certificate.

## Tournament Analysis and assumption challenge

The computation tests an alternate vertex set explicitly.

- **Exact vertices:** atomic endpoint cells and the 156 residue-labelled speed
  options.  Their incidence hypergraph preserves the finite strict-witness
  predicate exactly.
- **Diagnostic tournament vertices:** the twelve residue-choice obligations.
  For a pair `r,s`, erase the possible coverage of the other ten groups and
  compare the cells exclusively ownable by `r` and by `s`.  The sign of this
  difference is the switch/gauge; residue order `1->2->...->12` resolves ties.
  The companion output records scores, directed triangles, SCCs, edge flips,
  and Hamiltonian-path count.
- **Destroyed information:** pairwise owner pressure does not record which
  single option in a residue group covers several obligations simultaneously.
  It is therefore telemetry, not the certificate.

Runners, complementary residue pairs, teeth, atomic cells, and proof
obligations were all considered as vertices.  Complementary pairs preserve
the equality-splice labels in (10) but lose overlap rank and third-runner
blocking.  Atomic cells plus option incidence are the smallest carrier used
here that preserves the exact bounded LRC predicate.  The challenged
assumption is that the right global object should still be a runner
tournament: the proof lives in endpoint-owned hypergraphic propagation.

The height-twelve fingerprint makes the loss concrete.  All `66` pairwise
owner-pressure comparisons tie.  The imposed tie path therefore produces the
transitive score histogram `{0:1,1:1,...,11:1}`, no directed triangle,
twelve singleton SCCs, zero edge flips, and one Hamiltonian path.  Meanwhile
the incidence hypergraph distinguishes thirteen solutions from more than
`2.3e13` candidates.  The bare tournament has collapsed exactly where the
owner hyperedges do the proof.

## Reproduction

```bash
c++ -O3 -std=c++17 \
  04-computation/lrc13_full_residue_endpoint_owner_h12_codex_S3.cpp \
  -o /tmp/lrc13_full_residue_endpoint_owner_h12
/tmp/lrc13_full_residue_endpoint_owner_h12
```

The stored output is
`05-knowledge/results/lrc13_full_residue_endpoint_owner_h12_codex_S3.out`.
The final reference rerun took `132.59` seconds wall time and reproduced the
stored output byte-for-byte.  SHA-256 checksums are

```text
64e9cbdfedd5b5c2ed782e4ac7628b3ba2fe3804c086ecffefe18bb10acbb93d  script
f4a53c2b0c22e8abb61cb3cc72a5bfebff164686a80d151df75f27ea0195be07  output
```
