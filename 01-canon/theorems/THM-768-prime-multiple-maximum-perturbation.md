---
id: THM-768
title: A unique prime-divisible maximum forces a strict lonely-runner witness
status: PROVED (elementary prime-residue perturbation)
source: codex-2026-07-14-S3 continuation (n=12 sporadic-branch audit)
depends_on: []
related:
  - THM-593   # unit-residue improvement / formal residue pinning
  - THM-763   # finite height for tight sets
  - THM-765   # hereditary primitivity
  - THM-766   # first-window tooth ladder
  - HYP-6820  # q<=25 and n=12 uniformity audit
formal_companion:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCResiduePinning.lean
---

# THM-768 — a unique prime-divisible maximum is strictly loose

## Statement

Let `p>=3` be prime and let `S` be a finite set of distinct positive integer
speeds with

```text
|S| <= p-1.
```

Write `V=max(S)`.  If `V` is the unique member of `S` divisible by `p`, then

```text
M(S) > 1/p.                                             (1)
```

In particular, let `A={a_1<...<a_12}` satisfy `M(A)=1/13`.  Then `A` cannot
have `a_12` as its unique multiple of `13`.

Combining this with the kernel-checked prime-residue pinning theorem gives the
following exact partition of any putative tight twelve-set:

1. **no multiple of 13:** the residues are exactly `{1,...,12}` modulo `13`;
2. **a multiple below the maximum:** `13|a_i` for some `i<12` and
   `13 does not divide a_12`; or
3. **a multiple at the maximum:** `13|a_12`, and at least one other speed is
   also divisible by `13`.

Thus the full-residue lift branch is forced on the zero-owner-free locus, and
the previously undifferentiated complementary locus reduces to two explicit
multiple-owner branches.

## Proof

Remove `V` and reduce the remaining speeds modulo `p`.  Every one of their
residues is nonzero, and there are at most `p-2` of them.  Hence some

```text
u in {1,...,p-1}
```

is absent.  Since `p` is prime, choose `a in {1,...,p-1}` such that

```text
a u = -1 (mod p).                                      (2)
```

Multiplication by `a` permutes the nonzero residue classes.  Therefore the
residues

```text
r_v = a v (mod p),       v in S minus {V},
```

belong to `{1,...,p-2}`: residue `p-1=-1` is absent by (2).

Take the explicit time

```text
t = a/p + (1+1/(2V))/(pV)       (mod 1).               (3)
```

Because `p|V`, the phase of the maximum speed is

```text
V t = integer + (1+1/(2V))/p.
```

It is strictly between `1/p` and `1-1/p`: indeed `V>=p>=3`, so
`1+1/(2V)<2<=p-1`.

For `v<V`, use the residue `r_v` above.  A real lift of its phase at (3) is

```text
x_v = (r_v + (v/V)(1+1/(2V)))/p.
```

Since `1<=r_v<=p-2` and `v<=V-1`,

```text
0 < (v/V)(1+1/(2V))
  <= ((V-1)/V)(1+1/(2V))
  = 1 - 1/(2V) - 1/(2V^2)
  < 1.
```

Consequently

```text
1/p < x_v < (p-1)/p = 1-1/p.                          (4)
```

At the explicit time (3), every speed therefore has circular distance
strictly greater than `1/p`.  This proves (1).
∎

## Why the maximum and uniqueness hypotheses are real

The perturbation fee paid by a speed `v` is proportional to `v/V`.  The
strict inequality `v<V` is exactly what keeps every occupied residue below
the upper danger boundary after the missing class is moved to `-1`.

If another speed `v<V` is also divisible by `p`, its unperturbed residue is
zero and its perturbed clearance is approximately `v/(pV)<1/p`; the argument
then fails for a genuine reason.  Likewise, if the unique `p`-multiple lies
below a larger nonmultiple speed, the phase displacement of that larger speed
can exceed one residue cell.  Those are precisely branches 2 and 3 above,
not omissions hidden by the proof.

## Relation to residue pinning and the endpoint object

`LRCResiduePinning.lean` proves that a `13`-tight-from-above twelve-family with
no `13`-multiple occupies every nonzero residue class exactly once.  THM-768
supplies the first uniform theorem on the complementary zero-owner locus.  In
endpoint language, a missing unit class is switched to the incoming boundary
`-1`; moving slightly to its other side clears every nonzero tooth, while the
unique zero owner at the maximum is moved just beyond its own boundary.

The proof is a one-wall version of the endpoint-owned incidence principle:
the residue class alone is insufficient; its owner height determines whether
the wall perturbation is legal.

## Tournament Analysis and assumption challenge

Runner vertices obscure the proof.  Use the nonzero residue classes as
vertices, mark occupied classes by their owner heights, and orient two classes
by which owner reaches the upper boundary first under positive perturbation.
Multiplication by `a` is the switch/gauge; the missing class is normalized to
`-1`, and the cyclic residue order supplies the tie Hamiltonian path.  The
unweighted tournament is transitive after this cut and carries no more than
the order.  The theorem needs the owner-height sidecar `v/V` and the zero-owner
bit.  Quotienting either away destroys (4).

Alternate vertices considered were runners, fixed circle sections, boundary
events, residue classes, and proof obligations.  The smallest carrier
preserving the strict-witness predicate here is the marked residue-boundary
obligation, not the runner tournament.

## Scope

THM-768 does not prove that the `n=12` sporadic branch is empty.  It makes its
residue split quantifier-honest:

- the zero-owner-free branch is the full-residue lift problem;
- the unique-maximum-zero-owner branch is empty; and
- the residual zero-owner problem splits according to whether the maximum is
  a nonmultiple; if the maximum is itself a zero owner, there are at least two
  zero owners.

The next task is to combine the component/splice defect with those two
multiple-owner branches, rather than assuming every tight set is already a
full-residue lift.
