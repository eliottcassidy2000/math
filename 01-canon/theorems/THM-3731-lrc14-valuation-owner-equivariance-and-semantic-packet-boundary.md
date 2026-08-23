---
id: THM-3731
title: "LRC14 valuation-owner equivariance and semantic-packet boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The 165 THM-2349 rows index
  sorted valuation-orbit representatives Q.  On the labelled valuation lift,
  exactly 15 repeated profile orbits obstruct a deterministic S3-equivariant
  minimum-depth owner selector depending only on those valuations.  Retaining pairwise-distinct positive
  coefficients repairs the owner label by argmin_i(nu_13(c_i),c_i).
  An exact abstract-interface THM-2478 tensor hostile shows that positive aggregate drift can
  still vanish on the canonical owner atom.  No scalar-cover row is excluded
  and no LRC(14) conclusion is proved.
source: root + lrc14-circulant-branch / 2026-08-22
audit: >
  PASS after deleting an overbroad signed-coefficient aside.  The independent
  audit rederived the stabilizer iff criterion, the 165=150+15, 945, and 990
  counts, S3 equivariance, positive-rational-dilation invariance, THM-2349
  positivity composition, and the exact owner-loop drift 21/742586.  Normal
  and optimized runs byte-match the frozen transcript.
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
related:
  - THM-3718-lrc-complete-atom-orbit-defect-saturation-and-semantic-boundary
script: 04-computation/lrc14_valuation_owner_equivariance_thm3731.py
output: 05-knowledge/results/lrc14_valuation_owner_equivariance_thm3731.out
script_sha256: 98506113bdb967ca4cff94cc8e793c88dae3364f57c7fa491e33096731ce687b
output_sha256: 6820f630597ca8be9fe8e3d553e2b890fcd6c51e1465020e3160f0b96b4e5874
hash_basis: raw LF bytes
---

# THM-3731 -- the valuation quotient forgets the owner tie-break

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  This theorem separates
the static owner-label problem from the owner-supported semantic-arrival
problem.  The former is repaired by full positive coefficients; the latter
remains open.

## 1. The exact valuation quotient

Let

```text
Q={(1,b,c):5<=c<=19, 1<=b<c}.
```

The THM-2349 ledger consists of these `165` sorted valuation profiles:

```text
150 strict profiles (1,b,c) with 1<b<c,
15 repeated profiles (1,1,c).
```

It is not a ledger of `165` coefficient realizations.  Let

```text
X={lambda in {1,...,19}^3:sort(lambda) in Q},
G=S_3,
A(lambda)=argmin_i lambda_i.
```

Here `G` relabels the three blockers.  The labelled lift has

```text
150*6+15*3=945
```

points.

## 2. Stabilizer-fixed section criterion

More generally, let a group `G` act on states `X` and labels `J`, and
let `A(x) subset J` be a `G`-equivariant admissible-choice field.  On the
orbit of `x`, an equivariant deterministic selector

```text
f(gx) in A(gx),                 f(gx)=g f(x)           (1)
```

exists if and only if

```text
A(x)^(Stab_G(x)) is nonempty.                          (2)
```

Indeed, if `h x=x`, equivariance forces `f(x)=h f(x)`.  Conversely,
choose a stabilizer-fixed `j in A(x)` and define `f(gx)=gj`.  If
`gx=g'x`, then `g'^(-1)g` lies in the stabilizer and fixes `j`, so the
definition is well-defined.

For a strict profile, `A(lambda)` is its singleton minimum and (2) holds.
For every `(1,1,c)`, the stabilizer swaps the two admissible depth-one
labels, so (2) fails.  Thus exactly

```text
150 selector orbits pass,             15 fail.         (3)
```

The smallest obstruction is `(1,1,5)`.  The canonical quotient-level
repairs are the set-valued selector `A(lambda)` or the marked lift
`(lambda,j)`, `j in A(lambda)`.  The latter has

```text
150*6+15*3*2=990
```

points.  The obstruction concerns a single-valued section only.

## 3. Full positive coefficients repair the owner label

Let `c=(c_1,c_2,c_3)` be pairwise-distinct positive rational blocker
coefficients.  Define the lexicographic key

```text
j_can(c)=argmin_i (nu_13(c_i),c_i).                    (4)
```

It is unique.  For every blocker permutation `pi in S_3`,

```text
j_can(pi.c)=pi(j_can(c)).                              (5)
```

For every common positive rational dilation `d`,

```text
nu_13(dc_i)=nu_13(d)+nu_13(c_i),
dc_i<dc_j iff c_i<c_j,
```

so

```text
j_can(dc)=j_can(c).                                    (6)
```

On a strict profile, (4) chooses the unique depth-one blocker.  On a repeated
profile, it resolves the equal-valuation pair by the retained coefficient
order.  For a hypothetical scalar-cover realization on the 165-profile
frontier, THM-2349 proves that both repeated-profile depth-one exclusive-owner
sets are positive, so (4) selects a positive owner on that realization.
THM-2309 may then run its owner-aligned pivot construction with this label.

The key does not factor through the valuation quotient: reversing the two
coefficient magnitudes inside an equal-valuation pair reverses the chosen
owner without changing the profile.  The exact symmetry claim is (5) plus
(6), not an undefined broader cover-equivariance claim.

## 4. The first unpaid semantic interface

Let

```text
O=(0,0,0,0,0,1,0)
```

be the unique complete local semantic source atom: five guard/unit-safe bits,
the chosen owner dangerous, and the other nondeep blocker safe.  Retaining
the label (4) does not prove that THM-3718's adaptive drifting atom is `O`.
The first required implication would be

```text
positive aggregate complete-atom drift
  => positive drift on the canonical owner atom O.     (7)
```

It is false at the maintained abstract tensor interface.

Later interfaces are distinct and are not reached by (4):

- a terminal word is the temporal edge weight
  `W_k(tau,sigma)=mu(P_tau intersect E_j intersect T^(-k)Q_(j,sigma))`,
  not a state label;
- the predecessor owner root and the translated THM-2365 deep root retain
  different ancestry/sheet coordinates;
- THM-2531's excluded target is empty and is not a positive future arrival.

Thus the coefficient sidecar repairs only the static owner section, not an
owner-supported temporal/root/arrival packet.

## 5. Exact owner-loop hostile

On `F_13^3`, define

```text
H_O(r,s,t)=1/4 * 1_(r-t=1),
H_N(r,s,t)=1/4 * 1_((r,s,t)=(1,1,0)).                 (8)
```

Both tensors are nonnegative and diagonal-zero.  Use the literal THM-2365
projection and drift

```text
(PH)(r,s,t)=1/169 sum_(a,b) H(r+b,s+a,t+b),
D(H)=1/2197 sum_(r,s,t)(H-PH)^2.                       (9)
```

Exact arithmetic gives

```text
P H_O=H_O,
D(H_O)=0,
sum_r H_N(r,0,0)=0,
D(H_N)=D(H_O+H_N)=21/742586>0.                        (10)
```

Place `H_O` on the semantic owner atom, `H_N` on one nonowner atom, and
zero on the other complete atoms.  Every current nonnegativity,
diagonal-zero, and aggregate tensor identity survives, yet all positive
drift lies off `O`.  Therefore aggregate THM-2365 drift plus complete-atom
partition identities cannot prove (7).

This is an interface hostile, not a realized scalar-cover row.  A
cover-specific constraint could still exclude it.

## 6. Exact residual and next typed test

The `165` ledger stores no genuine cover coefficient realizations; the
maintained all-profile coefficient family is explicitly a noncover hostile.
The next lawful physical object must be formed before integration on one
ancestry base:

```text
C_(sigma,r_owner,r_deep;s,t)
 =mu(E_(j_can) intersect T^(-k)Q_(j_can,sigma)
     intersect owner-root(r_owner)
     intersect deep-probe(r_deep;s,t)
     intersect selected endpoint mask).               (11)
```

The next decisive test is whether one actual scalar-cover row forces a
nonzero THM-3713 edge/three-site mask on some table (11).  A positive result
would begin the owner/word/root/arrival bridge; a negative row would identify
the first destroyed coordinate.

No scalar-cover row is excluded here, none of the `165` frontier rows is
closed, and `LRC(14)` remains open.

## 7. Reproduction

```bash
python3 -B 04-computation/lrc14_valuation_owner_equivariance_thm3731.py
python3 -B -O 04-computation/lrc14_valuation_owner_equivariance_thm3731.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
