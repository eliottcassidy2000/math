---
id: THM-2810
title: "Factorial-Hankel faithfulness and bounded radial-carrier no-go"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT AUDIT.  In characteristic zero the factorial functional
  L(s^n)=n! has nondegenerate Hankel pairing in every finite degree, with
  determinant product_(j=0)^d(j!)^2.  Its algebra annihilator is therefore
  zero: it factors through no proper polynomial-algebra quotient and admits
  no fixed finite-dimensional multiplicative, linear-state, periodic, or
  confluent-Prony carrier for the whole tower.  The cyclic alias
  s^a(1-s^N) is explicit.  In characteristic p the annihilator is exactly
  (s^p), giving the sharp finite-carrier boundary.  This does not refute
  effective HYP-8765, replace THM-2022, globally de-factorialize Wick
  channels, or prove a new GMC result.  Until independent promotion, no
  proved result may depend on this candidate.
source: root/gmc-factorial-carrier-boundary-2026-07-28
depends_on: []
related:
  - THM-1620-the-pochhammer-bridge-toral-legendre-radial-hermite
  - THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
  - THM-2638-radial-height-graded-wick-decoder-and-laplace-forgetting-boundary
  - THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate
script: 04-computation/gmc_factorial_hankel_faithfulness_thm2810.py
output: 05-knowledge/results/gmc_factorial_hankel_faithfulness_thm2810.out
script_sha256: 6bb4be6f620149e11e49ad2bf34fbd0138623a2327028a191605175cde9922cf
output_sha256: fcf74518f1f2df4cd307ee055929a6b54fe7bbffa3e1ea7842ae0522432b6954
hash_basis: LF-normalized bytes
---

# THM-2810 -- the factorial tower has no bounded multiplicative carrier

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT AUDIT.**

THM-2638 identifies radial height as the coordinate forgotten by Gaussian
radial integration.  A natural response is to seek a finite quotient,
periodic address, or Prony state that remembers just enough height to recover
the factorial readout.  This theorem proves that no fixed tower-wide carrier
of those kinds can exist in characteristic zero.

The obstruction is not asymptotic.  The factorial functional is faithful as
an algebra-pairing functional: every nonzero polynomial is detected after
multiplication by another polynomial.  Finite cyclic carriers fail for the
most literal possible reason--they alias two exponents whose factorials are
different.  Characteristic `p` is the exact opposite boundary: the
annihilator becomes `(s^p)`.

## 1. The factorial Hankel factorization

Let `K` be a field of characteristic zero and define

```text
L:K[s] -> K,                         L(s^n)=n!.          (1)
```

For `d>=0`, let

```text
H_d=((i+j)!)_(0<=i,j<=d).                                (2)
```

Put

```text
D=diag(0!,1!,...,d!),
P_(i,k)=binom(i,k),                    0<=i,k<=d,         (3)
```

with `P_(i,k)=0` for `k>i`.  Vandermonde's identity gives

```text
(D P P^T D)_(i,j)
 =i!j! sum_k binom(i,k)binom(j,k)
 =i!j! binom(i+j,i)
 =(i+j)!.                                                 (4)
```

Therefore

```text
H_d=D P P^T D.                                           (5)
```

The Pascal matrix `P` is unit lower triangular, so

```text
det(H_d)=product_(j=0)^d (j!)^2 !=0.                     (6)
```

This is the determinant formula already adjacent to THM-1620's Pochhammer
bridge, now used as an algebra-faithfulness theorem rather than only a
positivity/Favard input.

## 2. The algebra annihilator is zero

Define the largest multiplication-stable part of the kernel by

```text
Ann(L)={q in K[s]: L(qg)=0 for every g in K[s]}.          (7)
```

Let `q!=0` have degree at most `d`.  Its coefficient vector is a nonzero row
in the bilinear pairing

```text
B_d(f,g)=L(fg),          f,g in K[s]_(<=d),              (8)
```

whose matrix in the monomial basis is `H_d`.  By `(6)`, `B_d` is
nondegenerate.  Hence some `g` of degree at most `d` satisfies

```text
L(qg)!=0.                                                 (9)
```

Thus

```text
Ann(L)=0.                                                (10)
```

This has an exact quotient interpretation.  The functional `L` factors
through `K[s]/I` precisely when

```text
I subset Ann(L).                                        (11)
```

Indeed `(11)` is necessary because `q in I` implies `qg in I`; it is
sufficient because `L` then has a well-defined linear readout on the
quotient.  Equation `(10)` proves:

> **Factorial-faithfulness theorem.**  In characteristic zero, `L` factors
> through no proper algebra quotient of `K[s]`.

More generally, suppose `A` is a finite-dimensional `K`-algebra,
`pi:K[s]->A` is an algebra map, and `lambda:A->K` is linear with

```text
L=lambda o pi.                                          (12)
```

Then `ker(pi) subset Ann(L)=0`, so `pi` is injective.  An infinite-dimensional
polynomial algebra cannot inject into finite-dimensional `A`.  Therefore no
fixed finite-dimensional multiplicative carrier represents the whole
factorial tower.

## 3. Quantitative finite-state and Prony obstruction

The same Hankel rank gives a representation-free invoice.  Suppose a
`D`-dimensional time-homogeneous linear state satisfies

```text
n!=ell(T^n v),                  0<=n<=2d.                (13)
```

Then

```text
(i+j)!=ell(T^i T^j v),          0<=i,j<=d,               (14)
```

factors `H_d` through the `D`-dimensional state space.  Hence

```text
rank(H_d)<=D.                                            (15)
```

But `(6)` gives `rank(H_d)=d+1`, so

```text
D>=d+1.                                                 (16)
```

In particular, no fixed finite-dimensional state represents `n!` for all
`n`.  Every finite-node Prony sequence and every finite confluent-Prony
sequence has a finite Jordan-state realization, so both are excluded.  The
same argument excludes every fixed-order constant-coefficient recurrence.

Equation `(16)` is stronger and more precise than saying that factorial
growth is too fast: even over an arbitrary characteristic-zero field, before
absolute values exist, the Hankel rank forces the carrier dimension to grow
with the requested height.

## 4. The sharp cyclic alias

For `N>=1`, let

```text
pi_N:K[s] -> K[C_N]=K[z]/(z^N-1),       s |-> z.         (17)
```

For every `a>=1`,

```text
H_(a,N)=s^a(1-s^N)                                      (18)
```

satisfies

```text
pi_N(H_(a,N))=0,                                        (19)
```

while

```text
L(H_(a,N))=a!-(a+N)!=0.                                 (20)
```

Thus every cyclic exponent carrier aliases a genuine factorial value.  At
the LRC-adjacent modulus `N=13` and `a=1`,

```text
L(s-s^14)=-87178291199.                                 (21)
```

The failure survives the standard two-charge Gaussian lift.  With
`s=ZW`, put

```text
P=W+Z s^(a-1)(1-s^N).                                   (22)
```

Under the complex Gaussian Wick functional, only equal total `Z,W` charge
survives, so

```text
E[P^2]=2L(H_(a,N)).                                     (23)
```

For `(a,N)=(1,13)`, this is

```text
-174356582398.                                          (24)
```

The cyclic quotient sees `(18)` as zero and cannot see `(24)`.  This is the
minimal explicit witness for the general quotient obstruction `(10)--(12)`.

## 5. Exact characteristic-`p` boundary

Now let `K` have characteristic `p>0` and define

```text
L_p(s^n)=n! in K.                                       (25)
```

Since `n!=0` for `n>=p`,

```text
(s^p) subset Ann(L_p).                                  (26)
```

On degrees `0,...,p-1`, the factorization `(5)` still holds modulo `p`.
Every diagonal entry `j!`, `j<p`, is a unit, and `P` remains unit lower
triangular.  Hence the `p x p` Hankel block is nondegenerate.

If `q` is not divisible by `s^p`, its residue in `K[s]/(s^p)` is nonzero.
Nondegeneracy supplies `g` of degree below `p` with

```text
L_p(qg)!=0.                                             (27)
```

Combining `(26)--(27)`,

```text
Ann(L_p)=(s^p).                                         (28)
```

Thus characteristic `p` has the exact finite multiplicative carrier

```text
K[s]/(s^p),                                             (29)
```

and the first singular Hankel block is size `p+1`: its last row and column
are zero.  This boundary clarifies why prime-local arguments can create
finite height carriers.

It does not identify `(29)` with THM-2022's proof mechanism.  That theorem
normalizes a lowest balanced face and uses Frobenius/Lucas whole-layer
preservation before returning to characteristic zero; it does not assert
that the original characteristic-zero factorial functional factors through
`(29)`.

## 6. Same support is not enough

As a sharp control, take

```text
H_-=6s-s^3,                         H_+=6s+s^3.          (30)
```

They have identical exponent support, but

```text
L(H_-)=0,              L(H_-^2)=504,              L(H_+)=12. (31)
```

Their two-charge lifts satisfy

```text
(E[P_-^2],E[P_-^4])=(0,3024),
E[P_+^2]=24.                                            (32)
```

Thus support, periodic residue, or a one-bit sign-blind carrier cannot
recover factorial nullity.  This is a control, not the headline novelty:
THM-2173 already gives the general sparse projective moment-floor mechanism.
The new content here is the zero annihilator, quotient classification,
finite-state rank invoice, and characteristic boundary.

## 7. Consequence for the effective GMC frontier

The connection contract is:

| item | exact content |
|---|---|
| source | radial-height polynomial before applying `L` |
| proposed shortcut | one fixed finite quotient/state/periodic carrier |
| first failure | `Ann(L)=0`, equivalently unbounded Hankel rank |
| minimal witness | `s^a(1-s^N)` for a cyclic carrier |
| sharp positive boundary | `Ann(L_p)=(s^p)` in characteristic `p` |
| surviving route | growing height, clock dependence, nonlinear ideals, or Frobenius whole-face preservation |

This does **not** refute HYP-8765's proposed `(k-1)R` cutoff or pair-radical
claim.  A carrier whose dimension grows with that finite cutoff is compatible
with `(16)`.  So is the clock-dependent recurrence

```text
(n+1)!= (n+1)n!,                                       (33)
```

which is not time-homogeneous.  THM-2639's nonlinear resultant/cumulant ideal
also lies outside the linear-state no-go, as does THM-2022's prime-local
whole-face mechanism.

Finally, `(10)` is faithfulness of the standalone radial functional.  It
does not separate Wick channels inside one scalar Gaussian moment and does
not license global termwise de-factorialization.  The MISTAKE-211 and
MISTAKE-215 boundaries remain intact.

## 8. Exact companion

Run

```bash
python 04-computation/gmc_factorial_hankel_faithfulness_thm2810.py
python -O 04-computation/gmc_factorial_hankel_faithfulness_thm2810.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_factorial_hankel_faithfulness_thm2810.out.
```

The companion verifies `D P P^T D` and the determinant formula through
degree ten, exact Bareiss divisibility, `96` cyclic aliases, the
`N=13,a=1` two-charge hostile, and the same-support control.  For
`p=2,3,5,7,11,13` it checks the full-rank `p x p` block, rank of the next
block, zero terminal row/column, and the determinant residue.  The universal
annihilator, quotient, and state-rank statements are proved above rather
than inferred from these finite controls.  The script has explicit exception
gates, no truth-bearing Python assertions, no floating point, and no scratch
dependency.

**Awaiting independent audit; not QED.**
