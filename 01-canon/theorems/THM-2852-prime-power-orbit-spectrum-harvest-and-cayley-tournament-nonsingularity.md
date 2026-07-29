---
id: THM-2852
title: "Prime-power orbit-spectrum harvest and Cayley-tournament nonsingularity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The normalized
  THM-2807 three-address mask has all 13^6 Fourier characters and full
  regular rank.  Each of THM-2829's six borrow-aware ancestry-arrow masks,
  including q11/h9, has all 13^5 characters and full regular rank.  More
  generally every Cayley tournament on an odd finite p-group has
  nonsingular adjacency.  These are convolutional algebra statements;
  they do not manufacture a physical translation action or a positive
  inverse.
source: root/prime-power-orbit-spectrum-harvest-2026-07-28
depends_on:
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
  - THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary
  - THM-2829-q11-semantic-reselection-and-fine-ancestry-phase-obstruction
related:
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
script: 04-computation/lrc14_prime_power_orbit_spectrum_harvest_thm2852.py
output: 05-knowledge/results/lrc14_prime_power_orbit_spectrum_harvest_thm2852.out
script_sha256: 6d2e8f56c474f3295bbfdd104e0dd25e10252d50ec01c9781e713e17dddfd31e
output_sha256: bebab55c4f9093d9b4ecb23cd8d514366ae90bd5ed49686eaed02ff8e7444d32
hash_basis: LF-normalized bytes
---

# THM-2852 -- prime-power orbit-spectrum harvest

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The group-algebra mechanism

Let `G` be a finite `p`-group and let

```text
A=sum_(g in G) a_g g in Z_(p)[G].
```

Write `epsilon(A)=sum_g a_g`.  THM-2839 proves

```text
A is a unit in Z_(p)[G] iff epsilon(A) is a p-unit.      (1)
```

Indeed, over `F_p` the augmentation ideal is nilpotent.  Thus
`A=s*1+n`, with `s=epsilon(A)!=0` and nilpotent `n`, is invertible.
For the regular multiplication matrix,

```text
det R_A=s^|G|=s                    in F_p.               (2)
```

The last equality is Frobenius because `|G|` is a power of `p`.
For cyclic `G`, `(1)` says equivalently that `A` has no zero at any
complex `|G|`-th root of unity and that its cyclic translates form a
basis of the rational regular module.

## 2. The THM-2807 three-address simplex

Put

```text
N=13^6=4,826,809,

n_0=3,454,614,
n_+=3,454,627,
n_a=4,143,978.
```

THM-2807 proves that its selected `tau=12` cells are three literal whole
cylinders with one common nonzero coefficient

```text
c=790,161,473,087,466,480.
```

Normalize away this characteristic-zero scalar and form

```text
H=Z^(n_0)+Z^(n_+)+Z^(n_a)
 =Z^(n_0)(1+Z^13+Z^689364)       in Z[C_N].             (3)
```

Since

```text
epsilon(H)=3!=0 mod13,                                  (4)
```

equations `(1)--(2)` give

```text
H is a unit in F_13[C_N] and Q[C_N],
H(xi)!=0 for every xi^N=1,
rank_Q Circ(H)=N,
det Circ(H)=3 mod13.                                    (5)
```

Consequently the `N` translates of `H` form a basis of `Q[C_N]`.

The normalization in `(3)` is load-bearing:

```text
c=0 mod13.                                              (6)
```

Thus the physical coefficient vector `cH` has the same characteristic-zero
zero set and rank, but it is not itself a modular unit.  Nor does THM-2807
supply a physical action realizing all `N` formal translations.

## 3. Every borrow-aware THM-2829 arrow is spectrally full

Let

```text
D=13^5=371,293.
```

For a source ancestry flag `S(a)` and target flag `T_r(a)`, THM-2829's
natural lift of the six translation arrows

```text
(q,h)=(0,0),(0,7),(3,10),(3,4),(11,2),(11,9)            (7)
```

has Boolean mask

```text
L_(q,h)(a)
 =S(a)T_(q+h mod13)(a+floor((q+h)/13)).                 (8)
```

The proved whole-cylinder counts, in the order `(7)`, are

```text
|L_(q,h)|
 =(66099,65612,65579,65612,65579,65098)
 =(7,1,7,1,7,7) mod13.                                  (9)
```

Therefore every indicator polynomial

```text
h_(q,h)(Z)=sum_(a in L_(q,h)) Z^a in Z[C_D]             (10)
```

is a unit in `F_13[C_D]` and `Q[C_D]`.  Individually,

```text
h_(q,h)(xi)!=0 for every xi^D=1,
rank_Q Circ(h_(q,h))=D,
det Circ(h_(q,h))=|L_(q,h)| mod13.                      (11)
```

In particular the wrapped `q11/h9` arrow retains all `371,293` ancestry
characters after its compulsory `+1` borrow.

The same conclusion holds separately for the auxiliary masks

```text
raw target:             (66099,65652)=(7,2) mod13,
same-section common:    (66099,65612)=(7,1) mod13,
reverse q11/h9 common:   65619=8 mod13.                 (12)
```

No sum of these ranks is asserted.  A Boolean mask is a diagonal physical
support predicate, whereas `(10)` uses it as a convolution kernel.  Full
convolutional spectrum therefore does **not** identify the required
target-active ancestry translation, contract THM-2847's `E3` realization
horn, or repair its semantic-word/current typing.

## 4. Cayley tournaments on odd `p`-groups

Let `p` be odd, `|G|=p^d`, and choose

```text
S subset G\{1},
S intersect S^(-1)=empty,
S union S^(-1)=G\{1}.                                  (13)
```

This is a Cayley-tournament connection set.  Its group-algebra adjacency
element

```text
A=sum_(s in S)s
```

has

```text
epsilon(A)=|S|=(p^d-1)/2=-1/2 modp.                    (14)
```

Hence `A` is a unit by `(1)`.  If rows are sources and columns targets,

```text
Adj_(g,h)=1_S(g^(-1)h),
[R_A]_(h,g)=1_S(g^(-1)h),
```

so the column matrix of right multiplication is

```text
[R_A]=Adj^T.                                           (15)
```

It follows that

```text
Adj is nonsingular over Z_(p), Q, and C,
det Adj=-1/2 modp.                                     (16)
```

Reversing every arc transposes the matrix and leaves `(16)` unchanged.
For nonabelian `G`, every Wedderburn block

```text
sum_(s in S)rho(s)
```

is invertible.  This matrix-block statement is the correct nonabelian
Fourier conclusion; it is not a scalar-character assertion.

The same argument applies to any integral Cayley multidigraph on a finite
`p`-group whose total connection weight is prime to `p`.

## 5. Sharp boundaries

1. A rational mask can have augmentation one and still vanish:

   ```text
   (1/p)sum_(r=0)^(p-1)Z^(r p^(d-1)).
   ```

2. Before scaling, this subgroup sum has mass `p` and rank `p^(d-1)`,
   not `p^d`.
3. The group order must be a `p`-power.  On `C_6`,
   `1+Z^2+Z^4` has odd mass three but rational rank two.
4. A nonnegative convolution mask with more than one support point has no
   nonnegative inverse.  Full spectrum is not positive transport.
5. Cayley-tournament nonsingularity gives no Hamiltonian-path, extremal,
   path-homology, or LRC row conclusion.

## 6. Exact evidence and independent audit

The exact companion checks the THM-2807 addresses, gaps, determinant
residue, and normalization boundary; all THM-2829 arrow and auxiliary
augmentation residues; all four tournaments on `C_5`; all sixteen on
`C_9`; and one hundred fixed orientations of the nonabelian order-27
Heisenberg group.  The latter controls verify tournament validity,
right-multiplication/adjacency transpose, reversal, and determinant
residue.  It also checks the `C_169` subgroup/rational-scaling hostile and
the `C_6` non-prime-power hostile.

An independent audit rederived `(1)--(2)`, checked the primitive-scalar
typing, reconstructed every cited residue, verified the nonabelian
transpose convention, and replayed all finite controls.  Normal,
optimized, and stored transcripts agree exactly; LF-normalized hashes
match the frontmatter.

**QED.**
