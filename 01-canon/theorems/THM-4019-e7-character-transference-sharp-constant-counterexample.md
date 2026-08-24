---
id: THM-4019
title: "E7 counterexample to sharp arbitrary-lattice character transference"
status: >
  PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For the standard
  E7 Cartan Gram matrix and its unique nonzero mod-two Gram-radical
  characteristic, the primal parity minimum is 6 and the odd dual minimum is
  3/2, so their product is 9>7.  An integral orthogonal-block family gives
  product d+2r>d in every rank d=7r+k with r>=1.  Thus THM-4015 equation (27),
  the sharp arbitrary-lattice extension with constant sqrt(d)/2, is false in
  every dimension at least seven.  THM-4015's first-kind theorem survives,
  and the weaker arbitrary-lattice d/2 bound remains OPEN.
source: character_transference_high_rank_probe / all-frontiers session, 2026-08-24
audit: >
  PASS.  A non-LDL exact verifier proves positive definiteness, finds the
  unique mod-two radical, enumerates all 126 E7 roots in a Cauchy-certified
  box, checks the complete 3^7 dual box, and verifies the block family through
  rank thirty at 2,515 optimization-stable gates.  An independent rational
  LDL sphere atlas checks all characteristics of A6,D6,E6,A7,D7,E7,A8,D8,E8,
  A9,D9 plus 1,152 seeded exact controls, with no cutoff among 3,509 instances
  and 518,300 recursive nodes.  Normal and optimized transcripts agree with
  the frozen LF outputs.
depends_on: []
related:
  - THM-4009-euclidean-covering-transference-short-relation-compression
  - THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression
  - THM-4015-first-kind-character-sensitive-foster-transference
script: 04-computation/e7_character_transference_counterexample_thm4019.py
output: 05-knowledge/results/e7_character_transference_counterexample_thm4019.out
independent_script: 04-computation/character_transference_rank6_9_exact_atlas_thm4019.py
independent_output: 05-knowledge/results/character_transference_rank6_9_exact_atlas_thm4019.out
script_sha256: f15fcd5628e2db1f1fd412a7dbf099c3bbd705d1d1e507a9815ac77896df437d
output_sha256: 6e4b5912a9695193548ec03b710c95c2c38ab7e53b2ad55ebe792d13ad403c49
independent_script_sha256: 70c30da55ede8c951f27661fc9f029e4e64b45ed6bebc2f1f25b5966c7014e91
independent_output_sha256: 2369d2e606e8eccf37d7245b18a638b46f19bed539714b627699330578cb4706
semantic_sha256: 5498c181601ed78f857d86f40f466fc606f6b87e51e2fe2d2a192c591c440d66
hash_basis: raw LF bytes
---

# THM-4019 -- the E7 odd-radical counterexample

**PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  This theorem
refutes only the sharp arbitrary-lattice candidate isolated in THM-4015.  It
does not weaken THM-4015 on lattices of Voronoi's first kind.

## 1. Character coordinates

For a positive-definite integral Gram matrix `G` of rank `d` and nonzero
`u in F_2^d`, put

```text
A(G,u)=min_(x=u mod 2) x^T G x,
B(G,u)=min_(u dot z=1 mod 2) z^T G^(-1)z.             (1)
```

If `p` is a shortest vector in parity class `u` and
`c=p/2 mod L`, then

```text
dist(c,L)lambda_odd=(1/2)sqrt(A(G,u)B(G,u)).           (2)
```

Thus THM-4015 equation (27), the proposed sharp extension from first-kind to
arbitrary lattices, is exactly

```text
A(G,u)B(G,u)<=d.                                      (3)
```

## 2. The intrinsic E7 characteristic

Take

```text
G = [[ 2,-1, 0, 0, 0, 0, 0],
     [-1, 2,-1, 0, 0, 0, 0],
     [ 0,-1, 2,-1, 0, 0,-1],
     [ 0, 0,-1, 2,-1, 0, 0],
     [ 0, 0, 0,-1, 2,-1, 0],
     [ 0, 0, 0, 0,-1, 2, 0],
     [ 0, 0,-1, 0, 0, 0, 2]],
u=(0,0,0,1,0,1,1).                                   (4)
```

The leading principal determinants are

```text
2,3,4,5,6,7,2,                                       (5)
```

so `G` is positive definite of determinant two.  Direct reduction modulo two
gives

```text
ker(G mod 2)={0,u}.                                   (6)
```

The hostile characteristic is therefore the unique nonzero Gram radical,
not a coordinate-selected accident.  In particular

```text
Gu/2=(0,0,-1,1,-1,1,1) in Z^7,
u^TGu=6,
u dot (Gu/2)=3.                                       (7)
```

## 3. The primal minimum

Every `x=u+2k` satisfies

```text
x^TGx=u^TGu=2 mod 4.                                  (8)
```

Positive values below six could therefore only equal two.  If
`x^TGx<=2`, Cauchy's inequality gives

```text
x_i^2<=(x^TGx)(G^(-1))_ii.                            (9)
```

The exact inverse diagonal puts every such vector in the coordinate box with
caps `(2,3,4,3,2,1,2)`.  Exact enumeration finds all 126 E7 roots of norm two
and none is congruent to `u` modulo two.  Hence (8) gives a lower bound of six,
while `u` itself attains six:

```text
A(G,u)=6.                                             (10)
```

## 4. The odd dual minimum

Let `H=G^(-1)`.  If `z^THz<3/2`, dual Cauchy with `H^(-1)=G` and
`G_ii=2` forces `z_i^2<3`; every possible shorter integral dual vector lies
in `{-1,0,1}^7`.  The complete exact box has no nonzero vector below `3/2`
and 56 vectors at `3/2`.  The explicit vector

```text
z=(0,0,1,-1,0,0,0),
u dot z=-1,
z^TG^(-1)z=3/2                                       (11)
```

is odd, as is `Gu/2` from (7).  Therefore

```text
B(G,u)=3/2.                                           (12)
```

Combining (10) and (12),

```text
A(G,u)B(G,u)=9>7,
dist(c,L)lambda_odd=3/2>sqrt(7)/2.                    (13)
```

This also proves, by THM-4015, that this E7 lattice is not of Voronoi's first
kind.

## 5. Integral amplification in every rank at least seven

For `r>=1`, `k>=0`, and `d=7r+k`, define

```text
G_(r,k)=diag((3G)^(direct sum r),2I_k),
u_(r,k)=(u,...,u,1,...,1).                            (14)
```

The primal parity minimum separates by orthogonal blocks.  Each scaled E7
block contributes `3*6=18`, and each trailing odd Gram-two coordinate
contributes two:

```text
A(G_(r,k),u_(r,k))=18r+2k.                            (15)
```

For the dual minimum, total odd parity is the sum of the block parities.
Every odd scaled-E7 block costs at least `(3/2)/3=1/2`, and every odd trailing
coordinate costs `1/2`.  At least one block must be odd, while choosing one
odd block and zero elsewhere attains the lower bound.  This remains valid
when `k=0`.  Hence

```text
B(G_(r,k),u_(r,k))=1/2,
AB=9r+k=d+2r>d.                                       (16)
```

Taking `r=1,k=d-7` supplies an integral counterexample in every rank
`d>=7`.  Taking `r=floor(d/7)` gives the displayed block-family maximum
`AB=d+2floor(d/7)`.

## 6. Exact boundary and surviving questions

The independent exact atlas gives the named-root maxima

```text
A6:36/7, D6:6, E6:16/3,
A7:7,    D7:7, E7:9,
A8:64/9, D8:8, E8:8,
A9:9,    D9:9.                                       (17)
```

It also checks 1,152 seeded integral controls without a skipped enumeration.
These are **FINITE-EXACT** data, not a classification.  The sharp inequality
is PROVED for every lattice through rank three by THM-4015, remains OPEN in
ranks four through six, and is REFUTED in every rank at least seven.

The weaker arbitrary-lattice candidate

```text
dist(c,L)lambda_odd<=d/2,
equivalently A(G,u)B(G,u)<=d^2,                        (18)
```

remains **OPEN**.  The family (16) does not approach its quadratic threshold.
For LRC(14), (13)--(16) close a method only: the first-kind hypothesis in
THM-4015 cannot simply be discarded.  They say nothing adverse about a
particular hypothetical LRC counterexample lattice.

## 7. Reproduction

Run

```text
python -B 04-computation/e7_character_transference_counterexample_thm4019.py
python -B -O 04-computation/e7_character_transference_counterexample_thm4019.py
python -B 04-computation/character_transference_rank6_9_exact_atlas_thm4019.py
python -B -O 04-computation/character_transference_rank6_9_exact_atlas_thm4019.py
```

and compare with the two frozen outputs named in the frontmatter.  Both
normal/optimized pairs agree after platform-newline normalization. **QED.**
