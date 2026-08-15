---
id: THM-3393
title: "Hadamard order 668: explicit normalized certificate and bordered Goethals--Seidel construction"
status: >
  PROVED + VERIFIED-EXACT + GENERATOR-INDEPENDENT CERTIFICATE + INDEPENDENTLY
  AUDITED.  An explicit
  normalized 668-by-668 sign matrix is embedded as compressed row bitsets.
  A standard-library verifier checks its byte hashes, shape, normalization,
  and all 222,778 off-diagonal row distances, each exactly 334, hence
  H H^T=668 I.  An independent direct block reconstruction identifies a
  4-(166;82,83,83,83;164) supplementary difference family and a four-row
  bordered Goethals--Seidel array.  This proves existence at order 668 and
  bypasses, but does not solve, the length-333 Legendre-pair route of
  THM-2833.  It does not prove the full Hadamard conjecture.
source: codex-2026-08-14 independent inert audit of the sign-puzzle certificate
audit: independent normalized-text reconstruction, all-pairs bitset Gram check, SDS/PAF audit, and block-construction line audit by root-2608-crouzeix-puzzle
depends_on: []
related:
  - THM-2833-legendre-333-asymmetric-multiplier-exclusion
  - THM-3394-twelve-formerly-missing-hadamard-orders-through-2000
script: 04-computation/hadamard668_explicit_certificate_thm3393.py
output: 05-knowledge/results/hadamard668_explicit_certificate_thm3393.out
script_sha256: f637f2017954ee949e83871470179d529c69bd7add404c0924b0b4dcf2be030f
output_sha256: 83870dbcbe81ded0c718be1fa3696fb0016dd0fcba3eba1a27d1300e06e148de
hash_basis: working-tree bytes (LF)
---

# THM-3393 -- an explicit Hadamard matrix of order 668

**PROVED + VERIFIED-EXACT + GENERATOR-INDEPENDENT CERTIFICATE + INDEPENDENTLY AUDITED.**

There exists a matrix

```text
H in {+1,-1}^{668 x 668}
```

such that

```text
H H^T = 668 I_668.                                           (1)
```

The supplied representative is normalized: its first row and first column
are all `+1`.  Before this certificate, THM-2833 and its cited July 2026
source treated `668` as the smallest unresolved Hadamard order.  Equation
(1) closes that finite existence question in the repository.  It does not
settle the Hadamard conjecture for arbitrary multiples of four.

## 1. Generator-independent explicit certificate

The companion script embeds the normalized matrix itself, without executing
or interpreting the source puzzle.  Each row is an 84-byte little-endian
integer bitset:

```text
bit j = 0  <=>  H[i,j] = +1,
bit j = 1  <=>  H[i,j] = -1.                                (2)
```

Only bits `0,...,667` are data; the top four padding bits of every row are
required to vanish.  The 56,112 raw bytes are compressed with `zlib` and
encoded with base85.  The immutable hashes are

```text
raw row bitsets: 47e5b8b061401c8cf18c2dee97f581a0c31b3ca7ad76e2930e43f1e4f18b50ca
zlib bytes:      6847f4d26b3e9ad284e57b9284f4297856a38d74065ba2b9c972f6fed70738d8
normalized text: 73f1de1539849e1dc7e6085cc69c563fd2965c44970263e8203384bd1a46aa63
```

Here the text form has one 668-character `+/-` row per LF-terminated line.
The verifier checks all three hashes, rejects truncated or trailing compressed
data, checks all padding bits, and verifies normalization explicitly.

Let `b_i` be the bitset of row `i`.  For two sign rows the exact identity is

```text
<H_i,H_j> = 668 - 2 popcount(b_i xor b_j).                    (3)
```

The verifier enumerates all

```text
binom(668,2) = 222,778
```

unordered distinct row pairs and obtains Hamming distance `334` for every
one.  Thus (3) is zero off the diagonal, while every diagonal inner product
is `668`.  This proves (1).  No floating-point arithmetic, solver, external
package, shell, or decoded program enters this argument.

## 2. Exact construction sidecar

The certificate also has a compact structural reconstruction.  Put `q=166`
and split the embedded 664-sign seed into sequences `a,b,c,d` of length `q`.
For the negative supports

```text
D_x = {j in Z/166Z : x_j=-1},
```

the verifier checks

```text
(|D_a|,|D_b|,|D_c|,|D_d|) = (82,83,83,83),                  (4)

sum_x |D_x intersect (D_x-t)| = 164
for every nonzero t in Z/166Z.                               (5)
```

Thus these supports form a supplementary difference family with parameters

```text
4-(166;82,83,83,83;164).                                    (6)
```

Equivalently, their row sums are `(2,0,0,0)` and their periodic
autocorrelations obey

```text
PAF_a(t)+PAF_b(t)+PAF_c(t)+PAF_d(t) = -4   (t != 0).         (7)
```

Indeed, (5) gives

```text
4q - 4(82+83+83+83) + 4(164) = -4.
```

For the exact block description, let `A,B,C,D` be the ordinary circulant
matrices with first rows `a,b,c,d`, and let `R` be the order-166
anti-identity.  Form the standard Goethals--Seidel core

```text
G = [ A    BR     CR      DR   ]
    [-BR    A    D^T R  -C^T R ]
    [-CR  -D^T R   A     B^T R ]
    [-DR   C^T R -B^T R    A   ].                            (8)
```

Use the following three sign matrices of order four:

```text
T = [- + + -]       P = [- - - +]       Q = [+ + + -]
    [+ - + -]           [- - + -]           [- - + -]
    [+ + - -]           [- + - -]           [- + - -]
    [- - - -],          [+ - - -],          [+ - - -].      (9)
```

With `1_q` the all-one column, the unnormalized matrix is

```text
H_0 = [       T          P tensor 1_q^T ]
      [ Q tensor 1_q             G       ].                  (10)
```

This is the precise construction, not an analogy.  The standard
Goethals--Seidel multiplication and (7) give

```text
G G^T = I_4 tensor ((4q+4)I_q - 4J_q).                       (11)
```

Directly from (9),

```text
T T^T=P P^T=Q Q^T=4I_4,             T Q^T=-2P.              (12)
```

The bottom-left border contributes `4 I_4 tensor J_q`, cancelling the
`-4J_q` term in (11).  The top-right border and `T` give
`4qI_4+4I_4`.  Finally, the seed row sums `(2,0,0,0)` make the core part of
the top-bottom product `2P tensor 1_q^T`, while (12) makes the border part
its negative.  Hence all three Gram blocks of (10) equal those required by

```text
H_0 H_0^T = 4(q+1) I = 668 I.                               (13)
```

The script independently constructs (10) from the seed, checks its canonical
text hash

```text
bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0,
```

normalizes it by row and column sign switches, and requires byte-for-byte
equality with the explicit bitset certificate.  Thus the structural and
matrix-level verification paths meet at a frozen object.

## 3. Relation to the length-333 frontier

THM-2833 excluded a broad multiplier-invariant class of **Legendre pairs of
length 333** and explicitly left unrestricted Hadamard order `668` open.  The
present construction uses four sequences of length `166` plus a four-row and
four-column border.  Therefore it bypasses the Legendre-pair route; it neither
constructs a length-333 Legendre pair nor contradicts THM-2833's exclusion.
There is no justified inverse implication from a Hadamard matrix of order
`668` to a Legendre pair of length `333`.

The external status statement being superseded is Arthur F. Ramos, David B.
Hulak, and Ruy J. G. B. de Queiroz, *Multiplier obstructions for Legendre
pairs of length 333*,
[arXiv:2607.20765](https://arxiv.org/abs/2607.20765), whose abstract called
`668` the smallest presently unresolved order on 2026-07-22.  This citation
sets the prior boundary; it is not a dependency of the certificate proof.

## 4. Failure boundary and reproducibility

- A one-bit corruption in the first row or column fails normalization.  A
  one-bit corruption elsewhere changes that row's distance from the all-plus
  first row from `334` to `333` or `335`.  Thus the cheapest hostile is caught
  before any subtle block audit.
- Hashes, decoded length, padding, alphabet, normalization, all row pairs,
  all 165 nonzero SDS shifts, and construction/certificate equality are
  checked with explicit exceptions rather than Python `assert`.
- Ordinary and `python -O` runs are byte-identical.  The frozen output ends in
  `ALL CHECKS PASSED`.
- Row and column sign switches and permutations produce equivalent Hadamard
  matrices; the theorem claims this explicit normalized representative, not
  an equivalence-class classification, uniqueness, or a priority result.
- The other formerly unresolved orders below 2000 are separate certificates
  routed through THM-3394.  They are not needed for (1).

Reproduce with

```bash
python3 04-computation/hadamard668_explicit_certificate_thm3393.py
python3 -O 04-computation/hadamard668_explicit_certificate_thm3393.py
```
