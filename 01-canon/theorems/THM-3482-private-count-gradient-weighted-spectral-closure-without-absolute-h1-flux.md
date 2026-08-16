---
id: THM-3482
title: "Private-count gradient weighted spectral closure without absolute H1 flux"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the canonical
  THM-3473 owner order, orient the thirteen coactivity edges increasingly
  and weight each edge by the difference of its endpoint private-sheet
  counts.  This edge cochain is a coboundary and pairs trivially with every
  graph cycle, yet its weighted reduced-Laplacian determinant is strictly
  negative for every k>=1, with the three displayed residue-class
  factorizations.  This is graph-level, orientation-gauged spectral closure,
  not an LRC current or LRC(14) bispectrum theorem.
source: codex-2026-08-15-private-gradient-spectral-closure
audit: >
  self-contained polynomial arithmetic; exact enumeration of all sixteen
  spanning trees in each K4; symbolic derivation in all three residue states;
  six symbolic cycle-pairing zero checks; direct signed 7x7 Bareiss
  determinants for k=1..3000; constant-potential and killed-bridge hostiles;
  independent Prüfer and 256-minor Cauchy-Binet reconstruction; independent
  six-triangle cycle basis, typed hub-gauge repair, k=1 orientation-reversal
  hostile, and determinant-does-not-descend-to-H1 witness; dependency,
  semantic, security, ID, docs, and normal/optimized/stored replay gates
depends_on:
  - THM-3473-three-times-p-eight-owner-private-sheet-partition-and-irredundancy
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-3472-odd-modulus-zero-half-conjugacy-and-global-zmc-rank-equality
script: 04-computation/lrc_private_count_gradient_weighted_spectral_thm3482.py
output: 05-knowledge/results/lrc_private_count_gradient_weighted_spectral_thm3482.out
script_sha256: 3a3f439a88abe5a14180e850e016db2a0b6b327fb979b5df0106ee7aaa2cbbac
output_sha256: 6632ead82d4dcdbac3919200aa8790403a4b2d61a2da5d29893f208a08b03db9
semantic_sha256: 98cba04620048b6c6f8fab03518ab39e0125623f033423eda29245c9b53a0162
hash_basis: LF-normalized bytes
---

# THM-3482 -- a coboundary can be spectrally nondegenerate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The proof, exact companion, independent signed matrix-tree reconstruction,
and hostile controls pass their stated gates.  The audit repaired one typed
pairing in `(7)`; MISTAKE-409 records the correction.

## 1. Source, target, map, and orientation gauge

Use the THM-3473 owner order and coactivity packets

```text
A={1,4,5,6},       B={2,3,5,8},       C={5,7}.       (1)
```

Their two-section is two `K4`s and a bridge, glued at owner five.  Orient
every edge from its smaller to its larger owner index.  Let `B` be the
reduced `7x13` incidence matrix obtained by deleting the owner-five row, and
let `B_full` be the full `8x13` incidence matrix.

For an owner potential `f in Q^8`, define

```text
Phi([f])=B_full^T f in Q^13,                         (2)
```

where `[f]` denotes the class modulo constant potentials.  Thus the source,
target, and map are

```text
Q^8/<(1,...,1)>  --Phi-->  Q^13,
relative owner potential          signed edge response.             (3)
```

The preserved predicate in this theorem is nonvanishing of

```text
S([f])=det(B diag(Phi([f])) B^T).                    (4)
```

The map loses sheet locations, support multiplicities beyond the selected
potential, and absolute graph flux.  Its required sidecars are the labelled
owner order and the ternary state `k mod 3`.

The orientation is genuine data.  Reversing one edge reverses its current
coordinate but not the rank-one matrix `b_e b_e^T`, so it can change `(4)`.
The theorem therefore concerns the canonical labelled orientation, not an
unoriented graph invariant.  Taking absolute values would remove this gauge
but would also discard the sign/phase information needed by a future current
transplant.

## 2. The canonical private-count potential

Put `e_r(k)=1` if `k==r (mod 3)` and zero otherwise.  THM-3473 proves that
the eight private-sheet counts are

```text
f_k=(4k,
     4k-2e_1,
     4k-2e_0,
     4k,
     8k-2e_2,
     4k,
     4k-2e_2,
     4k).                                            (5)
```

Set `w_k=Phi([f_k])`.  Write the deleted-hub representative as

```text
fbar_k=(f_(k,i)-f_(k,5))_(i!=5) in Q^7.             (5a)
```

Then `w_k=B^T fbar_k`.  In particular, it represents zero in

```text
H^1(G;Q)=Q^13/im(B^T).                               (6)
```

Equivalently, every cycle `c in ker(B)` satisfies

```text
<c,w_k>=<Bc,fbar_k>=0.                               (7)
```

This vanishing is an identity in `k`, not a numerical accident.

## 3. The weighted determinant factors as 3+3+1

For arbitrary commutative signed edge weights, the matrix-tree identity and
the cut-vertex decomposition give

```text
det(B diag(w) B^T)=w_(5,7) Tau_A(w) Tau_B(w),         (8)
```

where each `Tau` is the sixteen-term spanning-tree polynomial of its `K4`.
The degrees are

```text
3 + 3 + 1 = 7,                                      (9)
```

while the carrier still has thirteen edge coordinates.  In particular,
scaling the owner potential by `lambda` scales `(4)` by `lambda^7`.

On packet `A`, owners `1,4,6` all have potential `4k`, while owner `5` has
potential `8k-2e_2`.  The three edges among the equal outer potentials have
weight zero.  Exactly the hub star survives in the tree sum, with one edge
oppositely oriented, so

```text
Tau_A(w_k)=-(4k-2e_2)^3.                             (10)
```

The bridge is even simpler:

```text
w_(5,7)=(4k-2e_2)-(8k-2e_2)=-4k.                   (11)
```

For packet `B`, the potentials in owner order `(2,3,5,8)` and the exact
sixteen-tree sums are

| `k mod 3` | potentials | `Tau_B(w_k)` |
|---:|---|---|
| `0` | `(4k,4k-2,8k,4k)` | `-8(2k-1)^2(2k+1)` |
| `1` | `(4k-2,4k,8k,4k)` | `-8(8k^3+12k^2-2k-1)` |
| `2` | `(4k,4k,8k-2,4k)` | `-8(2k-1)^3` |

These identities follow by direct expansion of the sixteen spanning-tree
monomials.  The companion derives their coefficient vectors using its own
integer polynomial arithmetic; it does not fit them from sampled values.

## 4. Universal nonvanishing

Combining `(8)--(11)` gives

```text
det(B diag(w_k) B^T)

 = -2048 k^4 (2k-1)^2 (2k+1),
       k==0 (mod 3),

 = -2048 k^4 (8k^3+12k^2-2k-1),
       k==1 (mod 3),

 = -256 k (2k-1)^6,
       k==2 (mod 3).                                (12)
```

Every factor in the first and third lines is positive before the leading
minus sign.  In the middle line,

```text
g(k)=12k^2-2k-1
```

has `g(1)=9` and `g(k+1)-g(k)=24k+10>0`; hence the remaining cubic is
positive for every `k>=1`.  Therefore all three lines of `(12)` are strictly
negative on their lawful lanes.  This proves

```text
det(B diag(w_k) B^T) != 0        for every k>=1.     (13)
```

There are nevertheless four identically zero edge weights in states zero
and one, and six in state two.  Thus all-edge support is not necessary
either.  The first three exact determinants are

```text
k=1:       -34,816,
k=2:      -373,248,
k=3:   -29,030,400.                                 (14)
```

The ternary subsets retain THM-3473's natural and harmonic coefficient
`1/3`; `(12)` shows that no residue lane is a spectral zero lane.

## 5. The corrected logical boundary

Equations `(7)` and `(13)` hold simultaneously:

```text
absolute graph H1 class of w_k:       zero,
owner-order weighted determinant:     nonzero.       (15)
```

The reason is type-theoretic.  The first predicate is a quotient-linear
statement about an edge cochain modulo gradients.  The second inserts the
thirteen coordinates diagonally into a degree-seven matrix-tree polynomial.
That polynomial does not descend to `H^1(G;Q)`.

Consequently a phase/holonomy sidecar is necessary if the target is nonzero
absolute graph flux, but it is not necessary for graph-level weighted
spectral nonvanishing.  MISTAKE-409 records the former overstatement.

The sharp hostiles are retained.  A constant owner potential kills every
edge and makes `(4)` zero.  Equalizing owners five and seven kills the forced
bridge and therefore the whole determinant, regardless of both tetrahedra.

## 6. Scope and exact companion

Run from the repository root:

```bash
python 04-computation/lrc_private_count_gradient_weighted_spectral_thm3482.py
python -O 04-computation/lrc_private_count_gradient_weighted_spectral_thm3482.py
```

The standard-library companion enumerates the sixteen trees in each `K4`,
derives every symbolic polynomial in all three states, checks all six cycle
pairings identically, and compares the factorization with direct signed
`7x7` Bareiss determinants for `1<=k<=3000`.  It includes constant-potential
and killed-bridge controls, dependency and security gates, and a frozen
semantic digest.

This theorem does **not** map a THM-2334 relation current to the thirteen
edges, identify a physical phase, prove LRC `7 tensor 13` bispectrum
nonvanishing, exclude a scalar row, or prove LRC(14).  It supplies a concrete
spectral target and proves that absolute graph holonomy is not a prerequisite
for that target.
