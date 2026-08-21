---
id: THM-3658
title: "LRC mod-169 carry-Fourier block intertwiner"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  On the pinned field, the split character basis of two separate mod-13
  digits and the cyclic character basis of their assembled mod-169 residue
  are joined by thirteen invertible 13-by-13 blocks.  The zero high-frequency
  block is a permutation; every entry of every nonzero high-frequency block
  and its inverse is nonzero.  The transform exactly intertwines branch
  reversal.  The assembled-label bijection does not intertwine the two group
  laws or their natural convolutions, as witnessed by one carry.  This is a
  linear character-basis bridge, not a physical-current chain map or LRC(14).
source: kps-s191 / THM-3657 carry-sidecar continuation, 2026-08-21
depends_on:
  - THM-3657-lrc-two-current-quotient-address-atlas-and-reversal-gate
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
script: 04-computation/lrc_mod169_carry_fourier_block_intertwiner_thm3658.py
output: 05-knowledge/results/lrc_mod169_carry_fourier_block_intertwiner_thm3658.out
script_sha256: a298cb764b269d0519883d5965797592cd561b9b392823a4a454ca9fa8870f8d
output_sha256: 98028e209dbbed0022e6dadb1eb246be3f8b2f501094bb76094de3e7e9615a58
semantic_sha256: 557033432ed61af275e58c73d3015346594f5d0dca67c1df3f2643e963e7ab3e
hash_basis: raw LF bytes
---

# THM-3658 -- the two mod-13 digits have an exact carry-Fourier bridge

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  This theorem constructs the lawful linear change of character
basis suggested by THM-3657's carry sidecar.  It simultaneously records why
that change of basis is not yet a current or convolution intertwiner.

## 1. Two group structures on one 169-element set

Work over the pinned field

```text
F=F_p,                  p=755373809845391722745761.     (1)
```

The exact elements

```text
omega=298763986285447441216949,
zeta =123453826432109539554819=omega^13              (2)
```

have orders `169` and `13`, respectively.  Let

```text
S={(r0,r1):r0,r1 in F13}.                             (3)
```

There are two relevant structures on `S`:

```text
split:   F13 direct_sum F13,
cyclic:  a=r0+13 r1 in Z/169Z.                        (4)
```

The assembly in `(4)` is a set bijection, not a group isomorphism.  Both
structures nevertheless give character bases of the same 169-dimensional
function space `Fun(S,F)`.

The split and cyclic characters are

```text
phi_(u,v)(r0,r1)=zeta^(u r0+v r1),
psi_k(r0,r1)=omega^(k(r0+13 r1)),                     (5)
```

for `(u,v) in F13^2` and `k in Z/169Z`.

## 2. The transform has thirteen exact blocks

Write uniquely

```text
k=v+13h,                 0<=v,h<13.                   (6)
```

Then

```text
psi_k(r0,r1)=omega^(k r0) zeta^(v r1),                (7)
```

so only split characters with second frequency `v` occur.  Fourier inversion
in the low digit gives

```text
psi_(v+13h)=sum_(u=0)^12 C_v(u,h) phi_(u,v),          (8)

C_v(u,h)=1/13 sum_(r=0)^12 omega^((v+13h)r)zeta^(-ur). (9)
```

Thus the full 169-by-169 change matrix is the direct sum of thirteen
13-by-13 blocks `C_v`.

For `v=0`, equation `(9)` is ordinary character orthogonality:

```text
C_0(u,h)=1_(u=h).                                     (10)
```

For `v!=0`, put `q=omega^(v+13h)zeta^(-u)`.  Since

```text
q^13=zeta^v!=1,
```

the geometric sum gives

```text
C_v(u,h)=(1-zeta^v)/(13(1-q))!=0.                    (11)
```

Every one of the 169 entries in each nonzero block is therefore nonzero.
The columns of `C_v` are the ordinary 13-point Fourier basis multiplied
pointwise by the nowhere-zero function `r -> omega^(vr)`, so every `C_v` is
invertible.  Applying the same geometric-sum argument to the inverse change
shows that every inverse entry is also nonzero for `v!=0`.

The exact nonzero-count profile is consequently

```text
(13,169,169,169,169,169,169,169,169,169,169,169,169), (12)
```

and all thirteen block ranks equal `13`.  In particular, a split character
with nonzero high-digit frequency spreads across all thirteen cyclic
frequencies in its congruence block, and conversely.  Carry does not mix the
high frequency `v`; it completely mixes the low frequency inside that block.

## 3. Branch reversal is exactly intertwined

THM-3657's simultaneous branch reversal is

```text
j(r0,r1)=(12-r0,12-r1),
a -> -a-1 mod 169.                                    (13)
```

Pullback by `j` acts on the two bases as

```text
j^* phi_(u,v)=zeta^(-u-v) phi_(-u,-v),
j^* psi_k    =omega^(-k) psi_(-k).                    (14)
```

Substitution of `(9)` proves that the block transform intertwines these two
actions.  The companion verifies all 169 coefficient-vector identities,
not just evaluations or orbit counts.  This is the character-level lift of
the reversal that governs THM-3657's `37+124+8` address atlas.

## 4. The carry hostile and the exact boundary

The set bijection `(4)` does not preserve addition.  The minimal witness is

```text
(12,0)+(1,0)=(0,0)             in the split group,
12+1=13, whose cyclic digits are (0,1), in Z/169Z.     (15)
```

Hence the identity on assembled labels does not intertwine the two natural
convolution products.  More abstractly, `Z/169Z` has exponent `169`, while
`F13^2` has exponent `13`; every group homomorphism from the former to the
latter has image dimension at most one.

There is no contradiction with `(8)`: a linear change between two character
bases of one function space need not arise from a group isomorphism of their
labels.  Equation `(8)` is exactly the carry-sensitive replacement for that
false identification.

THM-2334's 169 target twists form the dual of a split `F13^2` quotient.  If a
typed identification of that split quotient with the digit-frequency plane
is supplied, `(8)` gives an exact finite change to cyclic mod-169 characters.
This theorem does **not** supply that identification, transport THM-2334's
orbit coefficients, or prove that either character basis is realized by a
physical two-current chronology.

## 5. Exact verification and scope

Reproduce with

```bash
python3 -B 04-computation/lrc_mod169_carry_fourier_block_intertwiner_thm3658.py
python3 -B -O 04-computation/lrc_mod169_carry_fourier_block_intertwiner_thm3658.py
```

The assertion-free companion source-pins THM-3657's exact script/output,
checks both root orders, constructs and inverts all thirteen blocks, verifies
`28,561` evaluations in each direction, all 169 reversal identities, every
density/rank claim, and the carry hostile `(15)`.  Normal and optimized
streams are byte-identical.  The block, inverse, and semantic digests are

```text
block:    7cebcaa457bf9d1ffa9cb43582bb3a07c02b8849c2b675500dca65ea9d0ac1e7,
inverse:  d3d6363bdbf39b41812e8b93010d66c52fd7c7056f29b187fdbbb26824ab8da5,
semantic: 557033432ed61af275e58c73d3015346594f5d0dca67c1df3f2643e963e7ab3e. (16)
```

This is an exact finite-field linear character-basis theorem.  It proves no
current-to-digit chain map, chronology, physical entry, admissible
coefficient law, characteristic-zero transport, row exclusion, or LRC(14).
**QED.**
