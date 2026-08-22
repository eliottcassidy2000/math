---
id: THM-3345
title: "Prime-XOR ancestry path groupoid and source-dependent Berggren cost"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the first
  spine-selected affine V4 fibre, the two edges of every prime-XOR matching
  have different Berggren ancestry costs and reduced path labels, so neither
  the prime nor its folded weight defines a source-independent ancestry word.
  The exact endpoint paths instead form a flat source-dependent groupoid;
  ambient-tree holonomy and graph H^1 vanish.  A prime-5 branch transplant
  proves simultaneous unbounded Boolean rank and ancestry dispersion and is a
  rational transduction on one unary subbranch.  Rooted automorphisms and
  fixed descendant words are excluded, while a uniform/global source-reading
  transducer across all Boolean fibres, LRC transport, and JC flux remain open.
source: codex-2026-08-12-prime-xor-ancestry-groupoid
audit: >
  An independent hostile audit rederived all six longest-common-prefix
  lengths, path distances, and depth jumps; checked the constant-label,
  rooted-automorphism, flat-square, ambient-holonomy, frozen-basepoint, and
  finite-state scope boundaries.  A second scale audit reconstructed the
  eight-parent c=99905 fibre and all 28 edges.  A third independent audit
  verified the infinite prime-5 matrix transplant, exact 5-adic scope, CRT
  quantifiers, and unary rational transduction.  Normal, optimized, and stored
  transcripts match after LF normalization, and both recorded hashes match.
depends_on:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
related:
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
script: 04-computation/prime_xor_ancestry_groupoid_thm3345.py
output: 05-knowledge/results/prime_xor_ancestry_groupoid_thm3345.out
script_sha256: 41157d468658f2836af3a22435959d515ab230866099b93311fe261139bf07f8
output_sha256: a1735c1cd9a503d9db74c35b04318ce56a770592b9f8429bdac8b6b0e78b17f2
hash_basis: working-tree bytes (LF)
---

# THM-3345 -- prime-XOR ancestry paths are source-dependent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

No literature-priority or global-novelty claim is made.  The payload is a
typed comparison of the external Gaussian collision torsor with paths in the
ambient Berggren tree.

## 1. Inheritance and the two carriers

[THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
proves that the first spine-selected affine `V4` collision fibre occurs at

```text
c=1105=5*13*17.
```

Choose the factor-allocation chart

```text
X={000,100,010,001} isomorphic to F_2^3/<111>,            (1)
```

where the three singleton directions are labelled by `5,13,17`.  The four
parents and their unique root-to-child Berggren addresses are

| bit | parent | word | depth |
|---|---|---|---:|
| `000` | `(943,576,1105)` | `DAUD` | 4 |
| `100` | `(817,744,1105)` | `UDUA` | 4 |
| `010` | `(1073,264,1105)` | `DADDD` | 5 |
| `001` | `(47,1104,1105)` | `U^22` | 22 |

The external `K4` joins equal-hypotenuse parents.  It is not a subgraph of
the Berggren tree: every nonempty Berggren child word strictly increases the
hypotenuse.  An external edge must therefore be transported by the unique
undirected path that rises to a common ancestor and descends again.

[THM-3336](THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md)
gives the intrinsic folded weights of the three XOR directions:

```text
K_1105(5) ={5,221},
K_1105(13)={13,85},
K_1105(17)={17,65}.                                      (2)
```

These are symmetric edge colours.  They do not yet specify an ancestry path.

## 2. Exact ancestry paths and costs

Let `w_x` be the word of `x`, and let `r(x,y)` be the length of its longest
common prefix with `w_y`.  Write `bar(v)` for the reverse of a word with every
letter inverted in the free group `F(U,A,D)`.  The unique tree path is

```text
P(x,y)=bar(w_x[r:]) w_y[r:],
L(x,y)=|P(x,y)|=|w_x|+|w_y|-2r(x,y).                     (3)
```

Applying (3) to both edges of each matching gives the complete table:

| prime | folded weight | edge | `r` | path cost `L` | absolute depth jump |
|---:|---|---|---:|---:|---:|
| 5 | `{5,221}` | `000--100` | 0 | 8 | 0 |
| 5 | `{5,221}` | `010--001` | 0 | 27 | 17 |
| 13 | `{13,85}` | `000--010` | 2 | 5 | 1 |
| 13 | `{13,85}` | `100--001` | 1 | 24 | 18 |
| 17 | `{17,65}` | `000--001` | 0 | 26 | 18 |
| 17 | `{17,65}` | `100--010` | 0 | 9 | 1 |

The first spine collision `c=85`, with words `U^5` and `AD`, is a positive
control: its unique ancestry distance is `7`.

Three consequences are immediate.

1. A prime-XOR direction, even together with its folded weight, does not
   determine ancestry distance: the two costs in every matching differ.
2. It does not determine a fixed free-group path label, even up to reversal,
   because equal or inverse reduced labels have equal length.
3. No rooted Berggren-tree automorphism realizes an entire prime matching.
   Such an automorphism preserves depth, while every matching above contains
   an edge with unequal endpoint depths.

Still less can a nonempty positive Berggren word act on a fibre edge: it would
send one endpoint to a strict descendant with larger hypotenuse.  These are
obstructions to a constant path label, descendant-word action, and rooted-tree
automorphism.  They do **not** exclude an arbitrary finite-state transducer
reading the source address; one finite fibre cannot decide that global
language question.

## 2A. Scale control at the first eight-parent record fibre

**FINITE-EXACT EXTENSION UNDER INDEPENDENT AUDIT.**

The next Boolean record in THM-3334 occurs at U-depth `222`:

```text
c=99905=5*13*29*53,
X_c isomorphic to F_2^4/<1111>,                 |X_c|=8.  (3a)
```

Use the chart whose representatives begin with zero.  The complete ancestry
census is

| bit | odd leg | word | depth |
|---|---:|---|---:|
| `0000` | 96033 | `DA U^8 D^3` | 13 |
| `0001` | 67137 | `DADDDUDA` | 8 |
| `0010` | 32193 | `A U^11 A U^2` | 15 |
| `0011` | 70623 | `AUAAAA` | 6 |
| `0100` | 48063 | `U^4 A U^4 D U` | 11 |
| `0101` | 99807 | `U^6 D^22` | 28 |
| `0110` | 89823 | `U^22 A^2 D` | 25 |
| `0111` | 447 | `U^222` | 222 |

Each of the seven nonzero quotient directions is now a perfect matching with
four edges.  In lexicographic source-pair order, their folded weights, path
costs, and absolute depth jumps are

| direction | folded weight | path costs | depth jumps |
|---|---|---|---|
| `0001` | `{53,1885}` | `(17,17,31,203)` | `(5,9,17,197)` |
| `0010` | `{29,3445}` | `(28,14,28,238)` | `(2,2,14,194)` |
| `0011` | `{65,1537}` | `(19,23,225,41)` | `(7,7,211,3)` |
| `0100` | `{13,7685}` | `(24,36,40,228)` | `(2,20,10,216)` |
| `0101` | `{145,689}` | `(41,19,237,31)` | `(15,3,207,19)` |
| `0110` | `{265,377}` | `(38,230,26,34)` | `(12,214,4,22)` |
| `0111` | `{5,19981}` | `(235,33,43,17)` | `(209,17,13,5)` |

Every direction remains source-dependent and every full matching contains a
cross-depth pair, excluding one rooted-tree automorphism per direction.  All
`56` directed involutions and `336` ordered affine squares still satisfy
(5)--(6).  Within this fibre, the maximum observed path cost has grown from
`27` to `238`, and the maximum absolute depth jump from `18` to `216`.

This is a second-scale finite-exact witness, not by itself an asymptotic
theorem.  Section 2A alone does not prove unbounded path dispersion and does
not decide whether one uniform source-reading finite-state transducer works
across all Boolean fibres; Section 2B proves unboundedness by a separate
infinite construction.

## 2B. A prime-5 toggle proves simultaneous unbounded rank and dispersion

**PROOF EXTENSION UNDER INDEPENDENT AUDIT.**

The scale growth in Section 2A belongs to an explicit infinite branch
transplant.  For every integer `s>=1`, put

```text
t=25s+1,
C_s=2t(t+1)+1=1250s^2+150s+5.                            (3b)
```

The distinguished U-spine parent and its Gaussian spinor are

```text
S_s=(50s+3, 1250s^2+150s+4, C_s),
z_s=(t+1)+it=(25s+2)+i(25s+1),
word(S_s)=U^(25s).                                       (3c)
```

There is a literal prime-5 factorization

```text
w_s=(15s+1)+i(5s),
z_s=(2+i)w_s,
z'_s=(2-i)w_s=(35s+2)-i(5s+1).                           (3d)
```

Thus replacing the unique `2+i` factor by `2-i` gives the primitive parameter
pair `(35s+2,5s+1)`.  Its parent is

```text
T_s=(1200s^2+130s+3,
     350s^2+90s+4,
     1250s^2+150s+5).                                    (3e)
```

Primitivity follows from

```text
gcd(35s+2,5s+1)=1,            (35s+2)-(5s+1) is odd.
```

Moreover

```text
C_s=5(250s^2+30s+1),          C_s/5=1 mod 5,             (3f)
```

so `v_5(C_s)=1`: (3d) really toggles the sole prime above `5`, and the
intrinsic folded edge weight is `{5,C_s/5}`.

The transplant has an exact Berggren address:

```text
word(T_s)=DD U^(s-1) ADD.                                 (3g)
```

Indeed, using `U^n=I+n(U-I)+binom(n,2)(U-I)^2`, direct matrix multiplication
gives

```text
D^2 A U^(s-1) D^2 (3,4,5)^T=T_s,                         (3h)
```

where (3h) is the matrix order corresponding to the root-to-child address in
(3g).  Since the Berggren tree is unique, the two depths are exactly

```text
depth(S_s)=25s,                   depth(T_s)=s+4.          (3i)
```

Their words begin with different letters, so their longest common prefix is
empty.  Consequently the prime-5 edge has

```text
absolute depth jump =24s-4,
ancestry path cost  =26s+4.                               (3j)
```

Both quantities are unbounded.

Equation (3g) is also a positive branch-transducer result.  On the regular
unary language `{U^(25s):s>=1}`, the prime-5 toggle is the rational
transduction

```text
U^(25s)  |-->  DD U^(s-1) ADD.                            (3j')
```

Finite control counts input letters modulo `25`, suppresses the first complete
block, emits one `U` for every later block, and supplies the fixed boundary
words.  This is a source-reading transducer on one distinguished subbranch,
not a constant Berggren word and not a transducer for every Gaussian fibre.

This can be made simultaneous with arbitrary Boolean rank.  Given `h`, choose
distinct primes `p_1,...,p_h`, all `1 mod 4` and different from `5`.  The
discriminant-`-4` calculation in THM-3334 gives two roots of
`2t^2+2t+1=0 mod p_i`.  Chinese remaindering those choices with

```text
t=1 mod 25                                                (3k)
```

produces an infinite arithmetic progression of indices satisfying (3b) and
`p_1...p_h | C_s`.  Hence `omega(C_s)>=h+1`, so the fixed-hypotenuse Boolean
fibre has affine dimension at least `h`; taking a sufficiently large member
also makes (3j) exceed any prescribed bound.

Therefore, for every `h,B>=1`, some spine-selected fibre of Boolean dimension
at least `h` contains a prime-5 matching edge with both ancestry depth jump
and tree-path cost greater than `B`.  This proves simultaneous unbounded rank
and ancestry dispersion.  It still does not classify all matching costs or
decide a uniform finite-state transducer beyond the sublanguage (3j').

## 3. The exact lift is a flat source-dependent path groupoid

For a direction `p in {5,13,17}` define the oriented arrow

```text
a_p(x)=P(x,x+e_p).                                       (4)
```

Reversing the external edge inverts the path, so

```text
a_p(x) a_p(x+e_p)=1.                                     (5)
```

For distinct directions `p,q`, uniqueness of paths in a tree gives the flat
square identity

```text
a_p(x)a_q(x+e_p)
 =a_q(x)a_p(x+e_q)
 =P(x,x+e_p+e_q).                                        (6)
```

The exact companion checks all `12` directed involutions and all `24` ordered
affine squares.  It also checks all `60` directed closed walks whose two,
three, or four nonterminal vertices are distinct in the external `K4`; after
substitution of the source-dependent paths, every one freely reduces to the
identity.

This is structural, not a finite accident.  The path groupoid of a tree is
thin: there is exactly one reduced path between any two vertices.  Therefore
every closed ancestry path is trivial and

```text
H^1_graph(Berggren tree; A)=0                            (7)
```

for every constant abelian coefficient group `A`.  The external `K4` itself
has first Betti number three, but its exact ancestry-path realization sends
all three graph cycles to null paths in the ambient tree.

## 4. Why freezing one arrow per colour creates false holonomy

If one forgets the source and freezes only the three arrows based at `000`,
the quotient relation `e_p+e_q=e_r` is no longer respected.  The three reduced
defect lengths are

```text
(p,q;r)=(5,13;17): 39,
        (5,17;13): 39,
        (13,17;5): 37.                                   (8)
```

These nonzero words are hostile controls for a colour-only representation of
`V4` in `F(U,A,D)`.  They are not intrinsic curvature or cohomology classes:
the products in (8) identify arrows with different sources before composing
them.  Restoring the source-dependent second arrow gives (6) and zero
holonomy exactly.

Thus the lawful structure is

```text
external XOR torsor --source-dependent flat functor--> Berggren path groupoid,
```

not a translation-invariant action by three fixed ancestry words.

## 5. Tournament, cohomology, and transfer boundary

The ancestry path orients only after an ordered source and target are chosen;
reversing them gives the inverse path.  It therefore adds no intrinsic
tournament orientation to the symmetric XOR matching.

Equation (7) also separates two cohomologies that cardinality can blur.
THM-3336's `tau` and `lambda_p` are genuine `H^1` classes on the primitive
Gaussian associate multiplication group.  They are not Berggren-tree classes.
Any proposed word-current or D5 flux class must retain a nontrivial quotient,
current, boundary, or other sidecar; ancestry in the tree alone cannot carry
it.

Nothing here proves an LRC current or owner, a JC flux, a uniform/global
finite-state prime-toggle transducer beyond the unary sublanguage (3j'), or a
global bridge between the external `K4` and the internal four-box of THM-3339.
The theorem identifies the exact missing coordinate: source-dependent
ancestry path, not another scalar edge colour.

## 6. Exact controls and reproduction

The companion independently rebuilds all `15,904` Berggren vertices through
hypotenuse `99,905`.  It verifies both Boolean record fibres, all `34` external
edges, the `c=85` positive control, the constant-label and depth hostiles,
equations (5)--(6), the first-fibre closed-walk reductions and frozen-source
defects, and the second-fibre `56` involutions and `336` ordered affine
squares.

Reproduce with

```bash
python3 04-computation/prime_xor_ancestry_groupoid_thm3345.py
python3 -O 04-computation/prime_xor_ancestry_groupoid_thm3345.py
```

Both modes must match
`05-knowledge/results/prime_xor_ancestry_groupoid_thm3345.out` after LF
normalization.
