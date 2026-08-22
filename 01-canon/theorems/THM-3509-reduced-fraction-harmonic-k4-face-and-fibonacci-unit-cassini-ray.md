---
id: THM-3509
title: "Reduced-fraction harmonic K4 face and Fibonacci unit-Cassini ray"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. Every reduced fraction in
  (0,1) is equivalent to one primitive positive four-term Fibonacci window
  and to one primitive positive harmonic triangle on a labelled K4. The
  harmonic face has a gcd decoder, its linear octahedral completion gives
  the Euclid Pythagorean current, and Stern--Brocot mediants add the four
  vertex values. The Berggren action splits the fractions into odd/even
  first-spinor trees rooted at 1/2 and 1/3; their normalized unordered
  primitive triples form exact two-element orbits under the leg-swap
  involution. Consecutive Fibonacci ratios are exactly the unit-Cassini
  slice. K4 and L(K4) are honest carriers of sizes four and six, but the
  Cassini sign lives on one oriented antipodal two-edge fibre and is erased
  by the three-perfect-matching quotient. No LRC or Jacobian transfer is
  asserted.
source: codex/fib-k4-bridge/2026-08-16
depends_on:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order
related:
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
  - THM-3379-fibonacci-ray-t4-mod3-projection-and-binary-median-fibre
  - THM-3382-fibonacci-dual-index-harmonic-bifurcation-and-heap-rank-rigidity
  - THM-3455-berggren-q-spine-cap-seven-fibonacci-spectrum-and-rank-obstruction
  - THM-3487-fibonacci-two-24-state-bundle-obstruction-and-three-h1-repairs
  - THM-3497-berggren-ancestry-language-rational-density-and-fixed-harmonic-split
  - THM-3510-binary-shortlex-equal-level-count-log-density-boundary
script: 04-computation/reduced_fraction_harmonic_k4_fibonacci_thm3509.py
output: 05-knowledge/results/reduced_fraction_harmonic_k4_fibonacci_thm3509.out
script_sha256: c12dc71df13e0c627740aaebf971c07342ce3302d674b276cc1600cf69266dd9
output_sha256: 3f31f6f80989b360227c792d525e013ecff9b06a8df09d252d31af801d8fdb9a
semantic_sha256: 7ca802f1f8706658b62fb330504544b9c46d473db6481429cb011b078695a46e
hash_basis: LF-normalized UTF-8 bytes; exact integer semantic ledger
---

# THM-3509 -- reduced fractions, a harmonic K4 face, and the golden ray

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.**

## 0. Inheritance and novelty gate

THM-3333 supplies the reduced-fraction/Euclid light-cone map; THM-2596
separates the binary Farey and ternary Berggren ancestries; THM-3339 supplies
the recurrence window, branch action, and Fibonacci Cassini slice; THM-3454
supplies the U-spine; and THM-2753 supplies the faithful six-edge action and
lossy three-matching quotient. None of those results states the primitive
harmonic-face converse (6)--(10), its gcd decoder, the exact two-tree
leg-swap quotient (24)--(27), or the location of Cassini sign on the oriented
antipodal fibre (46). Those are the new joint statements proved here.

THM-3510 correctly warns that six `K4` edges do not automatically become a
six-vertex tournament: fifteen comparisons are required. Here those
comparisons come from the displayed integer edge products and are strict
away from the classified root tie. MISTAKE-394 and MISTAKE-352 forbid
orienting that tie or promoting a modular/matching shadow into an action.

## 1. The three equivalent primitive objects

Let

```text
R={(m,n) in Z^2: 0<m<n and gcd(m,n)=1}.                (1)
```

For `(m,n) in R`, put

```text
x=n-m,                    y=m,
W(m,n)=(x,y,x+y,x+2y)=(n-m,m,n,n+m).                  (2)
```

Then `W=(w_0,w_1,w_2,w_3)` is positive, primitive, and satisfies

```text
w_2=w_0+w_1,             w_3=w_1+w_2.                 (3)
```

Conversely, every primitive positive integral four-tuple satisfying (3)
has the unique form (2), with `m=w_1` and `n=w_2`. Thus reduced fractions
in `(0,1)` are in bijection with primitive positive four-term windows of the
Fibonacci recurrence.

Label the vertices of a tetrahedron by `0,1,2,3`, put

```text
e_ij=w_i w_j,                                           (4)
u=e_01=xy,
v=e_02=x(x+y),
z=e_12=y(x+y).                                         (5)
```

The triangle on vertices `0,1,2` is harmonic:

```text
vz=u(v+z),             equivalently 1/u=1/v+1/z.       (6)
```

It is primitive and has the exact gcd decoder

```text
(gcd(u,v),gcd(u,z),gcd(v,z))=(x,y,x+y).                (7)
```

In particular, the original fraction and its complement are

```text
m/n=u/v,                 1-m/n=u/z.                    (8)
```

Conversely, every primitive positive integer solution of (6) occurs
uniquely in this way. Indeed, write

```text
v=da,       z=db,       gcd(a,b)=1.                    (9)
```

Equation (6) gives `dab=u(a+b)`. Since
`gcd(ab,a+b)=1`, there is a positive integer `k` with

```text
d=k(a+b),                (u,v,z)=k(ab,a(a+b),b(a+b)). (10)
```

Primitivity forces `k=1`; then `(x,y)=(a,b)`. Equations (2), (5),
and (7) therefore give three mutually inverse descriptions:

```text
reduced fractions in (0,1)
  <-> primitive positive Fibonacci windows on four labelled vertices
  <-> primitive positive harmonic edge faces of K4.     (11)
```

The word *primitive* is necessary. The harmonic triple `(2,4,4)` is twice
the root face `(1,2,2)` and would require the nonintegral vertex scale
`sqrt(2)`; it has no primitive integral decoder.

## 2. The octahedral completion and Pythagorean current

The remaining three edge products are forced linearly by the harmonic face:

```text
e_03=u+v,
e_13=2z-u,
e_23=v+2z=e_03+e_13.                                  (12)
```

The three opposite-edge products agree:

```text
e_01 e_23=e_02 e_13=e_03 e_12
             =xy(x+y)(x+2y).                          (13)
```

The Euclid triple attached to `m/n=y/(x+y)` is read linearly from the face:

```text
(a,b,c)=(u+v, 2z, v+2z-u)
       =(n^2-m^2,2mn,n^2+m^2).                        (14)
```

More generally, before imposing (6), the exact defect identity is

```text
c^2-a^2-b^2=4[vz-u(v+z)].                              (15)
```

Thus the Pythagorean light-cone equation is exactly the harmonic-face
equation in these coordinates; no tournament statistic is used to obtain
it.

Since `gcd(x,y)=1`, the raw triple (14) has content

```text
gcd(a,b,c)=1 if x is odd,
gcd(a,b,c)=2 if x is even.                              (16)
```

This is the required opposite-parity sidecar: `m=y` and `n=x+y` have
opposite parity exactly when `x` is odd.

## 3. Stern--Brocot addition and edge polarization

For arbitrary reduced fractions `m/n` and `p/q`, not necessarily adjacent,
the mediant has seed and window

```text
(n+q-(m+p),m+p)=(n-m,m)+(q-p,p),
W((m+p)/(n+q))=W(m/n)+W(p/q).                           (17)
```

When the two fractions are Farey neighbors, (17) is the Stern--Brocot child
operation. Hence the entire Stern--Brocot tree linearizes by componentwise
addition on the four vertex values.

The six edge products do not add. Their exact operation sidecar is the
polarization law

```text
e_ij(W+W')=e_ij(W)+e_ij(W')+w_i w'_j+w'_i w_j.         (18)
```

Dropping the two mixed terms is the first failed implication in trying to
transport Stern--Brocot ancestry directly on edge weights.

## 4. Two Berggren trees on all fractions

On primitive positive spinors define

```text
A(x,y)=(x,x+y),
B(x,y)=(x+2y,x+y),
C(x,y)=(x+2y,y).                                       (19)
```

For `r=y/(x+y)`, these act by

```text
A:r |-> 1/(2-r)       in (1/2,1),
B:r |-> 1/(2+r)       in (1/3,1/2),
C:r |-> r/(1+2r)     in (0,1/3).                      (20)
```

Every primitive seed outside the two seams has a unique parent:

```text
y>x:       A^(-1)(x,y)=(x,y-x),
y<x<2y:    B^(-1)(x,y)=(2y-x,x-y),
x>2y:      C^(-1)(x,y)=(x-2y,y).                      (21)
```

The coordinate sum strictly decreases in (21). Coprimality leaves only

```text
x=y:       (x,y)=(1,1),       r=1/2,
x=2y:      (x,y)=(2,1),       r=1/3.                  (22)
```

Every map in (19) preserves the parity of `x`. Therefore all reduced
fractions in `(0,1)` split into exactly two disjoint ternary Berggren trees:

```text
x odd:     root (1,1), fraction 1/2, raw primitive triples;
x even:    root (2,1), fraction 1/3, raw content-two triples.             (23)
```

The Stern--Brocot tree and the two Berggren trees have the same fraction
universe but different ancestry operations: binary mediants add neighboring
windows, whereas (19) is a ternary linear cross-section.

## 5. Exact parity pairing and the ancestry quotient

Define the parity-changing involution

```text
J(x,y)=(2y,x)       if x is odd,
J(x,y)=(y,x/2)     if x is even.                       (24)
```

On fractions it is the single Mobius map

```text
J(r)=(1-r)/(1+r).                                      (25)
```

It has no rational fixed point in `(0,1)`, satisfies `J^2=1`, swaps the two
trees in (23), and conjugates the branch alphabet by

```text
J A=C J,                 J B=B J,                 J C=A J.                (26)
```

If `(a,b,c)` is the raw Euclid triple of `(x,y)`, the normalized triple of
`J(x,y)` is `(b,a,c)`. For odd `x` the new raw triple is
`(2b,2a,2c)`; for even `x` it is `(b/2,a/2,c/2)`. Consequently

```text
{reduced fractions in (0,1)} / <J>
      <-> primitive Pythagorean triples with unordered legs.             (27)
```

Every fibre in (27) has exactly two elements. In particular,

```text
1/2 -> (3,4,5),
1/3 -> (8,6,10) -> (4,3,5),                            (28)
```

which is the minimal parity hostile.

The full fraction, window, or primitive harmonic face retains its unique
Berggren ancestry. The unordered-triple quotient retains ancestry only
modulo the root swap and the letterwise exchange `A<->C` in (26). A still
coarser branch-count summary loses word order: from the odd root,

```text
chronological AB -> (5,3) -> 3/8,
chronological BA -> (3,5) -> 5/8.                      (29)
```

Thus equal branch counts are not an ancestry decoder.

## 6. Honest carriers of sizes four, six, and three

The four entries of `W` are the labelled vertices of `K4`. Its six edge
products are the vertices of

```text
L(K4)=the octahedral graph.                             (30)
```

Adjacency in (30) means that two `K4` edges share an endpoint. This is an
intrinsic binary relation and gives an honest size-six carrier without any
orientation choice. Its three antipodal pairs are exactly the three perfect
matchings

```text
M_1={01,23},       M_2={02,13},       M_3={03,12}.      (31)
```

The quotient `E(K4)->{M_1,M_2,M_3}` has two-element fibres. As in THM-2753,
the induced `S4` action on the three matchings has kernel `V4`; the quotient
is a useful size-three carrier but loses the choice of edge within every
antipodal pair.

Numeric comparison supplies tournaments only when it is strict. For every
primitive seed,

```text
e_01<e_02<e_03<e_23,
e_01<e_12<e_13<e_23.                                  (32)
```

The only possible cross-chain equalities reduce to

```text
x=y,        x^2=2y^2,        x^2+xy-y^2=0.             (33)
```

The last two have no positive rational solution, and coprimality makes the
first equality exactly `(x,y)=(1,1)`. Hence, except at `1/2`, comparison of
the four vertex values gives a transitive `T4` and comparison of the six
edge products gives a transitive `T6`. At `1/2`, both carriers have ties;
there is no `T4` or `T6`, and no tie is oriented artificially.

Opposite-product equality (13) is not by itself a rational vertex decoder.
For example, assigning weight `2` to all six edges satisfies (13), but any
rational lift would require `w_0^2=e_01e_02/e_12=2`. The missing square
class is a necessary toric sidecar.

## 7. U-spine and the Fibonacci unit-Cassini ray

The Stern--Brocot parabolic boundary

```text
r_t=t/(t+1),                t>=1,                       (34)
```

has

```text
W_t=(1,t,t+1,2t+1),
(u,v,z)=(t,t+1,t(t+1)),
det((t,t+1),(1,1))=-1.                                 (35)
```

Thus its harmonic identity is the telescoping unit-fraction law

```text
1/t=1/(t+1)+1/[t(t+1)].                                (36)
```

This is the U-spine boundary, not the golden path.

Let `F_0=0,F_1=1`, and for `k>=2` put

```text
r_k=F_k/F_(k+1),
W_k=(F_(k-1),F_k,F_(k+1),F_(k+2)).                    (37)
```

Consecutive `r_k` are Farey neighbors and, for `k>=3`,

```text
r_(k+1)=mediant(r_(k-1),r_k),
W_(k+1)=W_(k-1)+W_k.                                  (38)
```

They form the alternating Stern--Brocot ray converging to `1/phi`.

On any seed define the labelled antipodal difference

```text
kappa=e_03-e_12=u+v-z=x^2+xy-y^2.                     (39)
```

Then

```text
|kappa|=1
  iff (x,y)=(F_(k-1),F_k) for one unique k>=2,          (40)

kappa(W_k)=(-1)^k.                                     (41)
```

For completeness, if `x>y`, then `kappa>=5`, so every nonroot unit solution
has `y>x`. The descent

```text
(x,y) |-> (y-x,x)                                     (42)
```

preserves coprimality, decreases the sum, and negates `kappa`. It ends at
`(1,1)`; reversing (42) gives precisely the consecutive Fibonacci pairs.
This proves (40), not merely a finite pattern.

The corresponding Pythagorean triples are

```text
(F_(k-1)F_(k+2), 2F_kF_(k+1), F_(2k+1)).              (43)
```

Their raw content is two exactly when

```text
k=1 (mod 3),                                           (44)
```

because `F_(k-1)` is even exactly when `3|(k-1)`.

For `k>=3`, `F_(k-1)<F_k`, and the complete edge order is

```text
e_01<e_02<min(e_03,e_12)<max(e_03,e_12)<e_13<e_23.    (45)
```

Moreover `M_3` is the unique unit-gap matching: the other two absolute gaps
are `c>1` and `2F_k^2-F_(k-1)^2>1`. Only the middle comparison flips, with
sign `(-1)^k`. The vertex `T4` order is unchanged. More importantly, both
`03` and `12` map to the single matching `M_3` in (31), so the bare
three-matching quotient cannot carry Cassini orientation. The sign lives on
the oriented two-sheet fibre

```text
{03,12} -> M_3,                                       (46)
```

in the recurrence labelling. An isolated swap `03<->12` fixing the other
four edges is not induced by a vertex permutation in `S4`. Thus (46) is a
genuine labelled orientation sidecar, not an unlabeled tetrahedral gauge.
At `k=2`, ties remain and no tournament is claimed.

## 8. THM-3506's conditional odd-face orbit in the even tree

This section is a related sidecar, not a dependency of Sections 1--7. Let
THM-3506's formal face state be `q=(e,m)^T`, with

```text
q'=[7 -2;3 -2]q.                                       (47)
```

Use the present theorem's seed chart

```text
(x,y)=(e-m,m),                 (e,m)=(x+y,y).          (48)
```

Conjugating (47) by (48) gives

```text
(x',y')=N(x,y),        N=[4 4;3 1],        det N=-8.   (49)
```

For any coprime odd `0<m<e`, the seed has coprime even `x` and odd `y`.
Its harmonic face and raw current are exactly

```text
(u,v,z)=(m(e-m),e(e-m),em),                            (50)

(a,b,c)=(e^2-m^2,2em,e^2+m^2).                        (51)
```

Equation (51) has content two, so division by two gives precisely
THM-3506's primitive triple

```text
((e^2-m^2)/2,em,(e^2+m^2)/2).                         (52)
```

The map `N` preserves the primitive even-`x` lane as an arithmetic map. Its
output has even first coordinate and odd second coordinate. Any common
divisor of the two outputs divides `8x` and `8y` by the adjugate identity;
it is also odd, so coprimality forces it to be one.

The currently lawful instances are typed as follows:

```text
(e,m)=(7,3)     -> (x,y)=(4,3)     -> (20,21,29),
(e,m)=(43,15)   -> (x,y)=(28,15)   -> (812,645,1037),
(e,m)=(271,99)  -> (x,y)=(172,99)  -> (31820,26829,41621).               (53)
```

The first two were THM-3506's verified full packets `H,J`.  THM-3513 and
THM-3522 subsequently promoted `G,R_5,...,R_8`, and THM-3528 now proves every
later raw cleared norm has the full packet.  Thus iteration of (49) is an
actual all-level fixed-map packet orbit.  THM-3529 additionally proves every
complete packet is a finite-sheet unit.  This supplies no image equation,
irreducibility, or Berggren ancestry word.

Nor is (49) one fixed Berggren branch word. In the seed chart the generators
are

```text
A=[1 0;1 1],       B=[1 2;1 1],       C=[1 2;0 1],
det(A),det(B),det(C)=(1,-1,1).                         (54)
```

Thus the proposed premise that every branch word has determinant `+1` is
false; the correct statement is that every branch word is unimodular with
determinant `+/-1`. If `N=lambda P` in `PGL_2(Q)` for a branch-word matrix
`P`, determinants would force

```text
lambda^2=-8/det(P) in {-8,8},                          (55)
```

and neither value is a rational square. Equivalently, `N` has projective
determinant square class `[-2]`, while branch words have class `[1]` or
`[-1]`. Orientation reversal alone is not the obstruction, since a word
with an odd number of `B` letters also reverses orientation; the square
class in (55) is.

There is already a fixed-instance ancestry hostile. The first verified seed
has root word `B`, while its image `(28,15)` has root word `A^6B`. The latter
does not have `B` as a prefix, so the verified update is not even a descendant
move from the first node. The matrix `N` is an arithmetic self-map of the
even tree's vertex set, not an edge or fixed-word action of that tree.

## 9. Connection contract and scope

```text
source:       a reduced fraction m/n, equivalently seed (n-m,m);
target:       primitive recurrence window W and harmonic K4 face (u,v,z);
map:          (x,y) |-> W |-> pairwise edge products;
preserved:    fraction, coprimality, recurrence, Farey determinant,
              Pythagorean nullity, and labelled Cassini difference;
destroyed by unordered triple quotient:
              leg order, parity-tree root, and A/C ancestry labels;
destroyed by matching quotient:
              the edge choice and sign inside each antipodal fibre;
needed sidecars:
              x parity, recurrence vertex labels, mixed polarization,
              and toric square class outside the primitive harmonic image;
cheapest hostiles:
              1/2 versus 1/3, AB versus BA, k=3 versus k=4,
              (2,4,4), all-two edges, and det(N)=-8 versus +/-1.          (56)
```

This is an arithmetic and finite-representation theorem. It constructs no
physical lonely-runner current, owner, phase, row exclusion, or LRC(14)
map. It constructs no polynomial map, Jacobian flux, or Jacobian-conjecture
transfer. Section 8 does not turn THM-3506's conditional renewal into an
unconditional face orbit. Similar words such as *current*, *Cassini*,
*tree*, and *matching* do not supply either missing cross-domain map.

## 10. Reproduction

Run

```bash
python 04-computation/reduced_fraction_harmonic_k4_fibonacci_thm3509.py
python -O 04-computation/reduced_fraction_harmonic_k4_fibonacci_thm3509.py
```

The standard-library companion pins LF-normalized hashes of the six proved
dependencies plus THM-3506's related hash and uses `require`, never `assert`.
Normal and optimized runs
agree. It checks all reduced fractions with `n<=240`, both Berggren trees
through depth eight, the Stern--Brocot tree through depth twelve, every Farey
edge in the denominator-64 box, the primitive harmonic converse through
`v,z<=600`, Fibonacci rows through `k=120`, all unit-Cassini solutions in the
`600` box, the U-spine through `t=240`, THM-3506's three fixed/exposed pairs,
the determinant-square-class no-go, and the ancestry, parity, orientation,
toric-lift, and operation hostiles above. The exact semantic digest is

```text
7ca802f1f8706658b62fb330504544b9c46d473db6481429cb011b078695a46e.
```

**QED.**
