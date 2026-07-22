---
id: THM-2128
title: "The mod-169 seven-comb obstruction empties the seven-plus-one pencil"
status: >
  PROVED. In THM-2124's guard-13-blocker branch, suppose exactly seven of
  the eight nonblocker terminal residues have one projective direction
  modulo thirteen and the eighth is transverse. The seven parallel strips
  must cover every line of each guard-safe 13-root fiber. Two singleton
  phases are therefore impossible. A radius-2 character lemma either
  collapses the seven terminal bands to at most four distinct bands, or
  puts all seven terminals and the divided guard on one rational character
  line. The latter case reduces to a one-dimensional seven-comb containment.
  An exact mod-169 certificate excludes guard valuation one, and a two-level
  toothpick descent lowers every higher 13-adic valuation by two. Hence the
  (7,1) direction pattern is empty. This proof does not treat the all-aligned
  (8) pattern; THM-2131 subsequently excludes it by a next-digit lift.
source: codex-2026-07-22-LRC-seven-comb-descent
depends_on:
  - THM-2124
related:
  - THM-2097
  - THM-2120
  - THM-2122
  - THM-2123
  - THM-2125
  - THM-2131
script: 04-computation/lrc14_mod169_seven_comb_pencil_descent_codex_20260722.py
output: 05-knowledge/results/lrc14_mod169_seven_comb_pencil_descent_codex_20260722.out
script_sha256: 5fd6cccbe299e415e9736cba6a4c0c9f5c3b024caa2036c2a941aea8de9d1edc
output_sha256: c8761bac0a8dea18f217dc071a3a7e7032b7c8fbe427940eaa559d582099a89f
hash_basis: working-tree bytes with LF line endings
---

# THM-2128 -- the mod-169 seven-comb obstruction

Put

```text
epsilon=1/14.                                          (1)
```

The proof has three parts. First we isolate the exact finite obstruction at
`13^2`. Then we prove a two-level descent for one-dimensional combs. Finally
we show that THM-2124's seven-plus-one projective pencil necessarily enters
that one-dimensional theorem.

## 1. The finite mod-169 lemma

Work in `Z/169Z`. Define the guard-safe universe

```text
U={z:z mod 13 is one of 2,3,...,11}.                  (2)
```

Thus `|U|=10*13=130`. For every `r in (Z/169Z)^*` (equivalently `13` does
not divide `r`) put

```text
S_r={z in U:min(rz mod 169,-rz mod 169)<=12}.          (3)
```

Here residues are represented in `{0,...,168}`. The cutoff is exact because

```text
12 < 169/14 < 13.                                     (4)
```

> **Mod-169 lemma.** No seven sets of the form `S_r` cover `U`.

The proof below is a finite exact certificate small enough to display.

### Double toothpicks in ten columns

Partition `U` by its ten non-guard residues modulo thirteen. Multiplication
by a unit permutes the twelve nonzero residue columns. The twenty-four
nonzero residues in the strict terminal danger interval are

```text
+/-1,+/-2,...,+/-12 mod 169.                          (5)
```

There are exactly two in each nonzero residue column. Restricting their
inverse image to the ten columns of `U` gives

```text
|S_r|=20,
|S_r intersect {z:z=a mod 13}|=2       for a=2,...,11. (6)
```

Thus every `S_r` is a self-similar double toothpick in each of ten columns.
Also

```text
S_r=S_(-r),                    S_r=-S_r.               (7)
```

There are consequently `156/2=78` sign classes.

Suppose seven of them covered `U`. Every column would receive fourteen
incidences on thirteen points. It would therefore have exactly one point of
multiplicity two and twelve points of multiplicity one. In particular,

```text
sum_(i<j)|S_(r_i) intersect S_(r_j)|=10.              (8)
```

The seven sign classes are necessarily distinct: a repeated class already
contributes the intersection `|S_r|=20`, contradicting (8).

Negation in (7) has no fixed point on `U`, so every pair intersection in
(8) has even size. Let `P` be the graph whose edges are the positive pair
intersections. Equation (8) gives

```text
|E(P)|<=5.                                             (9)
```

Hence `P` is disconnected.

### The six-row zero-neighbor certificate

Let `Z(r)` be the set of sign classes `s` for which `S_r` and `S_s` are
disjoint. Multiplication of the variable `z` by a unit congruent to `+/-1`
modulo thirteen preserves `U`. It transports the family (3) and is transitive
on each of the six classes

```text
r mod 13 up to sign =1,2,3,4,5,6.                    (10)
```

For representatives `r=1,...,6`, direct reduction of (3) gives the following
zero-neighbor lists, with every label taken modulo sign:

```text
r=1:  7,8,9,10,11,12,14
r=2:  7,9,11,15,47,62,77
r=3:  7,8,10,11,16,51
r=4:  7,9,17,21,38,76
r=5:  7,8,9,11,18,23,59,64
r=6:  7,19,25,31,44,50,69,75.                        (11)
```

Intersecting these lists after the same transport, and taking bitwise unions
of the explicit sets (3), gives the complete certificate table

```text
type   |Z(r)|   max |Z(r) cap Z(s)|   max |Z(r) cap Z(s) cap Z(t)|
                 max |S_r union six zero-neighbors|

  1       7              4                       3              110
  2       7              4                       3              102
  3       6              4                       3               96
  4       6              4                       3               96
  5       8              4                       3              112
  6       8              4                       3              104.       (12)
```

Every maximum in (12) ranges over all other sign classes, not merely the
displayed representatives. There are only 78 masks of 130 bits; the companion
constructs them directly from (2)--(3), checks all pair and triple common
neighborhoods, and checks at most `binom(8,6)=28` six-neighbor unions per
representative. Its complete pair-intersection histogram is

```text
size:   0    2     4    6   8  10
count: 273 1690  728  182  52  78.                   (13)
```

Now take a smallest connected component `A` of the disconnected graph `P`.
It has size at most three, and every vertex outside `A` is a common
zero-neighbor of every vertex in `A`. If `|A|=2`, its five-vertex complement
contradicts the pair cap four in (12). If `|A|=3`, its four-vertex complement
contradicts the triple cap three. Thus `|A|=1`. The other six labels are zero
neighbors of its unique vertex, but the last column of (12) says their whole
union has size at most `112<130`. This contradicts the assumed cover and
proves the mod-169 lemma.

As an independent check, the companion exhausts all seven-subsets of the 78
sign classes by branch and bound. At a node with current union `M`, it adds
the seven largest individual marginal gains

```text
|(M union S_r) minus M|                                (14)
```

needed to complete the subset. Their sum ignores future overlap and is
therefore a rigorous upper bound. The exact maximum is

```text
max |union_(r in R)S_r|=116,
R={1,3,5,7,8,9,11},                                   (15)
```

leaving fourteen points. This second audit is stronger than the graph
certificate but is not used in the paper proof.

The frozen transcript is reproduced, with runtime checks active in both
modes, by

```bash
python3 04-computation/lrc14_mod169_seven_comb_pencil_descent_codex_20260722.py
python3 -O 04-computation/lrc14_mod169_seven_comb_pencil_descent_codex_20260722.py \
  | cmp - 05-knowledge/results/lrc14_mod169_seven_comb_pencil_descent_codex_20260722.out
```

## 2. A radius-two character lemma

We will use the following variant of THM-2120's finite-kernel rigidity.

> **Lemma.** Let `a,b,u` be nonzero characters of a connected two-torus. If
>
> ```text
> {Y:||a.Y||<epsilon and ||b.Y||<epsilon}
>       subset {Y:||u.Y||<=2epsilon},                 (16)
> ```
>
> then:
>
> 1. if `a,b` are rationally independent, there are integers `m,n` with
>    `u=ma+nb` and `|m|+|n|<=2`;
> 2. if `a,b` are rationally dependent, then `a,b,u` lie on one rational
>    character line.

In the independent case, the map `(a,b)` is a finite surjective torus map.
Condition (16) confines the subgroup `u(ker(a,b))` to the closed radius-`1/7`
arc. Every nontrivial finite circle subgroup has an element at distance at
least `1/3` from zero, so `u` kills this kernel. Character exactness gives

```text
u=ma+nb,                    m,n in Z.                  (17)
```

If `S=|m|+|n|>=3`, choose

```text
1/(7S)<delta<min(1/14,1/(2S))                         (18)
```

and prescribe `a.Y=delta sign(m)`, `b.Y=delta sign(n)`, using zero for a
zero coefficient. Then both inputs lie in their strict epsilon bands while
`2epsilon<S delta<1/2`, contradicting (16). Hence `S<=2`.

In the dependent case write `a=A alpha`, `b=B alpha` for the primitive
generator of their saturated rational line. If `u` were nonzero on the
connected circle `ker(alpha)`, its restriction would be surjective. A point
with `a=b=0` and `u=1/2` would contradict (16). Thus `u` is also an integer
multiple of `alpha`. This proves the lemma.

There is a useful seven-character consequence. Suppose (16) holds for every
pair among `c_1,...,c_7`. If some `c_1` is off the rational `u`-line, then:

* every character on the `u`-line equals `+/-u`, or, when `u/2` is integral,
  `+/-u/2`; this follows by comparing the off-line coefficient in (17);
* any second off-line character is, up to sign, `c_1+u` or `c_1-u`;
* both possibilities cannot occur, since their sum and difference are
  `2c_1` and `2u`, whereas the lemma would require an integer relation of
  coefficient `l1`-norm at most two equaling `u`.

Two proportional off-line labels are also impossible by the dependent clause
of the lemma, which would put `u` on their line. Hence each of the two possible
off-line supports occurs at most once.

Thus there are at most two off-line labels and at most two on-line danger
bands. The seven bands have at most four distinct supports, so their union
has Haar measure at most

```text
4/7.                                                   (19)
```

## 3. The one-dimensional seven-comb theorem

For a nonzero integer `H` and integers `r_1,...,r_7` not divisible by
thirteen, put

```text
C_H={t:||Ht||>1/7},
D_r={t:||rt||<1/14}.                                  (20)
```

> **Seven-comb theorem.** One never has
>
> ```text
> C_H subset union_(i=1)^7 D_(r_i)                    (21)
> ```
>
> even up to a null set.

We induct on the thirteen-adic valuation of `H`.

### Valuation zero

If `13` does not divide `H`, choose `z mod 13` so that `Hz` is, say, `6`
modulo thirteen and put `t=z/13`. Then

```text
||Ht||=6/13>1/7,
||r_i t||>=1/13>1/14                 for every i.     (22)
```

All inequalities in (22) have margin, so a whole open neighborhood of `t` is
uncovered. This refutes the almost-everywhere containment (21), not merely a
pointwise version.

### Valuation one

Write `H=13J` with `J` a unit modulo thirteen. At a point `t=z/169`, set
`w=Jz mod 169`. The guard is strictly safe exactly when

```text
w mod 13 is one of 2,...,11,                          (23)
```

and terminal `i` is dangerous exactly when

```text
w in S_(r_i J^(-1)).                                  (24)
```

All coefficients in (24) are units. The mod-169 lemma supplies `w in U`
outside all seven sets. Its strict margins again leave an open uncovered
neighborhood, refuting (21) almost everywhere.

### Two-level descent

Suppose `169|H`, write

```text
H=169J,                                                (25)
```

and assume (21) toward a contradiction. Let `E` be its null failure set and
put `B=[13](E)`. The map `[13]` on the circle is a finite smooth cover, so
`B` is null. Fix, away from endpoints and every translate `B-l/13`, a point

```text
y in C_(13J).                                         (26)
```

Then every thirteenth root above every translate `y+l/13`, `l in F_13`,
satisfies the assumed cover.

Its thirteen roots under multiplication by thirteen are

```text
x_j=x_0+j/13,                   j in F_13.             (27)
```

They all lie in `C_H`, since `13x_j=y` gives

```text
H x_j=13J y.                                          (28)
```

For each `r_i`, its danger set on (27) has one or two points. The exact grid
count is

```text
one point iff ||r_i y||<=1/14,
two points otherwise.                                 (29)
```

If `s(y)` labels are singleton, seven combs have total incidence `14-s(y)`.
Covering all thirteen roots forces

```text
s(y)<=1.                                               (30)
```

The set `C_(13J)` is invariant under `y |-> y+l/13`. Choose `y` outside all
thirteen translated endpoint/null sets, so (30) holds at every translate.
For a fixed label, the number of translates satisfying the singleton
condition in (29) is again one or two; away from endpoints it is one exactly
when

```text
||13r_i y||<1/14.                                     (31)
```

Summing (30) over the thirteen translates gives total singleton incidence at
most thirteen. If every label contributed two, the total would be fourteen.
Therefore some label satisfies (31), and

```text
C_(13J) subset union_i D_(13r_i)                      (32)
```

almost everywhere. Under the surjective substitution `z=13y`, (32) becomes

```text
C_J subset union_i D_(r_i).                           (33)
```

Thus a putative cover at guard `H` produces one at guard `H/169`. Repeating
reduces the valuation by two until it reaches zero or one, both already
excluded. This proves the seven-comb theorem.

The strict-boundary convention is harmless but important. Equations (29) and
(31) exchange a closed singleton test with a strict next-scale test. All
equality phases form a finite set on the circle. Removing them together with
the finitely many translated exceptional sets preserves full measure, and
the base witnesses (22)--(24) have strict margin.

## 4. Excluding the seven-plus-one projective pencil

Let `Gamma` be a rank-two character lattice, `K=Hom(Gamma,R/Z)`, and suppose

```text
g=13u,
c_1,...,c_7 mod 13 lie on one nonzero projective line L,
c_8 mod 13 is transverse to L.                        (34)
```

Assume the eight terminal danger bands cover the guard-safe region almost
everywhere. Let `E` be the null failure set and put `B=[13](E)`. Multiplication
by thirteen is a 169-sheet smooth cover of `K`, so `B` is null. For every
`Y` outside `B` with `||u.Y||>1/7`, all roots of `13X=Y` are covered; they
form one full `K[13]`-translate plane.

On this root plane, the first seven danger sets are unions of at most two
parallel affine lines with common normal `L`. The eighth set is transverse.
On any affine `L`-kernel line not selected by the first seven strips, the
eighth strip meets at most two of its thirteen points. It cannot cover that
line. Hence the first seven strips themselves select all thirteen parallel
lines.

For terminal `c_i`, the strip has one line exactly when

```text
||c_i.Y||<=epsilon.                                   (35)
```

If two of the first seven strips were singleton, their total line capacity
would be

```text
1+1+5*2=12<13.                                        (36)
```

Consequently, for every pair `i!=j`,

```text
{Y:||c_i.Y||<epsilon and ||c_j.Y||<epsilon}
 subset {Y:||u.Y||<=2epsilon}.                        (37)
```

Although the root-plane cover was obtained almost everywhere, (37) is an
exact containment: any strict counterexample would have an open neighborhood
of counterexamples and could not lie in the null exceptional set.

Apply Section 2. If some `c_i` is off the rational `u`-line, the seven
terminal bands have union measure at most `4/7` by (19). But the root-line
argument above also says that the seven bands alone cover the guard-safe set,
whose Haar measure is

```text
measure{X:||13u.X||>1/7}=5/7.                         (38)
```

This is impossible.

It remains that all seven characters lie on the rational `u`-line. Write

```text
u=k alpha,                  c_i=r_i alpha,             (39)
```

where `alpha` is primitive. Since `c_i mod 13` is nonzero, every `r_i` is a
unit modulo thirteen. Transversality in (34) makes `c_8` independent of
`alpha`.

The root-line argument already showed that the first seven bands alone cover
the guard-safe set almost everywhere. Disintegrate Haar measure over the
primitive character `alpha`. A scalar phase is uncovered by those seven
bands exactly when

```text
||13k t||>1/7,
||r_i t||>=1/14                   for i=1,...,7,       (40)
```

holds. Fubini therefore turns the seven-band cover directly into the
almost-everywhere one-dimensional containment

```text
C_(13k) subset union_(i=1)^7 D_(r_i).                 (41)
```

The seven-comb theorem contradicts (41). This eliminates the `(7,1)` pattern
and proves THM-2128. QED.

## 5. Scope, self-similarity, and Tournament Analysis

This proof does not itself exclude THM-2124's all-aligned `(8)` pattern. Eight
double toothpicks have three units of incidence slack on a thirteen-point
fiber, so the singleton inequality (36), the mod-169 seven-set lemma, and the
two-level descent do not transfer unchanged. THM-2131 subsequently closes that
pattern with a different self-similar object: a 169-point next-digit section.
Ranks nine through twelve have still more slack.

The challenged assumption was that the projective conclusion of THM-2124
had discarded too much affine data to be useful. In the `(7,1)` case the
transverse eighth strip cannot repair even one missing thirteen-point pencil
line. This recovers a critical seven-edge cover on every line, and its
toothpick count reproduces itself after two powers of thirteen. The
self-similar object is the operation

```text
guard H -> H/169,        singleton phase -> next-scale danger,             (42)
```

not one static chord graph.

Candidate tournament vertices were terminal labels, the 78 unit sign
classes, the ten residue columns, the thirteen root sheets, and proof
obligations. The intrinsic binary relation in the finite certificate is
zero versus positive intersection of the sets `S_r`. Its positive graph and
zero-neighbor sidecar preserve the cover obstruction; orienting positive
edges by their intersection size, with label order for ties, destroys the
disconnected-component argument's zero cuts and adds no information. The
faithful carrier is therefore the undirected intersection graph together
with exact 130-bit masks and the continuous null-set sidecar. No tournament
fingerprint is claimed as an equivalent certificate.
