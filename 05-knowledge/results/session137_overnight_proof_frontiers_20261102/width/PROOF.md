# Sharp context divisors, exact residue states, and a branching-resource obstruction

**Status: PROVED (all-order statements); FINITE-EXACT (explicit small state
closures and independent literal checks).** The full odd Hamiltonian-path
spectrum, and whether its spectral substitution width is bounded, remain
**OPEN**. No new universal spectral omission is asserted. No external
priority claim is made.

## 0. Inheritance, question, and live board

The read-only baseline is the Session136 packet at published commit
`bb7cb2d1f61a14632389814469317131d897b640`:
`05-knowledge/results/session136_full_tournament_spectrum_20261102/`.
The closest mechanism is `gluing/PROOF.md`, Sections1–5: unordered path-cover
substitution, bounded-context tomography, and a finite **capped-output**
algebra. Its Section8–9 establishes the spectral resource `w(h)`, including
`w(13)=w(19)=5` and `w(31)=6`. These are inherited, not rederived here.
The source names THM-3121,
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`;
the sparse local baseline does not contain that canon file, so the explicit
kernel statement and proof in the Session136 packet are the inherited input.

The canonical hostile is transitive padding: all `L_n` have `H=1`, but future
homogeneous contexts distinguish them. The corrected near miss is using the
old cyclic detector's nonzero diagonal as though it were invertible modulo
an arbitrary modulus. The least-used sidecar is the **unordered collection of
return excursions** between repeated visits to the same module.

The live concepts are: spectral width (anchor); arithmetic context precision
(niche); simultaneous branching versus unary wrapping (wildcard); exact
finite-state composition; and the difference between a factor's H and its
context responses. The result does not settle the anchor, but gives a new
exact arithmetic state space and an infinite constructor-resource separation.

Rédei parity and the odd-cycle-collection formula (OCF) are inherited. The
latter is explicitly used in Session136 `ear_geometry/positive_cut_geometry.md`,
Section1, with THM-4114,
`01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md`,
as its named source. We give all required elementary consequences below.
The Stirling lift is not claimed as a new discovery here.

## 1. The exact universal coefficient divisor

For a nonempty tournament T, let `p_T(c)` count covers of its entire vertex set
by c **unordered nonempty directed paths**, and put `p_T(c)=0` for `c>|T|`.
Thus `p_T(1)=H(T)`. A one-module context is a tournament Q with distinguished
vertex x; `Q[T]` replaces x by T homogeneously. All other Q vertices remain
labelled singletons. The inherited maximal-run splitting gives

```
H(Q[T]) = sum_{c>=1} a_c(Q,x) p_T(c),
a_c = c! W_c,
```

where `W_c` counts directed quotient words using each exterior vertex once
and x exactly c times, with no consecutive x letters. Let

```
d_c = c!(c-1)! .
```

**Theorem1 (sharp coefficient ideal).** For every context, `d_c | a_c`.
For every c>=1 there is an actual context on `2c-1` vertices with `a_c=d_c`;
it has no coefficients above c. Therefore the gcd of the c-th coefficients
over **all** tournament contexts is exactly `d_c`.

**Proof of divisibility.** Split a word at its x occurrences into a prefix,
`c-1` internal return segments, and a suffix. Each return segment begins at
an out-neighbor of x and ends at an in-neighbor of x. Permuting the internal
segments preserves every required directed adjacency. Their exterior vertex
sets are disjoint and nonempty, so the permutation action is free. Thus
`(c-1)! | W_c`, including the `c=1` convention `0!=1`.

A return uses at least two exterior vertices: a one-vertex return would
require both `x->v` and `v->x`. Consequently a context with m exterior vertices
has `a_c=0` if `2c-2>m`.

For sharpness write `s=c-1`. Take vertices
`x,a_1,...,a_s,b_1,...,b_s`, orient `x->a_i`, `b_i->x`, and orient
`a_i->b_j` iff `i<=j`. Orient each of A and B transitively (these internal
orientations do not affect the leading coefficient). Call this context E_s;
E_0 is the identity context. A word with `s+1` x occurrences must use all its
`2s` exterior vertices in s two-vertex return segments, with empty prefix and
suffix. It therefore corresponds to a perfect matching `a_i->b_{pi(i)}`.
The inequalities `i<=pi(i)` force pi to be the identity, by summing them.
The s matched returns can appear in exactly `s!` orders. Hence
`W_{s+1}=s!` and `a_{s+1}=(s+1)!s!=d_{s+1}`. Higher terms vanish by the exterior
budget. QED.

The new detector is jointly optimal in exterior size and positive leading
coefficient. Its first rows in `(p(1),p(2),...)` coordinates are

```
E_0: (1)
E_1: (3,2)
E_2: (9,20,12)
E_3: (31,166,300,144)
E_4: (117,1368,5448,8064,2880).
```

All coefficients in column c are multiples of `d_c`. By contrast the inherited
cyclic detector `C3[T,L_s,L_s]` has leading coefficient
`(s+1)!(s!)^2=d_{s+1}s!`. That larger coefficient is adequate over integers,
but can discard additional modular information.

## 2. Exact modular tomography and the minimal operational state

Fix an integer M>=1. Say `T ==_M U` if every homogeneous one-module context
has `H(Q[T])=H(Q[U]) mod M`, including the identity context. Set

```
kappa(M) = min {c>=1 : M divides d_c}.
```

The sequence d_c is a divisibility chain, since `d_{c+1}/d_c=c(c+1)`.
Therefore all coordinates at and above kappa vanish modulo M.

**Theorem2 (exact equivalence, not just sufficiency).**

```
T ==_M U
iff d_c p_T(c) = d_c p_U(c) mod M for all c>=1
iff p_T(c) = p_U(c) mod [M/gcd(M,d_c)] for 1<=c<kappa(M).
```

For contexts with at most m exterior vertices, replace the upper coordinate
bound by `min(kappa(M)-1, floor(m/2)+1)`.

**Proof.** Sufficiency is Theorem1 plus the kernel. For necessity apply E_0 to
recover the first weighted coordinate. Suppose the previous weighted
coordinates agree and apply E_{c-1}. Its lower coefficients are integer
multiples of the corresponding d_j, its top coefficient is exactly d_c, and
higher coefficients vanish. Subtracting the already equal lower terms gives
the c-th weighted congruence. This induction does not divide by a nonunit
modulo M. The exterior budget proves the bounded-context statement. Finally
`M | d_c z` iff `M/gcd(M,d_c) | z`. QED.

In particular, if c is the first **weighted** modular coordinate at which two
actual tournaments differ, their minimum detecting exterior size is exactly
`2c-2`. This is conditional on the existence of that pair; no assertion that
every hypothetical formal first difference is realized is needed.

The canonical signature may store either weighted residues `d_c p_T(c) mod M`
or the reduced residues in Theorem2. Equality of signatures is **exactly**
operational equivalence. This is a minimal sufficient invariant in the sense
of identifying exactly the contextual equivalence classes, not a claim that
all formal tuples are realizable or that a particular bit encoding is optimal.
There are at most

```
product_{1<=c<kappa(M)} M/gcd(M,d_c)
```

classes. At M=1 there is one class and the signature is empty.
For a prime power,

```
v_p(d_c) = 2 sum_{j>=1} floor((c-1)/p^j) + v_p(c).
```

For composite M, `kappa(M)` is the maximum of the prime-power cutoffs.
CRT combines the congruences; it does not assert that every independent
choice of prime-power class has a tournament representative.

**Hostile to the old decoder.** At M=8, L_3 and L_4 have signatures
`(1,3,1)` and `(1,3,0)`, in the precision convention `(mod8,mod4,mod2)`.
Every context with at most three exterior vertices is blind to the difference.
The old four-exterior cyclic detector has row `(11,34,24)` and responses
137 and393, equal modulo8. The sharp E_2 has responses81 and221, different
modulo8. Nonzero integer triangular diagonals are not invertible modular
diagonals. The missing coordinate is restored by the sharp coefficient ideal.

## 3. A representative-free finite algebra for every quotient bank

Operational congruence is preserved by arbitrary tournament substitution:
replace one child at a time, viewing the entire rest of the construction as
its homogeneous context. Here is an explicit arithmetic transition, so this
is not merely a nonconstructive finite-index observation.

For a quotient Q of order q and positive vector `c=(c_1,...,c_q)`, form the
loopless digraph D_Q(c) by replacing vertex i by **c_i independent labelled
vertices**, with Q's orientations between blocks and no edges within a block.
This auxiliary digraph is not asserted to be a tournament. Let

```
K_Q(c,d) = p_{D_Q(c)}(d).
```

The general path-cover substitution identity is

```
p_{Q[T_1,...,T_q]}(d)
  = sum_{c_i>=1} K_Q(c,d) product_i p_{T_i}(c_i).        (1)
```

Indeed select the unordered maximal runs in each child, and regard those
actual paths as the distinguishable independent copies of its quotient
vertex. A cover of the independent-copy graph specifies how the runs are
concatenated. Missing same-color edges enforce maximality, so no cover is
counted twice. This is the inherited path-cover kernel in a useful integral
coordinate system.

**Lemma (coefficient precision).** For every i,c,d,

```
d_{c_i} divides d_d K_Q(c,d).                          (2)
```

**Proof.** Fix the auxiliary graph D=D_Q(c). For every sharp detector E_s,
`H(E_s[D])` is divisible by `d_{c_i}`: collapse color i and permute its
return excursions exactly as in Theorem1. The other labelled exterior
vertices may now have missing arcs, which does not affect that free action.
The color-i block is independent, so all its `c_i` vertices occur as separate
runs. Now induct on d in

```
H(E_{d-1}[D]) = sum_{r<=d} a_r(E_{d-1}) p_D(r).
```

The top coefficient is d_d. Each lower term is an integer multiple of
`d_r p_D(r)`, divisible by `d_{c_i}` by induction. Subtracting proves (2).
This proof explains why the weighted *output* coordinate is essential.

To compute an output signature, choose any integer representatives of each
input reduced residue. In (1), multiply the d-th row by d_d and retain only
`1<=c_i<kappa(M)`. Omitted monomials vanish modulo M by (2). Changing one
chosen input representative by its allowed modulus
`M/gcd(M,d_{c_i})` also changes that monomial by a multiple of M, again by
(2). Thus the weighted output is well-defined independently of all lifts.
Convert it back to the reduced-residue convention by solving
`d_d r = weighted_output mod M`; division is only performed after removing
`gcd(M,d_d)`, when the remaining factor is a unit.

Consequently every fixed finite quotient/seed grammar has an effective
finite-state algebra for **all-depth modular H and all homogeneous modular
responses**. Close the finite set of seed signatures under the finite
operations, retaining expression witnesses. Termination is finite-state;
soundness and completeness are induction on generated terms. This says
nothing about integer equality, the capped spectrum at unbounded M, or
tournaments outside that grammar. The algorithm is effective, not asserted
polynomial in the bit length of M or q.

### A quotient-aware refinement: inactive repeated letters

The divisibility proof also gives a useful refinement before the factorial
profile weights are inserted. For any directed quotient-word content vector
c, let `R={i:c_i>=2}` and let `W_Q(c)` count its directed words. If i is
**inactive in the induced quotient Q[R]** (not on a directed cycle there), then

```
(c_i-1)! divides W_Q(c).                             (10)
```

A return from i to itself cannot stay entirely inside R, since that would be
a closed directed walk through i in Q[R]. Each of its `c_i-1` return segments
therefore contains a letter outside R, and every such letter occurs only once
in the whole word. These unique letters make the return-segment permutation
action free even after the other repeated letters have been collapsed.
This proves(10). In particular, if Q[R] is transitive then
`(max_i c_i-1)! | W_Q(c)`. More generally use the maximum count among inactive
repeated letters. There is no assertion that all these factorial factors
multiply: the separately obtained divisibilities only imply their lcm.

The bound is sharp for each c>=2. Take the actual quotient
`C3[x,y,L_{c-1}]`, give x and y count c, and give each of the c-1 exterior
vertices count1. The repeated support is the transitive pair x->y. The
quotient-word pattern is forced, and only the exterior vertices can be
permuted, giving exactly `W=(c-1)!`.

The activity sidecar is essential: in Q=C3 at content `(3,3,3)`, the three
cyclic starting positions give `W=3`, not divisible by `2!`. Thus replacing
“inactive in the repeated support” by “any repeated letter” is false. This
is a quotient-aware arithmetic pruning rule, not a new scalar spectral
obstruction and not a claim that quotient H determines its word kernel.

## 4. Exact small universal state counts

These are all-tournament class counts, not just counts in a small grammar.
Their upper bounds follow from identities below; their lower bounds have
actual width-three witnesses. Write a state in reduced-residue coordinates.

| M | signature | exact nonempty contextual classes |
|---|---|---:|
| 3 | `(H mod3,p2 mod3)` | 5 |
| 4 | `(H mod4,p2 mod2)` | 3 |
| 8 | `(H mod8,p2 mod4,p3 mod2)` | 14 |
| 12 | `(H mod12,p2 mod6)` | 11 |

### Elementary identities used for the upper bounds

Let `Q_T(z)=sum_c p_T(c)(z)_c`, with falling factorials `(z)_c`.
Then

```
Q_T(-z)=(-1)^n Q_T(z), n=|T|.                         (3)
```

For clarity, one elementary verification uses permutations. If a permutation
of the n vertices has b backward consecutive adjacencies, cutting it into
c directed paths requires all b backward cuts and chooses the rest in
`binom(n-1-b,c-1-b)` ways. Since ordered covers count each unordered cover c!
times, Vandermonde gives

```
Q_T(z)=sum_{permutations pi} binom(z+n-1-b(pi),n).
```

Reversing a permutation changes b to `n-1-b`. The identity
`binom(-z+n-1-b,n)=(-1)^n binom(z+b,n)` proves (3).
At z=1 modulo3, terms with c>=3 vanish, so

```
p2=0 mod3 if n is odd; p2=H mod3 if n is even.         (4)
```

Rédei parity on each part of a set partition gives

```
p_T(c)=S(n,c) mod2.                                  (5)
```

In particular p2 is odd for n>=2; p3 is odd precisely when n>=3 is odd.
For the M8 upper bound we also need

```
p2=H+2 mod4 when n>=4 is even.                        (6)
```

Here is a direct inherited-OCF derivation. Modulo4,
`H=1+2t`, where t is the number of all directed odd cycles (up to cyclic
rotation), modulo2. In the set-partition formula for p2, expand H on the
two blocks using OCF. The empty cycle collection contributes
`S(n,2)=2^(n-1)-1`. Each individual odd cycle of length ell contributes
`2 S(n-ell+1,2)`; two or more cycles have weight divisible by4.
When n>=4 is even, no odd cycle is Hamiltonian and `n-ell+1>=2`, so every
single-cycle Stirling factor is odd. Therefore
`p2=3+2t=H+2 mod4`. This is exactly the needed special case of the inherited
Stirling/OCF lift; the general lift is not claimed new here.

### Classification and actual witnesses

At M3, (4) gives `(h,0)` or `(h,h)` for h in Z/3, a union of five states.
At M4, (5) and H odd give precisely `(1,0),(1,1),(3,1)`; the zero p2 state
only occurs at order1.

At M8 there are the singleton state `(1,0,0)`, the order-two state `(1,1,0)`,
four states for even order>=4 from (6), and eight possible odd-order>=3
states with H odd, p2 odd, and p3=1. These give at most14 distinct states.
Every one occurs in the following table. A word is read chronologically from
K1, with `P(T)=1->T` and `C(T)=C3[T,1,1]`.

| state modulo `(8,4,2)` | word | H |
|---|---|---:|
| (1,0,0) | empty | 1 |
| (1,1,0) | P | 1 |
| (1,1,1) | PPC | 9 |
| (1,3,0) | PPP | 1 |
| (1,3,1) | PP | 1 |
| (3,1,0) | CP | 3 |
| (3,1,1) | CPP | 3 |
| (3,3,1) | C | 3 |
| (5,1,1) | PPCC | 93 |
| (5,3,0) | PC | 5 |
| (5,3,1) | PCP | 5 |
| (7,1,0) | CCP | 15 |
| (7,1,1) | CCPP | 15 |
| (7,3,1) | CC | 15 |

All these graphs have at most seven vertices; the replay checks their full
profiles by both subset partitions and literal permutations/cuts. The same
words provide all five M3 classes and all three M4 classes.

At M12, aside from singleton `(1,0)`, equations(4)–(5) give
`(h,3)` or `(h,h mod6)` for `h in {1,3,5,7,9,11}`. The two sets overlap
at h=3,9, giving ten states plus the singleton. Put
`B=C3[L2,L2,1]`, with profile `(11,39,33,10,1)` and hence state `(11,3)`.
The replay supplies actual J/C expressions starting with K1 and B for all
11 states, with at most nine vertices. Thus all four upper bounds are sharp.

This does **not** mean the width-three grammar realizes every integer count:
it already misses13 by the inherited Session136 theorem. Even exact
realization of every contextual state at several moduli is not integer
spectrum completeness.

## 5. An infinite branching-resource separation inside width three

Define the restricted grammar G by

```
G ::= K1 | J(G,G) | C(G),
J(A,B)=A->B,       C(A)=C3[A,1,1].
```

Joins may have arbitrary depth and arbitrarily complicated children. Every
cyclic node has at most one nonsingleton child; singleton means order1, not
H=1. Transitive padding is therefore counted as genuine branching resource.

**Theorem3.** The exact all-depth M12 contextual state set of G is

```
(1,0),(1,1),(1,3),(3,3),(5,3),(5,5),(9,3).              (7)
```

Consequently no G expression has `H=7 or11 mod12`.

**Proof.** The exact state transitions are

```
J: (h,p),(g,q) -> (hg, hg+2hq+2pg+2pq),
C: (h,p)       -> (3h+2p, 3h+8p),                      (8)
```

with the first coordinate reduced modulo12 and the second modulo6.
For J, choose whether/how the one or two components on each side are
matched across the join. For C, the full coefficient law is

```
p_C(d) = p(d-2)+(2d-1)p(d-1)+d(d+2)p(d)+d(d+1)p(d+1), (9)
```

where nonexistent coordinates and p(0) are zero. To see (9), fix c input
runs. The two exterior vertices admit respectively
`1,2c+1,c^2+2c,c(c-1)` choices of zero,one,two,three joining edges;
three joining edges require distinct input runs to avoid a directed cycle.
Each edge reduces the output path count by one. At d=1 this is `3H+2p2`;
at d=2 it is `3H+8p2+6p3`, whose final term vanishes modulo6.

Label the seven states in (7), in their displayed order, by A,B,C,D,E,F,G.
The entire transition table is below; the last seven entries of a row are
its J-products with columns A,B,C,D,E,F,G.

| state | C(state) | J with A | B | C | D | E | F | G |
|---|---|---|---|---|---|---|---|---|
| A | D | B | C | B | D | F | E | G |
| B | F | C | B | C | D | E | F | G |
| C | G | B | C | B | D | F | E | G |
| D | D | D | D | D | G | D | D | D |
| E | G | F | E | F | D | B | C | G |
| F | B | E | F | E | D | C | B | G |
| G | G | G | G | G | D | G | G | G |

Every entry follows by substitution into(8), so closure needs no census of
tournaments. Actual witnesses for A through G are respectively
`1`, `P(1)`, `P(P(1))`, `C(1)`, `P(C(P(1)))`, `C(P(1))`, and
`J(C(1),C(1))`, where `P(T)=J(1,T)`. Induction on expressions proves that
these are exactly the all-depth states. QED.

The larger width-three grammar escapes: the actual B above has H=11.
Joining it with `C3[L2,1,1]` of count5 gives55, hence residue7 modulo12.
For every k>=0, joins realize

```
11*25^k and 55*25^k
```

at width3, while Theorem3 excludes every one from G. Width2 permits only
transitive tournaments, so all these values have spectral width exactly3.
This is an **infinite strict hierarchy between constructor resources at the
same spectral width**, not a proof that spectral width itself is unbounded.
In any width-three J/C3 construction of a value congruent7 or11 modulo12,
some cyclic node must have at least two nonsingleton children. Arbitrarily
many serial wraps and arbitrary ordinal joins cannot replace that node.

Notice that unary P/C alone already realizes all14 universal M8 states,
yet even adding arbitrary joins still misses two H residues modulo12.
The detecting modulus is a resource-sensitive observer; a successful small
modulus is not evidence that all moduli will be equally blind.

## 6. Reproduction, negative controls, and boundaries

Run in this directory with standard Python only:

```
python3 -B context_residue.py > normal.json
python3 -B -O context_residue.py > optimized.json
cmp normal.json optimized.json
```

The output includes the **consequence objects**: exact state lists, actual
adjacency witnesses, modular grammar omissions, decoder rows, and the old
cyclic-detector failure. It checks:

* every one of1,099 labelled contexts through order5 for coefficient divisibility,
  with exact gcds `(1,2,12)`;
* all1,099 profiles against independent literal permutations and cuts;
* six sharp detector rows, through eleven quotient vertices;
* 320 independent direct graph substitutions against the representative-free
  modular kernel at M=3,4,8,12, using all ordered factor tuples from
  `{K1,L2,L3,C3}` and both quotient operations J,C3;
* the universal integral coefficient divisibility used to justify truncation;
* all content vectors with entries1..3 in all labelled quotients through order4
  for the inactive-repeat-support divisor, its sharp balanced families and
  the cyclic `(3,3,3)` hostile;
* all small universal-class witnesses, with literal cover checks for M3/4/8;
* all seven restricted grammar states and all eleven universal M12 witnesses;
* actual H=11 and55 branch escapes;
* the L3/L4 old-decoder hostile and the deliberate erasure of the highest
  modular coordinate.

Normal and optimized outputs are byte-identical. There is no dependency on a
repository Python module, network service, native solver, or assertion that
would disappear under optimization. The independent Session137 reviewer
separately enumerates literal x-words and sharp detectors.

**Loss ledger.** Source: actual homogeneous modules with unordered path-cover
profiles. Target: H responses modulo one fixed M. Map: maximal-run splitting
and return-excursion orbits. Preserved: exactly all homogeneous modular
responses, and closure under quotient substitution. Lost: integer magnitude,
labels/isomorphism, nonhomogeneous attachment responses, and any information
in coordinates annihilated by d_c. Needed sidecar: precisely the reduced
path-cover residues. Changing M changes the state space. Missing residue
classes are grammar-relative, not universal exclusions. General bounded
versus unbounded spectral width, full odd-spectrum realization, and an
all-modulus classification of width-three states remain **OPEN**.
