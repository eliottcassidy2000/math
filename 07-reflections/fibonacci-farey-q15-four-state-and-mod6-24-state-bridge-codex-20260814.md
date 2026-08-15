# Fibonacci, Farey fractions, q=15 modes, and the exact four-state quotient

**Date:** 2026-08-14  
**Status:** VERIFIED-EXACT finite quotient and automaton + elementary analytic
density/harmonic argument; unnumbered.  The construction identifies a
pointed four-cycle and a 24-state frame-line quotient.  It does not transport
a physical LRC owner/current, recover a Berggren word, or create a canonical
tournament.

## 1. Inheritance and the typed question

THM-3333 proves that primitive rational slopes are the faithful Farey carrier
and that Gaussian squaring transports determinant squared to the
Pythagorean light cone.  THM-3339 puts consecutive Fibonacci triples on
three exact Berggren rays and exposes six orders of three `V4` matching
channels, with a missing affine owner.  THM-3379 turns one ray bit into a
period-three subset of the natural numbers with an exact harmonic constant.
The q=15 triphase probe independently produces four unit sign-classes

```text
{1,2,4,7}=(Z/15Z)^x/{+-1}                               (1)
```

cycled by multiplication by two.

The typed question is whether these “four” and “six” objects are connected
by a lawful map.  The answer is yes at the level of one finite congruence
quotient, and no at the level of physical owners or canonical tournaments.

| field | connection |
|---|---|
| source | primitive Farey slopes, Fibonacci recurrence, and Berggren branch words |
| target | `P^1(F_3)`, the q=15 unit-sign cycle, and a 24-state mod-six frame-line packet |
| map | reduce a slope modulo three; reduce an oriented Farey frame modulo two; point the q=15 identity class |
| preserved | finite recurrence action, branch permutation, mod-three direction, mod-two frame order |
| destroyed | determinant magnitude/sign lift, exact fraction, word, height, affine owner, phase and current |
| required sidecar | an integral `Gamma(3)`/Farey lift and a pointed affine owner |
| cheapest hostile | `0/1` and `2/5` are distinct in `P^1(F_3)` but have determinant two, so finite adjacency does not lift to a Farey edge |

## 2. Fibonacci becomes a four-cycle after projectivizing mod three

Use the four projective states

```text
x0=[1:0], x1=[0:1], x2=[1:1], x3=[1:2] in P^1(F_3).    (2)
```

The Fibonacci matrix

```text
G=[0 1;1 1]                                              (3)
```

acts by

```text
x0 -> x1 -> x2 -> x3 -> x0.                             (4)
```

Indeed `G^4=-I mod 3`.  The ordinary Fibonacci vector modulo three has
Pisano period eight, while projectivizing identifies a vector with its
negative and halves the period to four:

```text
[F_k:F_(k+1)] = x_(k+1 mod 4).                          (5)
```

All four states already occur among reduced fractions in `[0,1]`:

```text
0/1 -> x1,       1/3 -> x0,       1/2 -> x3,       1/1 -> x2. (6)
```

If two reduced fractions are Farey neighbors, their determinant is `+-1`,
so their reductions are distinct in `P^1(F_3)`.  The converse fails:

```text
det((0,1),(2,5))=-2,
[0:1] != [2:5]=[1:1] in P^1(F_3).                       (7)
```

Thus the four-point quotient preserves one implication and forgets the
integral distance.

## 3. The pointed q=15 cycle is equivariantly the same cycle

Choose the sign representatives in `(1)`.  Multiplication by two gives

```text
1 -> 2 -> 4 -> 7 -> 1.                                  (8)
```

Point the identity class and define

```text
1<->x0,       2<->x1,       4<->x2,       7<->x3.       (9)
```

Equations `(4)` and `(8)` make `(9)` an equivariant bijection.  It is
pointed, not canonical without the selected identity/cusp and generator.

There is also a CRT description.  Every sign class has a unique
representative congruent to one modulo three.  Reduction of that
representative modulo five identifies `(1)` with `(Z/5Z)^x`; after the sign
normalization, multiplication by two becomes multiplication by three, a
generator of that cyclic group of order four.

Inversion acts on `(1)` by

```text
1 -> 1,       2 <-> 7,       4 -> 4.                    (10)
```

On `P^1(F_3)` it is represented by the involution

```text
S=[1 1;0 2],      x0->x0, x1<->x3, x2->x2,             (11)
SGS=G^(-1).
```

Hence `<G,S>` is the order-eight dihedral group (the repo's `D4`).

This reflection is physically visible in the q=15 cover.  A unit speed `u`
blocks the sign-pair of `u^(-1)`:

```text
u=1: {1,14},     u=2: {7,8},
u=4: {4,11},     u=7: {2,13}.                           (12)
```

Thus owner classes run forward around `(8)`, while their blocked unit pairs
run backward by the reflection `(10)`.  This is an exact finite owner/block
duality, not yet a physical ancestry transport.

## 4. The full ternary Berggren tree is a four-state S4 automaton

Reduce THM-3339's parameter matrices modulo three.  On the state order
`(x0,x1,x2,x3)`, their permutations are

```text
A=(0 1 3),        fixing 2;
B=(0 1 3 2);
C=(0 3 2),        fixing 1.                             (13)
```

They generate all of `S4`.  In particular the ternary tree does not act as
one `C3`; its three letters have orders three, four, and three on this
quotient.

If a node has state `i`, let each of its three children contribute to its
new state.  The transition matrix is

```text
M=[0 2 0 1]
  [0 1 0 2]
  [2 0 1 0]
  [1 0 2 0].                                             (14)
```

It is three-regular on rows and columns, and

```text
chi_M(X)=(X-3)(X+1)(X^2+3),
M^4=2M^3+6M+9I.                                         (15)
```

Starting at `x0`, the first level vectors are

```text
n=0: (1,0,0,0)
n=1: (0,2,0,1)
n=2: (1,2,2,4)
n=3: (8,4,10,5)
n=4: (25,20,20,16)
n=5: (56,70,52,65).                                     (16)
```

The non-Perron eigenvalues have moduli one and `sqrt(3)`, so uniformly from
any start state

```text
#state_i at depth n = 3^n/4 + O(3^(n/2)).                (17)
```

This is a genuine recurrence law for the entire ternary tree, not only the
three Fibonacci rays.  It retains only the mod-three endpoint direction;
many branch words collide.

## 5. Tree states as subsets of the natural numbers and harmonic series

Encode a ternary word `w` of length `n` by the heap address

```text
H(w)=(3^n-1)/2 + value_base3(w).                         (18)
```

Depth `n` is then one consecutive interval of `3^n` natural numbers.  Any
initial heap segment decomposes into `O(log N)` complete ternary subtrees.
Applying `(17)` to each subtree gives, for every state fibre `S_i`,

```text
#(S_i intersect [1,N]) = N/4 + O(sqrt(N)).               (19)
```

Partial summation therefore yields

```text
sum_(h<=N,h in S_i) 1/h = (1/4)log N + C_i + O(N^(-1/2)). (20)
```

For any Boolean subset `J` of the four states, the pulled-back tree subset
has density and harmonic coefficient

```text
delta(S_J)=|J|/4,
sum_(h<=N,h in S_J)1/h=(|J|/4)log N+C_J+O(N^(-1/2)).     (21)
```

This is the precise sense in which every selected collection of the four
finite states becomes a subset of the natural numbers and hence a subset of
the harmonic series.  The constant depends on the chosen branch-letter heap
order; the coefficient does not.

There are two different reciprocal weights in play.  Index weights `1/h`
in `(20)` diverge.  Value weights such as `1/F_h` converge exponentially.
Calling both “the Fibonacci harmonic subset” would erase the carrier.

For comparison, the six speeds in the canonical q=15 trimode edge have the
finite root mass

```text
1+1/2+1/3+1/4+1/5+1/7=1019/420.                         (22)
```

That is a finite speed-set mass, not the logarithmic tree-index mass in
`(20)`.

## 6. Why four states and six states naturally make 24

An oriented Farey edge is an integral determinant-one frame.  Reduction
modulo two gives one of the six elements of

```text
SL_2(F_2)=S3.                                            (23)
```

Equivalently, it gives one of THM-3339's six orders of the three nonzero
`V4` channels.  Reduction of its first column modulo three gives one of the
four directions in `P^1(F_3)`.  The combined finite carrier is therefore

```text
(oriented frame mod 2, first-column line mod 3),
6*4=24 states.                                           (24)
```

This is not numerology.  There are exactly

```text
|SL_2(Z/6Z)|=144                                         (25)
```

matrices, generated by the two integral elementary transvections.  The map
to `(24)` is onto and every state has fibre six.  It is a literal congruence
quotient of the oriented Farey-frame carrier.

Along the single Fibonacci path, the six-state order has period six and the
four-state line has period four.  Their product visits twelve compatible
states, not all 24.  The full 24-state packet is the ambient finite quotient
suggested at the end of THM-3339; identifying its four-point coordinate with
THM-3339's affine owner would still require a calibration and a current test.

The six objects in `(23)` are frame/order states.  They are not the six
vertices of THM-3339's separate edge-product `T6`.

## 7. The exact size-four tournament boundary

For chosen nonzero representatives `v_i` of the four projective lines,
orient `i,j` by the sign of

```text
det(v_i,v_j) in F_3^x={+1,-1}.                           (26)
```

Changing `v_i` to `-v_i` reverses every edge incident to vertex `i`.
Therefore `(26)` defines a switching class, not one intrinsic tournament.
Modulo the global simultaneous sign there are exactly eight sections.  They
give:

```text
4 tournaments with score sequence (0,2,2,2),
4 tournaments with score sequence (1,1,1,3),
every one containing exactly one directed 3-cycle.       (27)
```

The two score types are dual.  No tournament on four vertices is invariant
under the four-cycle `(4)`: its square swaps the endpoints of an opposite
pair, forcing a supposed invariant arrow to point both ways.

Thus three related finite objects must remain distinct:

```text
P^1(F_3) adjacency: K4, symmetric;
Fibonacci action:   directed C4, with two missing diagonals;
determinant section: one of eight switch-equivalent tournaments. (28)
```

This gives a precise realization of tournaments with missing or both-way
edges.  Forcing all three into one tournament destroys either the recurrence
or the projective gauge.

## 8. What is genuinely transplanted

The exact pointed cospan is

```text
q15 unit sign-class --(pointed cycle)--> P^1(F_3)
                                           ^
                                           |
                         primitive slope mod 3.          (29)
```

It preserves the four-cycle and inversion reflection.  The Berggren letters
act on the target by the full `S4` permutation packet `(13)`.  To lift back
from `(29)` one must restore:

```text
an integral primitive slope,
the determinant-one neighboring frame,
the exact branch word,
and, for LRC use, the owner mode, source phase and cochain. (30)
```

THM-3339's branch affine image is also an abstract order-eight dihedral
group, but its vertices are affine `V4` owners and its cocycle carries a
signed Pythagorean current.  The q=15 vertices are cyclic unit sign-classes.
The pointed isomorphism of the two abstract dihedral groups does not identify
their predicates.

## 9. Reproduction and scope

Reproduce with

```text
python 04-computation/fibonacci_farey_mod3_q15_four_state_probe_20260814.py
python -O 04-computation/fibonacci_farey_mod3_q15_four_state_probe_20260814.py
```

The semantic digest is
`f111020b7aebaa3b4108800fe981923e77befc586c8521818cd8aa8600c87baf`.
The companion checks the `P^1(F_3)` actions, q=15 cycle/reflection, 555
Farey-edge reductions and a converse hostile, full `S4`/dihedral groups,
the ternary transition recurrence through depth 20, all eight tournament
sections, all 64 possible invariant-tournament controls, and the complete
144-to-24 mod-six quotient.

No Farey converse, word reconstruction, affine-owner choice, physical LRC
current/phase, or JC flux follows.
