# Distinct summands, diagonal forests, and the affine Collatz lift

**Status: PROVED elementary structure + FINITE-EXACT controls.** No Collatz
convergence theorem, cycle classification, Goldbach theorem, or LRC result is
claimed. This is the summand lane of the 2026-09-17 arithmetic-braids session.

## Inheritance and concept board

The closest proved mechanisms are
[THM-362, natural-operation-graph-shadows](../../01-canon/theorems/THM-362-natural-operation-graph-shadows.md),
[THM-2422, operation-fibres-summand-closure-and-twin-center-ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md),
and [THM-2433, operation-fibre-deletion-incidence-and-startup-scar](../../01-canon/theorems/THM-2433-operation-fibre-deletion-incidence-and-startup-scar.md).
THM-362's additive shadow permits equal summands. THM-2422/2433 retain the
swap-fixed diagonal and distinguish strict from weak parent fibres. That
distinction is essential here.

The historical source is
[summand-graph-fermat-zeckendorf.md](../../07-reflections/summand-graph-fermat-zeckendorf.md),
whose definition uses distinct positive summands, but whose Fibonacci-depth
and startup-restoration claims need the corrections in section 8 below.
[HYP-3003, summand-multiplicand-farey-basis-merge](../hypotheses/HYP-3003-summand-multiplicand-farey-basis-merge.md)
is a synthesis, not an additional proved dependency.

The live board is:

| Object | Predicate/invariant | Operation | Lost coordinate / hostile |
|---|---|---|---|
| Strict additive shadow | A missing ascending arc is a doubling | Delete equal parents | Allow equal parents: the missing set disappears |
| Dyadic forest | Odd core with integer height | Collapse a whole tower | Odd core alone loses the time spent halving |
| Collatz parent fibre | `a+b=z`, `a=2b-1` | Opposite affine shifts of parents | The source `1` is on the excluded diagonal |
| Prime-adic inverse braid | Full residue period of one inverse fibre | Scale the target by a power of two | `5n+1` has the same fullness and multiple cycles |
| Multiplicative shadow | Missing strict arcs are squarings | Exchange addition for multiplication | Square-root ancestry is not Pythagorean primitivity |

Anchor: make the user's summand/Collatz correspondence exact. Niche: recover
the operation-diagonal mechanism and historical correction. Wildcard: test
whether full inverse-fibre braiding distinguishes `3n+1` from `5n+1`.

## 1. The user's complement is exactly the doubling forest

Let `P={1,2,...}`, and let the ambient ascending digraph be

```text
A = {(x,z): x,z in P, x<z}.
```

Define the strict summand shadow by putting `x->z` whenever a positive
integer `y!=x` satisfies `x+y=z`. The companion is unique, namely `y=z-x`.
Consequently

```text
x->z belongs to G  iff  x<z and z!=2x,
A minus G = D = {x->2x:x>=1}.                              (1)
```

This proves the proposed complement statement in its precise universe.
The word "complement" means complement inside `A`, not reversal of arrows
and not complement inside every ordered pair of positive integers. If
equal parents are allowed, the graph is all of `A`, as in THM-362, so its
complement inside `A` is empty.

Every positive integer has a unique form

```text
n=2^k u,          u odd positive,          k>=0.            (2)
```

Thus `D` is the disjoint union of rays `u->2u->4u->...`, rooted at the odd
integers. Its edges are adjacent doublings; an arrow `u->8u` belongs to its
transitive closure, not to `D` itself. Reversing `D` gives precisely the
halving arrows used by Collatz.

The paired-parent information matters for a generation process: knowing
`x->z` does not certify that the companion `z-x` has been generated already.
On the fully labelled additive graph the companion can be reconstructed;
it is lost upon discarding labels or treating ordinary graph reachability
as two-parent generation. The faithful rule is the hyperedge `{x,y}->z`.

## 2. The multiplicative counterpart is a squaring forest

For any cancellative commutative operation `*`, the companion `y` in
`x*y=z`, if it exists, is unique. Hence deleting equal-parent witnesses
deletes precisely the diagonal arcs `x->x*x`, restricted to the chosen
ambient loopless shadow. This is the graph version of THM-2422/2433's
swap-fixed fibre correction.

For ordinary multiplication the weak loopless shadow is the proper-divisor
graph. Its strict version removes exactly

```text
x->x^2,              x>=2.                                 (3)
```

Every `n>=2` has a unique expression `n=b^(2^k)` with `b` not a square:
write `n=product p^e_p`, let `g=gcd{e_p:e_p>0}`, and set `k=v_2(g)`.
Thus the removed arcs form squaring rays rooted at nonsquares, with `1`
isolated in the loopless graph. Addition's doubling forest and
multiplication's squaring forest are exact instances of the same
equal-parent deletion mechanism.

Connection contract: operation fibres map to their diagonal orbit forests;
the predicate preserved is equality of the two parents; other legal
parent pairs are discarded; the sidecar is the full fibre. This is an
actual bridge between doubling and squares. It does not identify a
Collatz orbit with a primitive Pythagorean triple.

## 3. Exact reconstruction of shortcut Collatz

Define

```text
C(n) = n/2              if n is even,
       (3n+1)/2         if n is odd.
```

For an odd `n`, the selected additive witness is

```text
(a,b,z)=(n,(n+1)/2,(3n+1)/2),
a+b=z,                 a=2b-1.                            (4)
```

For `n>1` the parents are distinct and positive, so the odd arrow belongs
to `G`. The sole positive exception is `n=1`, for which `(a,b,z)=(1,1,2)`
lies on the removed diagonal. The shortcut cycle is `1->2->1`.

Write `3u+1=2^a v` with `u,v` odd and `a>=1`. In coordinates (2),

```text
(u,k) -> (u,k-1)        if k>=1,
(u,0) -> (v,a-1).                                         (5)
```

Collapsing the dyadic rays therefore gives the odd-only map

```text
T(u) = (3u+1)/2^v_2(3u+1).                                (6)
```

Keeping the exponent `a` reconstructs every omitted intermediate state
and the number of shortcut steps (`a`) between consecutive odd states.
For every positive initial integer, reaching the cycle `{1,2}` under `C`
is equivalent to its odd core reaching `1` under `T`. This is a lossless
reformulation when tower heights are retained, not a proof of convergence.

The formulation of `C` and its relation to the unhalved Collatz map are
also given in [Lagarias, *The 3x+1 Problem: An Overview*, section 1](https://arxiv.org/html/2111.02635v1).
All identities in this note are proved directly.

The graph-colour rule alone is insufficient. Define `B(1)=2`, let all
even nodes halve, and put `B(n)=n+2` for odd `n>1`. Every such odd arrow
is in `G`, yet `3->5->7->...` diverges. Changing just `B(5)` to `6`
creates the extra cycle `3->5->6->3`, still with all nontrivial odd
arrows in `G`. The needed information is the particular arithmetic
selection (4), followed by global control of its iterates.

## 4. The braid lifts to opposite affine shifts of the summands

The fixed-target-odd-core recurrence in the parent session is

```text
R(n)=4n+1,              F(n)=(3n+1)/2,
F(R(n))=4F(n).                                            (7)
```

On the paired witness (4), this becomes the exact map

```text
(a,b,z) -> (4a+1,4b-1,4z).                                (8)
```

It preserves `a+b=z` and `a=2b-1`; it scales the target by four and adds
two to its dyadic height. The two unit corrections cancel in the sum.
Thus the small affine rule generates an entire inverse fibre through a
precise conserved predicate, rather than a numerical resemblance.

The distinct-parent boundary is sharp: `(1,1,2)` maps to `(5,3,8)`.
Forward iteration preserves strictness once `n>1`, while inverse
iteration can land on the diagonal at `1`. This is exactly where the
trivial shortcut cycle attaches to the missing-edge forest.

For `t>=1`,

```text
R^t(n)-n = (4^t-1)(3n+1)/3,
v_3(R^t(n)-n)=v_3(t).                                     (9)
```

The second equality follows from `v_3(3n+1)=0` and the elementary odd-prime
valuation identity `v_p((1+p^r c)^t-1)=r+v_p(t)` when `p` does not divide
`c`. Hence `R` has one cycle of length `3^s` on odd residues modulo
`2*3^s`. It visits source classes `1->5->3->1 mod 6`.

This is an inverse-fibre orbit: every member has the same next odd target.
It is not a forward Collatz orbit. Source residue fullness says nothing
by itself about where that target's own outgoing edge lands.

## 5. Full prime-adic braiding also occurs with multiple cycles

The preceding mechanism generalizes. Let `p` be an odd prime, choose
`q=2^d` with `p | q-1`, and put

```text
c=(q-1)/p,              R_p(n)=qn+c,
F_p(n)=(pn+1)/2,        r=v_p(q-1).
```

For odd positive `n`,

```text
F_p(R_p(n))=qF_p(n),
R_p^t(n)-n=(q^t-1)(pn+1)/p,
v_p(R_p^t(n)-n)=r-1+v_p(t).                               (10)
```

Because `q` is even and `c` is odd, `R_p` preserves odd residues. Its
period modulo `2*p^s` is precisely

```text
p^max(0,s-r+1).                                           (11)
```

In particular, if `r=1`, it is a single cycle on every odd residue class
modulo `2*p^s`. If `r>1`, it fixes the initial `r-1` base-`p` digits and
fullness fails. The mechanism is the exact prime valuation of `q-1`.
Taking `q` to be the least power of two congruent to one modulo `p`
generates each nonempty inverse fibre of `F_p`; some target cores may
have no inverse if powers of two do not exhaust all units modulo `p`.

The summand lift also persists:

```text
(a,b,z)=(n,((p-2)n+1)/2,F_p(n))
    -> (qa+c,qb-c,qz).                                    (12)
```

For `p=5`, take `q=16`, `c=3`, `r=1`. The inverse fibre has the same full
prime-adic property, but the shortcut map has the explicit extra cycle

```text
13 -> 33 -> 83 -> 208 -> 104 -> 52 -> 26 -> 13.             (13)
```

Its odd-only cycle is `13->33->83->13`, with halving exponents `(1,1,5)`.
There is also the distinct cycle `1->3->8->4->2->1`. Thus full inverse
braiding, summand arrows, and dyadic return edges coexist with multiple
cycles. This is a hostile control against promoting local braid fullness
to a global convergence theorem.

For a cyclic odd word `n_0,...,n_(k-1)` satisfying
`pn_i+epsilon=2^a_i n_(i+1)`, define `A_j=sum_(i<j)a_i`. Direct composition
gives the exact global compatibility equation

```text
n_0 = epsilon * sum_(j=0..k-1) p^(k-1-j) 2^A_j
      / (2^A_k-p^k).                                     (14)
```

This identifies a missing coordinate: the ordered exponent word, not
merely local residue coverage. For (13), the numerator is `25+10+4=39`
and the denominator is `128-125=3`, yielding `n_0=13`.

## 6. What happens if even inputs also use `3n+1`?

There are two different extensions.

**Unhalved extension.** An even input `n=2^k u`, `k>=1`, supplies the new
arrow

```text
2^k u -> 3*2^k u+1,
```

whose target is odd and is `1 mod 6`. Every target `v>=7` with `v=1 mod 6`
has exactly one such even predecessor, `(v-1)/3`. These arrows leave
positive-height tower sites, unlike the usual shortcut's odd arrows.
If `n->3n+1` replaces the parity rule at every step, then

```text
n_t = 3^t(n_0+1/2)-1/2,
```

so every positive orbit grows without bound. If the arrow is added as
an optional move while retaining halving, the resulting graph is
nondeterministic. Existence of a path to one and convergence of every
allowed path are then different statements; the latter is already false.

**Half-step extension.** Applying `n->(3n+1)/2` at an even integer leaves
the integers immediately. With that affine rule at every step,

```text
n_t=(3/2)^t(n_0+1)-1.
```

When `n_0` is even, its denominator in lowest terms is exactly `2^t` for
every `t>=1`, since `3^t(n_0+1)-2^t` is odd. It never returns to an
integer. These two extensions must not be conflated.

## 7. Signed chains and three explicit minus cycles

Negation conjugates shortcut `3n+1` on negative integers to shortcut
`3n-1` on positive integers: `C_+(-n)=-C_-(n)`. Modulo six it fixes
class `3` and swaps classes `1` and `5`, exactly as proposed.

The odd-only minus map has the following three directly checked cycles:

```text
1 -> 1,
5 -> 7 -> 5,
17 -> 25 -> 37 -> 55 -> 41 -> 61 -> 91 -> 17.              (15)
```

Their existence is proved by substitution. The assertion that these are
the only cycles is not established here and must not be inferred from a
finite search. Under `3n+1` on signed integers they become negative
cycles. Both signs retain the dyadic-tower decomposition, so the number
of cycles is not determined by that decomposition alone.

## 8. Historical corrections and strongest survivors

The old summand reflection's claim that adjoining any one of `{1,4,6}`
restores all three is false: positive sums cannot create `1` from seeds
`{2,3,4}`. In fact, THM-2422 proves the original closure is
`P minus {1,4,6}`. Adjoining `1` fills all three; adjoining `4` fills `6`
but never `1`; adjoining `6` fills neither `1` nor `4`. The missing
coordinate was direction: dependency does not mean mutual reachability.

For the old table using seeds `{1,2}`, the exact synchronous closure at
depth `t` is `[1,2^t+1]`. Proof: from `[1,M]`, distinct sums produce
every integer from `3` through `2M-1`, and no larger integer; consequently
`M_(t+1)=2M_t-1`, `M_0=2`. Therefore

```text
d(1)=0,             d(n)=ceil(log_2(n-1)) for n>=2.        (16)
```

The apparent singleton `34` at depth six is a truncated table: the full
layer is `34,...,65`. In particular `35` and `55` also have depth six.
The Fibonacci-parent construction for `55=34+21` takes seven levels,
whereas `55=27+28` takes six. Thus the Fibonacci path is a valid sparse
construction but not the universal shortest-depth backbone. THM-2422
section 3 independently supplies the corrected synchronous law for
seeds `{2,3}`: after stage four its frontier is `27*2^(t-4)+1`.

The same reflection also incorrectly replaces all triangular numbers by
the finite seed set `{1,3,6}` in a three-summand claim: `11` is already
missing and every number above `18` is impossible. The classical
three-triangular-number statement uses the entire triangular sequence.
No analogy involving that finite seed substitution is inherited here.

## 9. Exact controls and remaining question

Reproduce with:

```text
python 04-computation/experiments/arithmetic_braids_20260917_summand.py
```

The script independently builds strict parent-pair shadows on `1..256`,
and finds `32640` ascending arcs, `32512` strict additive arcs, and exactly
`128` missing doublings. The multiplicative comparison removes exactly
`15` squaring arcs. It computes synchronous closure from `{1,2}` through
depth eight, with frontiers `3,5,9,17,33,65,129,257`; verifies all displayed
cycles; reconstructs shortcut Collatz from tower coordinates through
`1000`; and enumerates every odd residue for primes `3,5,7,17` through
modulus `2*p^4`. Periods at the fourth level are respectively
`81,625,2401,83521`. It also checks the parent lifts and half-integer
extension by exact rational arithmetic. All controls passed.

These are finite verification scopes supporting the algebraic proofs,
not trajectory evidence promoted to a theorem. The unresolved target is
to control the composition of successive target-core choices, retaining
the ordered exponent/height word. A complete inverse fibre, by itself,
has no information about that next composition.
