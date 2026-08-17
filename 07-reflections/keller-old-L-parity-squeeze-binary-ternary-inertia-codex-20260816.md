# The old Keller divisor is a binary fold carried through a ternary tower

**Status: proved synthesis around independently hostile-audited THM-3537 plus a sharply separated
VERIFIED-NUMERICAL wreath-state atlas through depth five.**  The exact
depth-two result is

```text
normalization exponent: 4,
inertia cycles:          4+2+1+1+1=9,
canonical x_2 order:     exponent 8, index length 2.       (1)
```

The first four numerical rows `1,4,14,46` suggested a two-exponential closed
form; the depth-five value `142` refutes its prediction `146`.  What survives
is the finer finite-state wreath description, not that closed form.

## 1. Inheritance pass

The closest proved mechanism is THM-3533's one-step `1+2` normalization:
one finite branch and one tame quadratic escaping branch.  The corrected
near miss is THM-2582's square class: it sees that old `L` is even at depth
two but cannot see whether its positive coefficient is `2,4,6,...`.  The
canonical hostile is the primitive rescaling `u` versus `Lu`, showing that
power-order multiplicity is not a field invariant.  The least-used sidecar
is THM-3529's finite-sheet unit, which says the distinguished finite branch
does not meet `L` again.

Those four inputs fit together more tightly than any does alone.

## 2. A three-value interval plus one parity bit

At depth one, `delta_1(L)=1`.  At depth two, discriminant transitivity gives

```text
delta_2(L)=3+c.                                           (2)
```

The finite child contributes zero by the finite-sheet unit.  The only live
relative algebra has degree three, so tame permutation ramification gives
`0<=c<=2`.  Meanwhile the intrinsic depth-two class is `[H]`, making the old
`L` coefficient even.  Therefore

```text
3+c is even and c in {0,1,2}  ==>  c=1.                  (3)
```

This is a reusable proof move: first trap an effective multiplicity in an
interval shorter than the parity modulus, then use the square class to select
the unique value.  It is the discriminant analogue of a bounded lift plus a
residue certificate.

## 3. The lawful binary/ternary connection

The level-one orbit partition is `1+2`.  On the next inverse step:

- the finite singleton has three unramified children, giving `1+1+1`;
- the escaping quadratic branch sees a relative `1+2` cubic type, whose
  ramification indices multiply by the inherited two, giving `2+4`.

Thus

```text
(1+1+1)+(2+4)=9.                                         (4)
```

This is the precise sense in which a binary fold lives inside the ternary
inverse tree.  It is not a claim that the entire ternary tree is binary, nor
that the affine point-count tree follows the same law.  The object here is
local inertia on normalized geometric sheets.

The literal coordinate sees an additional collision.  Both the quadratic
and quartic root packets reduce to `x_2=0`; their eight cross-pairs each have
difference valuation `1/4`.  Squaring adds four to the intrinsic different:

```text
field ramification 4 + coordinate collision 4 = order exponent 8.  (5)
```

So the normalization/order distinction is itself a two-layer `4+4`, but the
two fours have different meanings and must not be identified.

## 4. The two-exponential guess fails exactly where the child state changes

Chordal continuation around repeated shrinking transverse loops gives stable
permutation exponents

```text
d_1,d_2,d_3,d_4,d_5 = 1,4,14,46,142.                    (6)
```

The first four values fit

```text
d_n=2*3^(n-1)-2^(n-1)                 for n<=4.         (7)
```

But `(7)` predicts `d_5=146`, while the independently repeated depth-five
continuations give `142`.  Thus `(7)` is **REFUTED at depth five**.  Its
putative recurrence

```text
d_(n+1)=3d_n+2^(n-1)                                    (8)
```

has correction terms `1,2,4,4`, not `1,2,4,8`.  The generic ternary
replication and binary fold remain real ingredients, but their counts are
mediated by an `S_3` child-product state.  This is not Fibonacci—the visible
bases `3,2` were already different from the roots of `r^2-r-1`—and it does
not identify the Keller tree with a Berggren tree.

The full numerical cycle rows are

```text
n=1: 2,1
n=2: 4,2,1^3
n=3: 8,4^2,2,1^9
n=4: 12,8^3,6,4^3,1^27
n=5: 36,18,8^9,4^9,1^81.                               (9)
```

The appearances of `6` and `12` at depth four refute a naive story in which
every nontrivial cycle merely doubles forever.  At depth five they triple to
`18` and `36`, while the `4` and `8` cycles replicate ternarily.  The failed
formula discarded exactly this child-product distinction.

## 5. The next automaton

For a parent inertia cycle of length `ell`, multiply the three-sheet child
permutations encountered once around that cycle.  The product lies in `S_3`:

```text
child product identity       -> three cycles of length ell,
child product transposition  -> cycles of lengths ell and 2ell,
child product 3-cycle        -> one cycle of length 3ell.  (10)
```

This is the smallest faithful state space.  It keeps
the intrinsic binary relation—monodromy along a parent loop—without forcing a
tournament.  A tournament on four vertices would add orientations to six
pairs that the local cover does not supply.  The depth-two `4`-cycle is a
cyclic carrier with two unobserved diagonals, not a canonical `T4`.

The ancestry-local continuation now verifies the following transition
profiles, independently at each prescribed radius:

```text
1->2:  length 1 gets identity; length 2 gets transposition.
2->3:  fixed cycles get identity; lengths 2 and 4 get transpositions.
3->4:  fixed, one length-4, and length-8 get identity;
       length-2 and the other length-4 get 3-cycles.
4->5:  fixed, all length-4, and all length-8 get identity;
       lengths 6 and 12 get 3-cycles.                    (11)
```

Each profile reconstructs the complete next cycle multiset from `(10)` and
passes the block-projection check.  The cheapest decisive all-level test is
therefore no longer another scalar fit.  It is:

1. compute the wreath lift product in `S_3` for every parent-cycle state;
2. identify a finite closed state set under `(10)`;
3. prove its transition counts from the reciprocal inverse chart; and
4. only then solve the resulting state-vector recurrence for cycle counts.

## 6. Subsets and harmonic series

The failed formula has an exact failure set

```text
E={n>=1 : d_n != 2*3^(n-1)-2^(n-1)} subset N.            (12)
```

The faithful carrier is the indicator word or generating series
`sum_[n in E] z^n`.  The harmonic statistic `sum_[n in E] 1/n` is a useful
size probe but cannot recover `E`; distinct finite subsets already have the
same reciprocal sum.  This repeats the packet-return lesson: a subset of the
harmonic series is an indexed Boolean word first and a scalar sum second.

The verified window says `E` misses `{1,2,3,4}` and contains `5`.  Nothing is
yet proved about later membership.

## Connection contract

| field | exact answer |
|---|---|
| source | old-`L` tame inertia in the fixed Keller inverse tower |
| target | parent-cycle state plus accumulated child product in `S_3` |
| map | wreath restriction followed by product around a parent cycle |
| preserved | orbit length, tame different exponent, block ancestry |
| destroyed by exponent alone | cycle type and child-product state |
| exact depth-two output | `(4)(2)(1)^3`, exponent four |
| coordinate sidecar | canonical `x_2` index length two |
| numerical hostile | `(7)` fits depths one-four and is refuted at depth five |
| not transported | Berggren ancestry, signed `C4` class, LRC current, JC classification |

This lane advances the fixed Keller passport while keeping the attractive
tree/tournament analogies typed: the binary fold is real, the ternary carrier
is real, and their interaction is a wreath recursion rather than a cosmetic
tournament.
