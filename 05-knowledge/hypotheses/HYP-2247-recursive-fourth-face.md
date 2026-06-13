---
id: HYP-2247
status: OPEN carrier-method conjecture with finite scouts
source: codex-2026-06-05-S672
related:
  - HYP-2246
  - HYP-2245
  - HYP-2243
  - HYP-2240
  - HYP-2239
  - HYP-2238
  - HYP-2189
---

# HYP-2247: Recursion Is the Fourth Representation Face

## Claim

The repo's sum/product/fraction trichotomy should be split into a four-face
carrier:

```text
sum       = additive packets, moments, walls, counts
product   = local factors, shells, obstructions, automorphism orbits
fraction  = boundary/address state: owner, branch, lift, deleted card
recursion = transition law under iteration, descent, or outer extension
```

The important distinction is:

```text
fraction remembers a boundary state;
recursion proves how that state survives, shrinks, or explodes.
```

So the old loop

```text
sum -> product -> fraction -> sum
```

sharpens to

```text
sum -> product -> fraction -> recursion -> sum.
```

The final arrow means that a recurrence/outer-extension law can be unrolled
back into a new additive ledger.

## Source Anchors

Zhou-Markov's arXiv:0911.1933 uses recurrent integral sequences to prove
irrationality of `pi`, `tan(r)`, and `cos(r)` in the stated rational regimes.
For the `pi` proof, if

```text
I_n = integral_0^pi ((pi*x-x^2)^n/n!) sin(x) dx,
```

then

```text
I_0=2, I_1=4,
I_n = (4n-2) I_{n-1} - pi^2 I_{n-2}.
```

Under the false assumption `pi=a/b`, denominator clearing makes `b^n I_n` an
integer, the integral estimate makes it tend to `0`, and the recurrence keeps
the needed sequence alive.  The recurrence is not just a formula for a value; it
is a proof-producing machine.

Paris-Harrington gives the coloring extreme.  The finite statement is a Ramsey
coloring theorem with a relative-largeness boundary `|H| >= min(H)`.  Its
minimal witness function `sigma(n,k)` eventually dominates the PA-provably
total functions, and proof-length refinements tie the boundary to
`F_{epsilon_0}`.  Thus the recursive face can carry more proof strength than
the static coloring object visibly suggests.

An incoming S671 stub had already phrased this as an ultrafilter coloring
recursion: colorings are side choices over tuple atoms, bad colorings form an
outer-extension tree, and relative largeness prevents homogeneous traces from
being postponed forever into the tail.  That stub is preserved as
`HYP-2247-paris-harrington-ultrafilter-coloring-addendum.md`; the full
PH-specific S672 lane is `HYP-2247-paris-harrington-ultrafilter-coloring.md`.

## S672 Computation

The script `04-computation/recursive_fourth_face_coloring_s672.py` records:

1. Zhou-Markov recurrence polynomials `P_n(t)` for `I_n=P_n(pi^2)`:

```text
P_0(t)=2
P_1(t)=4
P_2(t)=24 - 2*t
P_3(t)=240 - 24*t
P_4(t)=3360 - 360*t + 2*t^2
...
```

2. A baby Paris-Harrington edge-coloring search: 2-color `K_N` while avoiding
   homogeneous `H` with `|H|>=3` and `|H|>=min(H)`.

| N | avoiding colorings | search nodes |
|---:|---:|---:|
| 3 | 6 | 13 |
| 4 | 18 | 63 |
| 5 | 12 | 265 |
| 6 | 0 | 987 |

This baby case collapses at the ordinary `R(3,3)=6` scale, but it displays the
right state-machine view: colorings are not merely counted; they are extended
while a boundary predicate is monitored.

3. A face tournament over:

```text
sum, product, fraction, recursion, raw_scalar
```

using the observable

```text
(static exactness,
 local factor leverage,
 boundary memory,
 outer-extension survival,
 growth witness,
 algorithmic use,
 coloring transfer).
```

The resulting tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1}
directed_3cycles=0
hamiltonian_paths=1
top_order:
  recursion > fraction > product > sum > raw_scalar.
```

This is the main methodological correction: after splitting recursion out of
fraction, recursion becomes a proof-strength apex rather than another symmetric
static representation.

## Repo Transfers

### S662-S665 Triune Cycle

HYP-2238 through HYP-2240 treated continued fractions, carry continuants,
deleted-card owners, branch sheets, and generics as the fraction face.  HYP-2247
separates two things that were blended:

```text
fraction  = retained state
recursion = transition law of retained state
```

This keeps the useful trinity while explaining why the "fraction" face kept
doing recursive work.

### LRC14

The current LRC14 no-leak target should be read as:

```text
sum       = odd-wall/pair-sum danger packets
product   = C=27 gcd shells and local residue factors
fraction  = carry-owner continuant v=r+27k
recursion = coherent +27 outer-extension theorem
```

The needed proof is not only "the owner coordinate separates sampled fibers."
It is:

```text
every allowed carry recursion either preserves the floor family or pays
positive strict tax.
```

### A000568

S671/HYP-2246 found:

```text
Phi(T)=multiset_v(iso_class(T-v), L/M/U side of outdeg(v))
```

injective through `n=8`.  In HYP-2247 language:

```text
product   = Aut(parent)-orbit incident cube
fraction  = deleted-card plus L/M/U owner address
recursion = n -> n+1 endpoint purity theorem
```

The next step is not just computing `n=9`; it is proving that the half-filter
address survives endpoint recursion.

### Paris-Harrington and Colorings

Ordinary Ramsey colorings are a sum/product object: count colors on tuples and
search for homogeneous blocks.  Paris-Harrington adds an observer boundary:

```text
|H| >= min(H).
```

That boundary turns a finite coloring statement into a recursion-strength
witness.  This matches the repo's owner/observer theme: adding a small
boundary condition can move a problem from static enumeration to proof-strength
growth.

### Outer Extension and Friedman Toys

S668/HYP-2243 was already the recursive face in finite form.  Coarse quotients
of colored rooted-tree extensions leaked, while embedding profiles were pure.
The extension operation is the recursion; the embedding/downset profile is the
fraction state that makes the recursion usable.

### Unit Distance and Cauldron Colorings

For unit distance:

```text
sum       = unit-edge spine / edge-gain counts
product   = Eisenstein direction and norm support
fraction  = point-deletion frontier owner
recursion = Moser slab or ear-extension ledger
```

For cauldrons/Schur:

```text
sum       = A+B=C additive witness
product   = residue/color-class sieve
fraction  = active cauldron/removal state
recursion = online placement/removal game
```

Both are coloring problems where the missing proof object is an extension law,
not a larger static count.

## Assumption Challenge

For LRC/Tournament Analysis, vertices need not be runners, arcs, or even
tournaments.  In S672 the meaningful vertices are representation faces, color
states, boundary owners, endpoint-extension states, proof obligations, and
outer-extension laws.

Preserved predicate:

```text
which retained state can be recursively extended without mixing the target
fiber?
```

Destroyed information:

```text
raw labels, exact geometric embeddings, and static scalar values when they do
not affect extension purity.
```

## Next Tests

1. Push HYP-2246 to `n=9`: does the half-filter trace stay pure under endpoint
   recursion?
2. Split every LRC14 carry table into `(fraction state, recursion law)`, then
   identify which coherent `+27` lifts are true transitions and which are
   accidental sampled rows.
3. Build a PH micro-lab for `3`-uniform, `3`-color relative-largeness states
   using a SAT/backtracking frontier rather than brute-force coloring.
4. For unit distance `n=21/22`, hold edge/direction product shadows fixed while
   varying point-deletion recursion owners.
5. For OCF/H(T), search for deletion/substitution continuants whose recurrence
   unrolls to Hamiltonian-path additive payloads.
