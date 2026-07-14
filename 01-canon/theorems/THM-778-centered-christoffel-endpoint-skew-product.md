---
id: THM-778
title: Centered-Christoffel endpoint words and the prime-sheet Euclidean skew product
status: PROVED (elementary midpoint-rank, Euclidean-cocycle, reconstruction, and skew-product theorem; finite-exact audit through 6,400 ordered pairs and five complete LRC movies)
source: codex-2026-07-14-S7
depends_on:
  - THM-773
related: [THM-536, THM-565, THM-637, THM-745, HYP-4078, HYP-6280]
verification: 04-computation/lrc14_centered_christoffel_endpoint_skew_product_codex_S7.py
  (+ 05-knowledge/results/lrc14_centered_christoffel_endpoint_skew_product_codex_S7.out and .json)
---

# THM-778 — Centered-Christoffel endpoint skew product

## 1. Endpoint clocks and the centered pair word

For a positive integer speed `u`, write its strict-danger endpoint clock as

```text
E_u={e_(u,i)=(2i+1)/(2u) : 0<=i<u} subset (0,1).
```

For two clocks `u,v`, the number of `v`-events strictly before the `i`-th
`u`-event is exactly

```text
C_(v|u)(i)
  = #{j : e_(v,j)<e_(u,i)}
  = ceil((v(2i+1)-u)/(2u)).                               (1)
```

Consequently the complete two-clock merge word is recovered without sorting:
before `u_i`, emit `v` until the emitted count reaches `C_(v|u)(i)`; then emit
the simultaneous block `uv` precisely when

```text
v(2i+1)=u(2C_(v|u)(i)+1),                                (2)
```

and otherwise emit `u`.  Emit any remaining `v` events at the end.  We call
this the **centered-Christoffel word** of `(u,v)`: it is the rational mechanical
cutting word sampled at cell centres rather than cell corners.

Let

```text
g=gcd(u,v),        p=u/g,        q=v/g.
```

Then:

1. the word for `(u,v)` is `g` repetitions of the word for `(p,q)`;
2. simultaneous blocks exist if and only if both `p` and `q` are odd;
3. in that case there is exactly one tie in each repetition, hence exactly
   `g` simultaneous blocks.

Thus a reduced continued fraction without the common-scale and parity data is
not a faithful endpoint address.

## 2. The Euclidean parity cocycle

For `s in {0,1}`, introduce the two centred lattice cosets

```text
F^s_(q,p)(i)=ceil(qi/p + q/(2p) - s/2).                   (3)
```

The increment word

```text
d^s_i=F^s_(q,p)(i+1)-F^s_(q,p)(i),       0<=i<p,         (4)
```

is cyclically balanced.  If `q=ap+r`, `0<=r<p`, then its letters are `a` and
`a+1`, with exactly `r` copies of `a+1`.  More precisely, put

```text
s'=(a-s) mod 2,       z=(a-s+s')/2.
```

Then the quotient-stripping identity is

```text
F^s_(q,p)(i)=ai+z+F^(s')_(r,p)(i).                        (5)
```

So an Euclidean shear removes the partial quotient `a` while updating one
parity bit.  Exchanging the two coordinate axes and iterating (5) is the
ordinary Euclidean algorithm for `q/p`; its quotient list is the simple
continued fraction of `q/p`.  Formula (5) is the additional midpoint cocycle
absent from the usual corner-based Christoffel word.

The pair word starts in the midpoint coset `s=1`; subsequent Euclidean
shears transport that bit by (5).  Thus the data

```text
(g, continued fraction of q/p, transported phase bit s)  (6)
```

therefore give a recursive decoder for the pair word.  The phase bit is not
decorative: it records which half-integer lattice coset survives each shear,
and the odd/odd case is exactly where the cutting line passes through cell
centres and produces ties.

## 3. A lossless global wall index

Let `W=(w_0,...,w_(r-1))` be any labelled family.  For owner `a` and local
event `i`, define the **centered Beatty rank**

```text
R_(a,i)=i + sum_(b != a) C_(w_b|w_a)(i).                  (7)
```

Then `R_(a,i)` is exactly the number of individual events of all owners that
occur strictly before `e_(w_a,i)`.  Hence:

- two indexed events have equal rank if and only if they lie on the same wall;
- sorting the distinct ranks and grouping equal ranks reconstructs the full
  simultaneous-wall schedule;
- the pairwise centered-Christoffel words are a redundant but lossless
  presentation of the global endpoint order.

This supplies an objective index for every wall event, analogous to the
tiling-fibre index in HYP-6825 but retaining the metric owner clock.

There is an exact involution.  If `N=sum_a w_a` and a wall of size `h` has
rank `R`, then reflection sends

```text
(a,i) -> (a,w_a-1-i),       x -> 1-x,
R -> N-h-R.                                                 (8)
```

In particular the owner-block word is palindromic.  Equation (8) is the rank
form of the global `t -> 1-t` symmetry, including the correction by the whole
simultaneous block rather than by one event.

## 4. Prime-seven fibre: an exact Euclidean skew product

Assume now `7` divides none of the `w_a`.  On the first chamber after `x=0`,
put every owner token at `k_a=0 in F_7`.  When a wall block `A` is crossed in
the positive direction, update

```text
k_a -> k_a-w_a^(-1)  (mod 7),       a in A,               (9)
```

leaving the other coordinates fixed.  At the wall itself, delete the tokens
of the owners in `A`.  Between walls the token state is constant.

The resulting skew product is exactly the strict prime-seven sheet movie of
THM-773:

```text
centred CF / Beatty wall schedule
          x
owner-labelled translations of F_7^r
          -> observation X^7-X | product_a(X-k_a).        (10)
```

After the complete word, owner `a` has crossed `w_a` events, so its net
translation is

```text
-w_a w_a^(-1)=-1  (mod 7).                               (11)
```

This recovers THM-773's global sheet carry rather than resetting it at
`x=1`.  At a wall, the polynomial observation in (10) is applied only to the
non-event tokens.

Under reflection, every non-event token transforms as

```text
k_a(1-x)=-1-k_a(x)  (mod 7).                              (12)
```

Since this affine map permutes `F_7`, covered walls occur in reflected pairs
with the same owner block.  Their owner word is therefore a palindrome.

## 5. Exact finite audit and the eight-owner wall word

The companion audit checks all `6,400` ordered pairs `1<=u,v<=80`:

```text
direct merge / formula mismatches             0
gcd-repetition mismatches                     0
odd/odd tie-count mismatches                  0
cyclic-balance mismatches                     0
Euclidean parity-cocycle mismatches           0
```

It then reconstructs five full LRC endpoint movies from (7), with zero rank
or fibre-transition failures.  The two seven-owner movies from THM-773 retain
their exact chamber counts `2` and `6`.

For HYP-6835's eight-owner survivor

```text
W=(108,169,143,213,206,197,30,162),
```

the audit finds `1,228` individual events grouped into `1,205` walls.  There
are `32` covered chambers and exactly `10` covered walls.  All ten walls are
simple, as THM-773 requires, and their owner-speed word is

```text
(162,108,108,206,197,197,206,108,108,162).                (13)
```

This word is palindromic by (8)/(12).  It includes the published wall
`x=19/216` owned by `108`; the other nine are new exact members of the same
absent-owner/heptagon transport skeleton.  Their exact least-path `n7-a267`
mask word is

```text
25773,32153,31115,14635,615,
30093,31115,615,14233,6035.                              (14)
```

Unlike (13), (14) is not palindromic.  This is not a failure of reflection;
it is a failure of the 25-mask quotient to retain enough labels.  An exhaustive
audit of all `5,040` owner-to-sheet permutations applies `k -> -1-k` before
recomputing the least-path mask.  Only `9` of the `25` source masks have one
possible reflected image; a source mask can have as many as `7` reflected
images.  Hence reflection on owner-labelled heptagon states does **not**
descend to a function on mask indices.

The degree-one redundancy sidecar is better behaved.  In either chamber
adjacent to a covered eight-owner wall, the event owner's token is the unique
duplicated sheet, so

```text
product_a(X-k_a)=(X^7-X)(X-d).
```

Across the ten walls its before/after root word is exactly

```text
((1->0),(4->6),(2->4),(0->2),(6->5))^2.                  (15)
```

Thus the mask quotient destroys the reflection lift, while the redundancy
root exposes a genuine period-five transport signal.  This identifies the
quotient polynomial, together with owner labels, as the sharper state for the
next Euclidean substitution automaton.

The open intervals between consecutive covered walls contain respectively

```text
57,301,3,24,329,24,3,301,57                           (16)
```

ordinary wall blocks, or

```text
58,306,3,25,336,25,3,306,58                           (17)
```

individual owner events.  Their full eight-coordinate owner-count vectors are
also palindromic.  Therefore the 1,205-wall movie reduces exactly to five
marked transport blocks and their reflection for the purpose of following its
covered-wall skeleton.

Thus continued fractions do not by themselves prove that the skeleton tears,
but they give its full recursive base word and reduce the next question to ten
owner-labelled mask walls in this example.  Any proposed substitution automaton
must retain the owner assignment above the mask and should act on the
redundancy root in (15).

## Proof

The inequality `e_(v,j)<e_(u,i)` is

```text
j < v(2i+1)/(2u)-1/2.
```

The number of nonnegative integers satisfying this strict inequality is the
ceiling in (1), including the integral boundary case.  This proves the merge
decoder.  Splitting `i=hp+i_0`, `0<=h<g`, places each reduced word in the
interval `(h/g,(h+1)/g)` and proves repetition.  A tie obeys

```text
q(2i+1)=p(2j+1).
```

Coprimality forces `2i+1=p ell` and `2j+1=q ell`; this is possible exactly
when `p,q` are odd, with `ell=1,3,...,2g-1`.  Hence there are `g` ties.

For (5), substitute `q=ap+r` into (3):

```text
F^s_(q,p)(i)
 = ai + ceil(ri/p+r/(2p)+(a-s)/2).
```

The definitions of `s'` and `z` give `(a-s)/2=z-s'/2`, proving (5).
Taking consecutive differences gives the two-letter alphabet and the count
`r`.  A length-`L` factor sum is
`F(t+L)-F(t)`; as `t` varies, two ceilings translated by the same real number
differ by at most one.  This proves cyclic balance.  Integer shears and axis
exchange are exactly the two operations of the Euclidean lattice algorithm,
so iterating (5) yields the continued-fraction decoder with its phase cocycle.

For (7), the first term counts the earlier events of owner `a`, and each
summand counts the earlier events of owner `b` by (1).  Thus the rank is the
strict predecessor count.  Equal-time events have the same predecessor set;
if `x<y`, the nonempty wall at `x` occurs before `y`, so the rank strictly
increases.  This proves global reconstruction.  Reflection sends every local
index as in (8); the events strictly before `1-x` correspond to the `N-R-h`
events strictly after the original wall, proving the rank identity.

Finally

```text
floor(w_a x+1/2)
```

increases by one precisely at the events of owner `a`.  THM-773's token
formula therefore gives (9), deletion at equality, and (11).  Replacing `x`
by `1-x` gives (12).  Polynomial divisibility is THM-773's exact coverage
criterion, completing the skew-product proof. ∎

## Tournament Analysis and information audit

- **Vertices:** owner clocks in the displayed tournament; endpoint events,
  Euclidean cells, simultaneous blocks, and token obligations were also
  considered and are the richer theorem carriers.
- **Pairwise observable:** which owner has the earlier next endpoint.
- **Switch/gauge:** chronological next-event order versus the Euclidean
  shear/axis-exchange presentation of the same pair word.
- **Tie Hamiltonian path:** owner index order inside a simultaneous next-wall
  block.
- **Fingerprint:** at every chamber the tournament is transitive, with score
  histogram `{0:1,...,r-1:1}`, no directed triangle, singleton SCCs, and one
  Hamiltonian path.
- **Destroyed by the isomorphism node:** the labelled Hamiltonian-path movie,
  all continued-fraction digits and phase bits, simultaneous-block sizes,
  centered ranks, inverse steps, token assignment, and metric wall positions.

On the eight-owner movie a single transitive node hides `948` distinct labelled
next-event paths and `6,620` pair-edge flips.  The isomorphism class is therefore
a constant base label, while the labelled Euclidean path is the dynamic stalk.
Even the 25-mask refinement is not closed under reflection without its
owner-labelled lift.
