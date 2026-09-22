# Catalan gaps, signed Collatz cycles, and the carry that survives

**Status: PROVED elementary statements with explicit scope; CITED general
Catalan theorem; FINITE-EXACT independent cycle replay for least odd period
at most ten. OPEN: classification of all signed cycles and positive Collatz
convergence.** No novelty or prize claim is made. Catalan's theorem does not
turn a cycle equation into a convergence proof.

## 1. Inheritance and the connection contract

The ordered cycle gate, minimal supporting parameter, content, and signed
conjugacy are inherited from
[signed cycles, §§1–3](arithmetic_braids2_20260917_signed_cycles.md).
The unit-gap classification is already in
[the mod-six synthesis, §3](collatz_mod6_20260917_synthesis.md), with the
necklace construction in the experiment
`collatz_mod6_20260917_three_n_plus_k_catalan.py`.
[Wild typing, §2](collatz_mod6_20260917_wild_typing.md) already records an
independently audited height-census completeness bound through period 23 at
parameter one. Our period-ten enumeration is a different exact certificate,
not an extension of that record. The present note rederives the needed
mechanisms so that every retained and discarded coordinate is visible.

Closest proved mechanism: an ordered word determines an affine map and
its rational fixed point. Canonical hostile: the known negative period-seven
cycle has raw gap `-139`, yet is integral. Corrected near miss: an integer
cycle requires cancellation by its carry, not a unit raw gap. Least-used
sidecar: the distinction between a primitive cyclic word and its repetitions.

Anchor: the precise Catalan consequence. Niche: complete bounded-period
enumeration without a height cutoff. Wildcard: type the proposed gluing as
a two-sheet dynamical system, and test what happens when the sheet is lost.

| Live concept | Exact map or predicate | Information needed after projection |
|---|---|---|
| Signed state | `n=s m`, with `m>0` odd | Sheet `s` selects the additive parameter |
| Ordered word | `w -> (K,L,B)` | Totals `(K,L)` alone lose the carry |
| Integer realization | `q=abs(Delta)/gcd(B,abs(Delta))` | Supporting parameter must be divisible by `q` |
| Repetition | `(B,Delta) -> S_r(B,Delta)` | Least period; the reduced denominator survives |
| Unit power gap | `abs(2^K-3^L)=1` | Only a small sufficient subclass of `q=1` |
| Positive orbit | Actual prefix satisfies a descent inequality | Adaptive valuations and ordered carry |

The source of the proposed connection is a Collatz exponent word; the
target is a difference of perfect powers. The map forgets the order and
retains `(K,L)`. It preserves the affine multiplier `3^L/2^K` and its gap,
but destroys the translation `B/2^K`. Restoring `B` is necessary and
sufficient for the integer cycle gate below. The cheapest decisive test is
the pair of words `(1,3)` and `(2,2)` at the same clock `(K,L)=(4,2)`.

## 2. What signed gluing actually means

For odd nonzero `b`, define on nonzero odd integers, whenever `3n+b!=0`,

```text
U_b(n)=(3n+b)/2^v2(3n+b),          v2 uses absolute values.
```

Odd dilation gives `U_(db)(dn)=d U_b(n)` for every nonzero odd `d`.
In particular,

```text
U_b(-n)=-U_(-b)(n).                                      (C1)
```

This is a forward conjugacy after changing the parameter. It preserves
arrows and exponents; it does not reverse arrows. The inverse relation is
instead `n=(2^k u-b)/3` with `k>=1` and `2^k u=b mod3`, and generally has
infinitely many predecessors.

For `b=1`, neither half-line changes sign. Writing `n=s m` gives exactly

```text
(s,m) -> (s,U_s(m)),       s in {1,-1}, m>0 odd.           (C2)
```

Thus signed `3n+1` is the disjoint union of positive `3n+1` dynamics and
positive `3n-1` dynamics with a sheet marker. Forgetting the sign is not
a quotient dynamical map: at the smallest hostile magnitude `m=3`,
`abs(U_1(3))=5` while `abs(U_1(-3))=1`. At `m=1` both magnitudes are one.
Consequently the negative cycles cannot be reached from positive starts
of `U_1`, and identifying equal magnitudes erases a necessary coordinate.

The three **known** negative cycles of `U_1` correspond to the three known
positive cycles of `U_(-1)`. Neither phrase asserts that the list is complete
at arbitrary period.

## 3. The exact cycle gate, including the guards

For a nonempty positive word `w=(k_1,...,k_L)`, put `K_0=0`,
`K_i=k_1+...+k_i`, `K=K_L`, and

```text
B(w)=sum_(i=0)^(L-1) 3^(L-1-i) 2^K_i > 0,
Delta(w)=2^K-3^L != 0,
q(w)=abs(Delta)/gcd(B,abs(Delta)).                         (C3)
```

The composed affine relation is `2^K n_L=3^L n_0+bB`.
Its fixed point is `n_0=bB/Delta`. Both `B` and `Delta` are odd, and
`Delta` is coprime to three.

**Gate theorem, inherited and rederived.** The word is an exact signed odd
integer cycle word at parameter `b` if and only if `q|b`. Repetitions are
allowed; the least period can be shorter than `L`.

Necessity is the fixed-point equation. For sufficiency rotate the word one
place and call its carry `B'`. Direct expansion gives

```text
2^k_1 B'=3B+Delta.                                      (C4)
```

If `Delta|bB`, oddness of `Delta` implies `Delta|bB'`. All rotated fixed
points are therefore integers; they are odd because each rotated carry,
`b`, and `Delta` are odd. Equation (C4) says that successive points satisfy
`3n_i+b=2^k_i n_(i+1)`. The right-hand successor is odd, so `k_i` is the
actual valuation, not merely a formal division exponent. Rotation closes
the cycle. Every node has the sign of `b/Delta`.

Writing `B/Delta` in lowest terms also gives the inherited content identity
`gcd(n_i,b)=abs(b)/q`. Thus a primitive integer cycle has minimal parameter
`abs(b)=q`; when `b=+/-1`, integrality is exactly `q=1`, equivalently
`Delta|B`.

Here is the minimal same-clock carry obstruction:

| Word | `K,L` | `Delta` | `B` | Source at `b=1` |
|---|---|---|---|---|
| `(1,3)` | `4,2` | 7 | 5 | `5/7` |
| `(2,2)` | `4,2` | 7 | 7 | 1 |
| `(3,1)` | `4,2` | 7 | 11 | `11/7` |

The first and last are two markings of the rational cycle
`5/7 -> 11/7 -> 5/7`. Rational valuations are exact here because both
denominators are odd. This is a hostile to losing integrality, not a
counterexample to integer Collatz.

## 4. The complete unit-gap subclass needs only elementary arithmetic

Mihăilescu's general theorem says that the only consecutive positive perfect
powers with both exponents greater than one are `8` and `9`. Its original
source is [*Primary cyclotomic units and a proof of Catalan's conjecture*,
J. reine angew. Math. 572 (2004), 167–195](https://doi.org/10.1515/crll.2004.048).
The publisher's bibliographic record was checked; its full proof was not
audited here. The precise statement is also checked in
[Bilu's 2004 Bourbaki exposition](https://numdam.org/item/SB_2002-2003__45__1_0/).
For the fixed bases two and three, the following argument is independent
of that deep theorem and handles the exponent-one boundary.

**Unit-gap classification.** For `K>=L>=1`,

```text
abs(2^K-3^L)=1 iff (K,L) in {(1,1),(2,1),(3,2)}.          (C5)
```

For gap `+1`, `K=1` is impossible and `K=2` forces `L=1`.
If `K>=3`, the equation forces `3^L=7 mod8`, whereas a power of three is
one or three modulo eight. For gap `-1`, `K=1` gives `L=1` and `K=2`
is impossible. For `K>=3`, `3^L=1 mod8` forces `L=2m`; then
`(3^m-1)(3^m+1)=2^K`. The two positive powers of two differ by two,
so they must be two and four. Hence `m=1` and `(K,L)=(3,2)`.

Every word in (C5) has `q=1`, giving the complete unit-gap cycle table at
every odd nonzero parameter `b`:

| Clock `(K,L)` | Word(s) | Cycle | Raw gap |
|---|---|---|---|
| `(1,1)` | `(1)` | `(-b)` | -1 |
| `(2,1)` | `(2)` | `(b)` | 1 |
| `(3,2)` | `(1,2)`, `(2,1)` | `(-5b,-7b)` | -1 |

The two length-two words mark one cycle. All these cycles have content
`abs(b)`, so their primitive parameter versions occur at `b=+/-1`.
The empty word has `K=L=0` and gap zero, and is outside the gate theorem.

**Exact limit.** The third known negative cycle at `b=1` is

```text
-17 -> -25 -> -37 -> -55 -> -41 -> -61 -> -91 -> -17,
w=(1,1,1,2,1,1,4),  L=7, K=11,
Delta=-139, B=2363=17*139, q=1.                         (C6)
```

This is an integral primitive cycle outside the unit-gap class. At `b=-1`
its negation is a positive cycle. Catalan classifies (C5) but does not
prevent the cancellation in (C6).

## 5. Raw power gaps change even when the underlying cycle does not

For the concatenation `w^r`, direct geometric-sum expansion gives

```text
S_r=sum_(j=0)^(r-1) 2^(K(r-1-j)) 3^(Lj),
B(w^r)=S_r B(w),   Delta(w^r)=S_r Delta(w),
q(w^r)=q(w).                                           (C7)
```

The same factor cancels from numerator and denominator. Thus `q` is a
cycle-realization invariant under repetition, while the raw gap is not.
The smallest hostile is the fixed point one: word `(2)` has `B=Delta=1`,
but its repeated word `(2,2)` has `B=Delta=7`. Repeating the word `(1)`
of the fixed point minus one produces arbitrarily large negative raw gaps.
Using least-period words removes this artificial source of large gaps,
but does not remove the genuine period-seven cancellation (C6).

## 6. A transparent all-height certificate through period ten

**Clock bound.** Every signed integer cycle at `b=+/-1` satisfies

```text
L <= K <= 2L.                                          (C8)
```

Indeed `2 abs(n_i) <= abs(3n_i+b) <= 4 abs(n_i)` at every odd node.
Multiplying ratios around a cycle gives `2^L <= 2^K <= 4^L`.
Equality on the right forces every node to be `b`; equality on the left
forces every node to be `-b`. For a nonconstant least-period cycle both
bounds are strict. A cycle with the sign of `b` has `2^K>3^L`; the
opposite-sign case has `2^K<3^L`, also immediate from (C3).

It follows that every cycle of least odd period at most ten occurs among

```text
sum_(L=1)^10 sum_(K=L)^(2L) binom(K-1,L-1)
 = sum_(L=1)^10 binom(2L,L) = 250952                    (C9)
```

ordered words. The supplied script enumerates all of them, applies (C3),
reconstructs every accepted node using the actual map, checks all actual
valuations, reduces to the first return, and identifies rotations only.
No reversal identification and no bound on node magnitude is used.

**FINITE-EXACT outcome:** 37 accepted words for each of `b=1,-1`, reducing
at `b=1` to precisely `(1)`, `(-1)`, `(-5,-7)`, and (C6). The `b=-1`
list is their negation, also directly reconstructed. This certifies only
least periods at most ten. The 37 count includes repeated markings:
ten copies each of the two fixed points, ten markings of the two-cycle,
and seven markings of the seven-cycle.

## 7. What still has to be proved

For a positive start `n>1`, a finite prefix of its **actual** halving word
strictly descends exactly when

```text
(2^K_L-3^L)n > B_L.                                    (C10)
```

Obtaining such a prefix for every positive odd `n>1` would prove Collatz
by induction. Catalan gives no rule selecting that prefix and no bound
on its length. The negative-cycle enumeration does not supply one either.
Even classifying every finite cycle would leave the possibility of an
infinite orbit escaping all finite sets. This is the precise adaptive
orbit obligation beyond the algebraic cycle gate.

The concept board therefore closes with three usable mechanisms: sign must
remain a sheet coordinate; carry must remain beside the power gap; least
period must remain beside the repeated word. The next meaningful target
is a restriction on the **actual** joint evolution of carry and exponents,
strong enough to force (C10). Another classification of raw unit gaps
cannot supply that missing implication.

## Replay and audit scope

[Experiment](../../04-computation/experiments/catalan_elliptic_20260921_catalan.py)
and [deterministic JSON](../../04-computation/experiments/catalan_elliptic_20260921_catalan.json):

```text
python 04-computation/experiments/catalan_elliptic_20260921_catalan.py
python -O 04-computation/experiments/catalan_elliptic_20260921_catalan.py
```

Besides (C9), controls cover all 1,274 words with `1<=L<=6`, `L<=K<=2L`
under rotation, 20,384 parameter tests for odd `-15<=b<=15`, 6,096 signed
dilation tests, unit gaps for `1<=L<=K<=64`, repetition factors for five
representative words and `2<=r<=5`, and the same-clock rational hostile.
Every check uses explicit exceptions and survives optimization. The
all-exponent and all-height assertions depend on the proofs above; the
finite boxes are reproducible controls, not replacements for those proofs.

Independent peer audit: **PASS** on the two-sheet conjugacy, valuation
guards in the rotation argument, all exponent-one boundaries, repetition
cancellation, the clock bound, and the complete bounded-period universe.
A separate read-only `main()` replay reproduced all JSON values, including
the 37 accepted words and four cycles for each parameter. Normal and
optimized stdout are byte-identical; saved JSON agrees after platform
newline normalization. Source SHA256:
`186495baec1707b2456363c25ce54641d69f0730864ecf05b652092a61a4672f`.
JSON SHA256, with the master replay's canonical LF serialization:
`80c97213a7c9ed3259c9593e90d3db2ddcc62119ede830022a5935a176927aff`.
