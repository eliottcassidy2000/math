# The all-length frontier for even prime-odd-part aliquot cycles

Research frontier, 14 October 2026. The [length-three dyadic classification](aliquot_dyadic_three_chamber_classification_20261014.md) is proved. The all-length conjectures below are not. All prior packets remain unchanged, and no global novelty is asserted for the elementary restrictions derived here.

## 1. An exact deterministic carrier

For every positive even integer `n`, write uniquely

`n=t p`, with `t=2^v2(n)` and `p` odd and positive.

Define

`F(n)=(t-1)p+2t-1=n-p+2t-1=(t-1)(p+1)+t`.

This map preserves positive even integers. Its exact discrepancy from the true aliquot map is

`s(tp)-F(tp)=(2t-1)(sigma(p)-p-1)`.

It therefore agrees exactly when `p` is an odd prime. For composite `p>1` it undercounts the true divisor sum; for `p=1` it overcounts the divisor contribution that the prime formula would have supplied.

The fixed points of `F` are exactly

`n=t(2t-1)`, for `t=2^a`, `a>=1`.

Only those with prime `2t-1` are actual perfect numbers in the prime-cofactor class. For instance, `120=8*15` is fixed by `F`, but is not perfect.

There is an exact equivalence, allowing repetition of a primitive cycle:

- A positive integral cyclic closure of dyadic prime-cofactor affine templates is a closed orbit of `F`.
- Every `F` cycle gives such an integral template closure by reading its exact dyadic cores.

Indeed, integrality of the affine equations forces every cofactor to be odd, so each prescribed power of two is the actual dyadic core of its integer. The affine system has a unique positive rational closure for each finite core word because its composed slope is strictly less than one. There is no extra chamber-consistency assumption to guess.

## 2. Two separate conjectures

**C1, cofactor-one obstruction.** Every nonconstant finite cycle of `F` contains a power of two, equivalently a state with odd cofactor one.

If C1 holds, every finite aliquot cycle confined to even integers with prime odd part is constant. This implication uses the prime hypothesis; it is not a statement about general aliquot cycles.

**C2, complete nonconstant surrogate-cycle classification.** Up to cyclic rotation, the only nonconstant primitive cycles of `F` are

`8 -> 22 -> 14 -> 10 -> 8`,

`16 -> 46 -> 26 -> 16`.

C2 implies C1. C1 alone leaves room for additional cycles through other powers of two. Neither conjecture follows from the proved length-three theorem. Neither makes a claim that every trajectory of `F` is bounded or eventually periodic.

Both displayed cycles are checked exactly. Their odd-cofactor words are respectively `(1,11,7,5)` and `(1,23,13)`. They are not actual aliquot cycles: `s(8)=7`, and `s(16)=15`.

## 3. Global identities retained at every cycle length

For a length-`k` cyclic closure `n_i=t_i p_i`, summing the edge differences gives

`sum_i p_i=2 sum_i t_i-k`.                                (3.1)

Let `x=min_i t_i`, and write `t_i=x u_i`. Reducing the outgoing equations modulo `x` gives `p_i=-1 mod x`. Thus

`p_i=x H_i-1`, with positive integers `H_i`.

The exact identities become

`sum_i H_i=2 sum_i u_i`,                                  (3.2)

`H_i=u_i+u_(i+1) mod x`.                                 (3.3)

Equation (3.3) is equivalent to the second-order congruence

`p_i=x(u_i+u_(i+1))-1 mod x^2`.

The position and successor of every core must be retained. A multiset of exponents loses this condition.

One immediate boundary consequence is that any cycle containing cofactor one has minimum dyadic core exactly two: substituting `p_i=1` in `p_i=-1 mod x` gives `x|2`. Thus C1 would in particular force every nonconstant surrogate cycle to visit valuation one.

## 4. A proved local tradeoff between size and dyadic core

Write `a=v2(t)` and `b=v2(p+1)`. From `F(n)=(t-1)(p+1)+t`,

- if `b<a`, then `v2(F(n))=b<a`;
- if `b>a`, then `v2(F(n))=a`;
- if `b=a`, write `p+1=t h` with `h` odd; then
  `F(n)=t[(t-1)h+1]`, and its dyadic core strictly increases.

In the last case,

`F(n)-n=t(2-h)`.

If `h=1`, the transition is `t(t-1) -> t^2`, so its successor has cofactor one. If `h>=3`, the integer strictly decreases even though its dyadic core increases.

When the dyadic core stays constant, `p+1` is a positive multiple of `2t`. The map is fixed if that multiple is one, and otherwise strictly decreases the integer.

Consequently:

**Local tradeoff lemma.** Every increasing edge whose successor does not have cofactor one strictly decreases the dyadic core. Every nonfixed edge whose core stays constant or increases strictly decreases the integer, unless it is the exceptional transition `t(t-1) -> t^2`.

In a nonconstant cycle avoiding cofactor one, increases in the integer and increases in its dyadic core therefore cannot occur on the same edge. This is a proved restriction, not a global Lyapunov function.

For example, the artificial loop `12 -> 14 -> 12` obeys this local anti-order rule using the actual cores `4,2`, while neither edge is the appropriate `F` edge. The rule alone cannot exclude cycles; additional arithmetic information is essential. A candidate proof must control the compensating blocks of decreasing cores and increasing integers, rather than treating the local statement as a global decrease.

## 5. A prime-specific odd-length restriction

There is one further unbounded restriction on the actual prime-cofactor problem:

**Proposition.** An odd-length aliquot cycle in even prime-odd-part chambers whose dyadic exponents are all odd can only be the fixed point `6`.

For odd `a`, `t=2^a=-1 mod 3`, and direct substitution gives `F(n)=-n mod 3`. Going around an odd cycle forces every member to be divisible by three. Because its odd part is prime, that odd part must equal three.

For odd `a>=3`, however,

`F(3*2^a)=5*2^a-4=4(5*2^(a-2)-1)`

has 2-adic valuation exactly two, contradicting the stipulated odd-exponent class. The remaining case is `a=1`, namely the fixed point `6`.

Primality is essential here. The same congruence only forces an arbitrary affine cofactor to be divisible by three, not to equal three.

## 6. Finite evidence and its exact limits

Two small falsification passes were run, with no counterexample to C1 or C2:

1. The functional graph restricted to even integers `2<=n<=10^6`, retaining an edge only when its true `F` successor also lies in that set, has exactly nine fixed points and the two displayed nonconstant cycles. The fixed points are
   `6,28,120,496,2016,8128,32640,130816,523776`.
   This classifies cycles **wholly contained in this window**. It does not establish convergence of every start below `10^6`, and it does not rule out a cycle that leaves the window.

2. For template lengths `4,5,6,7` with all exponents in `1,...,6`, the numbers of nonconstant cyclic representatives checked were respectively `330,1554,7820,39990`. The only integral closures found were the four-cycle template `(1,1,1,3)` with cofactors `(11,7,5,1)`, and, at word length six, a double traversal of the already-known three-cycle template `(1,1,4)`. There were no integral nonconstant closures at word lengths five or seven in this window. Template word length need not equal primitive period.

The dependency-free companion script reproduces these checks, all displayed cycles, the global identities on them, and the local tradeoff for every even input in the finite graph's vertex range. Its checks remain active under `-O`.

## 7. The remaining mathematical task

The carrier and arithmetic constraints are now exact: a deterministic map, a sum identity, a successor-indexed congruence, and a local exchange between size and valuation. The next useful step is a genuinely global invariant or a return-block argument. It must either force a visit to cofactor one in every nonconstant periodic orbit, or produce an explicit counterexample to C1.

A separate argument would then be needed to classify cycles through powers of two and establish C2. Repeating the finite census at a larger cutoff would not fill either gap.
