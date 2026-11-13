# Eight, nine, twenty-seven: sharp content, physical selectors, and retained clocks

2026-11-12 (America/Denver). **PROVED IN NOTES**, with independent reviews
and separately scoped **FINITE-EXACT** controls. **Chan–Pak remains OPEN.**
The numerals led to several different concrete objects; the useful synthesis
is through their operations and retained coordinates, not an identification
of unrelated state spaces.

## 1. The asymmetric9/27 gap is extremal, and repeated growth is classified

For ordered nonempty path-cover counts C_k(T), write

    h(T)=gcd_k C_k(T),
    g(T)=gcd_k (k−1)! C_k(T).

The inherited Session137/293 context theorem identifies g as the exact divisor of
all uniform-context Hamilton counts. The extra factorial permutes exterior
excursions; it is present even when the core has no automorphisms.

The new [content theorem](content/PROOF.md) proves

    h | g | h²/rad(h).

Thus excursion symmetry can thicken an existing prime but cannot introduce
one, and squarefree raw content cannot amplify. The asymmetric eight-vertex
example h=9,g=27 **attains the upper bound**.

More precisely, for each p|h let kappa_p be the first cover index where
v_p(C_k)=v_p(h), and put E(T)=product p^v_p((kappa_p−1)!). Then

    g(T^[r])=h(T)^r E(T)

for all r beyond an explicit finite threshold, at most the largest
v_p(h)−1 or1. Here T^[r] is the actual ordinal join of r copies.
For our core C this holds from r=1: h=9^r, g=3*9^r. Arbitrary singleton
padding preserves those values. At total padding m there are exactly
binom(m+r,r) distinct asymmetric members in this constructed family;
[their rational generating function](slot_polynomial/PROOF.md) is explicit.

The complete inherited bank through order7 has no amplification, so this
eight-vertex example is minimal in that finite classification. But8 is
not forced by a numerical factorial threshold alone: a formal degree7
polynomial passes the arithmetic tests and fails a physical triangle/cycle
realizability condition. The sharp general arithmetic lower bound is2p+1
for amplification at an odd prime p.

Incoming Session294 supplies a complete fourteen-class first-birth atlas;
Session295 independently supplies the power law with an exact least
threshold and the selective-port theorem below. These are credited and
combined in [INCOMING.md](INCOMING.md). The sharp content bound above is
an additional result.

## 2. Two actual arc flips form a universal mod-eight selector

For C with outgoing rows `(56,117,129,22,100,76,141,59)`, reverse1→6,
reverse3→2, or both. For **every fixed uniform exterior E**, put

    h=H(E[C]),  A=H(E[TT4]),  B=H(E[TT2]).

The four actual same-order counts are exactly

    h,  h+2A,  h+4B,  h−6A−12B.

They give every odd residue modulo8 exactly once, and all lie in the sharp
relative band

    (495/513)h <= output <= (517/513)h.

The mixed response is −8H(E[C3→TT1])=−8(A+2B); it is retained in the
formula. The individual positive responses hold for all uniform contexts
and for marked-cover coordinates whose marks stay exterior. With empty E,
the four counts are513,515,517,495 and their context divisors are27,5,1,9.
These are full response identities, not conclusions drawn from four totals.
See [constructors/PROOF.md](constructors/PROOF.md).

[Multiplication by9](BINARY_LIFT.md) then lifts this physical selector to
every odd residue modulo2^s. Powers of9 run through precisely1+8Z modulo
2^s; a bitwise algorithm finds the unique address. Ordinal joins with an
H9 graph and neutral padding put every residue at one common vertex order.
This is all-precision congruence coverage, not an exact-target theorem.
General scalar modular attainability was already known in Session137;
the new result is the explicit retained-exterior selector and lift.
An exterior cyclic p-module prevents a universal analogous selector modulo
any odd prime p, explaining the binary scope.

## 3. The27 obstruction also has exact neutral and reusable exits

A native switch on the Hamilton spine `(1,0,4,5,6,2,7,3)` preserves H513
but changes the context divisor **27→3**. This is optimal among every
count-preserving replacement: polynomial parity forces3|g(T) whenever
3|H(T). The reset response factors through the older five-vertex native
reset, and an explicit cyclic-wrapper family witnesses the reduced
valuation at every length.

Combining this physical reset with the incoming weighted-evaluation
identity gives an [exact repair ladder](content/REPAIR_LADDER.md): among
m ordinal blocks, replacing C by N one at a time keeps order8m and
count513^m fixed, while context-divisor exponents follow
`2m+1,2m,...,m+2,m`. The final two-step drop has a precise cause: the
first minimizing evaluator moves from4 to2 and loses its factorial3.
The2^m ordered block words are distinct isomorphism classes; each fixed
block multiplicity has the same entire uniform response.

To escape the last3, a named root-child program on the converse core
returns balanced primitive states for every k>=1:

    b_k=1746*2^k+241,       t_k=582*2^k−17,
    state=(b_k,t_k,t_k,b_k),       gcd(b_k,t_k)=1.

These strong graphs retain the original core as an induced subtournament
and enter the inherited repeatable C3/C5 controller with its guard intact.
The all-k coprimality proof is a genuine8/9 split: modulo73=8*9+1,
the powers of2 form a nine-element subgroup, and the forbidden residue28
has a nontrivial phase under x↦x^9. This is a finite character calculation
supporting an infinite physical continuation, not a primality claim about
the returned Hamilton counts.

A further [incoming six-return connection](constructors/SIX_RETURN.md)
keeps H17 and the same physical internal root while changing its state
`(17,8,9,18)→(17,2,9,44)`. Before this preparation, the specified return
program traps a factor9; afterwards it gives primitive balanced states
`b=105*2^k+7, t=35*2^k−11` for every k>=1. Every positive even k enters
the inherited repeatable controller. The root's identity is retained
during preparation; its indegree is allowed to change. This uses internal
marked response that the complete unrooted polynomial forgets.

## 4. Selective ports have an exact factorial precision budget

If at most r exterior vertices see individual core vertices nonuniformly,
the number h_k of Hamilton paths using k core runs is always divisible by

    ((k−r)_+)! ((k−r−1)_+)!.

The [proof](history/PORT_FACTORIAL.md) contracts exceptional excursions
into whole compound blocks, preserving their constrained seams, then
permutes the remaining free blocks. Modulo27 this makes every layer
k>=r+6 vanish. One selective vertex leaves only k<=6 relevant, regardless
of exterior size. The fully split eight-run layer needs at least three
selective vertices before it can break27; lower layers can escape sooner.

## 5. What the older square-sum and Collatz work actually contributes

| Object | Exact role of the small numbers | Retained coordinate |
|---|---|---|
| Square-sum path on1..15 | The9/16/25 reflections have endpoints8 and9: a deleted fixed point16/2 and an unmatched interval boundary. Their alternate orbit rotates by9 modulo16. | Square labels, endpoint type, host interval |
| Powers-of2/3 sumgraph | Consecutive targets8=2³ and9=3² compose to translation by1. | The two reflections, not a Collatz trajectory |
| Normalized slot polynomial | Degree<=8 needs at most modulus9 to observe its values mod3; degree9 can require27 in the integer-valued category. | Binomial coefficients and predivision precision |
| Actual asymmetric core | F(x)=x⁸+34x⁶+232x⁴+246x² is18 mod27 exactly at x=±4 mod9, and0 otherwise. | The second ternary digit of a unit square class |
| Collatz inverse observer | Roots7,25,43 have the same first menu mod9 but three different next menus after (16n−1)/3. A two-step chart needs27. | One additional input ternary digit and the chosen branch |

The [history note](history/HISTORY.md) recovers the September reflection
theorem, July graceful sign gauge, October indexed shared-clock work, and
the older hexagonal Collatz chart. The [slot proof](slot_polynomial/PROOF.md)
transfers the denominator-safe binomial-period theorem from the unrelated
Sun2/4/6/8 work. These transfers produced the new compound-block filtration,
precision observer, and physical selector; they do not assert Collatz
convergence or identify the orbit of the integer27 with the context ideal.

## Remaining target and reproduction

The unclosed target is a constructor selecting a prescribed exact odd count
while retaining a reusable source. Primitive interfaces, modular selectors,
and sharp positive responses now supply more structure, but no general
positive target selector or closure of the473 current repair families is
claimed. The large spectrum theorem remains open.

There is a sharp new application of an older obstruction:
[all uniform-context counts](constructors/LOW_CONTEXT.md) of a fixed
module M below `H(M)+F_M(2)` are exactly H(M) times actual tournament
counts. Thus even the primitive H517 source cannot reach3619 by any
uniform exterior. Selected-incidence changes or source replacement are
necessary to cross that boundary.

The stronger source-normalization hypothesis
`min_{H(T)=N} g(T)=gcd(N,3)` is explicitly OPEN. Among77 attained fibres
through order7,73 reach that lower bound and four need larger-order
constructions if it is true. Their exact witnesses and useful common-tail
responses are saved in the repair-ladder note; a bounded miss is not an
all-order counterexample.

The subsequent [bounded physical probe](content/NORMALIZATION_PROBE.md)
tested43,328 operations after source/sink padding through order10. All166
distinct neutral outputs retain the old divisor. This closes the declared
local operation menu, with the513 reset as a successful positive control;
it supplies no all-order negative conclusion.

Principal commands, from the repository root:

```sh
python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/content/check.py
python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/constructors/escape.py
python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/history/check.py
python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/slot_polynomial/check.py
python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/binary_lift.py
python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/constructors/six_return.py
python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/content/normalization_probe.py
```

The scripts write scoped receipts and can emit literal graph rows. See
[REVIEW.md](REVIEW.md) for independent proof review and finite universes.
No numbered canon ID or literature-priority claim is made.
