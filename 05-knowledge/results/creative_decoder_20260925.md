# Replacing the tournament decoder by a guarded join graph

**Date:** 2026-09-25. **Status:** PROVED elementary encoding, certificate transport,
scoped closure obstruction, and infinite certificate schemas; FINITE-EXACT
compiler/coverage experiments. Universal generation of a root certificate is
**OPEN**. This is an integrated construction from inherited arithmetic mechanisms,
with no claim that binary transducers or Collatz inverse cylinders are new.

## Inheritance, portfolio, and concept board

Closest proved mechanisms: the exact inverse fibres and sibling ladder in
[arithmetic braids, B1--B3](arithmetic_braids_20260917_collatz.md), the
[guarded inverse-completion theorem](arithmetic_braids2_20260917_inverse_completion.md),
and the [first-reset carry boundary](collatz_guards_20260921_valves.md).
The current tournament obstruction is the absence of pair modules in the
three-regular-block substitution, proved in
[the halving note](decoder_halving_20260925.md).
Canonical hostile: the positive minus cycles at 5 and 17. Corrected near misses:
relaxed zero-exponent inverse arrows need not be Collatz arrows, and local
residue coverage need not enlarge an inductively root-certified family.
Least-used sidecar: a certificate can use an intermediate sibling height;
quotienting to the smallest representative can discard that option.

Anchor: construct a sound structure-to-root certificate compiler. Niche:
transpose binary forward guards into ternary inverse guards. Wildcard:
complete a rule library by its suffixes and retain whole sibling ladders.

| Live concept | What the new object retains | Decisive boundary |
|---|---|---|
| Finite binary tape | actual integer, finite end marker, carry | scan termination is not orbit termination |
| Sibling ladder | all eligible intermediate representatives | canonical-only reduction loses 483's certificate |
| Affine port | exact word, guard, ordered carry, size threshold | a legal port need not decrease |
| Join graph | two routes to a common future, integer rank | equality without a smaller certified endpoint is insufficient |
| Inverse dual | ternary integrality, positive odd vertices | forward-closed families gain no new members |
| Root target | absorbing odd root 1; ordinary root 4 on expansion | minus cycle minima stay separate |

The selected research moves are “type every analogy and every implication” and
“compute the repair quotient before testing the residual defect.” Here the
second move has a concrete counterindication: keep intermediate representatives
when proof availability is not invariant under canonical reduction.

## 1. The object and its arithmetic semantics

Use positive odd vertices and

    U_sigma(n)=(3n+sigma)/2^v2(3n+sigma),   sigma in {+1,-1}.

A **guarded join** is a record `(x,y,a,b)` satisfying

    0<y<x,  x,y odd,  a,b nonnegative integers,
    U_sigma^a(x)=U_sigma^b(y).                              (J1)

It is an inference from a supplied certificate for y to one for x. It does not
assert that x follows the proof arrow to y. Proof arrows decrease ordinary
integer size; their associated arithmetic paths meet at a common future.
A finite collection of such arrows ending at 1 is a directed acyclic proof
graph. No complete pairwise relation or tournament orientation is required.

For the plus map, let h(y) be a certified number of odd steps taking y to 1.
Because U_+(1)=1, the exact transport rule is

    h(x)=a+max(h(y)-b,0).                                  (J2)

If b<=h(y), append h(y)-b steps after the common future. If b>h(y), that
common future is already 1. Thus h is a valid hitting-time bound, not
necessarily the first hitting time. The proof is well founded on the numeric
labels because y<x. The implementation rejects negative clocks, parity loss,
wrong sign, cyclic rank edges, and corrupt word metadata.

For each odd input greater than 1, expanding a root-1 certificate into the
ordinary map `n/2` or `3n+1` passes through 4 before 1. Input 1 reaches 4 in one
ordinary step; even inputs first halve to their odd part. This separates the
odd accelerated clock from the ordinary root-4 clock.

The [finite-word carry machine](creative_transducer_20260925.md) supplies local
arithmetic semantics: finite least-significant-bit-first words have an end
marker; halving deletes an initial zero; an odd step uses six carry transitions
and three terminal rules. Each individual scan terminates and exactly evaluates
the required integer map. Repetition of scans has no proved global rank.

## 2. A quotient must retain certificate alternatives

Put S(n)=4n+1. The inherited identity is U(S(n))=U(n). Consequently, if

    U(x)=S^j(y),  y<x,

then `(x,y,2,1)` is a legal join. Keep **every** eligible y along the sibling
ladder, as well as the direct forward rule `(x,U(x),1,0)` when U(x)<x.
Keeping only the smallest sibling is weaker.

A decisive example from the [sibling audit](creative_sibling_20260925.md) is

    U(483)=725,
    725 -> 181 -> 45 -> 11       (successive sibling cancellations),
    U(725)=U(181)=U(45)=U(11)=17.

The intermediate 181 has the decreasing certificate

    181 -> 17 -> 13 -> 5 -> 1,

where 17->13 is a direct odd step and every arrow is legal. The smallest
representative 11 is absent from the canonical-only inductive family. Using
181 gives the concrete join `(483,181,2,1)` and the actual trajectory

    483 -> 725 -> 17 -> 13 -> 5 -> 1.

In the exact universe of odd inputs at most 100000, direct decreasing forward
rules certify 252 inputs, retaining direct plus canonical sibling rules certifies
541, and retaining all sibling heights certifies 640. These are sizes of these
specific inductive families, not counts of Collatz-convergent integers. The
sibling note proves the mechanisms and gives generated infinite descendants.

## 3. Forward binary ports and inverse ternary ports

For an exact positive exponent word w=(k1,...,kL), put K=sum ki and

    C=sum_(i=0)^(L-1) 3^(L-1-i) 2^(k1+...+ki).

The empty prefix sum is zero. The actual odd-map identity is

    2^K z=3^L x+C.                                         (P1)

For positive odd endpoints, total integrality plus oddness verifies all
intermediate guards, by backward reduction modulo 3. No exponent zero is
permitted. This is inherited from the inverse-completion theorem.

**Forward port.** The source is one residue class

    x=(2^K-C)*3^(-L) mod 2^(K+1),
    z=(3^L x+C)/2^K.

When 2^K>3^L it strictly decreases exactly for
`x>C/(2^K-3^L)`. It supplies join `(x,z,L,0)`. Every verified decreasing
seed therefore produces an infinite arithmetic progression with a reusable
certificate of descent, rather than just one recorded trajectory.

**Inverse port.** The current vertex z is one odd class modulo 2*3^L, specified by

    z=C*2^(-K) mod 3^L,
    x=(2^K z-C)/3^L.

If 2^K<3^L, every positive odd integral endpoint has 0<x<z; reverse peeling
of each integral branch establishes positivity and exact exponents. This gives
join `(z,x,0,L)`. Its guard is ternary and its size contraction comes from an
expanding forward word. On the minus sheet replace C by -C and check the actual
size inequality; do not import the plus inequality without its sign.

For example w=(1,2) has K=3,C=5. Exactly the odd z=13 mod18 satisfy

    x=(8z-5)/9 < z,
    x -> (4z-1)/3 -> z.

At z=31 this reads `27 -> 41 -> 31`, giving rank reduction 31 to 27.
Combined with the sibling and one-step inverse rules, this removes classes
31,103,139 modulo144 from the former residual and gives local descent on
19/24 of odd residue classes. The residual classes modulo144 are

    7,15,27,39,43,55,63,75,79,87,91,111,123,127,135.

The finite bank of all words with L<=8 and 2^K<3^L contains 953 words. Exact
residue-union counting gives inverse coverage 1013/2187 and combined local
coverage 2329/2916 among odd integers. These fractions are FINITE-EXACT
computations of finite unions of arithmetic progressions, not orbit frequencies
or a density-to-convergence inference.

## 4. Why the inverse improvement cannot by itself finish the proof

**PROVED closure obstruction.** Let A be the least root-certified family
containing 1 and closed under the direct decreasing forward rule and the
eligible sibling joins in section 2. This includes the canonical-only variant.
Then U(A) is contained in A.

Induct on a certificate derivation. For the direct rule the successor is already
certified. For a sibling join x to y, write U(x)=S^j(y). By induction U(y) is
certified. If j=0, U(x)=y is certified. If j>=1, then

    U(U(x))=U(y)<S^j(y)=U(x),

since U(y)<=(3y+1)/2<4y+1. Therefore the direct decreasing forward rule
certifies U(x). A source-sibling cancellation obeys the same argument, or is
redundant in this least family.

It follows that adding **any finite inverse-orbit join**, of any length, cannot
enlarge A: if a certified y reaches x, forward closure already certifies x.
This is stronger than a failed finite search. The extra ternary rules may shorten
proofs, change normal forms, or help an incompletely saturated library, but cannot
supply the missing root entry to this closed family.

The exact experiment confirms the distinction. Among odd inputs <=10000 the
full sibling family has 154 certified members. All 953 inverse ports still give
154, despite improving the local descent coverage. A vertex reducing to a smaller
uncertified vertex is not yet a root certificate.

## 5. Completing a forward rule library by its suffixes

The implemented adaptive compiler processes odd inputs in increasing order.
It first tries all available guarded joins. If none decreases, it searches the
actual odd orbit until the first value below the input, with an explicit finite
experiment cap. That successful word is generalized using (P1).

Retain all its suffix words. Before the first descent, every intermediate value
is at least the starting value; therefore every suffix also ends at a smaller
value. Its positive carry forces its slope below 1. Each suffix is thus another
proved infinite descent port with its own exact guard and carry threshold.
This is proof-library completion, not a claim that the seed search always ends.

**FINITE-EXACT:** for odd inputs through 10000, all 5000 root certificates were
constructed. Without suffix completion the run learned 442 new seed obligations.
With suffix completion it learned 406, producing 3016 distinct infinite forward
ports. The same 406 seeds suffice when the inverse bank is omitted. The maximum
compiled bound is 96 odd steps. An independent ordinary-map implementation
verifies every input 1..10000 reaches 4, with maximum 259 ordinary steps.

This experiment tests the compiler and proof reuse; it is not advertised as a
new numerical verification bound for Collatz.

**PROVED finite-library obstruction.** Every finite collection of contracting
forward words has a maximum odd length L. Choose any even a>L+2 and n=2^a-1.
Its first L odd exponents are all 1, so no contracting forward word in that
bank applies. Also n is divisible by 3, hence is not the endpoint of any
nonempty odd Collatz word and has no inverse port of any length. For a>=4,
neither n nor U(n) admits a sibling cancellation, and U(n)>n. Thus these
local rules leave n unresolved.

The constructed bank has maximum word length 51; 2^54-1 is an explicit
verified unresolved input for this bank, and all larger even exponents give
further examples. This is a limitation of the bank, not a divergent orbit.
It identifies why an unbounded rule schema is necessary.

## 6. An unbounded schema that actually produces root certificates

The [descent lane](creative_descent_20260925.md) solves a reset port backwards.
For any integers a>=2,t>=0, set

    b=3^(a-1)(2t+1),
    N=2^a(2^b+1)/3^a-1.                                   (R1)

Here 3^a divides 2^b+1: the order of 2 modulo 3^a is 2*3^(a-1).
The quotient is positive odd. With C_+ the shortcut map, the exact itinerary is

    a odd steps, then b even steps, ending at 1.

Indeed `C_+^a(N)=2^b`. Since b>=3, it reaches 4 after exactly a+b-2 steps
along this certificate. For t=0 the first examples are

    a=2: N=3,          b=3;
    a=3: N=151,        b=9;
    a=4: N=26512143,   b=27.

This constructs certificates before tracing the starting orbit. Given any fixed
lookahead D, choose a>D: the source initially rises for more than D shortcut
steps, yet (R1) provides its complete root certificate. It is a specialized
explicit use of the inherited inverse-completion theorem, not a new universal
convergence theorem. The lane proves the exact congruence criterion, a target-v
version, and signed controls; it also checks symbolic certificates with millions
of binary digits without expanding millions of individual steps.

## 7. Typed connection and the remaining obligation

Source: finite binary tapes, guarded affine ports, and sibling-ladder joins.
Target: finite ordinary Collatz routes to 4. Map: evaluate the tape, verify each
word/guard, compose the common-future records by (J2), and expand the odd clock.
Preserved predicate: exact positive-integer trajectory and membership in the
specified root basin. Lost information: raw graph geometry and ordering of
irrelevant pairs; compressed words also omit the displayed intermediate nodes.
Required sidecar: sign, integer domain, ordered carry, actual exponent guard,
clock, size rank, and terminal marker. Cheapest decisive tests: the 483
intermediate sibling, the 31 inverse port, and the minus minima 5 and 17.

The remaining claim is specific: find an unbounded family of guarded forward
joins which, together with the retained sibling alternatives, gives every odd
integer greater than 1 a finite proof ending at a smaller already-certified
integer. A richer topology alone cannot supply that claim. Inverse saturation
and finite word libraries have the explicit obstructions above. The reset
schema proves such entry for infinitely many unbounded-run inputs, while its
uncovered cofactor cases remain OPEN.

No implication to graceful-tree or square-sum Hamiltonicity is established by
this construction; their separate label and adjacency predicates remain necessary.

## Reproduction and audit

```
python3 04-computation/experiments/creative_decoder_20260925.py
python3 -O 04-computation/experiments/creative_decoder_20260925.py
python3 04-computation/experiments/creative_decoder_audit_20260925.py
```

The [frozen output](creative_decoder_20260925.out) records the exact universes,
coverage fractions, learned-port counts, 11907 independent port-guard samples,
separate ordinary-map verification, and hostile certificate rejections.
The transducer, sibling, and descent lanes have separate scripts and outputs.
The [independent audit](../../04-computation/experiments/creative_decoder_audit_20260925.py)
checks 87040 word/domain/sign cases with separate repeated-division arithmetic
and 19991 signed joins. Its [output](creative_decoder_audit_20260925.out) records
eight rejected malformed certificates. The audit caught and repaired an even
inverse endpoint and a negative-clock acceptance hole before publication.
All checks use explicit exceptions and survive `python -O`.
