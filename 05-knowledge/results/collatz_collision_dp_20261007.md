# Exhaustive finite collision search at adaptive depth

**Status: PROVED elementary finite grammar and certificate transport;
FINITE-EXACT controls. Universal first-reset-2 coverage OPEN.**

The new selector searches the entire partner language of each inspected
Collatz prefix, rather than comparing that prefix only with a fixed rule
bank. It finds the new 1459 Mersenne collision automatically. A frontier
status remains explicit when its letter or cost budget runs out.

## Inheritance and research portfolio

The closest proved mechanism is
[THM-4555, uniform switches are collisions at minus one](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md),
especially its equal-time trailing-ones deletion rule. Its effective
fixed-length collision-core enumeration is prior work; this note uses a
different finite coordinate, the exact common denominator, to produce a
direct dynamic program. No new general collision principle is claimed.

The canonical hostile is the positive source `2^1459-1`, missed by the two
rules in [terminal lifts](collatz_terminal_lifts_20261007.md). The corrected
near miss is an apparent join obtained only by continuing the fixed point
1: exported certificates stop at its first hit. The least-used relevant
coordinate is the set of *all attainable partner lengths*, not only the
longest word: the available source run limits the legal deletion.

The concept board is: (1) rational collision anchors; (2) exact word-cost
gradings; (3) integer source guards; (4) first-hit terminal receipts;
(5) frieze coordinate changes; (6) Mersenne exponent phase cylinders.
The anchor is the residual branch; the niche is a complete finite language
parser; the wildcard is frieze ear insertion. The companion
[frieze analysis](frieze_collatz_reset2_20261007.md) tests that wildcard
against actual endpoint congruences. We use the meta-pattern
"Turn certificate failure into an address, then change sidecars": the
added denominator/length data eliminates a named rule-bank failure, but
does not imply every address will eventually be eliminated.

## A complete grammar for one anchor

Write `f_a(x)=(3x+1)/2^a`. A reduced word has first letter at least 2 and
all later letters positive. Its value at -1 has the unique form

    f_w(-1) = N/2^E,     N odd, E=sum(w)-1 >=1.

Let `L(N,E)` be the set of reduced words having this value. The following
recursion is exact:

* The one-letter word `(E+1)` belongs exactly when `N=-1`.
* Every longer word is uniquely a prefix `v` followed by a letter `E-f`,
  for some `1<=f<E`, where `3` divides `N-2^f` and

      v belongs to L((N-2^f)/3, f).

Indeed, if the predecessor is `M/2^f`, appending `E-f` gives
`(3M+2^f)/2^E`. Conversely each displayed predecessor and final letter
gives exactly the desired fraction. Once the first reduced letter has
occurred, every numerator stays odd and every denominator exponent stays
positive; there are no omitted integral intermediate states. Since `f<E`,
the recursion terminates. This proves completeness, including negative
numerators, without an individual upper bound on the valuation letters.

There are at most `2^(E-1)` reduced compositions of total cost `E+1`.
Memoization shares repeated states, but no polynomial-time or uniform
small-memory claim is made. The implementation retains one lexicographically
least witness for **each attainable length**. This is a lossless quotient
for the existence predicate "a partner has one of these lengths"; it
forgets other witnesses, their guards as separately presented words, and
any optimization objective other than the retained length and tie-break.

## From a parsed word to a paid receipt

Suppose the actual source prefix is `1^r u`, with `u` reduced. For every
attainable partner length `p` satisfying

    len(u) < p <= len(u)+r,

take the stored partner `v`, put `D=p-len(u)`, and set

    child = (source+1)/2^D - 1.

THM-4555 gives equal endpoints for `1^r u` and `1^(r-D) v`. The source's
actual prefix, not rational equality alone, supplies the native guard.
The child is positive, odd, and strictly smaller. The checker independently
replays both paths and compares their affine carriers.

ROOT requires a small distinction. The algebraic words may append actual
2-loops at 1, as for source 3 and child 1. The implementation verifies that
every discarded letter is exactly 2, removes that terminal suffix, and
exports first-hit paths. When a positive endpoint exceeds 1 no such
trimming occurs. If the source itself reaches 1 before any usable collision,
the selector exports that complete ROOT word directly.

`probe(source,max_letters,max_cost)` consumes successive actual reduced
letters. Each inspected prefix receives a complete finite partner search.
The explicit outcomes are `PAID`, `ROOT_CERTIFIED`, `ROOT`,
`NO_TRAILING_ONES`, `LETTER_FRONTIER`, and `COST_FRONTIER`.
`PAID` is a dependency on the smaller child, not a claim that the child
already reaches 1. The direct `ROOT_CERTIFIED` outcome includes its word.

For any fixed source and any finite guarded collision prefix with
`1<=D<=r`, sufficiently large budgets will inspect it and find a partner.
This is completeness **within this collision language**. It does not
establish existence of such a prefix for every source or termination of
unbounded budget escalation. For example, the bounded search certifies
127 directly at ROOT without finding a collision before that first hit.

## The new 1459 rule is found without a hard-coded pattern

At source `2^1459-1`, the first successful reduced prefix is

    u=(2,3,2,2,2,2,1,5,1,1,2,6).

The complete parser chooses

    v=(2,2,2,1,2,2,2,1,1,3,1,1,1,8).

Their common anchor is `22777543/2^28`; their lengths are 12 and 14.
Thus `D=2`, and the child is `2^1457-1`. At 11 reduced letters the selector
reports a letter frontier; at cost 24 it reports a cost frontier. The
independently discovered rule, its infinite exponent guard, and explicit
child grounding are in [reset-two rules](collatz_reset2_rules_20261007.md).

This is genuinely beyond the two inherited dispatch rules: its first
reduced letters are `2,3`, whereas the ordinary reset needs first letter
at least 3 and the old sporadic rule needs `2,6`.

## Arbitrary ternary depth after the binary collision

The incoming A2 correction to
[THM-4594, maximal class-decided sieve](../../01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md)
restricts its sign barrier to unrefined binary classes. Its negative-cycle
inverse branch gives a productive additional operation, not a no-go:

    G(m)=(8m-5)/9,   m=4 mod9,
    U_(1,2)(G(m))=m,
    G(m)+5=(8/9)(m+5).

For every positive odd m and integer r>=1, **r consecutive inverse turns
are legal iff `9^r` divides `m+5`**. Necessity and integrality follow from

    G^r(m)=8^r(m+5)/9^r-5.

At each intermediate step the child is an odd integer. For positivity,
`m+5` is an even positive multiple of `9^r`, so the final child is at
least `2*8^r-5>0`; the same argument applies at earlier turns. Each turn
strictly decreases a positive input since `G(m)-m=-(m+5)/9<0`.
The actual word from the final child to m is `(1,2)^r` and has no ROOT
padding: its fixed final letter 2 could reach 1 only from a nonintegral
predecessor for the complete block. Direct substitution also verifies
the exact first and second valuations.

Thus any checked common-future receipt from original n to child m can
replace its child by `G^r(m)` and prepend `(1,2)^r` to the old child word.
The original n and its actual route remain fixed. `negative_five_lift`
implements and independently replays this operation. It strengthens a
smaller-child dependency; it does not by itself ground the new child.

For the newly covered Mersenne phase, put `m=2^(K-2)-1`. Elementary
3-adic lifting gives

    v3(m+5)=v3(4(2^(K-4)+1))=1+v3(K-4),
    maximum r=floor((1+v3(K-4))/2).               (mixed guard)

This depth is genuinely unbounded on the repaired binary phase. For any
r>=1 and u in {1,2}, impose

    K=1459 mod2^21,
    K=4+u*3^(2r-1) mod3^(2r),     K>=1459.

CRT supplies an infinite arithmetic progression with exactly r available
inverse turns. `mixed_phase` computes its first member and period without
expanding a Mersenne integer. The resulting all-height paid child is

    h_r=-5+8^r(2^(K-2)+4)/9^r < 2^(K-2)-1 < 2^K-1.

At K=1459, exactly one turn is allowed, giving `h_1=(2n-11)/9`.
The same is true on the previous B8-disjoint subray K=1mod1458: there
`v3(K-4)=1`. Allowing different ternary phases is what makes r unbounded;
we do not falsely assign those depths to that single old subray.

This connects the local collision, negative-cycle affine anchor, and
two independent prime-adic guards while retaining their shared integer K.
It is stronger than a fixed core at arbitrary initial-run length, but
still covers only the displayed mixed progressions, not all reset-two
sources. The next child generally leaves the Mersenne form; recursive
closure is a separate obligation.

## Finite controls and unresolved families

The independent comparison enumerates all 16,383 reduced compositions of
total cost 2 through 15 and all their 7,602 distinct cost-indexed anchor
values. Exact Fraction evaluation agrees with the integer recurrence, and
the complete length-to-witness dictionaries agree with brute force.

For 128 sources `n=3 mod4`, `3<=n<=511`, with eight reduced letters and
cost cap 22, the outcomes are:

| Outcome | All 128 sources | Their 64 first-reset-2 sources |
|---|---:|---:|
| Paid smaller child | 66 | 2 |
| Complete ROOT word | 8 | 8 |
| Letter frontier | 53 | 53 |
| Cost frontier | 1 | 1 |

This is a bounded parser/control census, not new verification of these
small integers and not a universal coverage claim. The prior finite atlas
already grounds every one of them. Forged source letters, incorrect anchor
rewrites, valid abstract identities applied to the wrong source, boolean
and float aliases, and ROOT padding are all tested.

For pure-two reduced prefixes of lengths 1 through 16 the longest
partner has the same length. This is FINITE-EXACT only; no all-length
noncollision theorem is inferred. Odd Mersenne exponents close to 1 in
the 2-adic metric create arbitrarily long pure-two prefixes. These remain
a useful depth stress test even after the 1459 phase is repaired.

## Reproduction

Run `python -B 04-computation/experiments/collatz_collision_dp_20261007.py`
and repeat with `-O`. The [script](../../04-computation/experiments/collatz_collision_dp_20261007.py)
uses exact integers, an independent Fraction path, literal Collatz replay,
and the previously checked receipt auditor. Its [saved output](collatz_collision_dp_20261007.out)
records 119,313 explicit checks. This includes four literal mixed-depth
controls and 72 symbolic mixed-phase instances at depths 1 through 12,
both nonzero ternary digits, and three parameter lifts. No `assert` statement supplies a
load-bearing test, so optimized execution retains all checks.
