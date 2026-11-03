# Scope, exact universes, and proof review

Status labels: **PROVED scoped** means the all-order elementary proofs N1--N4, C1--C3, B1--B5 and D1--D3 in their named files. **FINITE-EXACT** means the bounded checks below. The proofs and code were reviewed by the same agent; separate algorithms are not an independent-person audit. No literature-priority claim is made.

## Reproducible universes

| Script | Universe | Separate path / decisive control |
|---|---|---|
| normal_forms.py | All8008 ten-symbol multisets of length<=6, all22880 applicable competing reductions | Canonical prime-exponent formula versus iterative rewriting |
| normal_forms.py |729 triples(c,d,e) in0,...,8; every integer1,...,10000 | Actual first-applicable fractions versus exact trace formula and strict height |
| normal_forms.py |3969^k, k=1,...,30; all reduced rational elasticities with denominator<=12 | Complete prime-exponent enumeration versus length/elasticity formulas |
| bank.py | Ten bank sources plus strong25, orders3,...,9 | Endpoint subset DP, unordered path-partition enumeration, independent ordered two-bin convolution |
| bank.py |45 pairs of private polynomial factors | Explicit integer Bezout identities |
| bank.py | Every ten-letter word of length0,...,4:11111 words | Direct gcd/product arithmetic versus minimal DFA; signed and unsigned counts through length12 |
| bank.py | Seven-state primitive DFA and thirteen-state refinement | Reachability, one-letter distinctions, full partition refinement |
| bank.py | All14 bit words of lengths1,...,3 | Literal strong wrapper, unique port degrees, SCC recovery after reversing numeric labels, explicit self-converse permutation |
| bank.py | Strong25 triangle;35 joint clock phases; eight7-adic lifts | Literal H329 DP and trial division47; modular powers and six test indices per lift |
| controllers.py | All16 four-bit guards;16 marked charts;49 hub routes;70 primitive selector transitions | Lcm criterion and literal native pair reversals; exact black-color and score checks |
| controllers.py | Eight arithmetic inputs including inert101 factor | Full integer/leaf/phase/native-clock trace |
| neutral_carries.py | All128 odd P from1 through255 | Independent full response formula versus closed carry law, exact falling-basis reconstruction, separately implemented native replay |
| neutral_carries.py | Literal P=1,3,5,7,11 sources | Endpoint DP and two-bin convolution; full path covers at1,3,5 and literal G at1,3,7 |
| neutral_carries.py | Five G iterations of the H497 pair | Nonnegative falling coefficients and exact scalar improvements |

`reproduce.py` runs all four scripts normally and under `-O`, checks byte-identical JSON receipts, and verifies pinned dependency hashes. It records interpreter/SymPy versions and receipt hashes. No network is needed.

Source hashes normalize CRLF to LF, since some inherited JSON paths lack an explicit checkout newline rule. This packet writes its own receipts with LF and pins that checkout convention locally.

## Positive and hostile boundaries

- Positive: the nine-vertex H497 pair is joined by the named two-move reset; its triangle continuations differ by48. The order-ten P5/P6 pair preserves the entire profile instead.
- Positive: one strong H25 source gives literal H329=7*47, and supplies an all-order recurrence with a uniquely lifting7-adic clock.
- Positive: every rational elasticity in[1,5/2] occurs on actual counts with prime support{3,5,7}.
- Hostile: the original ten-source G family never gives the required7p values; the repair changes the source library, not an abstract factorization.
- Hostile: equal divisor guards at25 and625 have different next guards. The finite selector does not replace unbounded arithmetic counters.
- Hostile: full-profile equality can hide an arbitrarily long ordered source register. Its native/rooted behavior must not be inferred from that profile alone.
- Inherited, pinned, **not rerun here**:81/27, common H462915 and differing named-port continuations;768 one-birth H55 controls with no H111. These remain source-identity boundary controls.

## Logical audit

1. **Spectrum:** the monoid M is an explicitly defined arithmetic object. Its global identification with actual H remains conditional. Local prime support{3,5,7} is fully witnessed, so N4 is unconditional.
2. **Irreducibles versus sources:** a numerical composite can be a useful new attachment source. Strong25 is not replaced by5 ordinal5 in a composition proof.
3. **Classes:** the rational series count injectively decoded subfamilies, not all tournament classes or all H fibres. SCC uniqueness is used on actual strong blocks.
4. **Controller:** black-line compilation changes the finite logical state; the integer and the source word remain explicit external coordinates. The hub's own H is not the represented integer.
5. **Clock:** the native controller has two moves per arithmetic step; the source reset also uses a birth/deletion, with the entire vertex clock recorded. H-neutrality refers to endpoints.
6. **Ports:** source row matrices and named apexes remain in reset receipts. Homogeneous dominance does not imply arbitrary marked-port dominance; no claimed F49/Q4 quotient discards the literal source.
7. **Ten and thirteen:** no theorem is inferred from equal cardinalities of unrelated state sets. K2,3 and the prior Petersen object are distinguished explicitly.
8. **All-size:** exact formulas and proofs provide the quantifiers; finite tests do not extrapolate them. Prime quotient47 is certified; no infinitude of prime terms is claimed.
9. **No open-problem overreach:** this packet proves no general Collatz convergence, Sha finiteness, LRC(14), exact full spectrum or full A000568 recurrence.

The only implementation repair during construction was ordering SCCs by arcs against a frozen copy of the component list: Python's in-place sorting hides its list during key evaluation. Relabelled recovery tests were added before any theorem relied on that implementation. This was an unpublished code issue, not a retraction of a canon result.
