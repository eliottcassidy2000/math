# Independent frozen audit: PASS

Subject: session203-binary-lift-obstruction/PROOF.md and accompanying verifier, counting/transfer code, and inherited positive-defects dependency. SOURCE_SHA256SUMS is copied from the frozen candidate; every entry was verified against that candidate. HASH_CHECK.txt records success. Candidate bytes were preserved.

Proof SHA256: 1174078417aa18fdab12e146d9b0d4781b35455e593dfbc491311bcf71b1233b.
Verifier SHA256: b8a8a1a9158d9eefa8f839522b252b3465ec394bb6e4ea7710ae0c05129f3e4e.

## Proof review

The stated obstruction is valid in the physical class G and for fixed positive integer c.

- The L<=n! injection is valid: concatenating the designated root-ending path and its ordinary companion preserves enough information to recover the cut from the unique root position. No concatenation edge is required because the target is arbitrary permutations.
- Weak-tail normalization is exact. A positive-start root lies in the initial strong component A. All later components must be in the ordinary bin of an L-cover, so h, s, L, and consequently B acquire the same tail factor q. Arbitrarily many singleton tail components cause no problem. With n=|A|, Camion's theorem (and the singleton boundary) gives qn<=h. Hence h<=c implies B<=2q n!<=2c!, as claimed.
- For weak inputs e=0. For strong inputs, n cycle cuts are distinct spanning paths, exactly one ending at the root, yielding h-e>=n-1. This proves the uniform bounded-J-defect implication. The gate implications correctly follow from their inherited defect lower bounds by h.
- The O-prefix identities are exact. R+b>B follows from the inherited strict Z doubling and e>=0. Therefore a positive O-prefix before Z produces >3B. Every two-nonneutral word produces >4B. The factorial exclusions apply even after an O-prefix because that prefix leaves B unchanged and stays in G.
- Thus every large bounded lift has exactly the literal form Z O^t. A nonempty terminal O run prevents the next bounded lift, while plain Z resets p0 and fixes delta. The requested indexing clarification is now explicit: B_0 is AFTER the first large lift, with p0=0, and the displayed ray formula holds for n>=0.
- The finite-seed argument does not require finitely many physical states below threshold. First-crossing scalar values lie in a finite interval. From each such scalar, the first large transition has finitely many scalar outputs; after it, the offset is an integer in {1,...,c}. These finitely many (scalar,offset) pairs determine all continuing rays. Infinite small-state variation and terminal O variants do not create extra scalar families. Consequently O(log M) counting and zero density follow.
- The F-certificate reduction has the correct offset 4j+d-2j_next, bounded in absolute value by D+32. Bounded shifts preserve O(log M). The conclusion remains a format obstruction, not a counterexample to F or to the complete finite-register presentation.

No substantive mathematical gap found.

## Verifier review and replay

Reviewed verify.py. Its require checks remain active under python3 -O. The finite state family is generated physically, not sampled from a relaxed cone. The word tests check arbitrary letters within their advertised length bound. Ray initialization is correctly after the first lift. Nonzero hostile/control counts prevent vacuous success.

Normal and optimized independent reruns are byte-identical to each other and to the author's NORMAL.json: 1,099 literal inputs, 772 positive inputs, 2,248 factorial checks, 1,542 gate checks, 432 physical states, 67,392 candidate word tests, 162 eligible large lifts, and 108 dead-terminal controls.

## Additional independent controls

- independent.py explores all distinct states from every five-letter physical word of length at most five, including repeated gates: 2,640 states. For c=1..4 it checks 1,628,172 candidate words of length at most three above threshold; all 171 bounded lifts have the asserted Z O* form. It also checks bounded non-Z defects and post-first-lift Z-ray formulas. Receipt: INDEPENDENT.json.
- weak_factorization.py reconstructs strong components independently on all positive-start labeled root-zero tournaments through order five: 772 inputs, including 201 weak inputs. It checks h/s/L/B factorization, q|A|<=h, the factorial bound, and weak e=0. Receipt: WEAK_FACTORIZATION.json.

These finite controls test implementation and boundary cases. The all-order obstruction rests on the proof reviewed above. No publication or repository edit performed by this reviewer.
