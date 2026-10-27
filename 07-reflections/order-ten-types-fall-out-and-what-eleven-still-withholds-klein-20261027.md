# The order-10 types fall out, and what eleven still withholds

klein, 2026-10-27 (overnight session; packet `05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/`, THM-4814, HYP-9422).

**The census was cheap once the registry was gone.** The order-10 prime cycle-index types took a
minute per part: nine types, and every row is forced by two one-line lemmas (a prime-order
automorphism of a prime tournament has at least two cycles of that length, because a lone cycle is a
module; an anti-automorphism has no odd cycle beyond a single fixed point). The anti-automorphism
data were then checked in the strongest way available: a self-converse decomposition identity that
returns the self-converse census `2, 2, 8, 12, 88, 176, 2752, 8784` at every order through 10.

**Eleven withholds exactly its primes.** The species identity carries the order-`n` primes as their
own `x^n` term, so "predict A000568(11) from the order-10 types" (HYP-9413 (d)) was a wish, not a
computation. What the identity does determine is the whole non-prime side, and with the known
A000568(11) and A002785(11) that turns into the first values of several sequences at 11: 796,876,429
prime tournaments, 259,307 of them self-converse, 447,114,461 intransitivity classes, and the exact
labelled counts with no OEIS input at all. The lesson generalises: a substitution species predicts
the composite layer of the next order; the atoms of that order must be counted or known.

**Two obstructions with different scopes.** The trace-norm separation (`(n−1)√n` for Paley against
`O(n log n)` for bounded banks) is real and I pinned the extremal side: the nuclear norm of a
tournament is maximised exactly by skew-Hadamard tournaments at even orders and by their one-vertex
deletions at odd orders; a source over `P_7` is one of the four maximizers at order 8. But for the
exception count the spectral bound is weaker than THM-4810's score bound on linear parses by a
`q^{3/2} log q` term, and it only bites for tree grammars in the thousands. Its real content is scope:
the score bound dies on regular tree parses (C3 towers), the spectral one does not. Neither
dominates; the dyadic-module session's uniform all-tree bound is the piece that makes the spectral
route meaningful.

**The black-core claim is a parity statement with a pure-blue exception.** From a regular
tournament, one native move changes a vertex's distance from the middle score by at most one, so a
singleton needs `m` moves; witnesses of length exactly `m`, all black, exist for every regular class
at orders 7 and 9 and for `P_11`. The claim's threshold `2m+1 >= 7` is not an accident of search:
the regular 5-tournament is the pure-blue atom `T5`, which has no black move at all.

**The Eisenstein remark generalises without effort.** `1+2ω = √−3` acting on the triangular lattice
is the `q = 3` case of `(I+S)/2` making the root lattice `A_{q−1}` a module over `Z[(1+√−q)/2]` for
every doubly regular tournament. The owner's three ROOT rows live in the `q = 3` case; the Paley
clocks of S63 and session74 have their own `O_{−q}`-modules of rank `(q−1)/2`.

**What I did not do.** No order-11 census (903,753,248 classes; the eight-part method would need
hours, and it would also deliver A280772(12)'s composite side). No reproduction of the `433 / 79`
numbers, which are not on main. No all-order proof of the black-core claim; the search is a
witness, the lower bound is the theorem.
