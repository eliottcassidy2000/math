# The summand alternation is an involution parity law

**kind-pasteur-2026-07-26-S131.** Provenance note, not truth source.

The owner's dispatch asked three things at once: understand which
previous terms combine to make each new A014574 term; understand why
the summand graph "as it goes to infinity is an alternating object";
and understand how the merged metagraph maps sets of tilings to
iso-class nodes with blue/black lines between tiling pairs. Working
them side by side exposed one mechanism appearing in all three, and
it is the repo's oldest one: **pair the objects by their natural
involution and read the fixed points**.

## 1. The alternation, three ways

**Integers.** The full summand fibre of `z` has
`r(z) = floor((z-1)/2)` unordered distinct pairs (THM-2422): the
period-two boundary quasipolynomial. The parity flicker between
`z = 2m` and `z = 2m+1` is exactly the presence or absence of the
diagonal pair `(m, m)` that the distinct-sum rule deletes. The
involution is the swap `(a,b) -> (b,a)`; even `z` is where it has a
fixed point.

**Twin ranks.** After thinning `N` to `K = A002822`, almost all of
that structure dies -- but the parity bit survives verbatim
(THM-2443): the ordered parent count `R(k)` is odd iff `k` is even
and `k/2 in K`. The dense alternation degenerates to a sparse
indicator, the doubling subsequence `{k : k, 2k in K}` (2,547 members
below 1.67e7). Center `12 = 6+6` is its startup avatar: the unique
center whose *only* parent pair is the diagonal one -- the same
special role the owner noticed for even nodes needing their half.

**Tilings.** The line space `L_n = X_n / <tau>` is paired by the
grid reflection `sigma`; black lines are the free orbits (hence even,
defect-by-defect), blue lines are the fixed points (hence a T-join on
the SC nodes): "blues odd, blacks even" (THM-643/THM-796). And one
level up, the same law decides pure-blueness: an SC class with a
one-element tiling fibre has `sigma` acting on a singleton, so the
tiling is fixed, so the class is pure-blue (THM-2444's rigidity
lemma). The n=9 refutation of the interleave formula came from
exactly this: the interleave was a shadow of the *count of fixed
singletons*, and that count jumped by two at n=9 (`H = |Aut| in
{1,3,9,9,27}`, all powers of three).

So: the "alternating object" the owner sees in the growing summand
graph, the odd/even blue/black law in the metagraph, and the parity
of twin-rank parent fibres are one lemma instantiated three times.
The productive move in each case was not the parity count itself but
**promoting the fixed-point locus to the next object of study**:
diagonal pairs -> doubling subsequence; sigma-fixed tilings -> blue
sub-cube `Fix(sigma) = Q_e`; singleton fibres -> rigid-SC classes.

## 2. Which previous terms combine (the quantitative answer)

The owner's recurrence reading of A014574 is the self-gap condition
`g_i in K` (THM-2422 eq 30), true on ~42% of transitions and
window-decaying. The genuine law is two-layered:

- **Existence** is carried by the full parent fibre: every rank
  `k >= 3` through 1.67e7 has distinct parents, with monotone-growing
  dyadic margins (min ordered count 5,176 in the top window). Under
  the open all-n statement (37) this is exactly what converts twin
  infinitude into Seymour's boundary-crossing form -- (37) makes the
  new OEIS conjecture *equivalent* to the Twin Prime Conjecture
  (THM-2443).
- **Selection** is singular-series-ranked, not size- or
  recency-ranked: the favoured partners are the gaps `g` maximizing
  `prod_{p|g} (p-2)/(p-4) * prod_{p|9g^2-1} (p-3)/(p-4)`, i.e. the
  primorial-adjacent centers 30, 42, 60, 12, 180 (HYP-9025; the raw
  gap histogram's wild oscillation flattens to a smooth decay under
  this normalization). "Which previous terms combine" has the same
  answer as every twin correlation question: the smooth ones.

The owner's `{2,3} -> 5 -> Fibonacci` intuition sits at the
degenerate end of this: in the dense integer graph the closure from
`{2,3}` is cofinite with hole `{1,4,6}` (THM-2422), and the
Fibonacci chain is the minimal-growth spine of the "must use recent
terms" restriction; in the thin twin graph the spine dissolves into
the singular-series selection above, but the startup holes `{4, 6}`
(the p=3-anomalous centers) mirror the integer holes `{4, 6}`
exactly -- THM-2433's deletion calculus is the bridge.

## 3. How the metagraph maps tiling sets to nodes (the owner's Q3)

Compressed answer, all proved: a tiling is a tournament with a marked
Hamiltonian path; the node map forgets the path and folds converse,
`pi : X_n -> M_n`. The inverse image of a merged node `u` is
`(2 - [SC]) * H(c)/|Aut(c)|` tilings (LEM-003 fibration; global
checksum `sum H/|Aut| = 2^m`). Lines are tau-orbits `{t, t+1}`;
colour is sigma-fixedness of the orbit; a line's boundary is the
unordered node pair `{pi(t), pi(tau t)}`, loops legal. Everything
else canon knows -- pendant law, T-join/even-graph split, parity of
cross-degrees, the failure of node-level Markov compositionality
(THM-796 SS8: nodes are base addresses, not stalks) -- follows the
same involution grammar. The open edges after S131: closed form for
rigid-SC(n) (data 2,1,2,2,3,3,5), whether (15,5,3) stays the unique
nonrigid pure-blue at all odd n, and n=10.

## 4. Stopping note

The parity mechanism is classical (Burnside/involution pairing); what
is repo-new today is the two fixed-point loci it opened (doubling
subsequence; rigid-SC classes) and the refutation it forced. Not
proposing a META-PATTERNS card yet: two of the three instances came
from the same session and the card would need counterindications
(where involution parity says nothing -- e.g. the black self-line
count at n=8 broke the naive SC/2 law, THM-849).
