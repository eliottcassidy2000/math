# Source switching at 55 and105, and singular golden words

2026-11-10. **PROVED scoped mechanisms, independently reviewed in this
packet; FINITE-EXACT target evaluations; full Chan–Pak conjecture OPEN.**

The target remains that the tournament Hamiltonian-path counts are exactly
the positive odd integers other than7 and21. This is Conjecture6.2 in
[Chan–Pak, section6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
which credits the2016 conjecture to bof and Gordon Royle. This session does
not claim the complete proof. It identifies a sharper obstruction to one
proposed induction and supplies concrete ways to move beyond restricted
source and reset languages.

## The decisive difference between55 and105

The [complete source classification](targets/PROOF.md) gives:

| Parent count | Sources supporting3h+2 | Sources supporting3h+4 | Common actual source? |
|---:|---|---|---|
|55|Two strong order7 cores, with arbitrary exterior padding|The same two cores|Yes|
|105|Two strong order7 cores and one order8 core, with exterior padding|One different strong order8 core, with exterior padding|**No, at any order**|

The two other strong H55 cores are exactly the inherited golden H55
examples. Neither supports111,167 or169, even after arbitrary exterior
transitive padding. The two active H55 cores each support all three.
This upgrades the earlier golden-language obstruction to an all-order
source classification using the inherited complete strong-child inventory.

The mechanism retains the actual source. For a coprime parent/child pair,
every SCC outside the new vertex's strong component is a singleton.
Delete the new vertex and remove only the source's *exterior* singleton
chains; retain every internal singleton gap. The resulting core and a pair
of padding requirements form an exact all-order attachment certificate.
Two branches have a common source precisely when their core sets intersect;
coordinatewise maximum padding then supplies that source. At105 the
intersection is empty. Independent canonicalizations and path counters
confirm every relevant deletion in the declared complete inventory.

The [previous threshold95 rule](../decomposition16_odd_spectrum_20261110/induction/PROOF.md)
explicitly allowed different sources for the two branches. It remains
OPEN and is not refuted by105. The stronger rule demanding both branches
from one current tournament is refuted. The missing coordinate is source
identity and attachment compatibility, not the scalar count105.

## A concrete reset that changes the source

Inherited [THM-1900(B)](../../../01-canon/theorems/THM-1900-the-insertion-response-calculus.md)
says a count-neutral one-vertex birth merely inserts a singleton SCC.
Consequently any sequence of neutral births/deletions and root changes
preserves the ordered nonsingleton SCC word. Those moves cannot switch
the105 branch families.

An actual two-replacement macro does switch them:

    105 -> 19 -> 105 -> 9 -> 105.

The starting source supports317 and the final source supports319. Six
original vertices survive, with all arcs between surviving vertices kept
at every step. The [reset proof and certificate](resets/PROOF.md) give all
rows, cuts and name maps. Any pair of those six names can be retained.
The four simultaneous-deletion counts for one retained pair change as

    (105,31,33,9) -> (105,25,9,1) -> (105,55,33,19).

Thus even a count-preserving macro with retained ports must update joint
data. It cannot reuse the old observer value. In the complete replacement
graphs at orders8 and9, the best possible minimum intermediate count is9;
this optimization is FINITE-EXACT, not an all-order barrier. The construction
is an explicit template at105, not a universal reset rule.

## The next incoming obstruction is173

The [fixed-core compiler](h5/PROOF.md) extends the previous H3 grammar.
For a named core K, a named cut O, and binary flank integers A,B, it gives

    H = qAB + eA + sB + h(g(A)+g(B)-2)
        + t + u(A mod2) + v(B mod2).

Here h=H(K), q is its ordered two-path-cover count allowing empty paths,
e,s retain the cut's cover endpoints, u,v its Hamilton endpoints, and
t=H(K+x). All coefficients come from the actual named core. The inherited
one-apex function g has an exact binary recurrence. Positivity supplies
a terminating inverse for every requested target and every fixed K.

Every H5 source is TT_a -> Q4 -> TT_b. The complete sixteen-cut compiler
therefore decides its full one-birth image, at every order. Combining
parents1,3,5 misses exactly7,21,173,335,663 among the odds through1001.
At173, no parent1,3 or5 exists;7 is forbidden, so every vertex deletion
has count at least9. The explicit rows

    (126,196,72,210,230,206,128,5)

have H173 and deletion counts(27,71,51,27,51,71,77,9), proving sharpness.
This is the next step after the sharp parent floor5 at55. It is not a
new forbidden value. At105 the parent3 lane works but the parent5 lane
fails at every order: more internal paths in a source do not imply a
larger attachment image.

The general formula also repairs uniform module collapse exactly:

    H(child) - H(K) H(collapsed child) = P_2(K) A(B+1)

for one uniform direction, with A,B exchanged in the other. P_2 counts
ordered covers by two nonempty paths. It is6 for C3 and22 for Q4.
This generalizes the prior MISTAKE-756 correction and retains module
reentry rather than assuming a contiguous path visit.

## A structural role for105 in the golden model

The [golden operator theorem G105](golden/NOTE.md) concerns the inherited
seven-letter one-reversed-bridge language and its exact three-coordinate
observer. A word can make **every matrix entry divisible by m** if and
only if every prime factor of m belongs to{3,5,7}. Thus105 is exactly the
squarefree support of uniform divisibility in this operation language.
The factor11 prevents such a reset modulo55.

Explicit mixed words W3,W5,W7 satisfy M(Wp)=p Up, with Up invertible
modulo every power of p. The7-word is

    W7 = q5 q0 q2 q5 q2 q5.

Its cost is24 vertices. Repeating it twice or three times gives raw zero
residues modulo49 or343. **The normalized phase survives bijectively**:
after dividing by7^k, the state is v U7^k. This is controlled divisibility,
not phase erasure or a count-neutral source reset. U5 is even integer
unimodular, so its word multiplies the gcd of all three coordinates by
exactly5. Small complete subspace walks find and audit these words;
no large Hamilton-count census is required.

The old golden/square-packing work supplies a precise comparison. Its
Fibonacci matrix has a finite golden quotient on which a fifth power is
the identity, while the integer lift retains nonzero lattice displacement.
The bridge contraction similarly hides surviving information in its
normalized lift. Its arithmetic norm is not a Euclidean packing metric.
The guarded independent-set operator and the square-sum graph supply
separate examples in which incidence and legality survive beyond a total.
Their detailed source/target maps and limitations are in the golden note;
no packing-to-tournament equivalence or physical packing optimum is claimed.

## What a complete constructive proof still needs

The [prime/7p reduction THM-4094](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
and the previous threshold95 reduction remain useful sufficient routes.
This session sharpens their procedural obligation: construct branch-ready
sources recursively, and prove an actual reset rule when readiness is
incompatible. Such a rule must carry named cuts, joint attachments and
normalized phase through each legal operation. The105 replacement and
the golden singular words solve different pieces of this problem; an
intertwining map or a common closed state space is still absent.

The81/27 control refers to the inherited pair's **automorphism orders**.
It is not the congruence ring modulo81 or27. The present results do not
collapse that pair, close the seventeen-coordinate observer, establish
unlimited ternary returns, or provide a new complete A000568 recurrence.

## Reproduction

    python3 -B 05-knowledge/results/chanpak_55_105_golden_20261110/reproduce.py

The default verifies the manifest, source pins, local links and saved
normal/optimized receipts. `--rerun` freshly executes the six checker
pairs. Proof reviews use independent algorithms: separate core
canonicalizations, finite-field vector sets versus row-reduced subspaces,
and literal Hamilton counts versus subset recurrences. The strong inventory
is an explicit inherited dependency. See [session audit and continuation](SESSION.md).
