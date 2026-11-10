# Independent audit: fixed-port arithmetic content and symmetry

Status: **PASS** for the all-order arguments and explicit witnesses. The finite
minimality statements are independently verified by literal enumeration through
order 5. The source's order-6 DP rerun is separately recorded when complete.
No Chan–Pak completeness or density conclusion follows.

## Mathematical audit

The inherited state and chronological S/B conventions agree with
`01-canon/theorems/THM-5085-fixed-port-tails-in-every-base-and-the-ternary-bank-to-ten-million.md`:
its x=2 dictionary is e=-Delta and d=h_u. This is an inherited transfer, not
new here. The oddness of h is the usual all-tournament Hamilton-path parity.

**Infinite gcd certificate.** Every descendant coordinate is an integer linear
combination of the four original coordinates. Hence c divides every count.
Conversely the four words empty, B, SB, BB give h, h+R-e, h+2R,
h+3R-e-2d. Their gcd g is odd, because it divides odd h. The third value minus
the first gives 2R, so g divides R; the second then gives e; the fourth gives
2d and hence d. Thus g divides c. This proves the infinite equality without
truncating the word language. Content divisibility, both new-prime criteria,
and the finite-bank prime pruning follow exactly as stated. In the prime
claim a composite content cannot divide a prime, while a prime content can
only supply that same prime; no sufficiency assertion is made.

**Automorphisms.** The A=Aut(C)_u actions on ordered spanning paths, paths ending
at u, paths spanning C-u, and ordered pairs (root-started P, complementary
spanning Q) are free. Fixing the pair fixes every vertex, including for empty
Q. This proves |A| divides all four coordinates. Restriction to C-u is faithful,
so the inherited actual-subgroup gap theorem applies to h/|A| and d/|A|.
For the singleton core d=|A|=1 handles the empty deletion separately. The
underlying all-order gap dependency is the expressly cited
`01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`,
including its inherited THM-4812 certificate dependency; this audit does not
freshly replay those larger certificates. Formal division by c is not licensed
by the group theorem, and the source correctly says so.

**Physical witnesses and information loss.** T=[24,1,19,22,2] is rigid: its
outdegree partition has one singleton and two pairs, each pair fixed because
swapping the endpoints of its internal arc is impossible. Its displayed
spanning directed cycle proves strong connectivity. Literal enumeration gives
z(T,1)=(9,9,3,3), z(T,0)=(9,11,1,1), and z(T^op,1)=(9,21,1,3).
The re-rooted and converse one-B counts are 19 and 29 respectively.
U=C3=>C3 with root in the second component has the same (9,9,3,3) state.
Rotating the other C3 is a rooted order-3 automorphism at every root; U has
no content-one re-root. At the corresponding fixed root, its converse has
state (9,36,0,3), still content 3. Thus the common state also fails to determine
whether converse releases this particular trap. Literal full groups have
orders 1 and 9. Equality after every common word follows by deterministic
transfer induction, not by finite testing. No same-order hypothesis is needed
for the nonexistence of a function of these four coordinates.

## Reproducible finite checks

Run from repository root:

    python3 05-knowledge/results/session237_chanpak_synthesis_20261109/symmetry_grafts_audit/check.py

The checker is independent of both source scripts: it enumerates every
ordered simple directed path directly, then counts complementary ordered
covers from the resulting subset counts. It exhausts every labelled tournament
of orders 1 through 5, all 5,405 rooted states, with no inherited filtering.
It checks actual rooted groups, all coordinate divisibilities, the h/d subgroup
gaps, the four-test gcd, content monotonicity, and onset for primes 3,5,7,11.
First occurrences independently agree: B increase at order4 (mask4/root1),
S increase at order5 (mask10/root4), and rooted rigidity with nontrivial content
at order5 (mask9/root1). Thus all claimed strict lower-order exclusions are
covered independently, without needing the order-6 computation.

The original literal witness script was also rerun successfully. Its seven
word comparisons through length2 are finite controls only. The source DP
loop's declared state count is
sum(n*2^(n*(n-1)/2), n=1..6)=202,013. Its literal-tail comparison count is
(1+4+24+256)*(1+2+4+8)=4,275. These are correctly described as rooted labelled
universes, not isomorphism classes. The source DP has a harmless unused `rev`
variable/comment; the actual second DP correctly initializes only root u and
counts paths starting at u in the original orientation.

Scope: the audit verifies the new algebra and finite witnesses independently;
it preserves inherited theorem provenance and does not claim literature novelty,
all-odd coverage, or positive-density misses. No substantive correction required.

## Completed reruns and quality gate

The original full DP probe finished successfully: `probe_rerun.out` records
all202,013 rooted labelled states and4,275 literal tail comparisons, matching
the declared universe and all first occurrences. This was a rerun of the
producer, not an independent order-6 implementation. The separate literal
order-5 checker was converted to explicit RuntimeError checks and rerun in
normal and optimized Python. Outputs `check.out` and `check_optimized.out`
are byte-identical; --fail exits1 in both modes, as recorded in `fail.out`
and `fail_optimized.out`. No assertion is disabled by python -O in that checker.

## Equal-order strengthening

**PASS.** The added witness X=[62,48,2,38,44,4] is a universal source
(vertex0) followed by T on vertices1..5. Both X and U have order6.
Literal enumeration independently gives z(X,2)=z(U,3)=(9,9,3,3),
z(X,1)=(9,11,1,1), and |Aut(X)|=1. Thus X has a primitive re-root
while U does not. The contrast survives retaining order together with the
four state coordinates, strengthening the original unequal-order example.

Converse states at the original roots are z(X^op,2)=(9,42,0,3) and
z(U^op,3)=(9,36,0,3), so even order plus the original four coordinates
fails to determine the converse state. This pair demonstrates differing
converse coordinates, not primitive converse escape: both displayed
converse states still have content3. The earlier T/U example separately
supplies the stronger difference in primitive converse behavior.

The source-padding explanation agrees with the literal counts: paths
starting inside T cannot reach the preceding source, and that source adds
no multiplicative path choices to either complementary or full paths.
After converse, it becomes a sink, doubling R and setting e=0 as in S.
The updated independent checker includes these X controls; normal/-O
outputs again match and both --fail controls still exit1.
