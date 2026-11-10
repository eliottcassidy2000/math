# Integral engines for the reversal-palette family

**PROVED:** factorial gauge extends to finite reversal palettes.
**FINITE-EXACT:** the attached exhaustive and partition checks.
**INHERITED:** the palette cycle-index mathematics from THM-4845,
`01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md`.
No new classification, literature priority, or speed record is claimed.

## Objects and conventions

Let A be a finite set of q colors with an involution rho fixing exactly t
colors. Thus q>=t>=0 and q-t is even. For each unordered vertex pair {i,j},
record a color with respect to the ordering i<j. Relabeling the vertices
applies rho whenever it reverses that ordering. Count these colorings up to
all vertex permutations. There are no loops. Graphs have (q,t)=(2,2),
loopless directed graphs have (4,2), oriented graphs have (3,1), and
tournaments have (2,0). For directed graphs the fixed states are neither arc
and both arcs; the two single-arc states are exchanged. For oriented graphs
only the no-arc state is fixed.

The zero-color case is valid: there is one empty assignment for n=0,1 and
none for n>=2. We consistently use 0^0=1, as does Python integer arithmetic.

## Pair-orbit formula and recurrence

For a permutation, an orbit of unordered pairs either returns each endpoint
to itself after one full orbit, or exchanges the endpoints. In the former
case it allows q choices; in the latter case it allows exactly t choices.
Between different vertex cycles of lengths a,b there are gcd(a,b) pair
orbits, all nonreversing. Within one ell-cycle there are floor((ell-1)/2)
nonreversing pair orbits, and one reversing orbit precisely when ell is even
(the antipodal pairs). These exhaust all pair orbits.

Use the same finite-horizon state (d,R) as in PROOF.md, now processing *all*
cycle lengths in descending order. Its integer weight is d! times the sum of
q^u t^v/z over the histories in the state. Insert k cycles of length ell.
Then D=d+ell*k and

    Delta_u = k floor((ell-1)/2) + ell*k*(k-1)/2 + k R_ell,
    Delta_v = k if ell is even, otherwise 0,
    C = (d+ell*k)! / (d! ell^k k!),
    W(D,R') += W(d,R) C q^Delta_u t^Delta_v,
    R'_j = R_j + k gcd(ell,j).

The factor C is the same integer conjugacy-class insertion factor proved in
PROOF.md. All continuations depend only on the retained future responses,
not on the full cycle history. Degree pruning and merging therefore remain
exact. Terminal division by n! gives the orbit count by Burnside. All weights
are nonnegative integers. When t=0, every even-cycle term vanishes, so omitting
even lengths is valid; for (q,t)=(2,0) the recurrence is exactly the tournament
recurrence. For t>0 all lengths must be retained: deleting the even lengths
would incorrectly count graphs as tournaments, already failing at n=2.

This is the inherited palette formula with a factorial-gauge implementation,
not an additional generating-function identity. It supplies efficient reusable
functions for A000088, A000273, A001174 and A000568 using one state engine.
The state count is reported, but no new asymptotic or performance claim for
the palette family is made.

## Independent checks and scope

Reproduction from repository root:

    python3 05-knowledge/results/session239_sequence_engines_20261110/burnside/verify_palette.py

Seven palettes are checked through n=12 against a separate rational partition
sum over all integer partitions, without the response quotient or integer
gauge. For every n=0..4, the verifier also enumerates **every** labelled
pair-color assignment and takes the minimum of its full S_n orbit as its
canonical representative. This is an independent group-action check, not a
second evaluation of the pair-orbit formula. The palette involution fixes
colors 0..t-1 and swaps the remaining adjacent pairs. At q=4,n=4 the universe
has 4^6=4096 assignments, each with 24 relabelings.

Hostile/boundary controls include the empty palette, a one-color palette,
two exchanged color pairs, and rejection of impossible involution parameters.
The tournament specialization agrees coefficientwise through n=60 with the
completed tournament engine. Historical sequence prefixes are supplementary
controls only, explicitly not newly checked against an OEIS web source.
The exhaustive check independently verifies initial terms of the three named
related sequences without relying on a copied table.

Saved outputs: `palette.json` and `palette.out`. Source files:
`palette_engines.py` and `verify_palette.py`. The completed tournament engine
and all of its prior benchmark files were left unchanged.
