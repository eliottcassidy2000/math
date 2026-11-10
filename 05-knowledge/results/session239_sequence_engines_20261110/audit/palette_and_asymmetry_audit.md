# Supplementary audit: reversal palettes and two asymmetric strong fibres

Status: PROVED identity audit; FINITE-EXACT independent implementation checks.

## Palette recurrence

The claimed recurrence is correct for a finite pair-color set of size q with
an involution having t fixed colors. Its necessary realizability conditions
are integers q>=t>=0 with q-t even; the implementation checks these. Vertices
have no loops. A color on an ordered pair reverses by the involution when
its endpoints swap.

Within one ell-cycle, unordered pairs of non-antipodal displacement form
floor((ell-1)/2) pair orbits. Returning to the initial unordered pair does
not reverse endpoints in these orbits, so each admits q assignments. If ell
is even, the antipodal pairs form exactly one additional orbit whose return
reverses endpoints; precisely t colors are allowed. Between different cycles
of lengths a,b there are gcd(a,b) pair orbits; endpoints cannot reverse
because they remain in different vertex cycles. Each has q assignments.

On inserting k cycles of length ell, the numbers of newly introduced ordinary
and reversing pair orbits are therefore

    du=k floor((ell-1)/2)+ell binom(k,2)+k R_ell,
    dv=k if ell even, otherwise0.

The factorial multiplier and future-response coordinate update are unchanged
from the previously audited tournament engine. Every transition is multiplied
by q^du t^dv. Equal degree/response states thus have identical continuation
weights; factorial-scaled states and final division remain integral. For
t=0, every even-cycle profile vanishes, so omitting even lengths is exact.
This also handles q=t=0: the empty/single-vertex structures are unique, and
all larger counts vanish. The code correctly uses the combinatorial empty
product convention0^0=1. No division by a color count is performed.

Graph(2,2) means two reversal-fixed states for absence/presence. Digraph(4,2)
means none/both fixed and two directed singleton states swapped. Oriented
graph(3,1) means absence fixed and the two orientations swapped. Tournament
(2,0) has only the two orientations. These interpretations are correctly
specified; digraphs here are loopless and permit reciprocal arcs.

## Independent checks

`palette_audit.py` imports only the target DP. It constructs full directed
color matrices for every labelled assignment through order4, tests fixedness
under every permutation directly by matrix entries, and averages the fixed
counts. It does not reuse the worker's direct partition formula or orbit
canonicalization. All nine realizable palettes with0<=q<=4 pass, including
(2,2),(4,2),(3,1), the empty palette, and two exchanged pairs. Horizon tests
at0,1,2,7 agree with the corresponding prefix computed at15. Six invalid
parameter tuples are rejected.

Normal and optimized runs pass with identical output:

    python3 05-knowledge/results/session239_sequence_engines_20261110/audit/palette_audit.py
    python3 -O 05-knowledge/results/session239_sequence_engines_20261110/audit/palette_audit.py

This audits the integral implementation and its finite evidence. Colored
Burnside counting and pair-reversal palettes are inherited mathematical
machinery, not an asserted new formula or literature-priority claim.

## Strong H203 and H217 are asymmetric

The implication is correct. Let T be strong with H(T) in{203,217}.
THM-4744,
`01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`,
uses the CITED Busch minimum

    f(3q)=3*5^(q-1), f(3q+1)=5^q, f(3q+2)=9*5^(q-1).

Its nondecreasing values give f(11)=225, so T has at most10 vertices. This is
an inherited all-order lower bound, not an inference solely from a finite
order11 census. The automorphism group embeds in S_n and acts freely on
ordered Hamilton paths. Thus its order divides both n! and H(T).

Since203=7*29 and217=7*31 and n<=10, the only possible nontrivial prime divisor
of the group order is7. Cauchy's theorem would give an order7 automorphism.
The newly audited strong sevenfold gap forces either n=7 with H175 or189, or
H>=2163, both impossible here. Hence the group has order1.

This is a restriction on strong realizers, not a nonexistence assertion.
The strongness assumption is essential to both the Busch bound and symmetry
gap. Targeted canon searches found no earlier203/217-asymmetry statement, but
no literature novelty or exhaustive priority conclusion is justified.
