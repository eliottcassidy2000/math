# Focused independent review: Session276

Status: PASS on the final statements reviewed below.
This reviews new session claims, not the repository's
historical theorem corpus. Chan–Pak remains OPEN.

## Selector lane: PASS

The unordered path-cover substitution identity is correct: each Hamiltonian
path uniquely decomposes the inserted module into maximal runs, and a j-path
cover has j! possible orders. The coefficient a1 is exactly H of the outer
tournament, even for the singleton outer context. The concrete C5 full cover
difference therefore proves delta=2H(B)+4a2. Positivity forces H(B)=1 if delta=2;
transitivity conversely prohibits separated visits, so the iff claim is sound.

`independent.py` imports no producer code. It enumerates every permutation of
each six/eight-vertex source. Invalid adjacent arcs force path breaks; valid
ones admit optional breaks. Dividing the resulting k-cover coefficient by k!
recovers the unordered cover profile. Separate split-permutation enumeration
recovers the entire start/end path-length response, confirming equality for
both source pairs. A separate subset DP counts the literal triangle
substitutions and confirms the C7 ordering reversal−14. Results are in
`independent.json`.

No assertion that an arbitrary fusion-central switch is substitution-central
survives; that is the explicit hostile control, not a failure of the theorem.
The general fusion and cover substitution mechanisms are inherited. The
concrete source/context comparisons are the scoped contribution.

## Corridor lane: PASS on current revised proof

The island theorem is correctly restricted to one/two apex-incidence changes,
neutral tails, and apex-preserving core automorphisms. It does not prohibit
the displayed internal-core escape or arbitrary reframing. The20-cut middle
Johnson layers are connected, and every departing core-only move costs at
least8. The revised distant-right-tail formula is
`H_core+2^d(21k+15)` for d>=1, distinct from the d0 expression `H_core+24k`.
`tails.py` constructs literal adjacency matrices independently of the
producer's quadratic cut evaluator. It checks all32 core cuts, both anomaly
sides and d0..3 (256 cases), then the entire128-cut nearest-two-tail cube;
minimum overall exit8 and minimum tail exit14 are both attained.

The ternary charge decoder correctly combines two facts: the first differing
coefficient must be even by its distinct2-adic valuation, and after that
its magnitude2 dominates every possible later contribution. Dominance alone
would not suffice for coefficient differences in{-2,...,2}. The resulting
fiber dimension/count and n>=2 endpoint-charge boundary are correct. The
width-at-most2 observation is explicitly finite, not promoted to a theorem.

## Wildcard lane: PASS on full proof and finite controls

Read the completed unique SCC grammar and checked all boundary cases,
including the two-vertex empty word, singleton apex cores A1, and absent C1.
The finite core bounds and product-bound total nonneutral length are sound.
The displayed ordinary generating function counts named word presentations,
not tournament isomorphism classes. Its coefficients are eventually
polynomial with degree at most2; no claim of positive numerator for all
allowed targets follows.

The productive minimal-DFA SCC claim also follows. A mixed-letter recurrent
walk repeated in an accepting context would require unbounded nonneutral
core transitions, contradicting the normal form. A unary cycle of length>1
cannot distinguish its states: for each fixed suffix, repeated insertion
makes H nondecreasing and equality/threshold acceptance an interval; periodic
distinction is impossible. Dead states are explicitly excluded from this
productive argument.

`languages.py` independently constructs literal graphs and uses its own
subset DP for all1,365 words of lengths0..5. Every output agrees with the
complete count-Moore21 artifact. All192 coefficients for odd targets1..63
and lengths0..5 agree with the finite-core generating-function artifact.
This is an independent control of the reported artifact, not an independent
re-enumeration of every strong core through the full Busch bound63. The
all-length conclusion rests on the normal-form proof and exact closure,
not those bounded controls.

## Digits lane: PASS through exact inverse bounds and lattice sieve

The vertex-transitive supplier proof retains the root-alone term h and
correctly obtains alternating coefficients by complement cancellation at
odd deletion order. The singleton deletion boundary is valid. The bilateral
sandwich formula is a direct exact multiplication of the physical carrier;
its central difference, unrestricted integer inputs, and parity terms are
correct. The nonnegative hyperbolic target bounds apply in the stated
interior-cardinality symmetric-core setting.

The gcd formula is exact as a lattice-generation statement only: values at
0,1,2,3 and the mixed11 difference isolate the displayed generators. The
C7 missing residues3,4mod7 and C9size4 modular sweep are sound and are not
promoted to absolute target coverage. C9size6 and size4 are now distinguished;
only the latter is used in the primitive example and final receipt.

`sandwich.py` independently constructs both rooted fusion and transitive
integer-input graphs without producer imports. Its104 literal checks include
all A0..63 one-apex counts, all C5 sandwich inputs A,B0..3, and the final C9
15381/15383 sources and three nontrivial sandwiches each. All pass.

## Selector extension: PASS

The all-length outer consumer C3[G,TT_m,1] for m>=1 has one-run coefficient
2^m+1 and two-run coefficient2^m-1. The nonempty subset before the unique
return port accounts for the latter. No third run is possible. The m0 ordinal
join is explicitly separate. `consumers.py` independently constructs the
literal graphs for both source pairs and m0..3 and confirms both exponential
response laws. These are consumer-relative laws, not proof of a universal
corridor.

## Provenance boundary

The balanced C5 counts41/43 are inherited from THM-4123,
`THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice.md`.
Generic central transport is inherited from Session229's rich-module pair
and the subsequent five-module bridge. General full-cover contextual
substitution, SCC multiplicativity, positive finite-semiring truncation and
standard automaton minimization are not new general theorems. The present
contributions are their precise physical suppliers, scoped consumers, inverse
formulas, all-order fixed-core barrier, and target-fibre normal form.

## Late constructive additions: PASS

The positive radix return map in digits section7 follows exactly by substituting
A=2a+r,B=2b+s into g(2a+r)=g(a)+(1-r)(a+a mod2). All seven updated
coefficients have the stated signs; u,v at leastdouble on every branch.
Initial-target budget therefore gives the stated finite depth bound, including
zero-budget acceptance and immediate over-budget rejection. The intermediate
state is explicitly formal, with physical meaning supplied by reconstructed
integer flanks, not misidentified as a tournament. A failed tree rejects only
the bank consumer, not Chan–Pak. No missing parity term or coefficient change
was found.

The corridor43→45 core-edge escape and its common frame are independently
checked in `escape.py`, together with all source-only lengths0..3 and positive
bifan lengths1..2 on each side. The all-order difference formulas follow
from the inherited packets and match these literal controls. Its+2step at
source length1 and failure of positive-bifan transport are correctly separated.

## Late all-order extensions: PASS on structural proofs

**One-sided escape transport.** Corridor section6 correctly computes the
43→45 gate's fusion difference. Its left multiplier equals2 exactly when
`U+d-H=1`, and the one-apex identity turns this into a positive selected-pair
sum. Hence at most one outgoing spine incidence is necessary and sufficient.
The right-context inequality forces H=V=e=1, equivalently a source root;
empty complements use H(empty)=1. The resulting complete shift language
`{0,48} union {57*2^t:t>=1}` is genuinely lacunary. Common frames extend as
claimed, but no path joining successive scales is asserted.
`transport.py` independently builds all63 unit-deletion contexts with spine
lengths0..5 and validates both iff classifications and every retained shift.

**Interior-insertion obstruction.** `insertion.py` is an independent
alternative to the producer's product-DFA proof. It enumerates the six H13
strong cores and25 H27 strong cores using literal adjacency matrices and DP
under the proved Busch cutoff. For subsequences of length3, every neutral
reservoir can be saturated at length3. The resulting35 maximal-clock
H27 templates avoid all six H13 cores (210 comparisons). By the normal-form
theorem this proves the obstruction at ALL lengths, not merely the tested
lengths. The same script independently checks the physical repair counts
13,25,13,15,27 and final reframe2010. The repair intentionally changes an old
arc and is not a pure insertion or a full0/2 count corridor.

**Fixed-apex hierarchy.** The compatible ordered apex partition correctly
retains every SCC's named apex block. Single-apex empty cores preserve all
reservoir positions. Multiplication of the strong-core bounds gives
`j+k-s<=B(N)`, so after passing to denominator `(1-z)^(k+1)` the numerator
degree is at mostB(N), independent ofk. Its possibly negative coefficients
are not identified with a physical count bank. The finite-difference and
count-one boundaries are sound. This proves a two-budget structural hierarchy,
not a uniform nonzero target numerator.

## Positive carry overlaps: PASS

**Compensated C9 diagonal.** The exact carry increment D_k(d), the product
change q(B-A-1), parity constant K, and coefficients a_k,C in the carry-stratum
Diophantine equation are correct. For q>h both coefficients are negative, so
an integer solution can be shifted into the nonnegative quadrant by the
explicit positive homogeneous direction. This proves sufficiency, not just
necessity of a congruence. The displayed C9 parameters retain their carry
lengths and parities and give delta4 for every j>=0; the central alternatives
then supply four consecutive odd counts. No increasing interval width,
arbitrary starting point or connected return path is implied.
`carries.py` uses a separately written iterative binary evaluator and checks
6,855 exact arithmetic identities, including parameter j=10^30. It does not
claim direct Hamiltonian-path enumeration at those huge vertex orders.

**C7 even-to-odd overlap.** Equation(9) and its gcd criterion follow because
`a+(a mod2)` runs through precisely the even nonnegative integers. The positive
solution direction makes both integer operands arbitrarily large; the explicit
C7 progression and its alternate even A are correct. The next-unit-jump lower
bound retains the parity sum and the rectangle mixed cost is exactlyQ.
`triples.py` independently counts the two eight-vertex cores, reconstructs all
four nineteen-vertex literal input graphs from independent fusion code, and
checks their reported literal counts against the independently evaluated
compiler. Huge-parameter controls include t=10^100. The producer's separate
nineteen-vertex DP is not redundantly claimed as a second full DP here.
The code0 C7 core orientation is now explicitly steps4,5,6, avoiding ambiguity
with the converse orientation used earlier in the note.

## Complete readiness replay: PASS

`readiness.py` rebuilds the inherited capped positive transfer and performs
explicit DFS reachability independently of the producer's reverse-topological
bitset computation. It recovers28,016 live states and1,080 count-Moore states,
checks every branch-mask row through source63, and confirms common-parent
readiness exactly at5,9 and no depth-two readiness through31. This is an
independent graph algorithm over the inherited transfer, not an independently
derived transfer. The separate literal-word tests audit that carrier locally.
The cap127 completeness argument is correct for both stated scalar ranges.

The full unit-context nonclosure theorem in corridor section6 also follows
exactly from DeltaV=2A-16qB with0<=A<q and B>=0. Zero requires A=B=0,
which is neutral padding. Association reduces any same-mark finite program
to such a sandwich. Context-dependent programs, different marks, and h>1
contexts are explicitly not ruled out.

## Final scoped repairs and boundaries: PASS

The carry lane's complete eight-neighbor parameter-overlap graph is a
matching. Same-sign diagonal increments exceed4 by the stated positive
product bound. The remaining diagonal's residue is
`8-4(A mod2)-5(B mod2)` modulo9, admitting among{-4,-2,0,2,4} only+4
from odd/even input. Its opposite-parity endpoint has only the reverse edge.
This is a classification of that parameter graph, not of the tournament
metagraph or arbitrary nonlocal parameter moves.

`repair.py` independently reconstructs the finite H159 receiver repair,
checks its entire(H,h,Q,U,V,d,e) packet by separate endpoint/cover enumeration,
and checks both supplied common frames. The direct157→R159 link changes
exactly three old arcs, whereas159→R159 changes two. The proof explicitly
enlarges the move language rather than silently using a forbidden move.
No all-order rearming law is claimed from this finite witness.

A minor ambiguous hostile-control sentence in the wildcard lane was repaired:
for01120 the position-sensitive analogue uses old arc(2,3) and new arc(0,3),
giving71 rather than desired51. Keeping the earlier numeric indices(1,2)
and(0,2) instead gives43. The final proof and receipt specify the intended
indices; the intended71 was independently confirmed. The obstruction and
positive H13fork are unchanged.

All reviewed positive results keep their absolute costs and actual consumer
or move language. No scoped negative result is promoted to a global omission,
and no unbounded family of bounded-width intervals is promoted to an eventual
threshold. Chan–Pak remains OPEN.
