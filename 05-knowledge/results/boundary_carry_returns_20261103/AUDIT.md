# Audit and exact scope

**Same-agent proof review plus independent implementations.** No independent
agent was commissioned and none is credited. The user did not authorize new
agents for this continuation.

## Proof checks

* BC1: run splitting is a bijection; each inter-run gap yields a disjoint
  crossing arc. The matching lower bound uses insertion into a directed
  word, where repeated v occurrences cause no failure. Positivity prevents
  cancellation of the largest exponential. At nu=0 the starting word is v.
* BC2: births are specified as literal cuts. The two words are both guarded;
  their spine difference changes exactly one old exterior pair. The extra
  sink prevents an endpoint clone from losing its common frame. All helpers
  are deleted by name; no old vertex is deleted. Different arc changes may
  reuse temporary names only after deletion with a freshly specified cut.
* BC3: exact coefficients, not numerical interpolation without a degree
  bound, establish the all-m identities. All permanent target vertices,
  both old pins and every clone address are present. Existing source pins
  do not imply unchanged future readiness; that stronger claim is excluded.
* BC5: root-side orbit banks plus all cross masks cover the full declared
  universe. Converse preserves the count function, not a free physical
  transport. Gcd-one removes constant SCC factors before using the strong
  census. No complete order-ten census is asserted.

## Independent paths

Producer `scout.py` uses subset endpoint DP and a positive cut formula for
all one-vertex births. It exhausts the stated two-complement/common-frame/
one-birth language, checking the ordered-run bound and one extra clone size.
It finds 14 first programs and 37 next programs after justified marked
deduplication. The two omitted isomorphic first targets have explicit maps
in `audit.json` fixing pins 1,5 and clone root 3. Those maps are used only
to avoid duplicate census work; they are not physical steps.

`audit.py` imports no producer module. For every permutation of the exterior
vertices it classifies gaps as forced root insertions, optional insertions,
or illegal, giving the full skeleton coefficients by binomial expansion.
It independently checks 420 rooted representatives at orders 2--6 against
the matching theorem, and independently covers all labelled tournaments
through six via relabelling orbits of the small banks. It checks every saved
positive program, including all 37 second-level short programs, not just
the chosen 97->293 witness. Boolean-matrix native moves replay the first
programs and three repaired continuations at m=1,2,4,7. This verifies
implementation; the proofs supply unbounded m.

For BC6, the producer evaluates the exact two-cover cut identity. The audit
enumerates every directed prefix literally to obtain its endpoint counts,
independently checking the four determining clone sizes for every exterior
pin of each new target. The fifth sample checks the degree bound; the proof
of that bound is necessary for the asserted constant defects 0,0,17.

BC7's exact moment-budget lower bound uses actual matching-star tournaments,
not arbitrary formal vectors. The producer uses subset convolution to count
ordered path covers. The audit checks the first three clone sizes of every
triangular witness by full vertex permutations and cuts, with no cover-DP
import. This independently validates the displayed determinants 1,6,144;
the general rank statement rests on BC7's positive top-coefficient proof.

BC8 uses the incoming general-child naturality theorem on this packet's
specific two-level programs. `child_audit.py` checks six literal steps for
each of C3 and the strong four-vertex H=5 child, including all helper paths.
Subset-DP counts of the substituted tournaments agree with the skeleton
formula whose child cover counts are independently enumerated by
permutations and cuts. The observed excesses are exactly c(H(child)-1).

`catalog.cpp` enumerates root-neighbourhood banks, counts H by subset DP,
and tests clone values through the proved maximum degree. Independently,
`audit_catalog.cpp` counts literal paths by recursive DFS, with cutoff only
after exceeding every target, then computes skeletons by exterior-word
permutations. The two algorithms agree on **every scalar and uniform hit
count** at orders 7--9. Their complete files are `catalog*.stats` and
`audit*.stats`. Order-ten cover-filtered outputs are existence witnesses;
each is independently audited by the skeleton algorithm before use.

`returns.py` selects target names at design time so the root star extends
the supplied source, and saves every birth cut and path word. Its chosen
repairs use 2,2,3 exterior flips. `audit.py` separately verifies their entire
old named relation and the target polynomials. The three programs are not
inferred merely from equal H values or isomorphic unmarked targets.

`arithmetic.py` computes multiplication in two explicit finite quotient
rings. Exhaustive power cycles prove the stated field/generator claims and
first-scalar return times. The rejected z+u generator has period 11 in the
projective quotient, whereas z^2+z+u has period 33; retaining this negative
control prevents a group-order/generator-order conflation.

## Corrections and stopping conditions

No inherited theorem is retracted. A fresh draft guess that z+u generated
the full 33-state projective clock failed immediately; its replacement and
mechanism are recorded in `CONNECTIONS.md` and the mistakes ledger.

The endpoint family exists at all m, but there is no claimed unbounded-depth
return induction. The next supplied sources are exactly the four BC3 target
relations, not arbitrary H=287/289/293/295 representatives. The field clocks
are not claimed as tournament automorphisms or native moves. No result here
proves Collatz, LRC, or the unrestricted {7,21} spectrum.

## Workspace preservation

Startup used an isolated worktree at cb6bbd43ec. The only apparent unrelated
change, `lrc14_removal_measure_kps_S38.out`, was a Git newline-filter mismatch:
the working bytes were verified identical to its HEAD blob (712 bytes).
They were left unchanged, with a local index assume-unchanged flag so it
could not be accidentally staged or block clean synchronization. The shared
checkout was not modified. Only this packet and explicit navigation/error
entries belong to this session.
