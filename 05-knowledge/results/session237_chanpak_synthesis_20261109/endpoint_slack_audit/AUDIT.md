# Independent audit: endpoint equality and pair energy

Status: **PASS**. The positive nonsink equality classification is proved at
all orders. The pair-energy identities and endpoint-minimum consequence are
correct. No Chan–Pak total-count completeness follows.

## Proof audit

1. **Strong-component reduction.** Hamilton paths traverse all SCCs in their
   total order and cannot return to an earlier one. Positive endpoint count
   puts a in the final SCC C. Both h and e factor by the product M of earlier
   H-counts. Equality gives M(H(C)-2e(C))=1 with positive integer M, so M=1.
   Each earlier SCC therefore has one Hamilton path and is a singleton, since
   a nontransitive tournament contains a cyclic triangle and insertion injects
   its three paths into distinct spanning paths. Their union is a transitive
   prefix. A singleton final SCC would make a a global sink, excluded here.

2. **Equality forces a unique exit.** The insertion map moving terminal a
   immediately before its first outneighbor is well defined and injective:
   its preceding vertex (if present) beats a, and deleting a then appending
   it recovers the original path. Its image ends in N-(a), disjoint from
   paths ending at a or N+(a). Thus 2e_a+e_O<=h. In a strong tournament,
   a Hamilton cycle cut at each possible endpoint proves every endpoint count
   positive. Hence h=2e_a+1 implies |O|=e_O=1, say O={b}, e_b=1.

3. **Transitive deletion.** In T-b, a is a sink. Every Hamilton path there
   ends a and extends with b. Conversely deleting terminal b from any path
   ending b leaves a path spanning T-b (necessarily ending a). Thus H(T-b)=1.
   This is exactly the transitive case. Its source x1 must be beaten by b,
   or x1 would be a source in T, contradicting strong connectivity.

4. **Partition count, including the adjacent correction.** In the apex family
   every path not ending a contains the edge a->b. Exactly one ends b, and
   the others are P,a,b,Q with nonempty Q. For min(Q)=i in S, all earlier
   positions are in P, and all later positions are independently assigned:
   m=sum 2^(l-i). Paths ending a are P,b,Q,a, excluding exactly those with
   nonempty P and max(P) in S. If min(Q)=i<j=max(P), the j-i-1 intermediate
   positions are free, giving 2^(j-i-1). If max(P)=i<j=min(Q), the gap must
   be empty, so j=i+1 and there is one additional partition. Therefore pair
   contributions are 2 at distance1 and 2^(d-1) at distance d>=2, exactly Q.
   So e=m-Q, h=(m-Q)+m+1, and h-2e=1+Q. This includes l=0 and S empty.
   No pair terms are negative, and every pair contributes positively;
   Q=0 iff |S|<=1. In the strong equality reduction 1 belongs to S, so
   S={1}. This completes the necessity argument and verifies the converse.

5. **Exceptional roots and consequences.** Without e>0, every nonsink of a
   transitive tournament is an equality root with (h,e)=(1,0). Excluding it
   is essential, not cosmetic. For n=3 the R_n representation can place any
   chosen triangle vertex in the penultimate position; for larger n that
   root is uniquely selected in the final strong component. The family gives
   (h,e)=(2^(n-2)+1,2^(n-3)). Endpoint universality in
   `01-canon/theorems/THM-5001-every-positive-endpoint-quotient-has-a-small-strong-realizer.md`
   explicitly supplies strong realizers, so f(N) is finite. Equality occurs
   exactly at powers of2, and oddness advances the off-power lower bound
   from 2N+1 to 2N+3. No unproved minimization formula is used.

## Provenance audit

The injection and the stronger bound 2e_a+e_O<=h are inherited from THM-4876
statement7 and its audit D, not new. The reversed-extreme-arc family and its
endpoint counts are also inherited. The audit D census through order9
supported but did not prove the classification in all orders.

`05-knowledge/results/session231_overnight_completion_20261109/descent/ROUND2.md`
Section4B treats endpoint equality b=0 together with cover defect delta=1,
as part of improved singleton-join equality delta+b=1. Its Section3 treats
delta=b=0 at h=3. It does not show that arbitrary endpoint equality forces
those defect values. The present argument removes the extra cover-defect
hypothesis and classifies endpoint equality directly, so its claimed advance
over that source is accurately delimited. The coefficient function Q is the
inherited Sessions197/198/205/206 and THM-4888/5085 pair-distance energy.
Identifying it with this physical endpoint slack is an exact synthesis, not
a new tail-transfer formula or external novelty claim.

## Independent finite controls

Run from repository root:

    python3 05-knowledge/results/session237_chanpak_synthesis_20261109/endpoint_slack_audit/check.py

This checker does not import the producer's DP or SCC shape recognizer. It
counts Hamilton paths by literal vertex permutations. Equality recognition
uses all relabellings of explicitly constructed transitive-prefix/R_n rooted
templates, rather than the proof's structural recognition algorithm.

Universe: every labelled tournament through order5, all 3,446 positive
nonsink rooted cases. Equality counts at orders1..5 are 0,0,6,48,360.
It separately records 566 nonsink zero-endpoint equality controls. All127
apex masks with l<=6 are checked by literal Hamilton paths against the three
identities, covering empty and singleton masks, adjacent pairs, separated
pairs, and higher interactions. These bounded tests supplement the proof;
they are not the source of its unbounded quantifiers.

The producer's larger DP census and apex masks through l=9 are separate
controls, not relabelled here as independently rerun coverage.

## Endpoint-observer automata: independent synthesis audit

**PASS.** Write E(m)=m-Q(m), read m in MSB-first binary, and output E(m)
modulo M. The exact minimal state count is 1 for M=1, 2M^2 for odd M>=3,
and M^2/2 for even M.

For odd M use state (r,t)=(m mod2M,Q(m) modM). The transitions are
(r,t)->(2r,t) on0 and (2r+1,t+r+(r mod2)) on1, with residues in their
respective moduli. All 2M^2 states are reachable by the independently proved
joint pumping of (m mod2M,q modM), because multiplication by2 is invertible
modulo odd M. Equal current and after-zero outputs force r=s modulo M
and t=u. Distinct states then differ by M in r, so have opposite parity.
After a one the output is the old output plus1 minus the old r parity;
this distinguishes the pair for M>=3.

For even M, the state reduces to (r,t)=(m modM,Q(m) modM), where t is an
even residue. Parity of r is well-defined modulo M, so the same transitions
are well-defined. Joint pumping gives every such pair, for M^2/2 reachable
states. Equal current and after-zero outputs force r=s modulo M and t=u;
thus the empty suffix or0 distinguishes every distinct pair. For M=2 the
two states already have distinct current outputs. M=1 is separately constant.

**Observer normalization matters.** The arithmetic packet's 2M^2 theorem at
all M>=2 concerns F(m)=(H(m)-1)/2 modM, not H(m) modM. For the latter,
H=1+2F, so its minimal count is that for F modulo K=M/gcd(M,2): it is1
when K=1 and 2K^2 otherwise. Therefore endpoint observation saves states
relative to normalized total observation at even M, but not relative to
actual H modM for even M>=4. At M=2 the endpoint needs two states while
actual H is constant. A statement merely comparing to "full H observer"
without the normalization would be false.

`endpoint_dfa.py` independently builds the reachable joint automaton and
runs Moore partition refinement for all M=1..40 on the three distinct output
functions E, F, H. Every count matches the formulas. This finite minimization
supports, but does not replace, the all-order reachability/distinguishability
arguments above. Normal and optimized runs are byte-identical; deliberate
--fail controls exit nonzero in both modes.

## Optimized-Python quality gate

The literal checker uses explicit RuntimeError checks, not Python assertions.
Normal and python -O runs are byte-identical (`check.out`,
`check_optimized.out`), and --fail exits1 in both modes (`fail.out`,
`fail_optimized.out`). Thus the recorded checks remain active under -O.

Final source cross-check: `../arithmetic_states/ENDPOINT_OBSERVER.md` explicitly
compares E modM with the **normalized full count** (H-1)/2 modM and therefore
matches the observer-normalization warning above. Its reachability,
distinguishability, M=1/M=2 boundary treatment, and fourfold even-modulus
saving are correct. The physical input-digit order is consistently MSB-first
for the arithmetic observer and reversed position assignment for the apex
family, with the retained word providing the construction sidecar.
