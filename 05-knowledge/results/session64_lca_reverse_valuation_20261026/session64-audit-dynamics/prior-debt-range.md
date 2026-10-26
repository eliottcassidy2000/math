# Guarded affine debt range summaries

Session56. PROVED: elementary finite-word monoid and range-query consequences. FINITE-EXACT: portable rational controls. OPEN: all-source Collatz coverage. No novelty claim, new source census, theorem identifier, or alteration of old canon.

## 1. Source pin and precise inheritance

SOURCES.json pins five raw-byte snapshots. The available math checkout reports HEAD 1b91742ce55a202b73c80d56561ba56d00b503a6; file hashes, rather than that commit alone, identify the bytes read. The old cumulative budget note (2026-10-05, sections 2–3) ALREADY proves the multiplicative summary (M,m) with concatenation (MN,min(m,Mu)), retaining checked source words and seed bank. Taking logarithms gives the additive monoid below; that mechanism is inherited, not rediscovered. Rational multiplicative comparison avoids numerical logarithmic errors in that old application.

THM-4735 reproduces the first-exit displacement ledger and its cone. Its direct charge_palettes_boundary_20261020g/charge.md dependency is absent in the accessible filesystem. We do not claim to have inspected it or rerun its census. The available theorem snapshot is sufficient to identify the inherited formulas, and the proofs used here are supplied independently. The 2026-09-26 clock note distinguishes immediate splitting from P5's EXTRA odd step before resetting: a split is not an arithmetic step. Session53's actual words distinguish odd U, shortcut T, and ordinary C clocks and require the extra final bit for an exact U word. Session33's clock-capacity theorem supplies the shared-port hostile below.

New connection: order-sensitive scalar debt summaries support disjoint range composition, segment trees, and binary lifting, but not the usual overlapping two-block sparse-table trick. Neither the algebra nor the data structure supplies legal native blocks or all-scale coverage.

## 2. Additive prefix minimum

For a finite list of real updates d_1,...,d_k, let Δ=sum d_i and m=min(0,d_1,d_1+d_2,...,Δ). Starting from B, every retained boundary including the input is nonnegative iff B+m>=0. The summary is (Δ,m), identity (0,0), and chronological composition is

(Δ_A,m_A) * (Δ_B,m_B) = (Δ_A+Δ_B,min(m_A,Δ_A+m_B)).

Proof: prefixes of AB either end in A or consist of all A followed by a prefix of B. Three factors yield total Δ_A+Δ_B+Δ_C and minimum min(m_A,Δ_A+m_B,Δ_A+Δ_B+m_C), independent of brackets. Therefore the operation is associative. On realizable summaries m<=min(0,Δ); that condition also ensures the displayed identity. The exact least real input is -m. For integer B it is ceil(-m), but integrality of intermediate states is a separate issue when updates are nonintegral.

Cheapest endpoint-only hostile: updates (-1,+1) return input 0 to 0 but fail at the first boundary. Reverse order (+1,-1) is feasible at 0. Equal totals do not determine feasibility or order.

## 3. Positive-affine prefix requirement

Fix a word of maps f_i(x)=λ_i x−δ_i with λ_i>0 and real δ_i. Every prefix j, including j=0, has f_[1:j](x)=Λ_j x−D_j, where Λ_0=1,D_0=0 and

Λ_j=λ_jΛ_(j−1),  D_j=λ_jD_(j−1)+δ_j.

Set Λ=Λ_k,D=D_k,R=max_(0<=j<=k) D_j/Λ_j. Exact real feasibility at every retained boundary is x>=R. This follows by dividing each inequality Λ_j x−D_j>=0 by its POSITIVE Λ_j. In particular R>=0 includes the initial-state requirement.

For A followed by B, define

(Λ_A,D_A,R_A) * (Λ_B,D_B,R_B)
 = (Λ_BΛ_A, Λ_BD_A+D_B, max(R_A,(D_A+R_B)/Λ_A)).

A must be feasible at x and its endpoint must be a feasible input to B. These two requirements prove the threshold formula exactly, not just sufficiently. The endpoint formula is ordinary affine substitution.

Identity: (1,0,0). The realizable summary space satisfies Λ>0 and R>=max(0,D/Λ); this condition ensures both identity laws. Closure follows directly: the new R is nonnegative and is at least D_new/Λ_new because R_B>=D_B/Λ_B. On this whole space the operation is associative: either parenthesization has endpoint parameters Λ_CΛ_BΛ_A and Λ_CΛ_BD_A+Λ_CD_B+D_C, and requirement

max(R_A, (D_A+R_B)/Λ_A,
    (Λ_BD_A+D_B+R_C)/(Λ_BΛ_A)).

The proof does not require δ_i>=0. With λ_i=1 and δ_i=−d_i, Λ=1,D=−Δ,R=−m, giving exactly section 2.

At integer inputs, feasible real-state trajectories are exactly x in Z with x>=ceil(R). A rational R is rounded by exact numerator/denominator arithmetic. This does NOT certify integer intermediates. Example f(x)=x/2 has R=0 and accepts integer x=1 over the reals, but its endpoint is not an integer. Integers on a specified congruence class x=r mod M require the additional guard. Their least member above R is r+M*ceil((R−r)/M), subject also to any explicit strict inequalities.

Strict positivity is different: for finitely many boundaries x>max D_j/Λ_j, hence integer threshold floor(R)+1. If only some boundaries are strict, retain which threshold is strict when maxima tie. No epsilon rounding is legitimate.

## 4. Ordered range queries, clocks, and failure of overlap

Append a nonnegative integer step count k; concatenate by k_A+k_B. If several clocks matter, append their entire vector and add it componentwise. The scalar summary alone forgets duration: one identity map and two identity maps share (1,0,0) but have different clocks. A coarsened itinerary must retain internal prefix minima or requirements; keeping only each block endpoint is invalid.

The operation is generally noncommutative and non-idempotent. One map x->x−1 has summary (1,1,1), whereas its duplicate has (1,2,2). Thus an overlapping two-block query can falsely count the overlap twice. For example, a three-letter word of x->x−1 has requirement 3, whereas combining its length-two prefix and suffix yields requirement 4. The input 3 is the cheapest boundary witness for that chosen overlap layout.

A segment tree combines disjoint covered intervals in left-to-right order in O(log n) monoid operations per query; each node stores O(1) rational coordinates and the clock/guard sidecar. A disjoint sparse table answers a nonempty interval by at most two disjoint summaries in the correct order after O(n log n) preprocessing; singleton queries are handled directly. Binary lifting composes consecutive disjoint dyadic blocks. For a deterministic successor map, the second block must start at the endpoint of the first, not at the original state. Query costs count algebra operations only: numerator and denominator bit lengths can grow linearly with itinerary length, so no constant bit-cost or universal speedup is claimed.

One may cache these summaries for fixed authenticated itineraries. The data structure cannot decide that an arbitrary state has the same itinerary, cannot discard source guards, and cannot manufacture available blocks at every scale.

## 5. Native Collatz interface: what is and is not the scalar state

For U(n)=(3n+1)/2^a with actual a=v2(3n+1), set z=n−1. Then

z'=(3/2^a)z + 4/2^a−1 = λz+p−d,
p=1 if a=1 else 0,
d=1−4/2^a if a>=3 else 0.

Thus a word has z_end=Λz+P−D, with nonnegative P,D and updates P_AB=Λ_BP_A+P_B and D_AB=Λ_BD_A+D_B. The inherited cone refers to D, NOT to the signed net intercept D−P. Applying section 3 with δ_i=d_i−p_i gives the exact displacement trajectory. Applying it with δ_i=d_i instead gives a lower envelope which has discarded positive charge; its feasibility implies nonnegative actual displacement, but is not equivalent to the native orbit's feasibility. For instance a=1 has actual z'=3z/2+1, not 3z/2. Never silently identify these two processes.

For any fixed exact valuation word a_1,...,a_k, put K_j=sum_(i<=j)a_i and C_0=0,C_j=3C_(j−1)+2^K_(j−1). Then

U^j(n)=(3^j n+C_j)/2^K_j.

The exact integer cylinder is n ≡ (2^K_k−C_k)*(3^k)^(-1) modulo 2^(K_k+1). Here inverse is modular. Final oddness and backward induction force every preceding state to be an odd integer and every valuation exact. Conversely a legal word clearly satisfies this congruence. Positive n gives positive states. A strict first-hit ROOT itinerary additionally requires every proper-prefix state greater than 1; if the word is a ROOT certificate, its endpoint must equal 1. First descent additionally compares each prefix to the ORIGINAL source. Those predicates are not the nonnegative-debt predicate.

The clock vector is (U,T,C)=(k,K_k,k+K_k). Composition adds these coordinates for the same actual path; any protocol that inserts another U step must add its valuation and its odd step. Arithmetic cylinders for AB are the domain of A intersected with the affine preimage of the domain of B, plus all stopping/source guards. This pullback operation is associative because function composition and set preimages are. Retaining the word itself provides an exact, if not size-compressed, guard sidecar. Endpoint congruence without the extra final bit is insufficient: the one-letter word (1) permits n=3 mod4, while n=1 mod2 also includes 1, whose actual valuation is 2.

## 6. Inherited debt cone and its exact boundary

For c>=0 let C_c consist of (Λ,D) with Λ>=1 and 0<=D<=c(Λ−1). It is closed under chronological composition: D_AB=Λ_BD_A+D_B <= c[Λ_B(Λ_A−1)+(Λ_B−1)]=c(Λ_AΛ_B−1), and nonnegativity and Λ>=1 persist. The old c=1 condition D<=Λ−1 therefore survives concatenation of blocks already certified in the cone. For a scalar leaf with δ>=0 this is exactly δ<=λ−1. If only the inequality δ<=λ−1 is assumed with signed δ, the upper inequality still composes for λ>0, but membership in the nonnegative-debt cone needs the separate signs and Λ>=1.

For scalar maps of such certified blocks, f(x)>=c+Λ(x−c)>=c whenever x>=c. This proves survival at BLOCK BOUNDARIES. It does not imply survival inside a block unless each internal prefix is certified or its R is retained. The two-step word x->x−2, then x->x+2 has endpoint (Λ,D)=(1,0) in every cone, but input c=1 fails internally. Its true R=2 detects this immediately.

For native Collatz, P must additionally be retained. The inherited cone closure provides a lawful account once actual guarded blocks exist. It does not prove every trajectory admits a useful decomposition into them, a last descending block, first exit, bounded reset debt, source-owned supply, or all-scale 153N coverage. THM-4735 itself labels universal coverage false for its threshold-two cone certificates and retains the direct-payment problem as open. We do not reproduce its large census or promote its pre-audit status.

## 7. Honest lost coordinates and hostile controls

The scalar monoid answers one predicate: nonnegativity at all specified cuts for a fixed positive-affine scalar itinerary, and its endpoint. It loses word identity, source cylinder, source endpoint, first-hit/first-descent information, and clock unless sidecars retain them. Word (1,3) versus (3,1) has the same multiplier but different affine carry and legal source class. Equality of even the entire scalar summary never supplies a ROOT certificate.

Independent vector debts can be handled coordinatewise only when each coordinate evolves independently and every coordinate must stay nonnegative. Summing coordinates is not exact: a vector update (−1,+1) at initial (0,0) preserves total 0 but violates the first account. General matrix-coupled updates require pullbacks of all prefix halfspaces, not one scalar R. Even entrywise nonnegative matrices do not make those halfspaces equivalent to one scalar threshold.

Shared ports are different from scalar credit. One unit eligible for either of two demands cannot satisfy both unit demands simultaneously, despite each singleton passing its own test. Session33's stronger three-label hostile keeps usable total and all singleton capacities equal: masks {1,2},{3} versus {1,3},{2} with unit capacities; demand (1,1,0) fails only the first bank. Required sidecars are joint eligibility, shared capacities, immutable source owners, and actual clock cells. Applying independent scalar debt summaries to both ports would double-spend an owner. The correct retained object is a network/feasible-region response, not a scalar cone.

## 8. Portable reproduction and stopping condition

Run python3 check.py > normal.json, python3 -O check.py > optimized.json, and cmp normal.json optimized.json from this directory. The checker uses only the Python standard library, exact Fraction arithmetic and explicit exceptions. It verifies pinned source bytes, exhaustive short rational affine words against direct trajectories, all word cuts, associative triples, exact integer thresholds, native cylinders and clocks, cone closure, and named hostiles. No old files or remote state are modified. Finite checks support the algebra and guard implementation; the proofs above establish the all-word theorem. Uniform supply and Collatz coverage remain outside the result.

Final finite receipt: 139688 explicit checks in each mode, byte-identical output. All 210 intervals (including empty intervals) of a 19-letter rational word are checked independently through a segment tree, disjoint sparse table, and binary lifting. The affine exhaustive universe is all 7,381 words over nine rational leaves of lengths zero through four; native controls are all 120 words over valuations 1,2,3 of lengths one through four, with four exact lifts each.
