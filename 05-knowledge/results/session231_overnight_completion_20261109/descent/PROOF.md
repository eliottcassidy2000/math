# Adaptive two-sided physical inversion for one-arc joins

Status: PROVED; independent audit PASS. See [audit](../descent_audit/AUDIT.md). Chan–Pak completeness remains OPEN. This is a finite exact inverse and descent theorem for a larger operation family, not an existence proof for a successful branch.

## Inheritance, choice of lane, and information ledger

The one-arc scalar law is KNOWN: THM-4875, `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, and `05-knowledge/results/session193_triangulated_proof_20261105/one_arc_joins/PROOF.md`, section 2. The nonsink bound is inherited THM-4876, statement 7. The rooted four-state oracle is Session228 `physical_state/PROOF.md`, independently audited in that packet. Session229 `inverse/PROOF.md` supplies the closest prior all-gadget descent result, for the narrower proper graft operation.

This lane recovers the older one-arc family because the finite coverage bank in THM-4875 belongs naturally to it, rather than forcing that bank into a fixed graft menu. No coverage computation or source-bank completeness is reclaimed. S86's MISTAKE-739 concerns invalid finite-to-uniform promotion; this packet makes no such promotion.

Concept board: (1) two-sided join law, (2) endpoint bounds, (3) physical-state oracle, (4) divisor inversion, (5) prime SCC normalization. Anchor: adaptive predecessors. Niche: exact physical bank inverse. Wildcard: whether normalization produces strong prime witnesses. Hostile: identical scalar counts do not preserve the required rooted interface.

Source objects are actual rooted tournaments. Their interface is (h,e,u), where h counts Hamilton paths, e counts those ending at the root, and u counts ordered two-path covers (P,R), P possibly empty and R starting at the root. The map to scalars preserves the join count but loses endpoints and response. Physical membership, not a numerical cone, restores the missing information. For membership only, the deletion count d is an existential sidecar; for literal output we retain a graph.

## 1. Proper joins and uniform two-input scalar descent

Take rooted tournaments (A,a) and (C,c). Let B be the converse of C. Start with the ordinal join A→B, then reverse its one arc a→c. Write x=(h,e,u) for A and y=(k,f,v) for C. Its Hamilton count is

    N = hk - ef + uv.                                      (1)

Paths avoiding the reversed arc are the hk ordinal paths, except the ef paths formerly using a→c. Paths using c→a correspond bijectively to one a-starting two-path cover of A and one c-ending two-path cover of B. This proves (1), including empty outside portions. This proof and formula are inherited, reproduced to fix conventions.

Call the join proper when at least one of a,c is a nonsink in its recorded input A,C. For any rooted tournament, u≥h. The inherited endpoint inequality gives e≤(h−1)/2 for a nonsink. For a sink, e=u=h. Hence a proper join satisfies

    if a is a nonsink:  2N ≥ (3h+1)k;
    if c is a nonsink:  2N ≥ (3k+1)h.                       (2)

Indeed uv≥hk and e≤(h−1)/2, f≤k prove the first; the other follows by symmetry. In particular

    h,k ≤ floor((2N−1)/3) < N,
    uv ≤ N−1.                                             (3)

For the first bound on h, use (2) directly if a is a nonsink. Otherwise c is a nonsink and (2) gives N≥2h, which implies 3h≤2N−1 because h≥1. The bound on k is symmetric. Finally properness implies ef<hk, proving the product bound. In particular each response is ≤N−1.

If both ports are sinks, (1) is exactly N=hk. This degenerate boundary must not be admitted as a strict-descent branch: taking h=1 can copy k=N unchanged. A proper branch is genuinely reducing in both input Hamilton counts, independently of gadget size, matrix invertibility, or a finite menu.

## 2. Exact all-input divisor inverse

For a target odd N≥3:

1. Enumerate odd h,k in [1,floor((2N−1)/3)].
2. Enumerate 0≤e≤h and 0≤f≤k, excluding e=h and f=k simultaneously. Optionally impose the nonsink endpoint bound whenever e<h or f<k.
3. Set t=N−hk+ef. Reject t≤0. For every positive divisor u of t, put v=t/u; retain u≥h, v≥k, and t≤N−1.
4. Test physical membership of (h,e,u) and (k,f,v). Return the literal one-arc join of witnesses for every successful pair.

Triple membership is a terminating subroutine, not an assumption. Enumerate positive odd d≤h and call the audited Session228 four-state oracle on (h,e,u,d). The deletion inequality d≤h is inherited. Any physical rooted triple has some such d; conversely a successful four-state test supplies its actual graph. Its explicit representative bounds are 1 when d=u=1; d+1 when u=d>1; floor(log2(u−d))+2 when d=1<u; and d+floor(log2(u−d))+2 when 1<d<u. Tuples with u<d are rejected. The algorithm is finite even without the optional Busch bound.

Necessity follows from (1)–(3) and actual physical witnesses. Sufficiency follows by constructing the two returned graphs and joining them. Thus the inverse accepts exactly the proper one-arc scalar image over ALL finite rooted inputs. It is not only a finite-bank inverse, and no cone point is accepted without a graph.

The two full triples are indispensable. TT2 rooted at its source has (h,e,u)=(1,0,2), while TT3 rooted at its source has (1,0,4). Joining either to the singleton triple (1,1,1) is proper, but the outputs are respectively 3 and 5. A recursive proof cannot replace the source by an arbitrary graph with the same scalar h=1. Decrease of h alone is therefore not a complete induction invariant.

## 3. Prime targets can be normalized to strong outputs

For any rooted A, remove the ordered SCC prefix P strictly before its root SCC, leaving A'. Then its entire rooted triple scales by H(P), exactly as in the Session228 four-state proof. The root reaches every vertex of A'. Apply this to both recorded inputs A,C. If the scaling factors are p,q, bilinearity gives N=pqN', where N' is the normalized join count.

Properness survives: a nonsink cannot have an out-neighbor in an earlier SCC. Thus N'≥3 by strict descent. If N is prime, p=q=1. Both normalized roots reach all their input vertices, without changing either triple.

A proper one-arc join of two such reachable-root inputs is strong. When each input has other vertices, paths from a reach all A vertices, then cross to B, whose vertices reach c, and c→a closes the route. When one input is a singleton, properness forces the other to have another vertex; the same route uses the direct cross arcs from the singleton to the nonroot vertices, or their converse. The case of two singleton inputs is excluded by properness. Conversely, if a cannot reach some A vertex, no path that enters A from B can reach it, and the output is not strong; the dual obstruction applies to C. Hence reachable-root normalization gives strong prime witnesses whenever a proper representation already exists.

This is a normalization theorem, not existence for every prime. It also implies primitive physical triples for prime N: gcd(h,e,u) and gcd(k,f,v) divide N by (1), and each is smaller than N by (3), so both equal 1.

## 4. Relation to prior coverage and precise missing step

Session229 made unrestricted rooted-graft inversion physically finite and decreasing. The present reduction does the same for the two-sided one-arc join interface used by THM-4875's large finite coverage banks. The two operations can be used as independent inverse branches; no graph-class containment claim is needed.

The missing claim is now explicit: for every admissible odd N, the divisor search must select two PHYSICAL triples satisfying (1), or another authorized constructive branch must succeed. The proof supplies neither that selection nor a finite universal bank. THM-4875's finite coverage is evidence about named ranges only. Even a successful scalar smaller-parent induction would not suffice unless it produces the required endpoint and cover coordinates. No new unrestricted interval, infinite prime family, or resolution of Chan–Pak is claimed.

## 5. Exact regression controls

`verify.py` uses the unchanged Session225 literal DP in `inherited.py`, copied from Session229's inverse packet. It takes all 1,099 labelled rooted-at-zero tournaments of orders 1–5, checks (2)–(3) on all ordered pairs, and checks (1) against independently constructed literal joins for all factors of orders≤4. It compares the divisor inverse with direct pair-image enumeration in the supplied physical triple bank for every odd target 3–41. This is an exact finite bank control, not a complete global image computation. The TT2/TT3 scalar-replacement hostile and the both-sink product boundary are included. All checks use explicit exceptions and survive python -O.

Run `python3 -I -B verify.py > replay.json`; compare with CHECK.json. Repeat with -O. No old source or repository file is edited.
