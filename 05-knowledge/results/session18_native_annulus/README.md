# Native annulus debt at every arrival cut

**Status:** PROVED exact finite-horizon identities and source-preserving endpoint recurrence; FINITE-EXACT computations in the declared universe; REFUTED unrestricted all-integer-scale contraction; OPEN the intended all-N>=4096 target, unbounded child grounding and universal Collatz entry.

This packet concerns the native primitive map P and first arrival in the inherited 627-point frontier V. It does not change to a fully odd U-word, replace the source measure, identify terminal fibres, or infer convergence from residues.

## 1. Inheritance, correction boundary, and live board

Closest mechanism: `native_cancellation_105_20261018g.md`, sections 2--3, fixes each named old charge under N -> 64N and proves the exact new-annulus source norm 153N. The inherited `root_eisenstein_first_arrival_20261018j.md`, section3, supplies the same arrived-polynomial consumer; its section5 already uses the equivalent three-endpoint scale sweep and checks separate-scale contraction through N<2^20 at all arrival prefixes. The sweep itself and that earlier contraction census are inherited, not new claims here. Our additional finite predicate is F(64N,d)<=F(N,d); two separate inequalities F(64N,d)<=0 and F(N,d)<=0 do not imply this ordering. The circuit completion in `eisenstein_circuit_completion_20261018g.md` removes additional complete zero packets but bounds residual mass **in terms of** energy, not energy from above. `guarded_children_105_20261018g.md` grounds 63 specified children; its all-height guarded half-descent is not universal ROOT entry. These descendants were recovered before this derivation. Proof review corrected the comparison notation: all new-source charges in the N ->64N ledger are A_(64N), not the zero extension of A_N; the implementations always used the enlarged-scale charge. The old-source equality is used only after this distinction.

Canonical hostile: a partial zero packet has a nonzero arrived vector. Corrected near miss: a zero packet containing both old and new sources preserves the enlarged energy but need not preserve the **difference** from the old energy. Least-used sidecar: an unordered cross term has a unique birth event, the larger of its two actual arrival depths.

Live board: named source / joint 105 register / P-depth and integer braid clock / terminal fibre / source-norm bank / mixed-packet old baseline. The exact endpoint sweep is a bridge from interval arithmetic to this original source measure, not a replacement dynamical model.

## 2. Original carrier and exact horizon bill

Let S_N be the positive odd integers s<4N with s not congruent to5 modulo8. Write B(s)=4s+1. The charge A_N(s)=(a,b) is the sum over B-chain points m<4N of

- L_(m mod3) below2N, where L=((1,2),(-2,-1),(1,-1));
- U_(m mod3) in[2N,4N), where U=((0,1),(-1,-1),(1,0)).

Let Q(a,b)=a^2-ab+b^2 and beta(v,w)=Q(v+w)-Q(v)-Q(w). Native P takes a=2 for s=1 mod8 and a=1 otherwise, forms(3s+1)/2^a, and strips each resulting5-mod8 braid by x ->(x-1)/4 while adding2 to its integer clock. First arrival in V has terminal u_s, full integer clock t_s and P-edge depth d_s. Sources in V arrive at depth and clock zero. If a source never reaches V, set its arrival depth to infinity; its terminal and clock are used only after an actual arrival. For the comparison N ->64N, put v_s=(-1)^t_s A_(64N)(s) for every s in S_(64N). On old sources this equals(-1)^t_s A_N(s) by fixed-source scale transport. At a general running scale M below, use the separately indexed v_s(M)=(-1)^t_s A_M(s), extended by zero before source admission.

All displayed first-arrival quantities at a finite horizon are defined by direct finite iteration even if some other sources have not arrived. For old sources s in S_N, the inherited identity A_(64N)(s)=A_N(s) preserves their routes and signed charges. Let O=S_N and J=S_(64N)\S_N; the source budget of J is exactly153N.

At horizon d, let W_u be the old arrived vector and H_u the new arrived vector in **actual** terminal u. Define E(N,d)=sum_u Q(W_u(N,d)), the Eisenstein norm energy of sources **already first-arrived in V by d**; it is not the single real-row energy. This is exactly the consumer of `native_cancellation_105_20261018g.md`, section3, equation(4), and the arrival-polynomial norm in `root_eisenstein_first_arrival_20261018j.md`, section3.

Moving sources retain their **original, unaggregated source norm** as collateral R_J(d). We do not coalesce moving endpoints or put their quadratic Gram energy into E. The distinct full-state identity in `root_sparse_boundary_module_20261017i.md`, R5 equation(18), contains both frozen and moving endpoint energies and additional merge terms. This packet does not pay that different full-state energy. No moving endpoint Gram is silently set to zero, and its control remains outside the present certificate. Let M_J(d) be the norm of new sources that have arrived, and R_J(d)=153N-M_J(d) their not-yet-arrived source norm. Define

C_JJ(d)=sum over unordered distinct arrived new sources in one terminal of beta(v_s,v_t),
C_OJ(d)=sum over arrived old/new pairs in one terminal of beta(v_s,v_t).

**Theorem A (exact cut criterion).** For every integer N>=1 and every finite d>=0,

Delta_N(d):=E(64N,d)-E(N,d)
 = M_J(d)+C_JJ(d)+C_OJ(d),

153N-Delta_N(d)=R_J(d)-C_JJ(d)-C_OJ(d).                 (A)

Consequently Delta_N(d)<=153N if and only if C_JJ(d)+C_OJ(d)<=R_J(d). The stronger arrived-budget statement is equivalent to C_JJ(d)+C_OJ(d)<=0. Neither condition follows from the source norm identity alone.

Proof: expand Q(H_u) into its named-source diagonals and unordered new/new cross terms; add beta(W_u,H_u), which expands into old/new pairs. Sum over u. Every source diagonal and every admissible unordered pair appears exactly once. No different-terminal pair occurs. Subtract from153N. This argument does not presume eventual arrival.

A useful equivalent form writes M(N,d) for the norm of **all** arrived sources at scale N and F(N,d)=E(N,d)-M(N,d). Then

C_JJ(d)+C_OJ(d)=F(64N,d)-F(N,d).                       (B)

Thus the stronger sufficient target is 64-scale nonincrease of the source-specific coherence F at each actual arrival cut. It is not an operator contraction on arbitrary weights.

### Unique event ownership and no double spending

Each new diagonal is booked at d_s. Each old/new or new/new cross term is booked at max(d_s,d_t). This gives disjoint event ownership, even when equal depths cause a batch arrival. The bank starts with153N and debits those diagonal and cross events exactly once. Its balance at cut d is the left side of(A). A negative cross term is a signed credit, not a second copy of a source's norm. Norms of moving sources remain in R_J(d).

If X_u is the vector of new sources first arriving at d and Y_u that of old sources first arriving at d, while W_u,H_u denote their previous-cut vectors, the exact event debit is

Delta_N(d)-Delta_N(d-1)
 = sum_u[Q(X_u)+beta(W_u+H_u,X_u)+beta(Y_u,H_u+X_u)].    (C)

Here Q(X_u) includes same-batch new/new pairs and beta(Y_u,X_u) includes same-batch mixed pairs. There is no arbitrary ordering within a simultaneous event.

## 3. Partial and mixed packets retain their boundary correction

Take any disjoint collection of source-authenticated zero packets, each with full signed sum zero in one actual terminal. Let O_u(d) be their aggregated arrived **old** subvector, J_u(d) their arrived new subvector, and B_u(d)=O_u(d)+J_u(d) their exposed boundary vector. Let R_u,T_u be the old and new arrived vectors outside these packets. This notation includes all cross terms **between** packet boundaries by aggregation.

**Theorem B (packet boundary identity).**

Delta_N(d)
 = sum_u[Q(R_u+T_u)-Q(R_u)]
 + sum_u[Q(B_u)-Q(O_u)+beta(R_u+T_u,B_u)-beta(R_u,O_u)]. (D)

Proof: W=R+O and W+H=R+T+B; expand their two squared norms.

For fully arrived zero packets, B_u=0, but the correction is still

-sum_u[Q(O_u)+beta(R_u,O_u)].                           (E)

It vanishes for packets entirely in the new annulus, because O=0. It does **not** generally vanish for packets mixing old and new sources. The scalar lattice example R=(1,0), one old packet member(1,0), and one new member(-1,0) gives actual increment1-4=-3; deleting the mixed pair from both ledgers would falsely report0. This is an algebraic hostile, not an asserted native-source witness.

Before complete arrival, B need not vanish, even for all-new packets. The receipt retains a literal native all-new opposite pair for each independently recomputed scale, including its distinct source IDs, charges, terminal, clocks and unequal depths. At N=4096, the all-new sources276145 and92051 arrive in terminal1 at depths4 and5, respectively, both with integer clock24. Their charges are(-2,-1) and(2,1), so the completed pair has norm zero but the depth4 boundary has norm3. Deleting it early loses a real boundary vector.

Formula(D), rather than the complete-packet slogan, is safe for every cut. It permits any packet sizes and does not require primitive-circuit saturation. Source disjointness is essential. Full-support integer zero-sum completion remains different from primitive-support circuit exhaustion.

## 4. The inherited source-preserving sweep and its coherence derivative

The three-endpoint mechanism is the inherited sweep in the Eisenstein note, section5, written with the step M ->M+1 rather than(M-1) ->M. We give its exact source-level proof to make this packet portable. Let rho(m) be the unique primitive B-ancestor obtained by repeatedly replacing m=5 mod8 with(m-1)/4. For every positive integer M, the only changes from A_M to A_(M+1) come from these original odd m:

1. m=2M+1 moves from upper to lower band, adding L_(m mod3)-U_(m mod3) to A(rho(m));
2. m=4M+1 enters the upper band, adding U_(m mod3);
3. m=4M+3 enters the upper band, adding U_(m mod3).

If ancestors coincide, their changes must be added before computing the source's norm change. No other odd m crosses a cutoff. This proves the recurrence exactly, including integer endpoints. A newly admitted primitive source previously has charge zero.

For a change delta A at named source s, put delta v=(-1)^t_s delta A. At each d>=d_s, update its actual terminal vector by delta v, its energy by Q(delta v)+beta(W_u,delta v), and arrived norm by Q(A+delta A)-Q(A). At d<d_s no arrived ledger changes. The integer clock and path depend on s, not on the scale M. This is an exact sweep over the original measures.

The map to a horizon vector loses source labels only at the explicit summation. The checker restores them with the endpoint ancestor, full route certificate and independently recomputed source rows. The mod15/mod35 registers agree modulo5 and decode to s modulo105; no independent525-state Cartesian replacement occurs.

### The residual analytic obligation in local coordinates

Let delta_s=A_(M+1)(s)-A_M(s), supported on the at-most-three primitive endpoint ancestors above. Set delta v_s=(-1)^t_s delta_s, and sum only arrived changes at the cut d. In each terminal, let D_u be their sum and W_u the original vector at scale M. Then the exact coherence increment is

g(M,d):=F(M+1,d)-F(M,d)
 = sum_u[Q(D_u)+beta(W_u,D_u)]
   -sum_s[Q(delta v_s)+beta(v_s(M),delta v_s)].           (G)

Equivalently, it is the sum of beta(W_(u_s)-v_s(M),delta v_s) over changed arrived sources, plus beta(delta v_s,delta v_t) for distinct changed sources in the same terminal. Grouping ancestors before this expansion prevents counting one physical source twice.

For G_N(d)=F(64N,d)-F(N,d), telescoping gives

G_(N+1)(d)-G_N(d)=sum_(j=0)^63 g(64N+j,d)-g(N,d).     (H)

This is an exact local recurrence for the outstanding cross debt, not a closed scalar recurrence: g still consumes actual terminal vectors W, source charges, arrival cuts and clocks. A uniform proof must control those native terms, or maintain the weaker bank in(A), at unbounded scales. Mod15/mod35 marginals or a finite list of already-paid cuts do not supply that control. There is not even simple nonincrease with the base N: at complete arrival, G_4096=-201674 but G_4097=-201616, an increase of58. The separate grouped-edit checker verifies(G)--(H) on all257 cuts for these65 elementary steps, independently rebuilds the final source universes, and retains this hostile.

## 5. Finite result and explicit failure boundary

The saved receipt declares its precise base-scale interval and all tested horizons. The intended production universe is every integer N=4096,...,16384, not only powers of two and not merely four selected examples. Every source at the enlarged scales is enumerated or introduced by the proved endpoint sweep. All routes are independently authenticated by the canonical edge and a strictly decreasing suffix depth. If a route repeats or exceeds the declared cap, the checker fails; no unresolved route is silently omitted.

For every base scale in this finite interval, the checker tests the **stronger** statement

E(64N,d)-E(N,d) <= M(64N,d)-M(N,d) <=153N              (F)

at every d=0,...,256. Every source is certified to arrive by that cap, so all later finite horizons are identical to the last one. Four widely separated base scales are also rebuilt from scratch by a separate source-enumeration path, including a separately summed H and old/new mixed term at each horizon. The sweep is not the only implementation path.

At N=4096 and d=1, new-only coherence is +204: Q(H)=2064 exceeds arrived-new norm1860. The mixed old/new term is -1129, giving total cross debt -925 and increment935. This rules out simply dropping the mixed term from the sufficient argument. Even the mixed bound is not fibrewise:204 actual terminals have positive local cross debt at that cut. Terminal9 has old vector(-1,2), new vector(2,2), new arrived norm4, and local increment6. Its old arrived sources are9 and6371; its new arrived sources are203889 and407779. Global signed credit, not a terminal-by-terminal arrived budget, pays the total. At N=5997,11002,16384 the peak increment occurs at depth76, strictly before complete arrival, so checking only the final ledger would miss the extremum.

The unrestricted all-N>=1 version of(F), even with the full153N allowance, is false. Its smallest N is7: at horizon8,

old energy=19, Q(H)=1086, mixed=-2,
arrived-new norm=1054, moving-new norm=17,
new/new cross=32, total cross=30,
Delta=1084 >153*7=1071.

The bank is -13. Every horizon for N=1,...,6 passes the full-budget inequality. The N=7 obstruction is outside the intended N>=4096 range and **does not refute that target**. It shows the precise failed implication: same-terminal positive coherence can exceed the remaining moving-source collateral, despite exact source budget and lawful clocks/registers.

The authenticated contiguous finite sweep finds the maximum ratio of full-budget use at N=5332, with peak increment581462 at depth150, against budget815796. This is a finite exact maximum, not an analytic bound for untested scales.

### Hostile 64-scale phase tests

Fixed-source transport does not transport an annulus violation: it fixes old labels, while the next annulus consists of new labels. Direct recomputation on the two rays beginning at the small failures gives these peak increments:

| Base N | Peak increment | First peak depth |153N budget |
|---:|---:|---:|---:|
|7|1097|22|1071|
|448|54514|92|68544|
|28672|2566744|156|4386816|
|8|1301|35|1224|
|512|62145|92|78336|
|32768|2858748|90|5013504|

Thus the initial full-budget failures disappear at the next two tested64-scales. However, the stronger arrived-only bound still fails at N=448,d=1 and N=512,d=1: the cross debts are47 and60, while the moving-source collateral is66707 and76477. This is an actual native example where retaining the moving ledger matters. At N=28672 and32768, the stronger arrived-only bound also passes every cut. The two largest tests lie above the contiguous interval and certify maximum first-arrival depth233. `verify_rays.py` directly rebuilds the six complete original source universes and their enlargements, retains source-row hashes and every horizon, and replays each full native route DAG. No cancellation or violation is lifted solely from scale transport.

No finite interval closes unbounded N, proves all 64-scale rays, or grounds all children. No universal Collatz conclusion follows. The remaining specific task is to prove(F), or the weaker criterion(A), for unbounded native source scales with the same source measure and actual cuts.

## 6. Reproduction, evidence and limits

Run from any directory, with the packet files adjacent:

    python3 -B verify_annulus.py --write
    python3 -B -O verify_annulus.py
    python3 -B verify_rays.py --write
    python3 -B -O verify_rays.py
    python3 -B verify_recurrence.py --write
    python3 -B -O verify_recurrence.py

The first writes the deterministic `annulus_receipt.json`; the second recomputes and requires byte-for-byte equality. No network, third-party package, repository import or source-tree layout is needed. All load-bearing checks are explicit exceptions, never Python assertions. Frontier data are copied from the inherited full receipt, whose original hash is preserved in `frontier.json`; the rank proof is replayed rather than trusted by hash alone.

The receipt records every scale's extremum and minimum bank, the digest of every exact horizon row, independently recomputed source-row hashes with full clock and joint-register schema, four detailed mixed/new ledgers, the literal partial packets, the small obstruction and the authenticated route-DAG hash. It is a reproducible finite computation, not a stored list of every expanded row: the checker regenerates the complete declared source universes. Hashes index the regenerated evidence and do not replace the computations.

All56,219 inherited source rows at N=4096,5997,11002,16384 have the same canonical JSONL hashes after independent rebuilding; those hashes are included in the source manifest. The source manifest records the inherited sources, exact file hashes, present packet hashes and commands. The final verification logs and numerical summary accompany the frozen deliverable.
