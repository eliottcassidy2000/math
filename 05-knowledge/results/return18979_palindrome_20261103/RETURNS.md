# Four literal returns and the size obstruction

2026-11-03. **FINITE-EXACT**, with an all-order consequence using inherited
proved classification and a cited order bound. Same-agent proof review;
independent computational paths, no external proof audit.

## R1. Source-preserving receipts, with explicit deletions

The sources are the actual two children in
[return6325_completion/README](../return6325_completion_20261103/README.md).
Their live root names are39 and43; their sole in-port is31 and spare is6.
`returns.json` records every relation, retained name, helper cut, legal
Hamilton guard, deletion and final two-birth continuation.

|Source H|Ready defect|Target order|Old names deleted|Native complements|Next H|
|---:|---:|---:|---|---:|---:|
|18977|1|16|none|20|56933|
|18977|2|14|5,13|18|56935|
|18979|1|14|2|24|56939|
|18979|2|13|2,5|22|56941|

These are actual operations, not a free replacement by an equal-H class.
Each completed arc repair and each final child is strong. Intermediate
helper-expanded native states are not asserted strong. Temporary helpers
number half the native moves, one live at a time; peak repair orders are
17,15,15,14. The final children have18,16,16,15 vertices.

Each arc repair inserts a helper, applies two off-spine complements with
their literal Hamilton guards, then deletes the helper. Matrix replay shows
that exactly the designated old arc reverses. The live root's entire
remaining star and its common entry/exit frame are preserved. The earlier
deletions are separately recorded; they must not be described as retaining
all old information. Root, in-port and spare survive every branch.

The exact returned profiles `(H,start_r,H(T-r),E_r(2))` are checked by an
independent forward endpoint DP. The inherited reset gives
`H_next=H+2E_r(2)=3H+2b`. This realizes all four displayed values but does
not return those four children to an input grammar.

## R2. Complete one-ready value fibres at18977/18979

By [READY sections1--3](../two_entry_funnels_20261103/READY.md), every strong
one-ready target is the indegree-one newborn of a strong unique-start
parent. The inherited [ready inverse O2](../ready_return_pentagon_20261103/OBSTRUCTION.md)
gives `H_parent <= (H_target+1)/2 <=9490`. Busch's sharp strong-order bound
has `f(19)=15625`, hence the parent has at most18 vertices.
[Published source, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

`funnel_bank.cpp` compiles the complete rooted unique-start grammar through
order18: **3524577 classes**. Its six positive parents are:

|Target H|Parent order|Canonical rank|Parent H|Parent E(2)|
|---:|---:|---:|---:|---:|
|18977|12|2368|729|18249|
|18977|12|3939|573|18405|
|18977|13|3724|1083|17895|
|18977|15|672|9441|9537|
|18979|13|3625|1543|17437|
|18979|13|4282|903|18077|

Thus strong one-ready targets at18977 have orders13,13,14,16; those at18979
have order14 only. These are rooted classes, not six labelled tournaments.
The saved18979 source has15 vertices. **Any H-preserving return from it to
a strong one-ready target must lose at least one old vertex**, if retained
names remain distinct vertices. Adding helpers cannot remove this endpoint
obstruction. Weak transitive padding lies outside that assertion; O1 of the
inherited inverse explains it explicitly.

The two-ready discovery is smaller in scope: all binary tails through
parent18 over the complete1387-object two-entry core bank through order9.
It yields three targets, including those used above. No complete all-order
two-ready inversion at these values is claimed.

## R3. Graded query compilation

The compiler inherits the exact P/E/S funnel recurrences from
[CALCULUS and FUNNELS](../polarized_path_calculus_20261103/FUNNELS.md).
At size n it stores bin coordinates only through `q=B-n+2`, B=18. A funnel
adds at least one vertex and requests at most one additional child bin,
so all dependencies are present. Shared subcores are evaluated once.
The three arrays retain starts and ends separately; the two-entry compiler
also retains the joint cross term C. Removing C changes the count.

The class order agrees with inherited rank/unrank. The final graded entries
satisfy `P_n(q)<=n!*binom(n+q-1,n)<=19!`, so signed64-bit storage is valid
at this bound; combinations are checked in128-bit temporaries. This is an
exact query algorithm on this grammar, not an A000568 extension or a
minimal carrier for arbitrary strong atoms.

Independent forward endpoint DP checks every class through9 and three
specified ranks at every order10--18,636 queries total, plus all nine
positive targets. The complete compiler search is replayed separately.

## R4. The new clone sidecars

For the four returned roots, the inherited CL1/CL2 formula gives:

|Returned branch|H after replacing root by TT_m|Defect m=1|Defect m>=2|
|---|---|---:|---:|
|18977,+2|19072*2^(m-1)-95|1|1|
|18977,+4|32974*2^(m-1)-13997|2|31|
|18979,+2|34820*2^(m-1)-15841|1|27|
|18979,+4|36070*2^(m-1)-17091|2|115|

The earlier defects3 and27 therefore lead to an exact larger sidecar, not
a closed two-symbol grammar. The first branch has reached the stable class
classified in [STABLE, S1](STABLE.md). The others expose deletion counts
31,27,115. Full parent interfaces are in `audit_returns.json`; merely
retaining the visible ready bit would discard them.
