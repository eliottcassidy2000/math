# A physical descending parent map, with an explicit choice of carry

A follow-up to the [synchronizing reset](SYNCHRONIZING_RESET.md) and [signed ternary parent overlap](PAIRED_AXIOMS.md). The arithmetic analogy with Collatz is useful here because a smaller parent can be read from an exact **state**, not because a known Collatz theorem gives new tournament counts.

## A fixed tournament can supply two next states

Let the root a have only one outneighbor v. Write (h,e,eta) for the ordinary count, paths ending at a, and paths of the deletion T−a ending at v, respectively. The marked starts read is t=h+eta; 0≤eta≤h because every path of the deletion extends by insertion and eta is just one subset of them.

Attach a transitive block of order k+1 by the one-back-arc seam from its sink b to a. The ordinary count is
```
N=(1+2^k)h+2^k eta.
```
Changing the root after the operation does not change N:

- **Sink (k≥1).** Root at b. It has a as its only outneighbor, and the literal rooted state (h,e,t,deletion) becomes (N,h,N,h). Thus the new carry is eta'=0.
- **Penultimate (k≥2).** Root at the vertex immediately before b. It has b as its only outneighbor, and the [previous reset theorem](SYNCHRONIZING_RESET.md) gives eta'=h and e'=(N−h)/2. Thus the new carry remembers the old absolute count. The two projective states come from the same literal graph.

The sink formula follows from the [one-back-arc transfer](PAIRED_RECURSION.md): with the sink of a nonsingleton transitive block, its full start read is zero, its full end read is one, and both remaining optional reads are 2^k. For an arbitrary old root (not necessarily sole exit), the same sink map is (h,e,t,deletion)↦(h+2^k t,h,h+2^k t,h). In particular it synchronizes a pair if just (h,t) scale by seven. Conversely, those are the two data that can be recovered from the sink output. If the upper graph already has an automorphism fixing its root, extension by fixing all tail vertices preserves it.

The penultimate choice preserves gcd(h,eta), since gcd(N,h)=gcd(2^k eta,h) and h is odd. The sink choice, in contrast, discards the carry: if one continues in this two-coordinate family from (N,0), the whole subsequent pair state stays divisible by N. Thus the two routes serve different purposes: **retain** the carry until a later prime needs it; **discard** it when a low-rank synchronization or the final count is enough. A zero-carry continuation is not a substitute for a fresh prime-capable port.

## Exact descent on the nonzero-carry states

For a target state (N,c) with an interior predecessor, the penultimate choice forces h=c. If q=(N−c)/c, take the unique k≥2 such that 2^k≤q<2^(k+1). Then compute
```
eta=(N−c)/2^k−c.
```
This branch is legal only if eta is an integer with 0<eta<c and the corresponding parent state is physically realizable. Endpoints eta=0,c can be handled separately; for primitive continuing states above h=1 neither can occur.

Equivalently the ratio eta/h lies in (0,1) and the next ratio c/N is 1/[1+2^k(1+eta/h)]. The intervals at successive k touch only at endpoints, so interior states have at most one candidate parent in this lane. The parent strictly decreases the ordinary count (N≥5h). Unlike an ordinary Collatz iteration, this rule cannot have a nontrivial cycle or diverge on states where it is defined. Unlike an arithmetic residue decoder, each reverse lift has a literal graph operation.

A concrete descent:
```
(941,165) → (165,29) → (29,5) → (5,1).
```
All three steps select k=2. The ordinary 5 state is constructed from the rooted transitive two-vertex tournament (H=1) by the two-vertex seam and marking the penultimate tail vertex. Its upper entrance is a different physical source with state (h,e,eta)=(35,12,7); the lower initially has (5,2,1). The first forward step forgets e and synchronizes the pair at 29 and 203. The next steps give 165 and 1155, then the prime **941 and its partner 6,587**, with nonzero carries throughout. Direct subset counting verifies the literal graphs. An alternate ordinary unit attachment gives the count-three triangle; no paired state at 1 or 3 is claimed.

There is a genuinely different branch when a case calls for it. Attach a directed triangle s→b→w→s by reversing just the seam b→a and mark s. The new ordinary and end reads are h'=6h+3eta−e and e'=2h+eta−e, and the new carry is zero. For the end count, old concatenations contribute h−e and the backward crossing contributes t via the unique triangle split with suffix w→s. Deleting s makes b a source of b→w, so no endpoint-at-b path is possible. Starting from the ordinary triangle (3,1,0) this gives full states (17,5,0), (97,29,0), (553,165,0). A separate apex lift from 25 gives 79. Thus the two axes for **79 and 553=7·79** can descend by different physical rules rather than require matching carries at the start. This also illustrates why the zero-carry trap above is scoped to the transitive-tail family, not all attachments.

## Keeping actual rotation while changing the root

The same sink projection splices the [rotational reader bank](PAIRED_ROTATION_BANK.md) into the rooted grammar. The replay supplies a rooted lower/upper pair from the (1295,−109) reader on TT2, initially with (h,t)=(5071,10842) and (35497,37947). Appending one ordinal sink to the upper, without changing its root, doubles its t and leaves h unchanged. Now the two reads align by seven. Apply a common TT2 sink migration and two common crossed singleton keeps to obtain **91,807 and 642,649**, literally verified. The upper graph retains the actual C7 action, fixing the old root and all new vertices. Its automorphism group is exactly C7: the action is free on paths, p=91807 is prime, and p exceeds the graph order.

## The signed Collatz parallel, precisely

On odd positive integers let C_±(x) be 3x±1 with all powers of two removed. Negating negative integers in ordinary 3x+1 gives C_−. Its three **known** cycles are (1), (5,7), and (17,25,37,55,41,61,91); the positive C_+ cycle (1) corresponds to 1→4→2→1 before acceleration. See [Barina et al., known cycles and signs](https://cmst.eu/wp-content/uploads/files/10.12921_cmst.2022.0000025_BARINA_ev.pdf). There is a literal contrast at the two tournament gaps:

- Under C_−, 1 and 3 reach the same fixed point 1, but multiplying by 7 splits them: 7 is on the (5,7) cycle while 21→31→23→17 enters the longer one.
- Under C_+, 21 goes straight to 1; 7 follows 11,17,13,5,1.

This does not imply a correspondence between the dynamical cycles and the forbidden counts. It illustrates how multiplying by seven can change the required branch, rather than commute with a descent. One is a unit, not a prime, but it naturally plays the seed role: the one-vertex tournament has one path, and ordinal multiplication by it leaves the count unchanged. The first two ordinary seeds 1 and 3 have no sevenfold partners; the rooted construction can carry them ordinarily and enter the paired regime later.

The remaining Chan–Pak step is to select a **physical carry for every scalar target**, or supply alternative parent moves for targets where this lane has no legal carry. The unique smaller parent theorem solves descent once that extra coordinate is supplied; it does not solve supply by assuming the full integer graph already forms a tree.

Replay: `python3 collatz_parent.py`; run `python3 collatz_parent.py --fixtures | ./rooted_pair_literal` to count the final rotational graphs independently.

## A genuinely total arithmetic selector, and the physical lift it asks for

There is also a way to make the **arithmetic** part of the analogy unconditional. For odd x choose the sign s with 3x+s divisible by 4 and take its odd part. For x>1 this strictly decreases, because 3x+1<4x. Its result is never a multiple of 3. On the ordinary axis the forbidden parent can therefore only be 7, not 21. On the upper quotient axis x denotes a request for 7x; its two forbidden inputs are 1 and 3, and only 1 can be the parent. This statement concerns the exact arithmetic map, not an existing graph lift.

Both defective parent rays have the form x=(y·2^a−s)/3, with a≥2, y=7 or 1, and s=(−1)^a. All but finitely many can be handled without an additive repair if a physical lift for a short signed itinerary can be constructed. For a≥6, take the opposite sign once and then two quarter steps. The three odd values are exactly
```
x_1=y·2^(a−1)−s,
x_2=3y·2^(a−3)−s,
x_3=9y·2^(a−5)−s.
```
The first transformation uses one factor of two; the other two use exactly two. Now x_3<x and all intermediate values are legal on the respective axes. Thus use the three-step itinerary as **one** descending arithmetic edge, recording its internal signs. The two intermediate values can exceed x, so a proof would need to reverse the entire certified itinerary from the smaller witness, not invoke induction on them. For the finite exceptions a≤5, the ordinary side needs the unit and the seeds 9,19,37,75; the quotient side needs physical upper seeds at quotients 5 and 11. The replay verifies all these base witnesses.

This gives a finite-seed reduction tailored to the holes: a physical inverse for each legal quarter step plus the three-letter bypass would imply the desired two-axis coverage. Neither family of physical inverse lifts is proved here. It is a different target than the transitive reset above, which is an actual physical move on a smaller set of states. The distinction prevents turning an unconditional arithmetic descent into an unjustified tournament theorem.

One infinite part of the upper defective ray has an independent physical bypass already. For a=2·2^j, j≥1, its quotient is
```
x=(4^(2^j)−1)/3 = product_(i=0 to j−1) (1+4^(2^i)).
```
Each factor is the exact count of a transitive tournament with the longest arc reversed (the first is 5). Replacing only that first block by the known count-35 tournament and joining the other blocks ordinally gives **7x** directly. Thus these bad-parent cases are not repaired by pretending there is a count-7 seed. They use a different factorization selected from the exponent. This does not treat the whole defective ray.
