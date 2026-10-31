# Exact source costs and a self-reproducing square-sum seam

**Status: PROVED scoped reductions and composition; FINITE-EXACT source-cost certificates; universal raw tournament spectrum OPEN.** Same-agent proof review with separate computational routes. No independent subagent was used or authorized. All tournament counts below are raw Hamiltonian-path counts, never H/Aut.

## Inheritance and board

The closest proved mechanisms are the [positive four-coordinate block action](../square_sum_spectrum_overnight_20261031b/README.md), [cut/interior coordinates](../square_sum_spectrum_overnight_20261031b/cut_interior_coordinates.py), [deficit-nine family](../square_sum_spectrum_overnight_20261031g/README.md), and [closed subcubic square-sum carriers](../square_sum_spectrum_overnight_20261031f/README.md). The repaired near miss is confusing the existence of a dyadic source with an affordable source for a specified output. The hostile is47: realizable globally, but absent from the dyadic-tail route. The least-used sidecar is the strong-component decomposition, now used on the parent rather than only on the final tournament.

| Lane | Object and operation | Retained information | Decisive result / next question |
|---|---|---|---|
| Anchor | All marked parents, bounded q | F,q,g,h; actual strong-block words | Finite exact automaton despite unbounded parent order |
| Anchor | Dyadic zero-tail source costs | Entry slope and positive integer tail length | Deficit9 has optimal slope512 |
| Niche | Two closed square-sum carriers | Actual cut edge, two seams, source bits | Exact product of Hamiltonian-cycle languages |
| Wildcard | Pythagorean ports under repeated gluing | Both endpoints are squares; exact scaling address | The new seam is itself a valid port |
| Reroute | Arbitrary integral slopes | Actual source, not merely an output integer | Test whether dyadic restrictions caused the early missing values |

## 1. A logarithmic-size representative for every bounded-cost source

Let A be a nonempty tournament, x a marked external vertex, and S its out-neighbors in A. Use the inherited state

    h=H(A), F=H(A+x),
    g=# Hamiltonian paths of A ending in S,
    q=# ordered path covers(P,Q) of A with Q nonempty and first(Q) in S.

The cover may have P empty. The inherited nonnegative cut/interior weights C,L give

    F=h+2C,   q=2C+g+2L,   c=q-F-g=2L-h.

**Injection bound.** q>=|S|h. For every v in S, cut each Hamiltonian path immediately before v. This injects the h paths into covers with first(Q)=v; these images are disjoint for different v. Thus q=0 iff S is empty, and for q>0,

    h<=q, 0<=g<=h, F<=q+h-g<=2q, -h<=c<q.

If c>0 then |S|>=2, since a singleton has L=0. Consequently h<=q/2 for a positive deficit.

**Kernel theorem.** For q>0, the same four coordinates have an actual parent representative with at most

    floor(log2(h^2)) + floor(log2(q)) + 1

vertices. The operation removes only unselected singleton strong components preceding the first component that meets S. The removed vertices and their original positions can be stored as a padding sidecar; this preserves the original integer clock and reconstructs the original tournament. It does not assert that the reduced graph is isomorphic to the original or that arbitrary re-marking ignores the sidecar. Named ports on removed vertices must be retained or explicitly reconstructed.

**Proof.** Write A as its ordered sum of strong components. Its h is the product of their raw path counts. Busch's strong-floor theorem implies H(B)^2>=2^|B| for every nontrivial strong component B. This follows directly at orders3,4,5, and increasing order by3 multiplies the floor by5 while 2^|B| increases by8. Hence all nontrivial components together contain at most floor(log2(h^2)) vertices.

For a block of m vertices, its number E of ordered two-path covers, allowing either path empty, is

    E=sum_{U subset V(B)} H(B[U])H(B[V(B)\U]) >= 2^m,

with H(empty)=1. The block action has q'=E q+q_B h. Once q is positive, each subsequent singleton at least doubles it; at most floor(log2 q) such singletons occur. Before the first positive q, an unselected singleton acts identically on the four coordinates and may be removed. There is at most one first selected singleton. Adding these bounds proves the claim. Strong components may also grow q, so this is a safe upper bound, not asserted sharp.

**CITED input.** Arthur H. Busch, [A Note on the Number of Hamiltonian Paths in Strong Tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/), EJC13(2006)#N3, Theorem1 and Corollary1. The exact floors at orders8,9,10,11 are45,75,125,225. The small floors through8 and the insertion argument are also inherited from [THM-1370-h-spectrum-omits-7-21-all-n.md](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md). Moon's [1966 Theorem2](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on_subtournaments_of_a_tournament.pdf) supplies the strong one-vertex deletion used in the second enumeration route.

## 2. Complete finite source automaton

For a strong block B with cut S_B retain its profile(H,E,p,q_B,s,t_B,F_B). Here p counts covers with P nonempty and last(P) outside S_B; s counts Hamilton paths starting in S_B; t_B counts those ending outside S_B. Appending the block to the parent acts by

    F'=HF+pq+(H-s)g+(F_B-H)h,
    q'=Eq+q_Bh,   g'=(H-t_B)h,   h'=Hh.

All coefficients are nonnegative, and F,q,h are nondecreasing. Start from(1,0,0,1), the marked vertex with an empty parent. Singletons are included as blocks. This start also has the same state as one unselected source-padding vertex if a nonempty parent is required.

For fixed Q,Hcap,Fcap, enumerate all strong blocks with H<=Hcap and every cut, then close the state under these actions while keeping q<=Q,h<=Hcap,F<=Fcap. This is an iff for the bounded state window: SCC decomposition supplies every actual source; monotonicity ensures that no allowed terminal source needs an intermediate state outside the window. Every accepted word constructs an actual source. The number of states is finite. The strong-floor theorem makes the alphabet finite as well.

Profiles with q_B>Q or F_B>Fcap may be discarded: every incoming h>=1 and F>=h, so q'>=q_B and F'>=F_B. This pruning is part of the declared universe, not a heuristic.

The q<=64,h<=31,F<=128 automaton has737 states and849 in-window transitions. Its2350 strong-cut profiles come from all strong tournaments through order7; the order8 floor45 excludes larger components. Every profile was independently checked by literal permutation covers, and every accepted state by adjacency subset-DP on its reconstructed graph. All witnesses have at most9 vertices in this census, though the theorem covers all parent orders.

For q<=256,h<=127,F<=512 there are16607 states and18721 in-window transitions. The profile bank is complete through order10. One route starts with every nauty strong representative through9, and augments the100 representatives with H<=127 to order10. A second route starts at order8 and uses two strong-ear augmentations. The two routes have identical profile sets:157058 before pruning,9499 after the declared q/F filter. The order11 floor225 closes the all-order boundary. See the verification receipts for the remaining computational audit details.

## 3. Exact dyadic entry costs

Appending a zero vertex z means A->z->x. It changes the state by

    (F,q,g,h) -> (F+q+g,2q,0,h).

After r>=1 such vertices the raw count is

    N=q*2^r+F+g-q=q*2^r-c.

The positive length condition is essential: the r=0 expression F+g is generally not F. A dyadic source has q=2^a. For fixed positive c its family begins at exponent a+1, not at an arbitrary exponent chosen after seeing the target.

**Exact small-output theorem.** Across all marked parents and all positive zero-tail lengths, the dyadic route's outputs in[3,127] are precisely

    3,5,9,17,33,41,65,67,69,71,73,91,95,101,111.

For N<=127, q<=64 since N>q. If c<=0, then q<=32 and h<=q. If c>0, h<=q/2<=32. Raw h is odd, so h<=31. Hence the complete small automaton covers every possible source, including arbitrarily large parents. Actual terminal graphs certify every listed output. In particular47 is a route obstruction, not a forbidden tournament count.

**Optimal deficit-nine entry theorem.** No actual source with c=9 has dyadic q<=256. Indeed h<=q/2<=128, and raw h is odd, so the larger complete automaton applies. None of its positive dyadic states has c=9. The inherited actual state(449,512,54,117) supplies c=9 at q=512. Therefore512 is the exact minimum slope, and exponent10 is the earliest positive-length dyadic-tail entry to this deficit line.

The conclusion is a finite exact certificate combined with an all-order reduction, not an extrapolation from a maximum parent order. The five earlier attained values23,55,119,247,503 need other certificates. The small table excludes23. For the other four, a positive-offset representation from the preceding power of two is impossible because its offset is at most q, at most one quarter of the next power; the required offset is larger. Their negative-offset representation would require the excluded c=9 source with q<=256.

The source-cost theorem does not prove universal spectrum completeness. It explains precisely which previously attained values cannot enter this particular useful recurrence. The retained actual words remain the objects on which later operations are allowed.

**Reroute, with actual sources.** Dropping the dyadic restriction recovers47 from(17,30,0,1), then one zero. Likewise23,55,119,247 have states(5,15,3,3), (9,41,5,5), (11,103,5,5), (15,223,9,9), each followed by one zero. Every state and resulting adjacency count is independently audited. The full bounded automaton finds a nonzero-slope one-zero certificate for every allowed odd output below256; its only missing odds there are7,21. This is a restricted-construction control, not a replacement for the already much stronger unrestricted finite coverage through99,999,999.

The next question is whether arbitrary-slope sources can be selected with a recursively controlled cost and actual terminal word. A scalar claim that every N is already an allowable source would merely restate the original problem. The source automaton supplies a precise object on which to test an induction.

## 4. Square-sum continuation

[PYTHAGOREAN-SEAMS.md](PYTHAGOREAN-SEAMS.md) proves the two-carrier connector, its exact cycle-language product and its recursive port. This advances the boundary-composition target on a sparse set of integer labels. It does not replace the stronger inherited existence theorem for initial intervals, or assert arbitrary endpoint/forest extension there.

The connection across the two lanes is operational: composing authenticated boundary languages gives exact product/transfer rules, but a numerical output alone loses the entry cost or port type. The tournament side needs q and the integer tail length. The square-sum side needs a common endpoint square class, two actual cut edges and the scaling address. The general seam is generated by Gaussian norm multiplication; classes1,2,5,10 have actual carrier examples, with one class10 bit explicitly paid. Neither result identifies tournaments with square-sum graphs.
