# A three-letter Chan–Pak construction, its prime arithmetic, and an exact module interface

2026-11-28. The [Chan–Pak conjecture](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf) asks whether the Hamiltonian-path counts of tournaments are exactly the positive odd integers other than 7 and 21. The new finite range below is a statement about a particularly restricted constructive grammar, not a claim to have exhausted the repo's earlier verification range or proved the infinite conjecture.

## 1. A root pointer gives a near-surjective integer grammar

For a rooted tournament (T,a) use four registers (h,e,u,g). They are the Hamiltonian-path count h, the number e ending at a, the deletion count g=H(T−a), and

$$
u=\sum_{S\ni a} H(T[V\setminus S])\,\#\{\text{Hamiltonian paths in }T[S]\text{ starting at }a\}.
$$

The empty tournament contributes one empty path. Start from a singleton with registers (1,1,1,1). Append a new site b, orienting all earlier sites to b except possibly the current root. Three legal edits are:

| Edit | Root/boundary choice | New registers |
|---|---|---|
| P | All earlier sites point to b; keep a as root | (h,0,2u,g) |
| K | Reverse just a→b; keep a as root | (h−e+u,g,2u−g,g) |
| M | Reverse just a→b; move the root to b | (h−e+u,h−e,h+u,h) |

For the flip, e old paths are lost, those that would have led directly from a to b as last site. Every gained path uses b→a: an arbitrary first path on the sites before b, then b, then a path starting at a on the remaining sites. This gives u. If the root stays at a, its new ending paths must use b→a at the very end; there are g. To compute the new u, placing b in the first part gives u, and placing it in the root-starting part gives u−g (append b unless that part would formerly have been the singleton a). If the root moves to b, the new u is h+u: either its part is the singleton b, or that part starts b→a. For P, b can occupy either part, without changing that part's count. These arguments apply to empty parts and the initial singleton.

This enlarges the [one-apex append machine](../tournament_computational_interfaces_fractran_20261127) by permitting a moving pointer. Every site still has at most one backward arc to an earlier site, and the root index never decreases. All edits use literal tournament arcs. The scalar h alone is not the whole state: the current observation and the three chronological suffix probes K, PK and KK give h, h−e+u, h+2u, h−e+3u−2g, and therefore recover all four coordinates.

**Finite result.** Exact state exploration gives **149,972 of the 150,001 odd heights through 300,001** from the singleton, with a maximum shortest word length of 20 in that bank. There are 29 misses in this grammar: 7, 21, and the 27 targets explicitly furnished by the fallback masks. The last grammar miss is **2,647**; every odd height from 2,649 to 300,001 has a literal witness from this one grammar. The fallback masks independently give the complete allowed interval through 300,001. This is data about an unusually narrow construction class, not cofinite coverage proved by the data.

To make the finite computation exact, h never decreases: u≥g≥e, since the singleton root suffix supplies g and truncating a path ending at a gives e≤g. The exploration discards a P when h+2u exceeds the cutoff, because all subsequent P's leave h unchanged and only double u; the next K or M would then exceed the cutoff. Trailing P's do not contribute other heights. With this pruning all registers are bounded, and deduplicating exact states is safe by the edit table. The C++ replay visited 56,438,310 states (it needs several GB of RAM). Python checks all generated words against the update rule, independently checks 1,093 complete rooted states with literal subset DP, and checks all fallback masks by subset DP.

A useful form for an exact-integer inverse is to subdivide by the last edit: P, K and M land respectively on e=0, e=g, and u=h+g−e. On the K face consecutive K's keep g fixed and send h+u−g expressions along dyadic progressions: after the first K on that face the difference u−h stays fixed, or equivalently (u−g,h−u)→(2(u−g),h−u). For M the predecessor h,e,u are recovered as (g',g'−e',h'−e'); only the old g has been discarded. The missing proof step is an integer descent that selects a smaller compatible state for every sufficiently large odd h, not just descent for a known word.

## 2. A modular strategy with an exact all-odd-prime theorem

Modulo every **odd prime p**, the reachable register vectors are *exactly* the three faces above. This includes every compatible vector, not just every scalar height residue. In particular there is no single odd-prime obstruction to that integer descent.

Proof. On a K face use x=u−g and y=h−u, keeping g fixed. K doubles x and fixes y. Follow PK by inverse powers of K to cancel its dilation; similarly for PPK. The resulting maps are

$$
E:(x,y)\mapsto(x+g/2,y-x),\qquad
F:(x,y)\mapsto(x+3g/4,y-3x-2g).
$$

These are invertible on the finite slice, so any inverse in this proof is implemented by a positive power of the corresponding forward word. If g≠0, E can set x arbitrarily. For p≠3 the commutator of E with KEK^−1 is a pure y translation by ±3g/4, so y is independently arbitrary. For p=3 use F, which fixes x and translates y by g. This proves transitivity in any nonzero-g slice once entered. The seed lies in g=1; an MK applied to a K-state with h=g' lands in slice g', with x=2u and y=0. This enters every slice.

At g'=0, x≠0 is already arbitrary by MK, and E then fixes x and translates y by −x. For x=0 and nonzero target y, take the nonzero slice g=−y and choose h=u=0 there; applying MPK gives the g'=0 state (x,y)=(0,−g). The zero state itself comes from MK. Thus the K face is full, and P and M on it surject onto their respective faces (solve the displayed linear maps). Inclusion in the three faces follows from the last edit, with the initial state already on K. Their union has 3p³−3p²+p vectors by inclusion–exclusion.

The same union was observed at moduli 9 and 25, but the proof given here is prime by prime, not a simultaneous CRT argument. Inverses used here are modular periodic words; they are not legal negative moves over the positive integers. This is a sharp example of full local reachability without an exact global descent, in the spirit of the earlier modular one-apex work.

## 3. A two-block construction and a finite negative test

Take independent rooted tournaments (A,a),(B,b); point every cross edge from A to B, then reverse a→b alone. Use the left port (h_A,e_A,u_A,g_A) and compute the right port on the **converse** of B, rooted at b. There is an exact bilinear law

$$
H=h_Ah_B-e_Ae_B+u_Au_B .
$$

The lost paths used a→b as their unique cross edge. A new path uses b→a, so consists of an A-prefix, a B-prefix ending at b, an A-suffix starting at a, and a B-suffix. The two splits are independent, with the A-prefix and B-suffix allowed to be empty. The formula gives a targeted search for exact n by pairing port values instead of searching orientations. It also admits a decisive pruning rule: both h_Ah_B and u_Au_B are at most the target, since e≤h and e≤u. As a test, none of the 27 remaining finite misses can be obtained from **two ports of the singleton grammar** by this one bridge (49,473 small states tested). This is a useful negative answer for that restricted composition, not an obstruction to two independently designed ports or multiple edited edges.

## 4. An exact module-substitution interface, and a no-go result for scalar gadgets

Suppose a template tournament has one site replaced by a module X: all outside sites see X uniformly. Write

$$
P_X(t)=\sum_{\phi:V(X)\to[t]}\prod_{i=1}^t H(X[\phi^{-1}(i)]).
$$

By the [Grinberg–Stanley odd-cycle collection expansion](https://arxiv.org/html/2307.05569), this is a polynomial:
$$
P_X(t)=\sum_{\mathcal C}2^{|\mathcal C|}
      t^{\,|X|-\sum_{C\in\mathcal C}(|C|-1)} .
$$
Here the collections are vertex-disjoint, nontrivial directed odd cycles. The identity holds because each cycle in the collection must be monochromatic, while unused vertices are free. Empty colors are allowed. The count B_j(X) of ordered j-part path covers with no empty parts is the binomial finite difference of P_X at 0. When a path through a substituted tournament visits X exactly j times, each legal outside skeleton has exactly B_j(X) fillings, because the module interface is uniform. Thus equality of P_X is sufficient for indistinguishability in **every** unmarked one-module substitution.

It is also necessary. Probe X at one vertex of a directed triangle and replace the other two sites by transitive r-chains. A visit away from X and back must use both chains, so there are at most r+1 visits to X. Exactly r+1 forces the skeleton X(YZX)^r, with all chain visits singleton. Its coefficient on B_{r+1}(X) is (r!)², while every other term is a known linear combination of B_j for j≤r. Starting with the zero-chain probe H(X) recovers successive moments. For X and X' on the same n vertices only r=0,...,n−1 are needed. This gives a complete test family for the stock of interchangeable modules; it is not a claim that every moment vector is realizable.

A new minimal collision control from the replay: on five sites, lexicographic forward masks 69 and 76 both give H=15, but P_X(2)=108 and 120 respectively. A directed-triangle substitution with two singleton sites gives H(X)+P_X(2), so their new counts are 123 and 135. Up to four sites, H and order do determine P_X because only the empty collection and single triangles contribute.

Finally **no fixed unmarked one-module gadget can universally send H to 2H±1**, or indeed produce any additive shift as a function of H alone, over all module sizes. If its template site lies on a directed triangle, substitute a transitive m-chain, with H=1 for all m. Every odd nonempty subset of the chain gives an odd directed cycle through the other two triangle vertices, so the substituted count is at least 1+2^m. If the template site lies on no directed triangle, it lies on no directed cycle (shorten a smallest cycle through the site using its chord to the third vertex). It is a singleton strong component, and ordered-join multiplication makes substitution simply H↦cH. A successful universal additive operation therefore needs marked incidence, type guards, or more than one unmarked module.

## Replay and next experiments

~~~sh
g++ -O3 -std=c++17 search.cpp -o /tmp/chanpak-shared-registers
/tmp/chanpak-shared-registers 300001 > /tmp/chanpak-register-bank
python3 replay.py /tmp/chanpak-register-bank
g++ -O3 -std=c++17 bridge.cpp -o /tmp/chanpak-bridge
/tmp/chanpak-bridge
~~~

For a light run use python3 replay.py alone (all independent small controls and masks). The complete word bank has SHA-256 706222c74a50b75986a6abafafac08fb2d3378e77c9963d02e70a0aea50cb147. The transcript is results.out.

Prioritize two next tests. First, search for bounded affine sections on the three integer faces that allow a certified descent of every sufficiently large h, using the prime theorem to rule out modular explanations for failure. Second, add a second marked root and measure its *joint* path-cover moments before trying a two-edge additive gadget; the one-bridge negative test and the unmarked module obstruction explain why merely making a wider isomorphism bank is unlikely to be the right interface.
