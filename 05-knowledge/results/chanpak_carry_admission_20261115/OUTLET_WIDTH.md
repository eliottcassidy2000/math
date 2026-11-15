# Outlet width: admitting a richer core by exact physical readouts

This is a different route from the binary carry. Instead of first choosing a rich tournament and guessing which microscopic coordinates will matter, start at what the final host can read. For uniform tournament substitution through at most (d) outlets, the complete Hamilton interface is **exactly the first (d+1) slot counts**. It has a short, physical test set. The first invisible difference can be read through a host having the minimum possible width and top multiplier.

The first-slot identity for one outlet is inherited from the [one-outlet theorem](../chanpak_chronological_core_changes_20261114/ONE_OUTLET.md); fixed-slot quantization and the example in §4 are inherited from the [slot-spectral note](../chanpak_crossing_scales_20261115/SLOT_SPECTRAL.md). The all-slot tridiagonal identity, arbitrary outlet bound, exact test set and sharp factorial multiplier are the extensions here. This packet is an interface theorem, not a proof of Chan–Pak.

## 1. Starting from the reader: an exact positive transfer

For a nonempty tournament (U), let (L_k(U)) count its ordered partitions into (k) **nonempty** directed paths. Let (Q_j(U)) count (j) ordered paths allowing empty slots, so (Q_1=H) and
```
Q_j(U) = sum_(1≤k≤j) binom(j,k) L_k(U).
```
The same formula gives finite triangular inversion from (Q_1,…,Q_j). For the empty tournament each (Q_j=1), including (Q_0(∅)=1); if (U) is nonempty then (Q_0(U)=0).

In a host tournament (T) with a marked vertex (v) of outdegree (d), replace (v) uniformly by (U); each vertex of (U) has (v)'s arcs to the exterior. In every literal ordered (j)-cover of (T[v←U]), compress each maximal consecutive (U)-run in a path to one placeholder. Keep the exterior vertices, the order of paths and the locations of their placeholders. This skeleton knows nothing else about (U). Each nonterminal run exits to a distinct exterior outneighbor of (v), while there are at most (j) terminal runs. Therefore the total (k) of distinct placeholders is at most (j+d).

Conversely a skeleton with (k) placeholders can be filled with the (k) nonempty paths of any ordered (k)-cover of (U); all contacts are permitted by the uniform arcs. These mutually inverse constructions prove
```
Q_j(T[v←U]) = sum_(1≤k≤j+d) a_(T,v,j,k) L_k(U),   a_(T,v,j,k) >= 0.
```
Here (a) counts only exterior skeletons, including path-slot labels. In particular two cores can be substituted interchangeably for Hamilton reads in **all** hosts with outdegree at most (d) if their first (d+1) slot counts agree. No parity or size match is needed. By taking converses the same conclusion holds for a bound on indegree instead.

For a *prescribed target* this gives an exact equation with a literal admissibility requirement: find a positive profile ((L_1,…,L_{d+1})) satisfying the reader's linear equation **and** realized by one physical tournament. No separate compatibility assumption for higher core slots is needed once these are fixed. Conversely the evidence below shows that in general the (d)-th truncation cannot replace this requirement whenever a pair with a hidden ((d+1))-st difference exists.

## 2. The first new slot has a sharp factorial tariff

Consider (j=1) and (k=d+1). All (d) exits are used, every exit leads to the beginning of one exterior excursion between two placeholders, and the last path position is a placeholder. Fix any exterior prefix before the first placeholder. The (d) nonempty excursions are freely permutable: each begins at an outneighbor of (v) and ends at an inneighbor of (v). Uniformity means every permuted list has the same validity. The first vertex of each excursion is distinct, so the orbit has exactly (d!) elements. Thus the top skeleton coefficient is always divisible by (d!).

A sharp reader needs an explicit host, not just a bound. Let (F(U)) be the tournament obtained by two new vertices (q,r) with all arcs (U→ q→ r→ U). Iterate it (d) times. The original (U) then faces precisely (d) exterior outneighbors (one per nested layer). Its top coefficient is exactly (d!), as the recurrence in §3 proves. Hence if the first differing slot of two cores is at (Q_{d+1}), then:

- every uniform host of smaller outlet width has identical Hamilton counts;
- at width (d), every Hamilton difference is a multiple of (d!ΔQ_{d+1});
- a literal width-(d) host attains exactly that difference.

The assertion is for a fixed pair and uniform contacts, not for two freely differing hosts or nonuniform access to the inside of (U).

## 3. An all-slot one-outlet identity and a complete calibration

This independent derivation explains why the iterated gates attain the tariff and gives the signs in the full one-outlet transfer. Let (v) have sole outneighbor (q), and (R=T-v). For each (j), write (A_j=Q_j(R)). Let (E_j,O_j) sum (Q_j) of the exterior cycle complement over the even and odd directed cycles of the original host through (v), respectively. A crossing cycle in the glued tournament includes (q), so there is at most one crossing cycle per vertex-disjoint odd packing and at most one (U)-run in it. Contraction of this run recovers the original exterior cycle. Colors of the crossing cycle contribute a factor (j). Cycles are not rooted multiple times; one recovers them by contracting the run and regarding directed cycles modulo rotation.

For the nonempty (U)-run and remaining (j)-cover of (U), the unsigned all-subset sum is (Q_{j+1}(U)), and the sum signed by the order of the run is (Q_{j-1}(U)). For the latter, use the [Rédei–Berge odd-packing expansion](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20). Every unused singleton or odd cycle of a glued packing independently contributes either (-1) if colored in the run slot, or (j) in the other slots. Thus the signed weight is the usual specialization at (j-1). Taking even/odd parts and removing the empty run in the even part proves
```
Q_j(T[v←U]) =
    (A_j − 2j E_j) Q_j(U)
  + j(E_j + O_j) Q_(j+1)(U)
  + j(E_j − O_j) Q_(j−1)(U).
```
For the gate (F), (A_j=j^2), (O_j=1) and (E_j=0), so
```
Q_j(F(U)) = j² Q_j(U) + j[Q_(j+1)(U) − Q_(j−1)(U)].
```
After (d) iterations, (H(F^d(U))) is a triangular integer combination of (Q_1(U),…,Q_{d+1}(U)) with leading coefficient (d!). The same leading coefficient applies to (L_{d+1}) since (Q_{d+1}) contains it with coefficient one. This proves the sharpness claim in §2 and gives an exact test set: **the (d+1) ordinary Hamilton counts (H(U),H(F(U)),…,H(F^d(U))) recover the entire interface**. They are necessary and sufficient for all width-(d) Hamilton consumers to be unable to tell two cores apart. This makes a richer core admission test independent of access to a dedicated slot counter.

If the core order is known to be (n), sufficiently many such observations also reconstruct its full slot polynomial: its form is (j^n+2sum_{s=1}^{⌊(n-1)/2⌋} b_s j^{n-2s}), by the odd-packing expansion. Knowing the same number of distinct nonzero (j)-evaluations as the number of unknowns is sufficient by a Vandermonde system in (j^2); transitive singleton and order-two cores have no unknowns. This is an exact reconstruction statement about the polynomial, not reconstruction of the tournament isomorphism class.

## 4. Backward use on a hidden physical switch

Take the inherited two order-eight cores with outgoing bit rows:
```
U = (92,141,32,228,206,19,38,101)
V = (52,229,208,231,10,148,49,81)
```
Their (H=345) and (Q_2=5640) match, but (Q_3=45225) and (44505). The earlier symmetric reader (C_3[U,TT_2,TT_2]) magnifies the split by four. The nested gates instead read the first new slot at **sharp** width and minimum multiplier:
```
H(F(U))   = H(F(V)) = 5985
H(F²(U))  = 118305
H(F²(V))  = 116865.
```
Their gap is (2·720=1440). Any width-two uniform reader for this pair has gap divisible by 1440. Thus to get an exact narrower step from these two cores, one must change something else, such as some exterior parameters along with the core, or admit nonuniform contacts. This does not forbid +2 switches in other families; it explains why just buying a third visit cannot make this particular hidden switch act as a unit.

## 5. Second round: a selective port bypasses the factorial tariff

The failure of the uniform read explains exactly what to vary. Allow each arc joining (U) to an exterior vertex to depend on its incident vertex of (U), and let (d) be the size of the union of exterior outneighbors of any (U)-vertex. The run-count bound (k≤j+d) still holds. A skeleton must now carry an entrance subset (vertices of (U) beaten by its exterior predecessor) and an exit subset (vertices beating its exterior successor) at each placeholder. With no predecessor or successor use all of (U), respectively. Counting ordered (k)-covers with each first and last vertex restricted to these listed subsets gives a positive skeleton law exactly as before. It does not collapse to the unmarked (L_k)'s. In particular free permutation of excursions no longer keeps their **ordered endpoint types** fixed.

This low-width escape is literal for the same two cores. Append one apex beating the three core vertices given by bit mask 49 on (U), and by bit mask 152 on (V). In other words the masks are `{0,4,5}` and `{3,4,7}`, and all other arcs go into the apex. These strong tournaments have (H=1939,1941), a genuine +2 instead of the minimum nonzero uniform gap 1440. Since the masks have the same size, one may first relabel (V) and use a common fixed marked set on both cores. Only one exterior vertex is required. The result is a localized constructive witness, **not** a switch with its whole selective continuation interface matched or an all-target renewal. It gives a precise next search target: constrain typed path covers instead of forcing a third unmarked slot to buy a unit correction.

Replay: `python replay_outlet_width.py` checks the two cores and nested gates by literal subset/end-point DP, compares skeletons with physical readouts, and independently verifies the tridiagonal identity from enumerated exterior cycles.
