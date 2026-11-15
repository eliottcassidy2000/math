# The fixed-order near-floor gap is at least n minus a logarithm

**Theorem.** Let f(n) be the Busch minimum for strong tournaments. Put
`s(n)=ceil(3 log_5(225n))+10`. For n≥100, every positive excess
`K=H(T)−f(n)` satisfies
```
K >= n − 2s(n) − 6.
```
In particular the liminf of the least positive excess divided by n is at least 1. This improves the half-order bound in [SLOTS.md](SLOTS.md). It does not assert sharpness, and it is not an exclusion from the unrestricted Chan–Pak spectrum.

The proof reuses the acute feature of small seams: an invalid arc in a Hamiltonian-cycle projection needs a supporting block. This time we also replace **other** blocks of the same cycle, not only the one supporting the invalid arc. This makes it impossible to hide a long chain wing inside a neutral slot.

## 1. The bad-cycle reduction works all the way up to n−2

Consider a split into two exact-floor chains A,B sharing v, with H(A)H(B)=f(n). Let L=H(T)−f(n), and assume 0<L<n−2. By [exact seam rigidity](../chanpak_exact_seams_douglas_20261115/README.md) the cross matrix follows a single interleaving of the factor cycle orders after v. There is exactly one cycle with both projections valid. The case of no other cycle has L≥2f(n)/45 by the Douglas bound, already above the range at issue for n≥100.

For any bad cycle, an invalid projected arc has a supporting segment along the actual cycle. Its endpoints lie in the factor in question, and all intervening vertices lie in the opposite factor away from v. There must be at least one intervening vertex since the projected arc is invalid. A rotation can repair that arc only by cutting an edge in its supporting segment.

In one factor, supporting segments of distinct invalid projected arcs have disjoint edge sets. In opposite factors, their intersection has **at most two edges**: every edge internal to the first segment's interior has both endpoints in the opposite factor off v, so cannot be in the second segment, which has no such edges. Only the first segment's two boundary edges can be shared. Thus with two invalid gaps at most two rotations can erase both, and at least n−2 paths are broken. This contradicts L<n−2. Therefore the bad cycle has exactly one invalid arc, on A say, and its B projection follows its unique Hamiltonian cycle. As before, a=|A|≥4. If D is its supporting B-block, t=|D| and b=|B|, the number of broken rotations is ℓ=n−1−t≤L. The number r=b−1−t of B−v vertices outside D satisfies r≤L+1−a, and b≥n−L.

For an invalid arc at the root, [the block replacement](CYCLE_BUDGET.md) gives H(R)−E_R(v) extra broken paths for R induced by D∪{v}. These paths start in D or at v, whereas original broken rotations never start in D; at most the one rotation from v overlaps. **If R is strong**, the number of replacements is at least t from its spanning cycle, so L≥ℓ+t−1=n−2. Thus strong R is impossible throughout this enlarged range. In particular the seam cannot be one-way: a one-way interleaving puts all of B−v in one cycle block, and a bad A arc can put this block only at v (an internal detour would need cross arcs in both directions). That would give R=B strong.

## 2. A large neutral wing cannot occur even outside the bad gap

Retain the convention of [SLOTS.md](SLOTS.md): a neutral slot is one with no Hamiltonian path of the factor making a backward step through that slot strictly before or strictly after v. Every such slot in a floor chain has a uniform cut X→Y, where X and Y precede and follow it in that factor's rooted cycle list. A positive slot has charge at least H/225 per opposite vertex placed in it. The charge counts valid-projection duplicate fibers and is disjoint from broken rotations.

Suppose a block J of B vertices all occupies the same **interior neutral** A slot, so X,Y are nonempty. In the bad cycle read from v, the projected off-root A sequence must put all X before all Y. Any consecutive Y→X is an invalid actual A arc, and it cannot be supported by a B segment traveling backwards in the Ferrers chronology. Every B vertex of J must consequently occur in the *same* gap of the A projection, from its last X vertex to its first Y vertex. Indeed, a B-block containing any of J cannot enter from Y or exit to X: its entries and exits in the B cycle order bracket that vertex in the Ferrers order. The two gaps incident to v also fail one of these constraints, since the first A vertex is in X and the last in Y.

The identified gap crosses a valid A arc and is not the invalid gap. Its B-block has uniform entrance and exit gates at its two off-root A endpoints. Replace that block by any of its Hamiltonian paths. Each replacement gives a bad cycle with the same faulty arc and the same ℓ broken rotations. Rotation sets from distinct cycles are disjoint. In particular
```
L >= ℓ H(J).
```
Here H of the whole block is at least H(J) by inserting extra vertices into a Hamilton path deterministically, with deletion recovering the original path. This estimate needed no bound on the number of bad cycles.

If v is interior to a B link and b≥s(n), the B slot charges forbid any A vertex within the B list (H(B)/225>n). Its single A slot must be neutral: otherwise the b−1 valid fiber extras and ℓ broken rotations give L≥b−1+ℓ=n+r−1. If the slot is exterior the seam is one-way. If it is interior, the preceding argument puts every B vertex in a valid A gap, leaving no support for the invalid one. Both cases contradict the bad-cycle reduction.

If v is a junction or extreme in B, discard the one adjacent link on each available side of v, leaving two exterior wings, either of which may be absent. Let p_i be the number of links' vertex costs in wing i, equivalently its order minus 1 (zero if absent). Then p_1+p_2≥b−1−7, since adjacent links cost at most 7 together. Consider any wing P with p_i≥s(n), attached at a boundary t0. By converse it suffices to take a left wing. The B-slot localization from SLOTS applies when b≥s(n): all A vertices can be only outside the B list or in its bounded central neutral band. In particular all P−t0 vertices are in a single A slot.

That slot cannot be neutral. If it is interior, apply the preceding block replacement with J=P−t0. Deleting the terminal mark of P still leaves the earlier links through the preceding junction, a floor chain of order at least |P|−4. Therefore H(J)>n by Busch's formula and the choice of s(n), contrary to L<n. If all A vertices precede P−t0, then they precede the entire B list (all slots at or before this wing are positive except the exterior slot), giving the prohibited one-way seam. If all A vertices follow P−t0, peel P off at t0: P−t0 uniformly beats the remaining vertices away from t0. The remainder T' is strong and its floor times H(P) equals f(n), by the floor leaf budget. If T' is nonfloor this costs at least 2H(P); if T' is floor, the [one-way surplus formula](../chanpak_oneway_seams_20261115/README.md) gives either zero or at least 5H(P)/9: the left chain is marked at its overall bottom, has no bad-prefix covers there, and its endpoint share is at most 4/9. Since H(P)>225n, all positive possibilities contradict L<n. Thus every wing with p_i≥s(n) occupies a positive A slot.

## 3. Sum the wing charges and lift through any decomposition

Suppose now L<n−2s(n)−6 in the exact-floor split. Then b≥n−L>2s(n)+6, so the localization applies and at least one wing has p_i≥s(n) in the junction case, since p_1+p_2≥b−8≥2s(n)−1. If a≥s(n), every positive A slot costs H(A)/225>n, contradiction. Thus a<s(n). Summing the slot charge over the off-boundary vertices of *all* large wings gives
```
L >= sum_(p_i>=s(n)) p_i >= b−1−7−(s(n)−1)
  = n−a−s(n)−6 >= n−2s(n)−5 .
```
This too is a contradiction. All comparisons with f use \(f(j)\ge5^{(j-1)/3}\) for strong j≥3, immediate from its three base values and the recurrence.

For an arbitrary strong tournament of order n≥100, outside the stability band the excess is already at least 2f(n)/25>n. Within the band choose a deepest nonfloor node U of order m, whose factors are exact. Its multiplier is at least n−m+1. If 48≤m<n, the [half-order bound](SLOTS.md) gives
```
K >= (m−5)(n−m+1)/2 >= n−6
```
by minimizing the quadratic at m=48,n−1. If m<48, parity gives K≥2(n−46)≥n−6. Hence an excess less than n−2s(n)−6 must occur at an exact-floor root split, which was already excluded.

## 4. What a seam just below the non-Ferrers threshold would need

For n≥500, if 0<K<n−2 then **every Busch split** is between exact floors. Indeed a proper deepest nonfloor node with m<48 has weighted excess at least 2(n−46)>n. For 48≤m≤n/2 its half-order bound gives at least (m−5)(n−m+1)/2>n. For n/2<m<n, the theorem just proved at m≥100 gives at least (n−m+1)(m−2s(n)−6)>n by minimizing the concave product at the endpoints m=n/2,n−1. Hence the root of any decomposition is the only nonfloor node. Stability applies since K<n<2f(n)/25.

In any such split select a bad cycle as in §1. Its **supporting factor B** necessarily has at most 2s(n)+6 vertices. Otherwise b≥2s(n)+7, and §2 already excludes an interior B root. In the junction case at least one wing is large and the slot argument of §3 again forces a<s(n) and at least n−2s(n)−5 valid-fiber extras. But the broken rotations are additional. In either replacement case, with k intact B links, L≥3^k−1 by [the cycle budget](CYCLE_BUDGET.md) and its proof, which applies to the single-invalid-gap situation of §1. This gives k≤log₃(n+1) and thus ℓ≥(n−2)/3−log₃(n+1)>2s(n)+5 for n≥500. In the last inequality substitute s(n)≤3 log₅(225n)+11, check n=500, and note that the resulting difference is increasing there and afterwards. Adding valid extras and broken rotations gives L>n, contradiction.

Now the faulty factor A has a=n+1−b≥s(n). Its positive slots cost H(A)/225>n, so every off-root B vertex lies in a neutral A slot. A nonroot invalid arc of A cannot be supported: any interior B vertex of the supporting segment lies between its two A endpoints in Ferrers chronology, but its neutral slot gives a uniform forward A arc between those endpoints. Thus the invalid projection arc in the large factor must touch v. These are additional, necessary structural conditions; they do not classify near-critical seams.

At L≥n−2 the cross matrix can cease to be Ferrers even on cycle rotations. Controlling and amplifying those cross configurations is a different problem.
