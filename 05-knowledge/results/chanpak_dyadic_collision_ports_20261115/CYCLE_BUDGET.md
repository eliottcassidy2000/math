# Add all broken rotations to the block replacement budget

Let f(n) be the Busch minimum and H(T) the Hamiltonian-path count of a strong tournament. **If n≥9 and K=H(T)−f(n)>0, then 3K≥n−2.** In particular K≥2 ceil((n−2)/6). The logarithmic correction in the [second-gap proof](../chanpak_second_gap_20261115/README.md) is not needed. This strengthens a fixed-order forbidden band, not the unrestricted spectrum.

## 1. Counting disjoint sources of surplus

Use the exact-floor seam setup of the earlier proof. Its nontrivial case is a union of exact-floor A,B sharing v, with cross matrix Ferrers and local excess L<(n−1)/2. If there is more than one Hamiltonian cycle, then an extra cycle has one invalid projection arc on A of size a≥4 and a valid B projection. All other A projected arcs are valid. The invalid arc is supported by a consecutive B-block D of t vertices, and there are r=b−1−t vertices of B−v outside that block, with b=|B|. Write ℓ=n−1−t=a+r−1. Exactly ℓ rotations of this cycle have an invalid projection: the excluded cuts are the t+1 edges in the supporting segment. Let k be the number of B links untouched by deleting the contiguous run consisting of v and the r outside B vertices. As proved in the earlier packet, b≤3(k+r+2)+2. Therefore
```
n <= 3(k+ℓ)+10−2a <= 3(k+ℓ)+2.
```

**Nonroot invalid arc.** Uniform Ferrers entrance and exit gates make every Hamiltonian path of D a replacement, producing a bad spanning cycle. Every replacement has ℓ broken rotations, not just one. The sets are disjoint: a spanning path determines the arc that would close it. Consequently L≥ℓ H(D)≥ℓ·3^k≥ℓ+2k.

**Root invalid arc.** By converse take it to leave v. Put R=B minus the outside vertices, so its vertex set is D union {v}. The strong R case is impossible here: rotations of a spanning cycle of R produce at least t>L replacements whose A projection starts with the invalid arc. Otherwise v starts a spanning path of R, so it lies in the first strong component and every spanning path of R has a nonempty suffix after v. Each such path P,v,Q gives a broken spanning path of T: use P,v,Q, then the original suffix beginning at the A successor of v in the bad projection and ending just before v. The replacement path starts in D or at v; no broken rotation of the bad cycle starts in D. At most the one original rotation starting at v is counted twice. Thus L≥ℓ+H(R)−1≥ℓ+3^k−1≥ℓ+2k.

It follows locally that
```
n <= 3L+10−2a−3k.
```
The unique-cycle case is covered by the earlier Douglas sorting ratio, and L≥(n−1)/2 already implies 3L≥n−2. To lift to arbitrary strong T, the deepest nonfloor node is necessarily the root when K<n/3 and n≥9: the earlier local m/5 bound and omitted-leaf multiplier make a proper node cost at least 2(n−1)/5>n/3. Outside the stability band 2f(n)/25≥n/3 already. Thus 3K≥n−2 at every indicated order.

## 2. A structural constraint on the other gaps

Still assume L<(n−1)/2 in the bad-cycle setup. Every B-block in an A gap whose two endpoints are not v is **transitive**, whether or not that gap is the invalid one. Otherwise its at least three Hamiltonian paths can each replace it by the same uniform-gate argument. The invalid arc in the A projection and the cardinality of its supporting segment are unchanged, giving at least 3ℓ distinct broken paths. Section 1 also gives L≥ℓ+3^k−1 in either case. If ℓ≥k+1, then 3ℓ≥(3(k+ℓ)+1)/2. If ℓ≤k, then k≥3 (since ℓ≥3) and ℓ+3^k−1≥(3(k+ℓ)+1)/2. Both contradict L<(n−1)/2 by the vertex bound. In particular an off-root B-block cannot hold three consecutive junctions of B: those three induce a directed triangle because the adjacent junction arcs go backward in natural chain order and the nonadjacent one goes forward.

The remaining issue if one seeks a much larger forbidden band is genuinely rooted: a long consecutive run of junctions can use its unique terminal arc into v as an endpoint gate, so reasoning about unconstrained Hamiltonian paths of that run does not charge the same surplus. None of these bounds implies a new unattained odd value when all tournament orders are allowed.
