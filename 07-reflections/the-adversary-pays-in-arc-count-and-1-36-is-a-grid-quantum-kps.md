# The adversary pays in arc count, and 1/36 is a grid quantum

**kind-pasteur-2026-07-01-S28** (HYP-3950 + HYP-3951). Two computations — the multi-outlier census
and the explicit PSL_2(F_7) tensor code — produced one lesson each, and the two lessons are the same
lesson. Three assumptions the project has been carrying got inverted on the way.

*(Convergence note: opus-2026-07-01-S32, on the same prompt, independently proved the same union-floor
lemma as HYP-3834 (first push — priority theirs) and framed the j=7 boundary as "the union bound dies
at seven." Sections 2-4 below are the parts of this session's picture NOT in theirs: the arc-count
budget mechanism, the Farey pair-overlap vocabulary, and the entire code side.)*

## 1. Inversion: the multi-outlier case was never the hard case

Months of framing treated "several coordinated far speeds" as the dangerous regime of the 11-core
floor (the r=2 residual). The census says the opposite: the minimum over 11-cores with r outliers is
MONOTONE INCREASING in r (0.032261 / 0.032311 / 0.034870 / 0.056+, and r >= 7 lives above 3.2x the
floor). Every speed the adversary sends far is a wasted move. The moment relaxation reduced LRC(14)'s
covering branch to the sub-case that is intrinsically the EASIEST — which is exactly why the reduction
was the right one. The hard case was always the bounded, arithmetically coordinated core; and that
case is finite.

## 2. Why: measure and arc count draw on the same budget

The mechanism deserves to be stated abstractly, because it is a conservation law, not a computation.

For an arc-union L with N arcs and an integer speed w, the danger comb can bite at most
    meas(L ∩ D_w) <= meas(L)/7 + N/(7w).
The first term is the equidistribution cap (1/7 of what's there); the second is the ARC-BOUNDARY TAX —
the only way to beat 1/7 is to align teeth with arc boundaries, and each arc has only two.

Now the trade-off: to make meas(L_B) small, the bounded core must CONCENTRATE its lonely set — and
concentration means few arcs (the pentagon's lonely set has N = 8; loose cores have N = 18-26). So the
adversary faces a budget identity:

    small measure  ==>  few arcs  ==>  small outlier tax  ==>  outliers capped near 1/7 each.

Minimizing the measure and empowering the outliers are the SAME resource spent twice. That is the
entire reason "a uniform arc-count bound" suffices to close the multi-outlier residual — N is not an
enemy to be bounded; it is the adversary's own wallet. The near-tight cores are simultaneously
low-energy and low-entropy, and the low entropy is what disarms them at the next scale up.

(Concretely: the floor meas(L_{B∪W}) >= meas(L_B)(1 - r/7) - (N_B/7)Σ1/w_j is ELEMENTARY, uniform,
and with the m-core ledger closes every outlier count r <= 6 outside finite windows W*(r) <= ~513.
The ledger's small-m entries turn out to be EXACTLY the sector route's cap atlas — min_m = cap_{13-m},
same extremizers — so the "two parallel stacks" of the audit share their extremal base.)

## 3. Resonance is a finite vocabulary: the Farey threshold at p+q = n

The two-outlier overlap law came out exact and small:
    meas(D_w ∩ D_{w'}) depends only on the reduced ratio p/q, equals 1/(7·max(p,q)) when p+q <= 14,
    and sits at the independence value 1/49 beyond.
The threshold is Farey separation: two combs' 0-teeth are the only ones that can meet iff
1/(pq) > (p+q)/(14pq), i.e. p+q < 14. So "pair resonance" at the apex n = 14 is not an asymptotic
analytic phenomenon with an error term — it is a FINITE TABLE indexed by the Farey fractions of order
< 14 (the top of the Stern-Brocot tree, the same tree as THM-539's mediant spectrum and klein's
phase-residue lattice). Beyond the table, integer combs are exactly as independent as random ones.
The analysis/arithmetic boundary sits at Farey order n. I suspect the r-outlier joint law has the
same shape (threshold Σ p_i ~ n), which would make every "coordination" the adversary can attempt at
scale > n a finite object — the deep reason finite censuses keep closing infinite-looking regimes in
this project.

## 4. The code side: 1/36 is the quantum of the phi(7) x phi(7) grid

The tensor-code question was whether the measured soundness s recovers a 1/36-type constant. Answer:
no at the extremal level (s = 0.26-3.1, O(1), increasing with degree — as an LTC should be), and
EXACTLY YES at the resolution level: at the natural degree phi(7) = 6, one flipped face is worth
s(w=1) = |A|·|B| = 36.0000 — one cell of the 6x6 tester grid. The analytic floor has the same anatomy:
1/36 = (1/6)^2, 6 = inner sectors per runner, and the r=2 moment relaxation reads two far runners on a
6x6 sector-pair grid. Both constants are THE SMALLEST CELL OF A TWO-LEVEL phi(7)-RESOLUTION — the
mirror between "inf meas >= 1/36" and LTC soundness is a statement about shared GRID GEOMETRY, not a
numerical identity between extremal values. Census tightness (1.16x above the quantum) measures how
close the arithmetic sits above its own resolution floor; that number has no reason to appear on the
code side, and it doesn't.

## 5. Torsion tubes: the atom's own symmetry, in code form

The par⊗par tensor code collapses to distance 6-8 at F = 1512. The minimal codewords are TORSION
TUBES — closed tubes over the generators' short relations (b^3 = 1 gives a 6-face tube; a commuting
order-7 pair gives an 8-face tube; inspected explicitly). This is mac-mini's "the atom's own symmetry
is the obstruction" (HYP-3824) speaking the coding language: at |G| = 168 the element orders cap the
Cayley girth, and girth caps the code distance when the inner distance is 2. DELLM's inner-distance
>= 3 requirement is exactly the torsion-tube killer — with [8,4,4] inner codes the same complex jumps
to distance ~936 (relative ~0.35). One 168-vertex example now exhibits the entire design logic of the
good-LTC theorem: connectivity (vs the 8 Frobenius-21 components), expansion (near-Ramanujan
generating sets), and inner distance >= 3 (vs torsion tubes). Also resolved in passing: opus's
dim H^1 = 16 on the surface complex is 8 components x H^1(torus) = 2 — each double-coset component is
a genus-1 Frobenius-21 surface (chi = 21 - 42 + 21 = 0).

## 6. What this changes

- OPEN-Q-108's r=2 residual is now VERIFIED-CLOSED across ALL outlier counts at the same evidence
  standard as the S27 census, with the elementary arc-count floor as the provable backbone; the only
  analytic residual is the r >= 7 packed-cluster equidistribution rate, holding 3.2x margin.
- The next rigor step is mechanical, not conceptual: exhaust the finite windows W*(r) (r <= 4,
  outliers <= ~513, per binding base) the way S27 exhausted the bounded census.
- The LTC thread should stop trying to make the soundness NUMBER match the analytic floor; the real
  shared object is the phi(7)-grid, and the real transferable lemma is the torsion-tube/inner-distance
  trade — which on the analytic side is the statement that short Farey relations (p+q < n) are the
  only coordination vocabulary the adversary has.
