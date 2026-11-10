# Scale-circle gluing: three narrow rails or one wide rail

Session247, 10 November 2026.

**PROVED-IN-NOTE:** the exact arithmetic scale-circle theorem, phase-count obstruction, tripling specialization, and affine seam-debt test below. **CONDITIONAL APPLICATION:** the reported finite physical deletion windows are taken as inputs, without auditing them. **OPEN:** any window-preserving all-order rail of tournaments. This note does not prove Chan–Pak, and it does not identify a numerical interval with an operation-ready tournament state.

## 1. Why this is a different metagraph question

The old half-tiling square-complex reflection suggests studying a quotient and its boundary obstruction, rather than merely adding more local exchanges. The pants-boundary work insists that cuffs require actual marked gluing data. The Gordian continuation theorem distinguishes equality of an observable from equality under future operations. Combining these ideas yields a much smaller global target:

> Do not transport each missing target through the metagraph. Transport a finite family of actual window-producing states around one multiplicative scale circle.

The quotient here is exact: positive real magnitudes modulo multiplication by an odd integer p. Its logarithmic coordinate is a circle of circumference log p. A physical interval window becomes an arc on that circle. The point of topology is only to organize overlaps and winding. Continuous points of the circle are not tournament witnesses.

The specific new signal, reported by the parent from fresh THM-5127 work, is a finite bank of actual deletion windows

    [3h,5h) for odd 237 <= h <= 1023,
    [3h,9h) for odd 365 <= h <= 1023.

IMPORTANT QUANTIFIER: the inherited window is a UNION over all realizers R with H(R)=h, not the response image of one fixed R. Different targets may require different deletion cores. Here and below a window means *every odd integer* in the stated half-open interval is the Hamiltonian-path count of a literal tournament furnished by a certified BANK of realizers at h. A state in this note means that bank or its exact target-witness relation, never one magic representative. No single-realizer window is assumed. The finite input is cited at frozen source commit `7d5cf47d10b1`, in `01-canon/theorems/THM-5127-heat-flow-reading-of-the-w-polynomial-tilted-ensemble-mean-and-the-mod-four-quadratic-form-law-of-a-vertex-cube.md`. It is inherited without re-verification. A finite bank does not establish its infinite continuation.

## 2. The scale-circle theorem

Let p>1 be an odd integer, and let alpha<beta be positive odd integers. Choose positive odd phase seeds

    a_1 < a_2 < ... < a_q < p a_1.

Suppose for every i=1,...,q and t>=0 there is an actual certified bank state S_(i,t), whose members all have base count p^t a_i, and whose union of witness constructors realizes every odd integer in

    I_(i,t) = [alpha p^t a_i, beta p^t a_i).

No equality, adjacency, isomorphism, or common deletion core between two states is required. The certificate must, however, concern actual integer counts, not just average, congruence support, or a real relaxation.

**Theorem (exact phase gluing).** Their union contains every odd integer at least alpha a_1 if and only if

    alpha a_(i+1) <= beta a_i           (1 <= i < q),
    alpha p a_1 <= beta a_q.                         (G)

### Proof

The lower endpoints occur in the displayed cyclic order: first all q phases at t=0, then all q at t=1, and so on. Their upper endpoints occur in the same order, because each upper endpoint is beta/alpha times its lower endpoint. Each endpoint is odd. The last odd integer in [L,U) is U-2. Thus the next interval [L',U') leaves no missing odd integer exactly when L'<=U. Applying that test between neighboring phases and at the wrap gives (G), including equality. Multiplying by p^t preserves the inequalities. This proves sufficiency.

If one inequality fails, its positive difference is an even integer. The odd integer U at the previous interval's upper endpoint is absent from that interval and lies below the next lower endpoint. Earlier intervals end no later, and later intervals begin no earlier, so U is missing. This proves necessity. QED.

The circle interpretation is now literal: logarithms of the seeds modulo log p mark q starting points, and each forward arc has length log(beta/alpha). Condition (G) says those arcs cover the circle. The witness choice is an arithmetic interval lookup, not an appeal to continuous intermediate values.

### Optimal phase-count obstruction

Multiplying all q inequalities gives

    p <= (beta/alpha)^q.                              (P)

Hence q>=ceil(log(p)/log(beta/alpha)). If the right side is an integer and equality is needed, rational/odd seed constraints still have to be checked separately. In the strict case, sufficiently large odd seeds approximating equally spaced logarithmic phases exist. This last existence statement concerns *numerical seeds*, not physically certified states.

The measure-of-a-circle lower bound and the multiplied arithmetic inequalities are the same obstruction. More metagraph connectivity cannot evade it for the specified rail model.

## 3. Three narrow rails: an explicit Chan–Pak reduction

Set p=3, alpha=3, beta=5 and use

    (a_1,a_2,a_3) = (237,343,495).

The exact seam inequalities are

    3*343 = 1029 <= 1185 = 5*237,
    3*495 = 1485 <= 1715 = 5*343,
    9*237 = 2133 <= 2475 = 5*495.

Consequently:

**Conditional narrow-rail corollary.** If the three physical window-producing banks at base counts 237,343,495 can each be continued to states at counts 3^t times their base count, retaining a full [3h,5h) window at every t, then every odd integer >=711 is a tournament Hamiltonian-path count.

Thus, after a separate finite spectrum certificate below 711 with the expected exceptions 7 and 21, these three rail-transport theorems would prove the full Chan–Pak spectrum statement.

All three initial seeds belong to the *reported* finite bank. Only their continuation is proposed. There is no claim that all tournaments with these counts have the property or that arbitrary states with equal H are interchangeable.

Three rails are optimal for this exact factor-three / [3h,5h) model:

    (5/3)^2 = 25/9 < 3,     (5/3)^3 = 125/27 > 3.

In particular two narrow rails cannot cover the scale circle, regardless of their initial phases. A scalar recurrence producing ever-larger h is therefore insufficient unless it also has enough logarithmic phases.

The first generation illustrates the glue:

    [711,1185), [1029,1715), [1485,2475),
    then the next generation begins at 2133.

Every interval here refers only to its odd integers.

## 4. One wide rail and the sharp critical seam

Set p=3, alpha=3, beta=9. One phase suffices, exactly at equality beta/alpha=p.

**Conditional wide-rail corollary.** If one actual bank with h_0=365 has window [3h_0,9h_0) and can be continued indefinitely by exact count tripling while retaining that window property, then every odd integer >=1095 is realized.

Indeed the intervals

    [3^(t+1) h_0, 3^(t+2) h_0)

tile the odd tail with no overlap and no hole. More generally, any odd h_0 works. This turns a universal extension theorem into a statement on only one orbit of bank-to-bank constructors. Members may depend on the target.

The price of requiring only one orbit is criticality: the boundary has no slack. A constructor giving count 3h+2 instead of 3h cannot be substituted silently.

## 5. Affine seam debt: what an operation must actually preserve

Suppose an actual transition takes state S of base count h to state S' of base count

    h' = a h + b,

with positive growth, parity preserved, and the windows

    [alpha h+c, beta h+d),
    [alpha' h'+c', beta' h'+d').

Assume endpoints are odd and the second lower endpoint is the next in increasing order. The seam is covered precisely when

    (alpha' a-beta)h + (alpha' b+c'-d) <= 0.           (D)

This is a complete integer test, not asymptotic notation. If alpha'a<beta, it holds beyond the explicit threshold

    h >= (alpha'b+c'-d)/(beta-alpha'a)

when the numerator is positive. If alpha'a=beta, the constant term decides forever. If alpha'a>beta, every sufficiently large seam fails. A finite-state cycle of constructors must satisfy this test at every phase, not just have a favorable average multiplier.

For the wide critical window and h'=3h+2e, the debt is 6e. Positive e leaves exactly 3e omitted odd integers at the seam: 9h,9h+2,...,9h+6e-2. Exact tripling has zero debt. This explains why precise constants in an interface recurrence can decide a global spectrum theorem.

The three narrow phases above have margins 156,230,342 before scaling. Those margins increase by 3^t, making bounded perturbations potentially absorbable. The exact test (D), together with a retained ordering of phases, determines how much perturbation is allowed. This is a reason to prefer a redundant multi-phase construction over a critically tight one-rail construction when physical interfaces introduce bounded corrections.

## 6. A constructive selector and certificate size

Under the rail hypotheses, a target odd N>=alpha a_1 can be selected by comparing N against the q lower endpoints at successive scales. Choose the last lower endpoint L<=N. Condition (G) guarantees N lies below its corresponding upper endpoint U. Then invoke that state's certified window constructor.

For fixed q, finding the scale requires O(log N) elementary multiplications/comparisons, or binary search in the scale if large-integer operations are counted. If the bank transport admits a target-specific lift adding at most r vertices per step, and all initial bank members used have order at most n_0, the selected core has at most n_0+r t vertices, so the rail portion of the witness has order O(log N). Any further cost and vertices required by the window selector remain explicit additional obligations. The theorem does not imply efficient enumeration of a window merely from its existence.

## 7. Why ordinary tripling is not the missing transport

There is a trivial count-tripling operation: take the ordinal sum with a directed triangle. Every Hamiltonian path concatenates the blocks, so H(C3 join T)=3H(T). It creates numerical rails at once.

But applying that operation to every witness in an old window merely multiplies all counts by three. It yields only multiples of three in the scaled interval. It does **not** produce a full odd window for the new base state. This is the exact failure of identifying a scalar rail with a window-producing rail.

The current session's marked-interface fixed-divisor theorem makes this warning more concrete: inherited module divisibility can survive bounded interfaces. The rail state must include an operation-ready deletion/endpoint profile capable of creating the missing residue classes at each new scale. One must prove that profile propagation, not infer it from H'=3H.

Similarly, overlapping numerical intervals need no actual metagraph edge between their endpoint witnesses. Requiring such an edge would impose an unnecessary problem. The only required marked transport is along each rail, where a new window is produced. This separates two tasks that are easily conflated:

1. **Horizontal gluing:** exact numeric seams, completely solved by (G)/(D).
2. **Vertical continuation:** physically transporting the window-producing profile, still open.

## 8. Provenance, inherited inputs, and next experiment

The conceptual sources are deliberately broader than direct spectrum searches:

- `07-reflections/the-half-tiling-is-an-abelian-square-complex-the-metagraph-is-an-anti-LTC-and-cut-cycle-is-alexander-duality-kps.md`: quotient/cell viewpoint. Historical analogies are not assumed as theorems here; no cut-cycle direct-sum claim is used.
- `05-knowledge/results/pants_boundary_holonomy_20261011b.md`, especially the marked gluing interface: boundary counts need a realizable matching. Here the realizable matching is a window-preserving rail constructor.
- `01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md`: equal observables need not have equal continuation profiles.
- `01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md`: order-join product geometry and multiplicative H, used only for the negative control in section 7.
- `01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md`: the inherited local unit-component theorem. The present result does not reprove local interval filling; it compiles finitely many recurring intervals into a global tail.
- `01-canon/theorems/THM-5127-heat-flow-reading-of-the-w-polynomial-tilted-ensemble-mean-and-the-mod-four-quadratic-form-law-of-a-vertex-cube.md`, frozen commit `7d5cf47d10b1`: the finite window bank, inherited without an audit in this lane.

**Next discriminating experiment:** select actual witness banks at 237,343,495 from the finite window data, retain their full owned deletion/endpoint certificates (allowing a different core for each target), and look for a bank-to-bank constructor (or a finite cycle of constructors) that triples the base count while preserving enough profile to certify the next full window. A successful step is evidence about this *profile transport*, not about the already-established arithmetic gluing. If exact tripling is obstructed, search bounded affine perturbations and use the three-phase slack via (D).

The genuine reduction is from all targets/all parents to three designated profile orbits (or one for the wide bank). The remaining obligation is stated explicitly rather than hidden inside an assumption of universal interval coverage.

## 9. The marked triangle interface makes the rail obligation concrete

Use the companion `../marked/PROOF.md`, Theorem 2, with a core R and one root x. Write h=H(R), let O be the vertices beaten by x, and let

    E = number of Hamilton paths of R ending outside O,
    S = number of Hamilton paths of R starting inside O,
    D = sum over nonempty proper L of e(L)s(R\L),
    B = sum over nonempty U of s(U)H(R\U).

Here e and s are the corresponding induced-subset endpoint counts, and H(empty)=1. The old rooted child has F=E+S+D. Append a directed triangle after R in the ordinal core, keeping the same root. The root beats t triangle vertices. Specializing the marked formula gives

    t=3:  3F,
    t=2:  3F+3h-E+3B,
    t=1:  3F+4h-2E+6B,
    t=0:  3F+3h-3E+9B.

For completeness, on the triangle h_2=3, E_2=3-t, S_2=t, and (D_2,A_2) for t=3,2,1,0 are (0,0),(2,3),(2,6),(0,9). Substitution into the exact two-block formula

    3S+hE_2+3D+hD_2+E S_2+B A_2

gives the four displayed formulas.

This is a scoped physical identity inherited from the companion's proved formula, not a claim that it preserves windows. The useful next target is now explicit: retain the JOINT attainable tuples (F,E,B) over the witness bank, and show that the four affine images cover the next window. Knowing only which F occur is insufficient. In particular the t=3 branch is only the divisible-by-three sublattice; t=2 and t=1 introduce residue information through E and h as well as range shifts through B. A bank-level transfer theorem for these joint supports would discharge the vertical transport required above without pretending that one fixed core supplies the whole window.
