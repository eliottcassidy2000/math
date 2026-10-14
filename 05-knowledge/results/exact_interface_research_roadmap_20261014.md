# Preserve the extension problem: a Collatz–colouring–Jacobian roadmap

Date: 2026-10-14. **PROVED** elementary interface and module lemmas; **OPEN / SPECULATIVE** two finite-word/finite-graph hypotheses. This is a research roadmap, not an equivalence between the three problems or a claim of novelty for CRT, constraint composition, or Brieskorn modules.

## 1. What the new results actually change

The common lesson is precise: retain the datum whose extension is being asked for.

- Collatz needs one fixed positive integer source. A carry count or total cylinder mass has already integrated that source away.
- Unit-distance colouring needs one compatible assignment on the entire seam. Separate four-colourability and independently chosen palettes have integrated that assignment away.
- A polynomial Jacobian mate needs the marked class of the constant `1` to vanish in a Hamiltonian cokernel. Critical-point counts alone have discarded other obstructions to that equation.

The new [polynomial carry-fibre theorem](collatz_polynomial_carry_fibres_20261014.md) closes the entropy shortcut: word and carry growth both have rate `beta`, and true native mass and minimum-cost representative mass both have rate `beta/3`. Fibres lose only `O(log L)` bits. No new weighted-entropy gap is proposed here.

The [common-centre seam theorem](independent_hinge_seam_four_colouring_20261014.md) closes independent-hinge growth as a route above four colours. The [affine test](affine_hinge_overlap_relation_20261014.md) produces real cross-edge conflicts, but all eight proper overlap patterns of `Y_1` still extend. Breaking a convenient formula is not yet a colouring obstruction.

The [critical-free toy clutch](toy_clutch_exception_20261014.md) has no finite critical point but still has no polynomial mate: the special-fibre equation requires a logarithmic primitive. This opens a replacement certificate, not a reason to repeat the failed critical-point extrapolation.

## 2. Collatz: keep the original source in the state

For a native valuation word `w` of length `L`, retain

    (A_w,B_w,r_w,U_w),
    F_w(n)=(3^L n+B_w)/2^A_w,
    n=r_w mod 2^(A_w+1),  n>=3 odd,  n<=U_w.

Here `U_w` is the minimum exact source threshold from all contracting prefixes, or infinity. This is the [source-fibre compiler](collatz_survivor_carry_frontier_20261014.md), not just its ternary projection.

The exact compositional state is `(n,x)`, original source and current value. A step is

    (n,x) -> (n,U_+(x)), retained only if U_+(x)>=n.

In particular, after composing a suffix, compare its states with the original `n`, not with the current endpoint `F_w(n)`. The latter would change the question.

An exact source-generating object is

    Z_(L,r,A)(T)=sum T^n,

where the sum runs over positive odd sources whose native length-L prefix has cost `A`, ternary carry `r`, and no descent below its original source. Each word contributes an explicit arithmetic-progression series, truncated at `U_w` when finite. Native cylinders of equal length are disjoint. Thus the coefficient of `T^n` retains the actual source; substituting a density or summing coefficients generally does not.

**Source-coherence lemma.** Let `m_L` be the least odd `n>=3` surviving `L` steps. Then `m_L` is nondecreasing, and `m_L -> infinity` is equivalent to eventual strict descent for every positive odd `n>1`. Proof: the finite-source survivor sets are nested; a bounded sequence of minima leaves a source in their finite nested intersection. Conversely, a fixed infinite survivor bounds every minimum. This is a precise target, not a new proof of descent. The all-ones Mersenne family remains the mandatory hostile: its finite witnesses grow with `L`, while its coherent 2-adic source is `-1`.

### Hypothesis C: cost–carry univalence

Among slope-admissible words of a fixed length `L` and fixed total cost `A`, the carry `B 2^(-A) mod3^L` determines the word.

This is now [proved through `L=21`](collatz_length21_shift_certificate_20261014.md): the mod-4 range theorem handles lengths through 20, and the compact shifted-numerator certificate below handles 21. Its all-length form remains a falsifiable hypothesis. A counterexample beyond that range would be two different authenticated words with the same triple `(L,A,r)`. Univalence makes every cost coefficient of the weighted carry polynomial `0` or `1` and allows the exact native source address to be reconstructed from `(L,A,r)`. It improves source recovery, not the exponential rate, and does not itself establish positive-integer descent.

**New range lemma.** Put `b_i=floor(log_2(3^i))`. At fixed admissible `(L,A)`, the exact carry numerator satisfies the attained bounds

    B_min=3^L-2^L,
    B_max=sum_(i=0)^(L-1) 3^(L-1-i) 2^min(b_i,A-L+i).

For the upper bound, every prefix cost is at most both `b_i` and `A-(L-i)`. Their minimum has positive successive increments and terminal cost `A`, so it is itself realized by an admissible word. The lower bound is attained by `L-1` initial valuations equal to one and the remaining cost in the last letter.

Two words with the same `(L,A,r)` have numerator difference divisible by `3^L`; both numerators are odd, so a nonzero difference has magnitude at least `2*3^L`. Therefore

    B_max-B_min < 2*3^L

proves univalence for that cost. This parity-only range test proves it for every cost through `L=12`: the maximum width is `1,041,348 < 1,062,882`. At `L=13` the maximum width is `3,644,236 > 3,188,646`, so this weaker test first stops there.

**Mod-4 refinement and theorem through 20.** Every slope-admissible word starts with `a_1=1`, so for `L>=2` its numerator satisfies `B=3^(L-1)+2 (mod4)`. A nonzero same-cost carry difference must therefore have magnitude at least `4*3^L`. If `D_L` is the largest width across costs, then `D_(L+1)=3D_L+2^b_L-2^L`, so `D_L/3^L` is nondecreasing. The single exact inequality `D_20=13,805,179,460 < 4*3^20=13,947,137,604` proves every length through 20. The independent frontier checks all `13,472,296` words at length 20. At length 21 this refined range test alone fails; the following certificate resolves that next case.

**Length-21 certificate.** For a candidate difference `B-B'=k*3^L`, a common prefix of cost `S` followed by first unequal letters `a,b` forces `v_2(k)=S+min(a,b)`. At length 21 only `k=+4` or `-4` is possible, with the high word beginning `(1,2,...)` and low word `(1,1,...)`. Exact completion bounds leave 11 high numerators. Ten shifted low candidates fail the positive affine decoder; the remaining one violates its fifth prefix budget. The explicit 55-subtree pruning certificate covers all 2,427,209 relevant length-20 prefixes while visiting only 168 states.

**All-length normalization and next target.** Since `B` omits the final valuation, every admissible length-`L-1` prefix extends to common total cost `b_L`. A same-cost carry collision at any cost exists iff `p -> 3B_p+2^A_p (mod3^L)` is noninjective on those prefixes. This removes cost stratification but does not prove injectivity. Continue with bounded numerator shifts and exact completion envelopes beyond 21; retain native decoding and every prefix constraint. Terminal slope alone already has a length-7 counterexample.

A separate temptation needs attribution: identifying all no-descent words with slope-admissible words leads to the classical coefficient-stopping-time question. [Garner (1981), §4](https://www.ams.org/journals/proc/1981-082-01/S0002-9939-1981-0603593-2/S0002-9939-1981-0603593-2.pdf) explicitly conjectures equality of coefficient and actual stopping times. This should not be renamed as a discovery of the current compiler.

## 3. Hadwiger–Nelson: the full seam is larger than the overlap

For a finite graph `G` and labelled boundary `S`, retain the Boolean relation

    R_G(S)={colourings of S that extend to a proper four-colouring of G}.

For a complete unit-distance union of two exact point sets, first identify equal points and compute every unit edge. Let `I` be the common vertices. Let `S_A` and `S_B` contain `I` and all endpoints, on their respective sides, of new cross edges between exclusive vertices.

**Exact seam-composition lemma.** The union is four-colourable iff there are assignments in `R_A(S_A)` and `R_B(S_B)` that agree on `I` and give different colours at every new cross edge. This is an existential join of two full relations with the cross-edge constraints. The proof is restriction and gluing.

Retaining only `R_A(I)` and `R_B(I)` is not sufficient when new unit edges join vertices outside `I`. Full seam relations, or equivalent exact constraint factors, are the appropriate interface. Independent palette permutations are not substitutes for shared labelled colours.

The new [type-M toolkit](../../04-computation/experiments/hn_records_20261014/README.md) already provides exact multiquadratic arithmetic, selector-based colour queries, mono-set enumeration, and core shrinking. Its 1,087-vertex Parts reproduction is a working baseline. Its 599-vertex intermediate is not a record. The reported UNSAT verdicts are solver verdicts, not DRAT-checked certificates. Reuse this machinery rather than duplicate a coordinate or SAT layer.

### Hypothesis H: a two-tag seam certificate

In the reproduced union `L727 union rho*S361`, with `rho=(7+i sqrt(15))/8`, tag a cross edge between unrotated points `u` and `v` by the ordered pair

    (|u|^2, |v|^2).

Retain every internal edge. Hypothesize that retaining cross edges with tag `(4,4)` together with **one** additional occurring tag `(s,t)` already gives a graph that is not four-colourable.

This is a bounded, completely falsifiable statement about one exact graph. It is motivated by the observed inadequacy of radius-two forcing alone and the successful combination of a rigid piece with a forcing piece. It is not established by those observations.

Test the `(4,4)` control, then each additional occurring tag with incremental selectors. For every SAT case, retain an explicit colour vector. For an UNSAT candidate, retain the full exact edge subset and obtain a checkable proof before promoting it beyond a solver verdict; then shrink the interior while preserving the seam obstruction. If every two-tag case is SAT, the hypothesis is refuted and the exact minimum number of interacting tags becomes the next measurement. Either outcome locates the needed interaction more precisely than adding another independent hinge. A certificate for this known five-chromatic graph would not increase the lower bound for the plane.

## 4. Jacobian: a marked Brieskorn class and generic-plus-jet detection

Let `k` be algebraically closed of characteristic zero, `A=k[x,z]`, and `P` nonconstant. Define

    D_P=P_x partial_z-P_z partial_x,
    M_P=A/D_P(A),  with t acting by multiplication by P.

This is a `k[t]`-module since `D_P(PH)=P D_P(H)`. It is a vector-space/module quotient, not a quotient by an ideal. A polynomial mate with nonzero constant Jacobian exists iff the marked class `[1]` vanishes; rescale the mate to make its Jacobian one.

This is the classical two-variable Brieskorn-module object under `f -> f dx wedge dz`, because `dP wedge dH=D_P(H) dx wedge dz`. See [Dimca–Saito (1999)](https://arxiv.org/abs/math/9906129). Finite generation must not be assumed: [Shabbir's homogeneous analysis](https://arxiv.org/abs/0705.0741) explicitly allows non-finitely generated torsion. The marked-class interpretation is the useful interface here, not a novelty claim about the module.

For every fibre jet, right exactness gives

    M_P/(t-c)^r M_P = coker(D_P on A/(P-c)^r).

A nonzero image of `[1]` in any such cokernel obstructs a mate. Reduced fibres alone may miss an element: in the model module `k[t]/t^2`, the nonzero marked element `t` vanishes in every reduced fibre, but is detected by the second jet. This is a module example, not a claimed polynomial counterexample.

**Generic-plus-bounded-jets lemma.** Let `M` be any `k[t]`-module whose entire torsion submodule is annihilated by a nonzero polynomial `g=product(t-c)^e_c`. An element `m` vanishes iff its image in `M tensor k(t)` and its images in every `M/(t-c)^e_c M` vanish.

Proof: generic vanishing makes `m` torsion. By CRT the jet vanishings imply `m in gM`, say `m=gz`. Since `m` is torsion, `z` is torsion too. The annihilator hypothesis gives `gz=0`. The converse is immediate. Finite generation is a sufficient condition for a torsion annihilator, but is stronger than needed. Without any such control, the divisible torsion module `k(t)/k[t]` shows why generic and all finite-jet tests can both be blind.

### The new time-form theorem already supplies generic detection

On a smooth generic-fibre component, the time form can be written

    omega=dz/P_x=-dx/P_z,  omega(D_P)=1.

A mate would give `dH=omega`. A nonzero logarithmic residue or period obstructs this, as in the toy clutch. More strongly, if the nonzero time form extends holomorphically to the complete smooth curve, it cannot be the differential of a rational function: a pole of that function would produce a pole of its differential, while a pole-free function on a complete curve is constant.

The accepted [all-degree accessory time-form theorem](all_degree_accessory_time_form_obstruction_20261014.md) now supplies exactly this certificate for `P_B=(Vz^2+Bz)^2+Az+x`, with `deg V=16`, `deg A=8`, the stated inherited multiplicities of the roots of `V`, and simple zeros of `A` at those roots. Every polynomial `B` is excluded as a Keller component: a common zero of `B,V` gives an explicit critical point; otherwise the complete puncture ledger makes the generic time form nonzero and holomorphic on every completed fibre component. This is **PROVED**, published in [commit cb61f25](https://github.com/eliottcassidy2000/math/commit/cb61f25f1ca39807b13dbacad429706204a46f37). It detects `[1]` generically and bypasses the stronger finite-critical-point hypothesis H-J, which remains open. No jet test is needed to re-exclude this completed architecture. The generic-plus-jets lemma is a tool for a new architecture in which the generic obstruction vanishes or needs a finer exceptional-fibre analysis.

The separate [balanced-degree extension](balanced_degree_clutch_time_form_obstruction_20261014.md) goes further: for `deg V=2a`, `deg A=a>=1`, every boundary-unit clutch `gcd(B,V)=1` has no mate, with no local-root assumptions; if `A` has a simple zero at every root of `V`, all polynomial clutches are excluded. Its proof is published in [commit 0ab22ce](https://github.com/eliottcassidy2000/math/commit/0ab22ce269b2b34887c17f22274a12a459070e63). Thus merely changing the balanced degree is not a fresh open sector.

## 5. Research order and stopping rules

1. **Collatz:** investigate cost-resolved source recovery and the exact positive-source coefficients. Stop treating either carry count or weighted mass as a route to a new exponential gap.
2. **HN:** use the full-seam join and test the finite two-tag hypothesis with exact geometry and retained certificates. Stop repeating independent common-centre hinge growth.
3. **JC:** the fixed `(16,8)` architecture and the stated balanced-degree sectors are now excluded. Map the sharp degree and local-order boundary where the time form first acquires a pole, rather than search more clutches inside that closed family. For a new architecture whose generic class vanishes, study a torsion annihilator and the corresponding finite fibre jets. Absence of finite critical points is never evidence for a mate.

The unifying research principle is conditional and operational: a compressed interface is useful only if its composition preserves the specified extension obstruction. The preserved objects above are different; their shared discipline is to keep the source, the seam assignment, or the marked class until the final existence question is decided.

The companion [range checker](../../04-computation/experiments/cost_carry_range_20261014.py) and [output](cost_carry_range_20261014.json) use exact integers. Normal and optimized Python outputs agree. No new large graph search is claimed here; the all-degree clutch theorem is inherited from its cited proof packet.
