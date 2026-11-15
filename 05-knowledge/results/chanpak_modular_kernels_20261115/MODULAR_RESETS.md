# Modular bottlenecks without rational ones

**Proved.** In the exact thirteen-unit ribbon grammar, avoiding rank-one transfer products over the rationals does not remove coherent phases at 3 and 5, and does not restrict the coherent sublanguages to these two odd primes. There are explicit strong periodic and exponentially branching families with rational transfer rank at least two at all depths, but eventual phase 1 at every fixed conductor (p^r), for each (p=3,5,7,11,13). On four of the five periodic examples, the hyperbolic kernel reads a second coordinate that does not collapse modulo (p).

## 1. Physical construction and transfer certificate

Use the thirteen units and initial letter 0 of the [adaptive cone](../chanpak_cone_collision_placement_20261115/SUFFIX_CONES.md), numbered in its order starting from 0. The initial vertex is never dropped. Here are five periods (w), their ordinary count transfer characteristic polynomials (x^6(x^2-tx+d)), and two equal-length macro choices. The choices are (Z_0=w^k) and (Z_1) obtained by replacing the initial units with the indicated tuple.

| (p) | (w) | (t) | (d) | (k) | replacement | vertices per macro |
|---:|:---|---:|---:|---:|:---|---:|
| 3 | (7,9) | 39 | −30 | 4 | (4,5) | 24 |
| 5 | (1,5) | 25 | 10 | 3 | (0) | 15 |
| 7 | (0,9,3) | 161 | 42 | 3 | (2) | 21 |
| 11 | (12,9,8) | 473 | 550 | 3 | (4) | 27 |
| 13 | (0,4,3) | 117 | −52 | 3 | (2) | 21 |

The replay builds the 119 nonnegative raw fragment states and checks the exact 8-coordinate intertwiner. Its reduced matrices have coefficients in (ℤ[1/2]); there is no division by an odd prime in any check or product. For each row, let (A) be the period matrix, (D=A^k) and (B) the matrix for (Z_1). The replay establishes all of the following with exact arithmetic:

* (p divides t), (v_p(d)=1), and (D∈p Mat₈(ℤ[1/2])). The ordinary rank of (D) is two. Its image (E) is precisely the nonzero primary factor of (A); on (E) the period is invertible over (ℚ).
* (B) also maps into (E) and is invertible on (E). Thus after the first macro, every product of macros has rank two. Since products cannot increase rank on either side, no factor of a concatenation can have rank zero or one. The same argument applies to every finite power of the original period, using the nonzero quadratic factor.
* Each of (0Z_0,0Z_1,0ww) is strong as a literal tournament, as verified by two directed reachability tests. Copies of these spanning intervals overlap in a vertex; the induced edge rule depends only on differences and on the later letter. Therefore every nonempty macro concatenation and every periodic word (0w^m) with (m≥2) are strong. The subsets needed for ordinary direct path-count checks are recounted by literal subset DP, separately from transfer.

Both macro choices have the same length but their first Hamiltonian counts differ in every row. Distinct symbol words produce (2^m) labelled ribbon words at macro depth (m); no claim of (2^m) different counts or isomorphism classes is required.

## 2. What the twist does, quantitatively

For a word (ω) of (m) independent uniformly chosen macros, write (J) for the number of (Z_0)'s. On the entire ambient (ℤ[1/2])-lattice its transfer is divisible by (p^J). In particular its integer (H(ω)) is divisible by (p^J). For every integer (a) and fixed (r ≥ 1),
$$\left|1-\mathbb E_\omega e^{2\pi i aH(\omega)/p^r}\right|\le 2\Pr(J<r)=2^{1-m}\sum_{j=0}^{r-1}\binom{m}{j}.$$
The normalized magnitude tends to **one**, not zero. At (p=5) both macros are (0 mod 5), so the mean is exactly one once (m ≥ r). The binomial bound also tends to zero exponentially in (m) for (r ≤ (1/2−ε)m). Since all the counts are odd, conductors (2p^r) give the corresponding constant sign. This extends the counterexample to a length-synchronized branching set, not just a single word at each length.

For the pure periods, the two nonzero eigenvalues both have (p)-adic valuation (1/2): the Newton polygon of (x^2-tx+d) is a single slope (-1/2). This is not rational rank loss. It is degeneration at a finite place, detectable from a characteristic polynomial of rank two. In particular the same phenomenon occurs at 7, 11 and 13 without invoking either Moon–Busch multiplier. This does **not** establish a global Fourier lower bound for the whole thirteen-unit language: the constructed words might have negligible weight under a specified full or tilted ensemble.

The underlying calculation is a caution for endpoint twists. If a walk operator is twisted only by the increment (H(y')-H(y)), the twist is a diagonal conjugation by (e^{2π i aH(y)/q}), not a new operator-norm gap. The meaningful estimate is decay of the endpoint character against the actual starting law and length conditioning (or decay on a specified centered space), with the absorbing residue strata separated. For pointwise coverage one needs a minor-arc integral and positive main term at the target, not a gap at just finitely many fixed rational phases.

## 3. The hyperbolic reader separates a number-theoretic coordinate

Let (T_m=0w^m), (h_m=H(T_m)), and (q_m=Π₂(T_m)), where (Π₂) counts two ordered disjoint possibly empty paths covering all the vertices. Substitute (T_m) into the first [THM-5133 hyperbolic kernel](../../../01-canon/theorems/THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md), fixing the other port to a singleton. The new ordinary count is
$$N_m=6q_m-h_m.$$
Thus (h_m) mod (p) alone does not decide which primes can divide the result.

To see this at all orders, allow up to **two** old heads in the fragment transfer and give each state with one or two final fragments weight two. One final fragment has two placements around an empty slot, and two distinct fragments can be ordered in two ways. There are 206 reachable raw states. Every completed two-cover restricts to pieces with at most two old heads (none of them can ever acquire an incoming arc from a future vertex), so the transfer computes (Π₂) exactly. Modulo the indicated prime the observed (Π₂) recurrence polynomials of the periodic sequences are:

| (p) | recurrence polynomial mod (p) | initial coordinates |
|---:|:---|:---|
| 3 | (x) | 2 |
| 5 | (x^3+3x^2+3x+4) | 2,4,2 |
| 7 | (x^3+5x^2+2x+1) | 2,2,4 |
| 11 | (x^3+x^2+3x+8) | 2,5,8 |
| 13 | (x^3+8x^2+7x+12) | 2,8,12 |

The verifier does not assume that a fit extrapolates: it verifies 206 consecutive zero recurrence residuals from the first possible residual. The residual itself is the readout of a vector propagated by the 206-dimensional **same** period transfer. Cayley–Hamilton makes the check sufficient at every depth. At (p ≥ 5) the companion mapping is invertible, so the initial triple recurs periodically and (q_m ≡ 2 mod  p) on an infinite progression of depths. For sufficiently large such (m), (h_m ≡ 0 mod  p) but (N_m ≡ 12≠0 mod  p). For example the mod-5 profile has period 8. This does not imply prime outputs.

At 3 the situation is different: [THM-4756](../../../01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md) gives (Π₂(T) ≡ (-1)^{|T|}H(T) mod 3) for every tournament. Hence if a module has (3 divides H), any unmarked quotient reading it through only these two coordinates retains divisibility by 3. In fact every ordered nonempty cover with at least three runs has a factor (3) from the free ordering of the runs. No unmarked finite substitution of that module can escape this factor. Selective contacts are a different operation and not covered by the obstruction. This is an algebraic reason why a single prime-blind correction cannot be assumed.

## Replay and provenance

Run `python verify_modular_resets.py` with NumPy and SymPy. The rank and matrix certificates redo the older adaptive cone's transfer reconstruction (not its cone inequalities). [Wave 28](../chanpak_waves_20261113/README.md#wave-28) supplies rational rank-one factorization, and [Wave 33](../chanpak_waves_20261113/README.md#wave-33-a-survey-of-tangential-prior-work-and-what-it-points-to) posed the restricted-twist and indefinite-reader tests. The new content is rank-two modular collapse, the coherent equal-length branch languages, and separation using the independently certified two-cover profiles. Next: characterize which reset blocks survive the *actual* tilt and test port-specific modular mixing after separating configurations by the valuation of the full visible profile. Chan–Pak remains open.
