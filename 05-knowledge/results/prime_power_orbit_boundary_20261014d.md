# Prime powers, signed orbit fibres, and a lower bound that survives refinement

**PROVED scoped:** elementary primitive-orbit inversion, the countable-tail countermodel, and the fixed-source first-hit refinement criterion. **FINITE-EXACT:** the declared functional-map census and four supplied signed cycles. **CITED:** the classical prime-race and Mertens results. The giant Mersenne obligation and universal Collatz coverage remain **OPEN**.

Artifacts: [checker](../../04-computation/experiments/prime_power_orbit_boundary_20261014d.py), [exact output](prime_power_orbit_boundary_20261014d.json), [prime-power inversion](prime_power_bias_20261014d.md), [signed inverse fibres](signed_collatz_pair_lifts_20261014d.md), and [power-map defects](residue_power_defects_20261014d.md).

## 1. Inheritance and the working board

The closest mechanisms are [THM-1926, tournament zeta over the strong core](../../01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md), the explicit signed Collatz cycle zeta in [the Lucas/discrepancy packet, section 5](collatz_lucas_monotile_discrepancy_20260930.md), and the [discounted Green weight](collatz_green_weight_20261005.md). These already separate periodic information from ROOT support. We recover that separation rather than claim a new cycle-zeta construction.

The canonical hostile is an escaping ray invisible to every periodic-point count. The corrected near miss is the [mod-six sandwich calculation](collatz_mod6_20260917_sandwich_bias.md): a square diagonal does not determine the sign of its retained off-diagonal residual. Its counterexample at 1453 is unrelated to the present prime count 1473. The least-used sidecar is a lower bound at an *unchanged source* under refinement.

The live board is **primitive object / repetition index / lost sign / native path / first-hit clock / nonvanishing weight**. The anchor is source-specific Collatz coverage; the niche is invertible prime-power bookkeeping; the wildcard is a countable limit with identical cycle data but different basin coverage. No tournament is imposed on objects without an intrinsic pairwise orientation.

The incoming [THM-4640 repair-algebra packet](../../01-canon/theorems/THM-4640-the-repair-algebra-at-minus-one-collisions-are-increasing-two-three-representations-ladders-and-gersonides-carries.md) suggested a concrete new hostile: squaring opposite anchor values creates false collisions. The companion power-defect packet classifies and checks those collisions while retaining the original source.

## 2. The square correction is exact; a permanent phase transition is not

Let `Delta=pi(x;4,1)-pi(x;4,3)` and let `J_a` count each prime power `p^k<=x`, in residue `a`, with weight `1/k`. The [prime-power packet](prime_power_bias_20261014d.md) proves

\[
J_1(x)-J_3(x)=\Delta(x)+\tfrac12\pi_{\rm odd}(\sqrt x)+R(x),
\qquad R(x)=o(\sqrt x/\log x).
\]

The square term is asymptotic to `sqrt(x)/log(x)`. The weighted race `J_1-J_3` is still present, so this does not assert a pointwise prime deficit of that size. The full residue-vector transform is invertible across lower cutoffs even though squaring collapses residues; the companion retains the exponent, residue pushforward and cutoff.

Exactly at **26861**, the two prime counts are **1473 and 1472**. At the next prime, **26863**, they tie; at **26879**, class 3 leads again. These are complete finite calculations, including independent trial division. Thus the first reversal is not a boundary after which the other class permanently leads.

[Granville and Martin, *Prime Number Races*](https://arxiv.org/abs/math/0408319), discuss both this race and the `pi(x)-Li(x)` race associated with Skewes. The relevant classical result is recurring sign change, not a single irreversible passage from a small-scale law to a different large-scale law. The useful research distinction is between a systematic offset and oscillations capable of overcoming it.

The Mertens function is a different observable, `M(x)=sum_(n<=x) mu(n)`. [Odlyzko and te Riele's original report](https://ir.cwi.nl/pub/6447/6447D.pdf) disproves its proposed square-root bound through zeta-zero computations and lattice reduction, without giving an explicit integer counterexample. Its positive smoothing kernel converts a finite oscillatory calculation into a rigorous limsup bound. This is an existence mechanism, not a rule locating every witness.

The analytic distinction is already visible locally: at a simple nonzero zeta zero `rho`, the Perron integrands for `psi` and `M` have residues respectively

\[
-\frac{x^\rho}{\rho},\qquad
\frac{x^\rho}{\rho\zeta'(\rho)}.
\]

These follow by expanding `-zeta'/zeta` and `1/zeta` at the zero. They share zero locations but have different weights. This local residue computation does not invoke convergence of an infinite zero sum, RH, or simplicity of all zeros. The mod-four race additionally involves the corresponding Dirichlet character. No prime-square correction is silently substituted for Mertens oscillations.

## 3. Primitive objects and repetitions: a genuine shared transform

For a map `T` with finitely many points fixed by each iterate, let `c_l` count primitive cycles of length `l`, and `F_m=#Fix(T^m)`. Then

\[
F_m=\sum_{\ell\mid m}\ell c_\ell,
\qquad
c_m=\frac1m\sum_{d\mid m}\mu(d)F_{m/d},
\]

\[
\zeta_T(z)=\exp\left(\sum_{m\ge1}\frac{F_m}{m}z^m\right)
=\prod_{\ell\ge1}(1-z^\ell)^{-c_\ell}.
\tag{1}
\]

The first identity counts points on cycles whose length divides `m`; divisor inversion proves the second; expanding each formal logarithm proves the third. These are formal series identities, so no analytic convergence region is needed here.

The analogous prime expansion is `log zeta(s)=sum_(p,k) p^(-ks)/k` for `Re(s)>1`. Both operations retain a primitive object and its repetition count. This is a structural correspondence, not a bijection from primes to Collatz cycles, nor a transfer of zeta-zero estimates to a Collatz source.

The four **supplied nonzero** cycles of the shortcut map

\[
T_+(n)=\begin{cases}n/2&n\text{ even},\\(3n+1)/2&n\text{ odd}\end{cases}
\]

on the integers have lengths `2,1,3,11`, through `1,-1,-5,-17`. Their catalogue contributes

\[
Z_{\rm known}(z)=\frac1{(1-z)(1-z^2)(1-z^3)(1-z^{11})}.
\tag{2}
\]

Including zero adds another factor `(1-z)^(-1)`. Under the odd-only accelerated clock, the same nonzero cycles have lengths `1,1,2,7`, giving `(1-z)^(-2)(1-z^2)^(-1)(1-z^7)^(-1)`. Changing clocks changes the trace; the number of cycles alone does not retain their lengths.

Equation (2) is the zeta of the supplied catalogue. Equality with the full integer map's zeta would exclude further cycles. It would **still not exclude escaping orbits**. The old Lucas packet explicitly called its zeta identity the *cycle half* of the conjecture; that qualification remains essential.

## 4. Identical cycle data can coexist with opposite coverage verdicts

Here is an exact infinite example, not a Collatz counterexample. Fix

\[
X=\{c\}\cup\{r_j:j\ge0\}.
\]

For each integer `N>=0`, set `T_N(c)=c` and

\[
T_N(r_j)=\begin{cases}
r_{j+1}&j<N,\\c&j=N,\\r_{j-1}&j>N.
\end{cases}
\tag{3}
\]

Every source reaches `c`, with exact first-hit time `|j-N|+1`. The strict potential `|j-N|` decreases until the final edge. Its only periodic point is `c`, so `F_m=1` for every `m` and `zeta_(T_N)=1/(1-z)`.

For each fixed `j`, however, `T_N(r_j)=r_(j+1)` for every `N>j`. The pointwise limit `T_infinity` fixes `c` and shifts every `r_j` to `r_(j+1)`. Every `r_j` escapes, while **every periodic-point count and the entire zeta remain unchanged**.

This strengthens a finite-horizon warning: for any prescribed depth, a completely rooted system can reproduce the escaping ray to that depth. What changes is the terminal's distance from the fixed source. On a finite functional graph, a complete cycle catalogue does imply eventual entry into that catalogue. Finiteness is the missing hypothesis in extending that inference to an infinite graph.

## 5. Exact criterion for a source-specific lower bound under refinement

Let `T_N,T` act on the same discrete set `X`. Assume that for every `x`, `T_N(x)=T(x)` for all sufficiently large `N`. Keep the same terminal set `C`, source `s`, and number `0<z<1`. Define

\[
H_N(s)=\begin{cases}z^{\tau_N(s)}&\tau_N(s)<\infty,\\0&\text{otherwise},\end{cases}
\]

where `tau_N` is the first entrance time into `C`, including time zero. Define `H` using `T`.

**Proved criterion.**

\[
H_N(s)\longrightarrow H(s),\qquad
s\text{ reaches }C\text{ under }T
\iff\liminf_{N\to\infty}H_N(s)>0.
\tag{4}
\]

**Proof.** If the actual first-hit time is `h<infinity`, the finitely many edges before that hit eventually stabilize. All sufficiently late approximants have exactly the same first-hit time `h`, hence the same weight `z^h`. If the actual orbit never hits, then for every fixed `B` its first `B` edges eventually stabilize and avoid `C`. Consequently `0<=H_N(s)<=z^(B+1)` eventually. Letting `B` grow proves convergence to zero. This proves both assertions. Finite `C` is unnecessary; it only needs to be unchanged.

Equivalently, an eventually uniform bound `H_N(s)>=epsilon_s>0` bounds the first-hit clocks by

\[
\tau_N(s)\le\left\lfloor\frac{\log\epsilon_s}{\log z}\right\rfloor.
\]

The constant can depend on `s`; it may not decay with refinement. In (3), every `H_N(r_0)>0`, but `H_N(r_0)=z^(N+1)->0`. This is precisely the distinction between *positivity at each finite stage* and *positivity that survives refinement*.

For a fixed map, the first-hit weight is the unique bounded solution of `H=1` on `C` and `H(x)=z H(Tx)` off `C`: iteration determines a finite-hit value, while boundedness forces zero on an orbit avoiding `C`. At `z=1`, that uniqueness fails and the constant function 1 solves the equations even on an escaping ray. The discounted quantity is a first-hit observable, not an unqualified stationary probability. This is a specialization of the repo's Green-weight mechanism, with the refinement quantifiers made explicit.

## 6. How to transport a positive bound without recycling a certificate

Positive averaging, as in the Mertens smoothing argument, can prove that *some* witness exists. For universal Collatz coverage, the target is the specific supplied source. The missing map must therefore be authenticated.

Suppose a finite native path sends `s` to `m` in `ell` steps with no earlier terminal. Then exactly

\[
H_z(s)=z^\ell H_z(m).
\tag{5}
\]

A collection of such verified paths can be averaged with nonnegative weights summing to 1: every transported summand already equals the same `H_z(s)`. Averaging values at unrelated sources lacks (5), so a positive average supplies no lower bound at `s`. The sign, actual path, terminal type and length are the required extra data.

The [signed-fibre packet](signed_collatz_pair_lifts_20261014d.md) provides that interface. For odd `q` coprime to 3, all positive odd incoming sources across the two maps are

\[
n_a=\frac{2^a q-b_a}{3},\qquad
b_a\in\{-1,+1\},\quad b_a\equiv2^a q\pmod3,\quad a\ge1.
\]

They satisfy `U_(b_a)(n_a)=q`, `b_(a+1)=-b_a`, and `n_(a+1)=2n_a+b_a`. This enumerates the full two-map inverse fibre. A terminal certificate for `q` under one map grounds the sources on that same map; it does not automatically ground the alternating sources under the other map. A numeric collision at `q=5` is the cheapest hostile: plus dynamics reaches 1, minus dynamics enters its `5,7` odd cycle.

There is a constructive positive control. At `q=26861`, the companion supplies and independently replays both ROOT words: 48 odd steps for plus, 45 for minus. Hence every

\[
n_a=\frac{26861\,2^a-(-1)^{a+1}}3,\qquad a\ge1,
\]

reaches 1 under the single autonomous map with sign `(-1)^(a+1)`. Its certified odd rank is 49 for odd `a`, 46 for even `a`. The first sources are 17907 and 35815. The prime-race meaning of 26861 and this finite Collatz verification are separate facts; the number's race role did not prove its ROOT certificates. This is a grounded infinite inverse family, not a coverage theorem for arbitrary inputs.

The companion also closes its exact intersection with the outstanding target. For an odd endpoint `q>1`, its combined inverse fibre contains a positive Mersenne number if and only if `(q+1)/3=2^r` for an integer `r>=1`. The only marked Mersenne ports are then `M_(r+1)` on plus with valuation 1 and `M_(r+2)` on minus with valuation 2. Directly expanding `3(2^E-1)+b` proves this classification, with `q=1` treated separately. Here `(26861+1)/3=8954` is not a power of two. Thus this new grounded fibre contains **no positive Mersenne number**, at any exponent. Its certificate interface transfers; its coverage does not resolve the old branch.

For raw odd branches, `(3n+1)/2` and `(3n-1)/2` are consecutive and lie in residues 2 and 1 modulo 3. Further halving needs its valuation parity. Likewise squaring modulo 3 sends `0->0`, `1->1`, `2->1`; the prime 3 must be excluded before saying that all prime squares land in class 1. The [power-defect packet](residue_power_defects_20261014d.md) makes these corrections executable and shows how dropping an anchor sign falsely licenses a Collatz switch.

## 7. Next proof obligation and reproducibility

The positive target is now explicit: construct a hierarchy of **source-preserving native rules**, with fixed terminal semantics, and prove a nondecaying lower floor at each uncovered source. Formula (4) explains exactly why that would suffice. Producing that floor remains the difficult part; the criterion does not manufacture it.

For the [uncovered giant Mersenne phase](mersenne_clocked_completion_20261014c.md), the previous all-word inverse search already rules out a smaller predecessor within 64 odd steps. The new sign and power transforms do not evade that result. A useful next candidate must supply a longer native path, a different grounded child family, or a first-hit estimate that remains positive without changing the source. Any proposed square/norm compression should retain the sign and clock before its coverage is assessed.

The checker covers every labelled functional map on 1 through 5 vertices (3413 maps), every period through 12, and formal series through degree 12. Independent iterate counts, cycle decomposition, divisor inversion and Euler products agree. It verifies the four supplied cycles under both clocks through series degree 132; it does not search for all integer cycles. Countable-family controls use `N=0..64`, `j=0..96`, and every prefix depth through 128. The all-depth assertions in sections 4 and 5 are proved above, not extrapolated from those controls.

```text
python -B 04-computation/experiments/prime_power_orbit_boundary_20261014d.py
python -B -O 04-computation/experiments/prime_power_orbit_boundary_20261014d.py
```

Both runs give **357651 explicit checks** and the same LF-normalized JSON SHA256:
`e2bb76647e9966561e3547f38f4c5ff1c57cfe81e772921289da31520f25b6e7`.

Independent proof/API review and optimized replay passed for this packet and all three companions. The four retained outputs total **1,053,815 explicit checks**; repeated elementary validations are included in that count. The scope of each universal statement comes from its proof, not the number of finite checks. No new theorem identifier or global coverage claim is promoted here.
