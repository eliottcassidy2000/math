# Strong-dephasing gain cycles and the planar-Jacobian wedge-rank sidecar

**Session date:** 2026-08-20
**Status:** internal derivations `PROVED`; symbolic controls `FINITE-EXACT`;
literature novelty `UNAUDITED`; planar Jacobian consequence `NECESSARY-DIAGNOSTIC
ONLY`; `JC(2)` remains `OPEN`.

Exact replay:
[`strong_dephasing_gain_cycle_expansion_codex_20260820.py`](../04-computation/strong_dephasing_gain_cycle_expansion_codex_20260820.py),
with stored output
[`strong_dephasing_gain_cycle_expansion_codex_20260820.out`](../05-knowledge/results/strong_dephasing_gain_cycle_expansion_codex_20260820.out).

## 1. Outcome first

The strong-dephasing bridge has a sharper form than its resistor limit.
Writing `epsilon=1/lambda`, the slow population generator on its invariant
graph is

```text
G_lambda=epsilon G1+epsilon^2 G2+epsilon^3 G3+O(epsilon^4),
G1=Pi K^2 Pi,
G2=Pi K^3 Pi,
G3=Pi K^4 Pi-2(Pi K^2 Pi)^2.                    (1.1)
```

Here `K=-i[H,-]` and `Pi` projects matrices to their diagonal.  The last term
in `G3` is a genuine retardation correction.  Replacing the fast coherence by
its zero-frequency stationary value would miss one copy of this term.

The three coefficients have distinct meanings.

1. `G1=2L_c`, with `c_xy=|H_xy|^2`: this is the resistor walk.
2. `G2` is the sum of triangle `sin(flux)` circulations.  It is real
   skew-symmetric and is the first possible detailed-balance defect.
3. `G3` is real symmetric.  Its phase-sensitive monomials are four-cycle
   `cos(flux)` terms.  A square contributes the local matrix with
   diagonal/adjacent/opposite coefficients `4/-8/12`.

For the Hermitian dilation of a planar-Jacobian response matrix, the graph is
bipartite, so every even-index coefficient `G2,G4,...` vanishes.  The first
information beyond conductance is therefore `G3`.  It recovers from the
dephased graph the pairwise column-wedge energies

```text
||A_col(u) wedge A_col(v)||^2
```

and hence the full scalar invariant `||wedge^2 A||_F^2`, the sum of the squared
`2x2` minors.  This exactly distinguishes the previous equal-conductance
rank-`1`/rank-`2` hostile.  It does **not** recover the oriented Wilson phases,
the target-range condition, higher compounds, or the coefficient Segre
sidecar.

A second boundary is equally important: `G_lambda` is the generator in the
raw population coordinates on the slow invariant graph.  Beyond leading
order it need not be a Markov generator.  A unit square with flux `pi` has a
negative virtual nonedge coefficient `(G3)_02=-16`.  The exact population
marginal remains positive because it includes an initial slip and memory;
one must not reinterpret every higher coefficient as a classical jump rate.

## 2. Setup and exact invariant-graph expansion

Let `V` be finite and let `H=H*` have `H_xx=0`.  Position dephasing with
off-diagonal decay rate exactly `lambda` is

```text
rho_dot=K rho+lambda(Pi rho-rho),
K rho=-i[H,rho].                                         (2.1)
```

Put `Q=1-Pi`, `p=Pi rho`, and `q=Q rho`.  Define

```text
A=Pi K Q,       B=Q K Pi,       C=Q K Q.                 (2.2)
```

Because `Pi K Pi=0`, the block equations are

```text
p_dot=Aq,
q_dot=Bp+Cq-lambda q.                                    (2.3)
```

For sufficiently large `lambda`, the slow spectral subspace is the graph
`q=R_lambda p`.  Invariance gives the matrix Riccati equation

```text
R A R=B+C R-lambda R.                                    (2.4)
```

With `epsilon=lambda^-1`, this is

```text
R=epsilon(B+C R-R A R).                                  (2.5)
```

At `(R,epsilon)=(0,0)`, the derivative of the left side minus the right side
with respect to `R` is the identity.  The finite-dimensional analytic
implicit-function theorem therefore gives a unique small analytic branch.
Writing `R=sum_(n>=1) epsilon^n R_n`, coefficient comparison gives

```text
R1=B,
R2=CB,
R_n=C R_(n-1)-sum_(i+j=n-1) R_i A R_j       (n>=3),      (2.6)
```

and `G_n=A R_n`.  In particular,

```text
G1=AB,
G2=ACB,
G3=AC^2B-ABAB,
G4=AC^3B-2ACBAB-ABACB.                                  (2.7)
```

Inserting `Pi+Q` between commutators gives

```text
AB=Pi K^2 Pi,
ACB=Pi K^3 Pi,
Pi K^4 Pi=AC^2B+ABAB,                                   (2.8)
```

which proves `(1.1)`.

For arbitrary initial coherences, variation of constants gives the exact
memory equation

```text
p_dot(t)=A exp((C-lambda)t)q(0)
        +integral_0^t A exp((C-lambda)(t-s))B p(s) ds.   (2.9)
```

Thus an arbitrary diagonal initial condition is not exactly on the slow
graph.  It acquires an `O(lambda^-1)` initial slip.  This is the precise
reason to distinguish the invariant-graph generator from a globally exact
population Markov semigroup.

On the resistor slow time `t=lambda tau/2`, `(1.1)` becomes

```text
d p/dtau=L_c p+(1/(2lambda))G2 p
                  +(1/(2lambda^2))G3 p+O(lambda^-3).     (2.10)
```

## 3. Entry formulas and conservation laws

Let

```text
c_xy=|H_xy|^2,       d_x=sum_y c_xy.
```

For `x!=y`, direct commutator expansion gives

```text
(G1)_xy=2c_xy,              (G1)_xx=-2d_x,               (3.1)

(G2)_xy=-6 Im(H_xy(H^2)_yx)
        =-6 sum_z Im(H_xy H_yz H_zx),
(G2)_xx=0,                                                     (3.2)

(G3)_xy=-8 Re(H_xy(H^3)_yx)+6|(H^2)_xy|^2
        -2(G1^2)_xy.                                    (3.3)
```

If desired, the phase-blind last term can be written

```text
(G1^2)_xy=4 sum_(z notin {x,y}) c_xz c_zy
           -4(d_x+d_y)c_xy.                              (3.4)
```

The diagonal of `G3` is fixed by conservation; an explicit form is

```text
(G3)_xx=2(H^4)_xx-2d_x^2-8 sum_y c_xy^2.                 (3.5)
```

The Hilbert--Schmidt adjoint parity of `K` yields

```text
G1^T=G1,       G2^T=-G2,       G3^T=G3.                 (3.6)
```

Trace preservation and unitality give, coefficient by coefficient,

```text
1^T G_j=0,       G_j 1=0.                                (3.7)
```

So the uniform population is stationary at every order, even when the raw
higher-order generator is not Metzler.

## 4. Gauge invariance and the cycle-length filtration

For a vertex gauge `U=diag(exp(i chi_x))`, replace

```text
H by U H U*,       rho by U rho U*.                      (4.1)
```

The diagonal populations are unchanged, so every `G_j` is gauge invariant.
The phase coordinates are therefore the Wilson products

```text
W(C)=product_(xy along C) H_xy,                           (4.2)
```

one for each independent cycle.  Reversing orientation conjugates `W(C)`.

The coefficient `G_j` is homogeneous of degree `j+1` in the hopping entries.
Every gauge-invariant phase-sensitive monomial in it is supported on a closed
walk of length at most `j+1`.  Hence:

- a cycle of length `ell` cannot affect populations before physical order
  `lambda^(-(ell-1))`;
- at its first irreducible appearance, an odd cycle contributes through
  `Im W(C)` and a skew circulation;
- an even cycle contributes through `Re W(C)` and a symmetric correction;
- algebraic sums of cycles with the same endpoints can still cancel.

The parity assertion follows because `Pi K^ell Pi` is skew for odd `ell` and
symmetric for even `ell`; at the first appearance, retardation products use
only shorter, phase-blind closed walks.

On a tree, every edge phase is a vertex coboundary, so all phases can be
gauged away.  This is an all-`lambda`, all-time statement about the full
Lindbladian, not merely an asymptotic cancellation.  The replay checks
`Pi K^n Pi` for phased and real paths through `n=8` as a hostile control.

For a bipartite graph, changing signs on one color class gauges `H` to `-H`.
But `G_j` is homogeneous of degree `j+1`; therefore

```text
G_(2r)=0 for every r>=1 on a bipartite hopping graph.     (4.3)
```

This places the square term at `lambda^-3`, the six-cycle term at
`lambda^-5`, and so on.

## 5. The triangle and square packets

### 5.1 Triangle: oriented circulation

For an oriented triangle `x->y->z->x`, put

```text
Phi_xyz=Im(H_xy H_yz H_zx).
```

Its contribution to `G2`, in cyclic vertex order, is

```text
6 Phi_xyz * [ 0 -1  1]
              [ 1  0 -1]
              [-1  1  0].                               (5.1)
```

It is a divergence-free circulation.  On a triangle, `Phi!=0` is equivalent
to `G2!=0` and, for sufficiently large `lambda`, the leading positive edge
rates dominate.  The truncated slow generator through `lambda^-2` is then a
valid nonreversible Markov generator with uniform stationary distribution.

On a larger graph, triangle circulations may cancel edgewise.  The correct
boundary is `G2=0`, not “every triangle has zero flux.”  If `G2=0`, the
generator is symmetric through `lambda^-3`; if `G2!=0` and the reduced matrix
is Metzler, uniform detailed balance fails at order `lambda^-2`.

### 5.2 Square: symmetric Wilson energy

For a chordless square `0->1->2->3->0`, put

```text
W=H_01 H_12 H_23 H_30.
```

The holonomy-bearing part of `G3` is

```text
Re(W) * [ 4 -8 12 -8]
         [-8  4 -8 12]
         [12 -8  4 -8]
         [-8 12 -8  4].                                 (5.2)
```

The same matrix is the coefficient of this square monomial inside a larger
graph; different squares then add.  Unlike `(5.1)`, `(5.2)` is symmetric and
sees `cos(flux)`, not the orientation of the flux.

The first higher-order Markov hostile already occurs here.  Give every square
edge unit magnitude and take `W=-1`.  Then the opposite pair `0,2`, which is
not an original edge, has

```text
(G3)_02=-16.                                             (5.3)
```

Thus no positive `G1/lambda` term dominates it on that nonedge.  The raw
slow-coordinate generator has a negative virtual transition for every finite
`lambda` at this asymptotic order.  This is consistent with the general
higher-order positivity warnings in adiabatic elimination, but the exact
example and coefficient here are internally derived; no novelty claim is
made.

## 6. The bipartite Gram/wedge theorem

Let `A` be any complex `m by n` matrix and form

```text
H_A=[0  A*]
    [A   0].                                             (6.1)
```

Label the first color class by input columns `u` and the second by output rows
`w`.  Its conductance shadow is `c_wu=|A_wu|^2`.  Define

```text
N_u=sum_w c_wu,
S_uv=sum_w c_wu c_wv,
g_uv=(A* A)_uv=sum_w conjugate(A_wu) A_wv.               (6.2)
```

Because two distinct input vertices are nonadjacent in `H_A`, `(3.3)` reduces
to the exact identity

```text
(G3)_uv=6|g_uv|^2-8S_uv,                   u!=v.          (6.3)
```

Therefore the conductances plus `G3` recover

```text
|g_uv|^2=((G3)_uv+8S_uv)/6,                             (6.4)

E_uv=N_u N_v-|g_uv|^2
    =sum_(w<z)|A_wu A_zv-A_wv A_zu|^2.                  (6.5)
```

Equation `(6.5)` is the complex Lagrange identity.  Summing it gives

```text
sum_(u<v) E_uv=||wedge^2 A||_F^2
 =1/2[(tr A*A)^2-tr((A*A)^2)].                           (6.6)
```

Consequences:

1. `conductance + G3` detects `rank(A)<=1` exactly, since `(6.6)` vanishes
   exactly when every `2x2` minor vanishes;
2. it recovers `tr((A*A)^2)`, hence the quantitative lower bound

   ```text
   rank(A) >= ceil((tr A*A)^2/tr((A*A)^2));              (6.7)
   ```

3. for a `2x2` response block, `(6.5)` is exactly `|det A|^2`.

The row-side analogue follows by exchanging `A*A` with `AA*`.

### Exact response hostile repaired by `G3`

The two actual Jacobian response blocks from the previous session are

```text
A_+=[-1  1],        A_-=[-1  1].
    [-1  1]             [ 1  1]                         (6.8)
```

They have identical unit conductances, but

```text
rank(A_+)=1,       rank(A_-)=2,
||wedge^2 A_+||^2=0,       ||wedge^2 A_-||^2=4.          (6.9)
```

Their `G1` matrices agree and `G2=0` for both; their `G3` matrices differ.
Thus the first finite-dephasing correction repairs exactly the rank blindness
exhibited by the resistor hostile.

## 7. Planar-Jacobian connection contract

For

```text
P=sum_a p_a x^a1 y^a2,       Q=sum_u q_u x^u1 y^u2,
```

the fixed-`P` response matrix is

```text
A_(w,u)=p_a det(a,u),       a+u=w+(1,1),                 (7.1)
```

and the Keller equation is `Aq=kappa e_00`.

| Contract field | Exact content |
|---|---|
| Source | finite coefficient response matrix `A_P` |
| Map | Hermitian dilation `A -> H_A`, then the invariant-graph coefficient `G3` |
| Preserved/recovered | support, `|A_wu|`, column norms, Gram-entry magnitudes, pairwise wedge energies, `||wedge^2 A||_F`, and the rank-`1` versus rank-`>=2` boundary |
| Destroyed | imaginary/oriented Wilson phases, individual square phases when several paths aggregate, phases of Gram entries, target-range orientation, higher compound minors |
| Required sidecars | a cycle-basis phase local system, higher compounds or full Gram phase, the target vector `e_00`, and the coefficient factorization `p_a det(a,u)` together with the global `p_a q_u` Segre constraints |
| Cheapest decisive test | the blocks `(6.8)`, for which conductance agrees but `G3` recovers wedge energies `0/4` |

The transfer is therefore sharper than

```text
response gain graph -> resistor shadow + arbitrary holonomy sidecar.
```

It is now

```text
G1: resistor shadow,
G3: aggregate four-cycle cosine / Gram / wedge-rank sidecar,
remaining data: oriented cycle phase + higher compounds + target + Segre.   (7.2)
```

This is a genuine pruning invariant, not a planar-Jacobian proof.  In
particular, positive `||wedge^2 A||` merely says the response has rank at least
two; it neither puts `e_00` in the range nor constructs a polynomial mate.

## 8. Detailed-balance and positivity boundary

The hierarchy should be read in three regimes.

1. **Singular limit.** `G1/lambda` is a symmetric Markov generator.  The
   resistor/electrical dictionary is exact after time rescaling.
2. **First circulation correction.** `G2/lambda^2` is supported on triangle
   edges.  For sufficiently large `lambda`, it is dominated by positive
   `G1` rates and yields a legitimate nonreversible perturbation whenever
   `G2!=0`.
3. **Virtual-transition regime.** `G3/lambda^3` may act on graph nonedges with
   either sign.  The invariant-graph population coordinate need not carry a
   Markov generator.  “Detailed balance” is then not the correct predicate
   until one either verifies the Metzler condition or chooses and audits a
   different reduced-state parametrization.

This boundary agrees with two pieces of literature, without importing their
theorems as proofs of the formulas here:

- [Bressanini--Benedetti--Paris](https://arxiv.org/abs/2204.01836) establish
  position-dephasing classicalization of continuous-time quantum walks and
  the large-dephasing localization on unrescaled time.
- [Azouit--Sarlette--Rouchon](https://arxiv.org/abs/1603.04630) develop
  invariant-manifold/adiabatic-elimination expansions for Lindblad systems.
- [Tokieda--Elouard--Sarlette--Rouchon](https://arxiv.org/abs/2303.04495)
  exhibit higher-order complete-positivity failures and stress the dependence
  on the reduced-state parametrization.

The bounded search used in this session did not locate the explicit
triangle/square matrices or the Jacobian-response wedge identity in this exact
setting.  That is **not** a novelty claim: a targeted literature and priority
audit remains mandatory before any such statement.

## 9. Exact controls and reproduction

The replay performs all of the following.

- exact block-superoperator verification of the Riccati coefficients;
- exact triangle `G2` and square `G3` matrices;
- phased-tree versus real-tree equality of `Pi K^n Pi` for `2<=n<=8`;
- the square `pi`-flux negative-nonedge hostile;
- exact real and Gaussian-integer complex Gram/Lagrange identities;
- the actual response rank hostile;
- an independent numerical solution of the Riccati equation, with triangle
  error scaling as `lambda^-4` and bipartite-square error scaling as
  `lambda^-5`;
- exact numerical invariance under a random vertex gauge.

Both ordinary and optimized Python runs return `VERDICT: PASS`.

## 10. Concept board and next frontier

| Object | New invariant | Lost coordinate | Cheapest next test |
|---|---|---|---|
| general Hermitian gain graph | `G2` triangle circulation and `G3` square Wilson energy | higher-cycle separation and a positivity-preserving reduced gauge | derive `G5` on `C6`, theta, and two-square graphs |
| bipartite response dilation | `||wedge^2 A||_F^2` from conductance plus `G3` | oriented Gram phase and higher compounds | find same `G1,G3` with different `rank>=2`, or prove the smallest impossible size |
| planar-Jacobian response | rank-`1` boundary and stable-rank lower bound | target range and coefficient Segre compatibility | apply to a finite `(72,108)` support passport after leaf peeling |
| sparse magnetic square | exact negative virtual transition at flux `pi` | a canonical Markov parametrization | compare raw-population Riccati gauge with CP-preserving reductions |
| shortest-cycle filtration | cycle length `ell` first possible at `lambda^(-(ell-1))` | cancellations among equal-endpoint cycles | isolate coefficients on chorded cycles and cycle bases |

The best immediate Jacobian experiment is not another conductance-only search.
For each leaf-peeled finite response passport, compute

```text
G1, G3, ||wedge^2 A||_F^2,
```

then compare the remaining target-range and Segre equations against a basis of
oriented cycle phases.  The first two layers are gauge-invariant and cheap;
the expensive phase solve is reserved for the survivors.

## 11. Scope boundary

- `PROVED`: formulas `(1.1)`--`(7.1)` in finite dimension, under uniform
  position dephasing and `H_xx=0`.
- `FINITE-EXACT`: the listed symbolic and numerical controls.
- `CITED`: only the general classicalization/adiabatic-elimination context in
  Section 8.
- `OPEN`: a positivity-preserving all-order population reduction; isolation
  of individual square phases in dense response graphs; recovery of higher
  compound ranks; any application that excludes the live planar-Jacobian
  degree passports.
- Not claimed: literature novelty, a counterexample to `JC(2)`, a proof of
  `JC(2)`, or an equivalence between finite-dephasing dynamics and the Keller
  response equation.
