---
id: THM-848
title: The Hamiltonian-path Walsh needle stalk gives the exact all-size H-drift and the first missing node colour
status: PROVED ALL n (needle boundary identity, Walsh/OCF drift, exact level-energy EGF, Krylov sufficiency/minimality, Mobius radial-dilation observable functional form, and n<=6 closed forms) + FINITE-EXACT (all converse-merged nodes n=4..7; hash-guarded directional-gradient classifier n=4..8; strong-lumpability boundary through n=5)
source: codex-2026-07-15-S15
depends_on: [THM-062, THM-076, THM-163, THM-203, THM-259, THM-833]
related: [THM-013, THM-201, THM-204, THM-217, THM-589, THM-785, THM-791, THM-810, THM-830, THM-840, HYP-2268, HYP-6900]
verification:
  - 04-computation/h_drift_walsh_needle_stalk_codex_S15.py
  - 05-knowledge/results/h_drift_walsh_needle_stalk_codex_S15.out
  - 05-knowledge/results/h_drift_walsh_needle_stalk_n8_codex_S15.out
  - 04-computation/h_drift_K_second_krylov_referee_codex_S14.py
  - 05-knowledge/results/h_drift_K_second_krylov_referee_codex_S14.out
  - 05-knowledge/results/h_drift_K_second_krylov_n8_codex_S14.out
  - 04-computation/h_drift_forward_polynomial_mobius_codex_S15.py
  - 05-knowledge/results/h_drift_forward_polynomial_mobius_codex_S15.out
  - 04-computation/h_drift_operation_congruence_and_farey_rewrite_codex_S16.py
  - 05-knowledge/results/h_drift_operation_congruence_and_farey_rewrite_codex_S16.out
---

# THM-848 - the Hamiltonian-path Walsh needle stalk

The conditional drift of the Hamiltonian-path count is linear, but not in
`H` itself.  Its exact state is the vector of even homogeneous Walsh layers
of `H`.  This vector is an all-size, isomorphism- and converse-invariant stalk
that can be attached to every node of THM-830's exact-defect groupoid.

At `n=5,6` there are only two nonconstant layers, so `(H,c3)` determines the
drift.  At `n=7` the degree-six layer enters and is genuinely independent of
`H`, `c3`, score sequence, blue/black phase, and local depth.  Retaining all
directional derivatives is stronger still: their unordered multiset is a
complete invariant of all 272 converse-merged size-seven nodes and, in a
separate hash-guarded audit, all 3,528 size-eight nodes.

These statements refine the node-coloured defect algebra.  They do not add
the LRC metric stalk and do not prove LRC(14).

## 1. Flip drift is the boundary census of discrete needles

Let `T` be a tournament on `n` vertices, put

```text
M=binom(n,2),
T^e = T with the unordered arc e reversed,
g_e(T)=H(T^e)-H(T),
b_H(T)=M^(-1) sum_e g_e(T).                              (1.1)
```

Regard a permutation of the vertices as a discrete Hamiltonian needle: its
`n-1` consecutive unordered pairs are the support of the needle, and it is a
directed Hamiltonian path exactly when all those pair constraints agree with
`T`.  Let `U_1(T)` be the number of permutations with exactly one failed
consecutive arc.  Then

```text
sum_e g_e(T)=U_1(T)-(n-1)H(T),
b_H(T)=[U_1(T)-(n-1)H(T)]/M.                             (1.2)
```

Indeed, flipping any of the `n-1` arcs of a present Hamiltonian path destroys
that path.  Conversely, a permutation with one failed adjacency is created
by flipping its unique bad arc.  Summing created minus destroyed needles over
all directions proves (1.2), with no isomorphism or small-size assumption.

Thus the apparently mysterious functional form of `H`-drift is already
exactly named: it is the codimension-one boundary-to-interior ratio of the
Hamiltonian-needle family.  The question "does drift depend only on H?" is
equivalent to "is U_1 constant on every H-fibre?"  It is not.

## 2. The even-path-forest formula

Use the full arc cube, not the fixed-path tiling cube.  For `i<j`, let
`x_(ij)(T)=+1` when `i->j` and `-1` otherwise, and extend antisymmetrically by
`x_(ji)=-x_(ij)`.  Write the homogeneous Walsh decomposition as

```text
H(T)=mu_n + sum_(r=1)^q H_(2r)(T),
mu_n=n!/2^(n-1),                 q=floor((n-1)/2).       (2.1)
```

THM-076's support and amplitude formula has the following evaluation form.
Let `E_(n,r)` be the edge sets `F` of size `2r` whose nontrivial connected
components are paths of positive even length.  Let `c(F)` be the number of
components.  Orient each path component arbitrarily and put

```text
omega_T(F)=product_(oriented path arcs i->j in F) x_(ij)(T). (2.2)
```

Reversing an even path changes `2a` signs, so (2.2) is independent of the
chosen component orientations.  For every `n` and `r`,

```text
H_(2r)(T)
 = (n-2r)!/2^(n-1)
   sum_(F in E_(n,r)) 2^c(F) omega_T(F).                 (2.3)
```

This is the exact surviving-subneedle expansion.  Odd components cancel by
orientation reversal; disjoint even subneedles survive, and the factor
`2^c(F)(n-2r)!` counts their two orientations and block completions.  Formula
(2.3) is also a direct evaluation version of THM-076, rather than a new proof
of the OCF.

For `r=1`, every support is a two-edge path.  Summing its three signed shapes
on each vertex triple gives `3` for a cyclic triple and `-1` for a transitive
triple.  Hence

```text
H_2(T)=a_n [4c3(T)-binom(n,3)],
a_n=(n-2)!/2^(n-2)
     =mu_n [8/(n(n-1))] / 4,                             (2.4)
```

or, with `c3*=binom(n,3)/4` and the THM-833 rate
`kappa_n=8/(n(n-1))`,

```text
H_2(T)=mu_n kappa_n [c3(T)-c3*].                         (2.5)
```

This is the exact reason the cyclic-triangle OU law is the linear shadow of
the `H`-flow.

## 3. Exact drift and the Krylov stalk

Let `P` choose one of the `M` arcs uniformly and flip it.  A Walsh character
of degree `2r` has hypercube eigenvalue `1-4r/M`.  Therefore

```text
P H_(2r)=(1-4r/M)H_(2r),                                (3.1)

b_H(T)=(PH-H)(T)
      =-(4/M) sum_(r=1)^q r H_(2r)(T),                  (3.2)

U_1(T)=(n-1)H(T)-4 sum_(r=1)^q r H_(2r)(T),             (3.2a)

E[H(T_t)|T_0=T]
      =mu_n+sum_(r=1)^q (1-4r/M)^t H_(2r)(T).           (3.3)
```

There is an exactly equivalent odd-cycle-forest coordinate.  In THM-062's
notation let

```text
K(T)=sum_I 2^(|I|+f_I)
    =2^(n-1) sum_I product_(C in I) 2^(2-|C|),           (3.3a)
```

where `I` ranges over vertex-disjoint directed odd-cycle collections.  The
forward-edge distribution gives `b_H=2(K-nH)/M`.  Comparing with (3.2)
proves the all-size dictionary

```text
K=nH-2 sum_r rH_(2r)
 =n mu_n+sum_r (n-2r)H_(2r),                             (3.3b)
U_1=2K-(n+1)H,                    K=nH+(M/2)b_H.         (3.3c)
```

Thus `K` is the Walsh degree-Euler image `(n-D)H`.  The pair `(H,K)` is an
alternative scalar carrier for the one-step averaged drift.  It is not a
replacement for `W_n` when one must propagate all future conditional means:
for `q>1`, (3.3b) is only one weighted moment of the layers.

Define the **Walsh needle stalk**

```text
W_n(T)=(H_2(T),H_4(T),...,H_(2q)(T)).                    (3.4)
```

It is invariant under relabelling and under converse, so it is well-defined
on ordinary and converse-merged metagraph nodes.  It retains `H` as the
zeroth layer moment and `b_H` as its first degree-weighted moment.

The stalk is not merely sufficient.  The `q` eigenvalues in (3.1) are
distinct, and every layer is nonzero by (2.3).  Thus the values

```text
H, PH, ..., P^(q-1)H                                    (3.5)
```

recover `W_n` by a Vandermonde inversion.  Equivalently, the smallest linear
`P`-invariant function space containing `H` has dimension `q+1`, including
the constant.  In this precise linear-dynamical sense, (3.4) is the minimal
all-size information needed to propagate every future conditional mean of
`H`.

For a fixed direction `e`, the stronger formula is

```text
g_e(T)=-2 sum_(S containing e) Hhat(S) chi_S(T).          (3.6)
```

So `W_n` closes the **averaged first moment**, not each directional flip.
The equivariant directional stalk

```text
Grad_H(T)={g_e(T):e in E(K_n)} as an unordered multiset  (3.7)
```

is also well-defined on converse-merged nodes.  Its first moment is (3.2),
its second moment is the conditional `H`-diffusion, and its full multiset
retains the direction-by-direction needle X-ray up to relabelling.

The multiset forgets which derivative reaches which neighbouring node.  Exact
edge transport requires the target-coupled orbit

```text
TargetGrad_H(T)={(g_e(T), pi_n(T^e)):e in E(K_n)}/Aut(T). (3.7a)
```

This distinction matters when the stalk is joined to a coherent relation
algebra: (3.7) preserves the response spectrum, while (3.7a) preserves its
incidence with projected node targets.

There is an all-size stationary fluctuation law as well.  Put

```text
E_(2r)=E_uniform[H_(2r)^2].
```

Walsh orthogonality and (3.6) give

```text
E_uniform[b_H^2]                 =(16/M^2) sum_r r^2 E_(2r),
E_uniform,T,e[g_e(T)^2]          =(8/M)    sum_r r E_(2r),
E_uniform[(H-mu_n)b_H]           =-(4/M)   sum_r r E_(2r),
Cov(H(T_0),H(T_t))               =sum_r (1-4r/M)^t E_(2r). (3.8)
```

In particular the middle two quantities differ by the exact factor `-2`,
the hypercube Dirichlet fluctuation-dissipation identity.  The even-path
forest formula also evaluates every level energy exactly.  Put

```text
Q(w)=(1+w)/(1-w),                 m=n-2r.
```

Then, for every admissible `n,r`,

```text
E_(2r)/mu_n^2
 = [w^r] Q(w)^m / (n)_(2r)                              (3.9)
 = 1/(n)_(2r) sum_(j=0)^min(m,r)
       binom(m,j) binom(m+r-j-1,r-j).                   (3.10)
```

Here is a direct proof, closing the enumeration problem left by THM-204.
Squaring (2.3) and using Walsh orthogonality assigns a forest `F` the weight
`4^c(F)`.  A labelled even path with `2a` edges has EGF contribution
`2x^(2a+1)w^a`: the factor four for its component weight cancels the factor
two from reversal in the number `(2a+1)!/2` of unoriented labelled paths.
Including isolated vertices, the weighted-forest EGF is

```text
exp(x + sum_(a>=1) 2x^(2a+1)w^a)
 = exp(x + 2x^3 w/(1-x^2 w)).                           (3.11)
```

After the substitution `y=x^2 w`, coefficient extraction gives

```text
sum_(F in E_(n,r)) 4^c(F)
 = n!/(n-2r)! [y^r] Q(y)^(n-2r).                        (3.12)
```

Multiplying (3.12) by the squared coefficient
`((n-2r)!/2^(n-1))^2` and dividing by
`mu_n^2=(n!/2^(n-1))^2` proves (3.9).  In THM-201's notation this says

```text
2g_r(m)=[w^r]Q(w)^m.                                   (3.13)
```

The first three coefficient polynomials are

```text
[w]Q(w)^m=2m,
[w^2]Q(w)^m=2m^2,
[w^3]Q(w)^m=2m(2m^2+1)/3.                              (3.14)
```

Consequently THM-204's historical replacement of the coefficient by `2m^r`
is exact for `r<=2` and also on the terminal boundary `m=1`, which explains
every successful check through `n=7`.  It first fails at `n=8,r=3,m=2`:

```text
[w^3]Q(w)^2=12 != 16=2m^3,
E_6/mu_8^2=1/1680 != 1/1260.                            (3.15)
```

Together with `E_2/mu_8^2=3/14` and `E_4/mu_8^2=2/105`,
(3.15) gives

```text
Var(H)/mu_8^2=131/560,                                  (3.16)
```

independently agreeing with THM-589's exact `W(8)=49752`, since
`49752/8!-1=131/560`, and refuting THM-204's prediction `59/252`.

Thus `H` is an exact finite mixture of `q` relaxation modes, not one scalar
OU mode.  The conditional diffusion can still vary inside a `W_n`-fibre
because (3.8) is a stationary average, not a pointwise closure.

## 4. Closed low-size laws and the first failure of pointwise reversion

For `n<=4`, only `H_2` is present, giving

```text
n=4:  b_H=(2/3)(3-H),             U_1=12-H.              (4.1)
```

For `n=5,6`, only `H_2,H_4` occur.  Eliminating
`H_4=H-mu_n-H_2` with (2.4) gives the exact all-tournament laws

```text
n=5:  b_H=(6c3-4H+15)/5,
      U_1=12c3-4H+30,                                   (4.2)

n=6:  b_H=(24c3-8H+60)/15,
      U_1=24c3-3H+60.                                   (4.3)
```

This explains HYP-6900's first split without fitting a nonlinear curve.  At
`n=5`, the two `H=15` nodes have respectively

```text
c3=4, b_H=-21/5;             c3=5, b_H=-3.              (4.4)
```

So `(H,c3)`, not `H`, is the exact two-layer coordinate.

Equivalently, THM-062's odd-cycle coordinates give

```text
b_H=(12-H-6t5)/5,
H=1+2c3+2t5,                                             (4.4a)
```

which is identically (4.2).  The OCF and Walsh descriptions retain the same
first-moment information but organize it by cycle length and cube degree,
respectively.

The common phrase "H drifts toward its random mean" is pointwise true at
`n=4,5`, but already false at `n=6`.  Here `mu_6=45/2`, while the exact merged
nodes

```text
n6-a10: H=23, c3=6, b_H=4/3,
n6-a26: H=25, c3=6, b_H=4/15                             (4.5)
```

both lie above the mean and drift farther upward.  Stationarity still forces
zero global mean drift, and the spectral components still contract according
to (3.1); neither fact implies pointwise scalar mean reversion after projecting
the vector `W_n` to `H`.

## 5. The degree-six obstruction at n=7

At `n=7`,

```text
M=21,       mu_7=315/4,
H_2=15c3-525/4.                                          (5.1)
```

The degree-six forests are exactly the `2520` unoriented Hamiltonian paths of
the complete graph, and (2.3) becomes

```text
H_6(T)=(1/32) sum_(unoriented Hamilton paths Q) omega_T(Q). (5.2)
```

Eliminating `H_4` from (3.2) gives

```text
b_H(T)=[105-8H(T)+60c3(T)-4H_6(T)]/21.                  (5.3)
```

Thus `H_6` is exactly the first new conditional-drift coordinate beyond
`(H,c3)`.

The verifier independently evaluates all forest sums and all 21 arc flips on
one representative of every one of the 272 converse-merged nodes.  The
successive partitions are

```text
colour                                      cells  max fibre  drift-split
H                                              77       8          50
(H,c3)                                        126       6          28
(H,c3,phase)                                  162       6          21
(H,score sequence)                            165       5          20
(H,c3,score,phase,local depth)                 211       5          11
W_7=(H2,H4,H6)                                156       6           0
Grad_H multiset                               272       1           0. (5.4)
```

A strong witness survives every scalar colour in the fifth row:

```text
node       H  c3  score                    phase       depth  (H2,H4,H6)              b_H
n7-a106   49   7  (1,2,2,3,3,5,5)         pure black    3    (-105/4,-21/4, 7/4)     6
n7-a118   49   7  (1,2,2,3,3,5,5)         pure black    3    (-105/4, -9/4,-5/4)    46/7. (5.5)
```

The new information is not a repackaging of score regularity or blue/black
phase.  It is the signed degree-six Hamiltonian-needle balance.

There is also a sharp boundary.  `W_7` has 156 cells and determines drift,
but 68 of its cells split by the conditional second moment
`M^(-1)sum_e g_e^2`.  For example four nodes share

```text
(H2,H4,H6)=(-345/4,63/4,3/4)                             (5.6)
```

and have diffusion values `9332/21, 908/3, 1340/3, 916/3`.  The Walsh stalk
is an exact first-moment closure, not a claim of Markov lumpability.  In
contrast, the gradient multiset (3.7) separates all 272 nodes.  The same
finite audit gives singleton gradient cells at `n=4,5,6` as well, with node
counts `3,10,34`.

There is a direct colour-refinement interpretation.  On the weighted
full-arc-flip quotient, `g_e=H(target)-H(source)`, so `(H,Grad_H)` is exactly
one multiplicity-aware neighbour-colour refinement round seeded by the node
colour `H`.  Its singleton census through `n=7` says that this one round is a
complete merged-node classifier in the audited range.  It does not say the
same refinement remains complete at all later sizes, or on tiling/line fibres
above a node.

The next size has now been tested.  A raw `2^21`-entry little-endian rank
atlas exported by THM-828/843, guarded by SHA-256

```text
30debad3387a4ea0ef51108ea132115efda2ac2fcdfcc2c5c1d4d23155095835,
```

supplies one representative of every size-eight merged node.  The optional
verifier lane computes

```text
n=8: merged nodes=3528, Grad_H cells=3528, max fibre=1.    (5.7)
```

This is a finite-exact classification conditional only on the checked atlas
artifact.  It strengthens the observed range of the one-round classifier; it
does not establish an all-size theorem.

## 6. Refining the node-coloured exact-defect algebra

Use THM-830's coordinates `(z,u,delta)` and node map

```text
c_n(z,u,delta)=pi_n(the tiling with halves u,u+delta and trace z). (6.1)
```

The raw defect algebra remembers the full additive arrow `delta`; its exact
relations satisfy `A_delta A_epsilon=A_(delta+epsilon)`.  The map `c_n`
forgets the marked path and colours those arrows nonlinearly by tournament
nodes.  Define node colours `col_H` and `col_grad`, then pull them back along
`c_n`.  The new, operation-typed seed colours are

```text
Col_H(z,u,delta)=col_H(c_n(z,u,delta))
  =(blue/mixed/black phase of c_n, W_n(c_n)),             (6.2)

Col_grad(z,u,delta)=col_grad(c_n(z,u,delta))
  =(blue/mixed/black phase of c_n, Grad_H(c_n)).          (6.3)
```

For either pulled-back colour `Col`, the exact disintegrated line tensor is
immediately computable:

```text
B_delta(alpha,beta)
 =1/2 #{(z,u):
        Col(z,u,delta)=alpha,
        Col(z+1,u+1,delta)=beta}.                        (6.4)
```

Equation (6.4), followed by colour refinement using all exact-defect
relations, is a finite coherent-configuration computation that can be run
now.  The theorem does **not** assert that (6.2) is already equitable.  What
it proves is the appropriate operation kernel statement:

- `W_n` is closed for every conditional mean under the uniform full-arc flip
  operator `P`;
- `Grad_H` retains every one-step directional `H` change;
- `TargetGrad_H` additionally couples those changes to projected edge targets;
- neither stalk is closed under a face, induced deletion, continued-fraction
  lift, or LRC transport merely because it is closed for `P`.

The defect and needle coordinates preserve complementary information.
`delta` remembers marked-path mirror failure and composes additively, but
forgets the invariant flip response after node projection.  `W_n` forgets
the marked path and exact defect, but retains the spectral response of `H`.
Their fibre product, rather than either scalar alone, is the correct next
node-coloured algebra.

## 7. The n=14 stalk

At `n=14`, the exact-defect word has 36 bits by THM-830, whereas the
nonconstant Hamiltonian needle stalk has only six coordinates:

```text
W_14=(H2,H4,H6,H8,H10,H12),
mu_14=42567525/4,
M=91,                                                     (7.1)

P-eigenvalues=(87,83,79,75,71,67)/91,
b_H=-(4/91)(H2+2H4+3H6+4H8+5H10+6H12).                 (7.2)
```

Thus six successive conditional expectations recover the complete mean-drift
stalk by Vandermonde inversion.  For individual arc transport the natural
upgrade is the 91-entry gradient orbit, or a compressed orbit fingerprint
whose next test is whether it remains complete beyond `n=8`.

The concrete preservation target for LRC14 is therefore not another scalar
node coordinate.  It is

```text
(boundary trace z, half object u, 36-bit defect delta,
 six-layer W_14, owner/carry/metric stalk).               (7.3)
```

The final fields in (7.3) remain external and indispensable.

## 8. Kakeya reading, Tournament Analysis, and boundaries

The Kakeya language here is precise but discrete.  Hamiltonian permutations
are needles in the tournament edge-coordinate cube; `H` counts fully admitted
needles, `U_1` counts their one-facet boundary, the even-path forests are the
reversal-stable subneedles, and `Grad_H` is the directional X-ray.  No
continuous Kakeya measure, dimension, maximal-function, or polynomial-method
theorem is being imported.

For Tournament Analysis, the seven vertices are the seven candidate node
colours in (5.4).  The pairwise observable is the number of unordered pairs
of the 272 nodes separated by the colour.  The first gauge uses raw separated
pairs; the switch divides by `ceil(log2(number of colour cells))` to measure
retention per address bit.  The displayed order supplies the tie Hamiltonian
path.  Both resulting tournaments are transitive, with score histogram
`{0:1,...,6:1}`, no directed triangle, singleton SCCs, and one Hamiltonian
path, but the switch flips 14 of 21 edges.  Raw completeness favours the
gradient; bit economy substantially reorders the middle refinements.

The challenged assumption is that Tournament Analysis vertices must be
runners, original tournament vertices, or arcs.  We considered those, path
gaps, Hamiltonian needles, Walsh layers, defect words, and proof obligations.
Candidate observers are the useful vertices here because the question is
exactly which quotient preserves flip response.

The quotient preserves tournament isomorphism, converse, `H`, every full-arc
Walsh layer, conditional `H` drift, and, for (3.7), the unordered directional
derivative field.  It destroys the labels of individual arc directions,
marked Hamiltonian path, exact mirror defect, face seams, runner speeds,
metric gaps, scale, owner, wall chronology, and loneliness.  No implication
for the fourteen-runner conjecture follows without those missing stalks. ∎

## 9. Live-pull referee addendum (codex-S14): `K` is the second Krylov needle

The first-moment law above has an exact next step.  Put

```text
S f(T)=sum_(e in E(K_n)) (f(T^e)-f(T)),    M=binom(n,2),
A_T(x)=sum_(j=0)^(n-1) a_j(T)x^j.
```

For a fixed vertex order, flipping a nonconsecutive pair leaves its number of
forward adjacencies unchanged.  Among its `n-1` consecutive pairs, each
forward pair lowers that number by one and each backward pair raises it by
one.  Double counting orders and flip directions therefore gives, with
`a_(-1)=a_n=0`,

```text
S a_j=(j+1)a_(j+1)+(n-j)a_(j-1)-(n-1)a_j,              (9.1)
S A_T=(1-x)[(1+x)A_T'-(n-1)A_T].                        (9.2)
```

Thus the top forward coefficients form a finite tridiagonal Krylov ladder.
Writing `b_r=a_(n-1-r)`, (3.3c) and (9.1) yield

```text
K=(b_1+(n+1)b_0)/2,
S K=b_2+b_1-Mb_0,                                       (9.3)

E[Delta K|T]
 =[a_(n-3)+a_(n-2)-M H]/M
 =-(4/M) sum_(r=1)^q r(n-2r)H_(2r).                    (9.4)
```

The last equality also follows immediately by applying the flip operator to
`K=n mu_n+sum_r(n-2r)H_(2r)`.  In particular `K` lies in the first Krylov
span generated by `H`, while its drift is the next degree-weighted moment.

There is also a direction-by-direction OCF current, stronger than (9.4).
For an unordered pair `e={u,v}`, let `C_e^-(T)` be the directed odd cycles of
`T` which use its present arc and let `C_e^+(T)` be the directed odd cycles of
`T^e` which use the reversed arc.  Then

```text
K(T^e)-K(T)
 =4[ sum_(C in C_e^+(T)) K(T-V(C))
    -sum_(C in C_e^-(T)) K(T-V(C)) ],                   (9.5)
```

with the natural convention `K(empty)=1/2`.  To prove (9.5), normalize
`Z(T)=2^(1-n)K(T)=sum_I product_(C in I)2^(2-|C|)`.  Every changed OCF cycle
contains both endpoints of `e`, so an independent OCF family contains at
most one changed cycle.  The remaining cycles form an arbitrary OCF family
on `T-V(C)`.  A changed cycle of length `l` consequently contributes

```text
2^(n-1) 2^(2-l) Z(T-V(C)) = 4K(T-V(C))
```

to `Delta K`, with the displayed sign.  This proves the formula for every
edge and every size; the empty convention handles a spanning odd cycle.

### 9.1 Exact fibre boundary through `n=7`

The hash-guarded referee independently replayed THM-848 and then evaluated
(9.1)--(9.5) on one representative of every converse-merged node through
`n=7`.  In the table, `split` is the number of fibres on which the indicated
conditional drift has more than one value.

| `n` | nodes | `H` fibres split by `b_H` | `H` fibres split by `b_K` | `K` fibres split by `b_K` | `(H,K)` fibres split by `b_K` |
|---:|---:|---:|---:|---:|---:|
| 4 | 3 | 0 | 0 | 0 | 0 |
| 5 | 10 | 1 | 1 | 0 | 0 |
| 6 | 34 | 5 | 0 | 0 | 0 |
| 7 | 272 | 50 | 50 | 10 | 0 |

So the first `H`-fibre split occurs at `n=5`, the first `K`-fibre split for
its own drift occurs at `n=7`, and `(H,K)` closes both averaged first steps
through `n=7`.  This last fact is unexpectedly rigid and breaks at the next
size; it is not an all-size two-scalar theorem.

Indeed, at `n=7` write the OCF level sums as

```text
alpha_6=1,
(A,B,C)=(alpha_4,alpha_2,alpha_0)
       =(2t3, 2t5+4 alpha_33, 2t7).
```

Then

```text
H=1+A+B+C,
K=64+16A+4B+C,
a_4=1191+15A-9B+15C.                                  (9.6)
```

The ambient kernel of `(H,K)` is `(1,-5,4)`: two level vectors with the same
pair would differ by `lambda(1,-5,4)`, changing both `a_4` and `S K` by
`120 lambda`.  The determinant after adjoining `a_4` is `-360`, so a third
moment is algebraically independent.  Nevertheless, the exact node census
has

```text
# (H,K) cells = # OCF-level cells = # Walsh-stalk cells = 156,
```

and realizes no nonzero step in that kernel.  Hence `(H,K)` actually retains
the whole three-layer stalk on the size-seven tournament locus.  The same
audit checks (9.5) on all `6,340` node-edge incidences through `n=7`; the
unordered `K`-gradient multiset separates all `272` size-seven merged nodes.

The optional hash-guarded `n=8` lane locates the exact break.  On all `3,528`
converse-merged nodes there are `1,727` `(H,K)` cells, and exactly **one** of
them splits by `b_K`.  Its two rows are

```text
rank   mask   H    K     a_6   a_5   S K    b_K       (alpha_5,alpha_3,alpha_1)
3338  51915  345  2820  2535  6165  -960   -240/7    (30,184,130)
2367   7767  345  2820  2535  6885  -240    -60/7    (32,174,138).
```

The level-vector difference is exactly `2(1,-5,4)`: the first ambient kernel
direction left unrealized at `n=7` appears at `n=8`.  Here

```text
a_5=4293+189 alpha_5-27 alpha_3+9 alpha_1,
```

so `Delta a_5=Delta(S K)=720` and the conditional-drift gap is `180/7`.
The scratch atlas is the same `2^21`-entry object guarded in §5 by SHA-256
`30debad3387a4ea0ef51108ea132115efda2ac2fcdfcc2c5c1d4d23155095835`.
This is a finite-exact statement conditional on that checked artifact.

For Tournament Analysis, the vertices were the observer quotients `H`,
`(H,K)`, `(H,K,a_4)`, the OCF level vector, and the `K`-gradient multiset,
not runners, original tournament vertices, or arcs.  Pairwise separation of
node pairs gives a transitive observer tournament (score histogram
`{0:1,...,4:1}`, no directed triangle, singleton SCCs, one Hamiltonian path).
Switching to separation per address bit flips seven of ten edges but remains
transitive.  This challenges the assumption that the relevant tournament
must live on the original object: here it lives on competing quotients of
the proof obligation.

The `n=8` witness confirms that the hierarchy does not stop.  The number of
nonconstant Walsh/OCF levels grows with `n`; `(H,K)` has only two moments.
Formula (9.2), or equivalently the full stalk of §3, is the honest all-size
closure. ∎

## 10. The exact functional form: Mobius radial dilation

Formula (9.2) diagonalizes without solving a new recurrence.  Put `m=n-1`
and

```text
Phi_k(x)=(1-x)^k(1+x)^(m-k),       0<=k<=m.
```

For `L A=(1-x)[(1+x)A'-mA]`, direct differentiation gives

```text
L Phi_k=-2k Phi_k.                                      (10.1)
```

The full forward-adjacency polynomial has the exact expansion

```text
A_T(x)=sum_(r=0)^floor(m/2)
       H_(2r)(T)(1-x)^(2r)(1+x)^(m-2r).                (10.2)
```

To prove it, sum `x^(number of forward adjacencies)` over all vertex orders
and expand each consecutive-edge factor into its constant and signed parts.
The degree-`k` term is a homogeneous degree-`k` Walsh polynomial times
`Phi_k`.  Reversal of the order cancels odd `k`; taking the coefficient of
`x^m` identifies the surviving coefficient with `H_(2r)`.  Equivalently, if
`A_T=sum_j a_jx^j`, the inverse Krawtchouk formula is

```text
H_(2r)(T)=2^(-m) sum_(j=0)^m K_(2r)(j;m)a_j(T),         (10.3)
```

and every odd transform is zero.  Thus the forward-count vector, the full
Walsh needle stalk, and the even Krawtchouk transform are three lossless
coordinate systems for the same averaged-flip state.

Now use the Mobius coordinate

```text
z=(1-x)/(1+x),
B_T(z)=sum_r H_(2r)(T)z^(2r).
```

Then (10.2) and (10.1) become

```text
A_T(x)=(1+x)^m B_T(z),
S A_T=(1+x)^m[-2z B_T'(z)],                             (10.4)
e^(tS)A_T=(1+x)^m B_T(e^(-2t)z).                       (10.5)
```

This is the requested functional form.  The continuous flip-sum flow is
literal radial dilation in `z`; the discrete uniform-flip chain multiplies
the `2r` mode by `(1-4r/M)^t`.  The apparently nonlinear scalar `H` profile
is evaluation at `x`'s top coefficient after projecting away the radial
modes.  The first split of an `H` fibre, and the later split of an `(H,K)`
fibre, are therefore projection collisions rather than failures of a hidden
one-variable mean-reversion law.

For Tournament Analysis, take the even modes as vertices and orient a pair
from slower to faster generator decay.  At `n=14` this is the transitive
seven-vertex tournament with score histogram `0,...,6`, no directed cycle,
seven singleton SCCs, and one Hamiltonian path.  This quotient preserves the
entire averaged flip semigroup and every `a_j`, but destroys arc labels,
directional currents, and all LRC metric data.  The useful vertices are thus
spectral proof obligations, not tournament arcs or runners. ∎

## 11. Koopman closure is not strong lumpability (codex-S16)

The radial formula closes the averaged observables; it does not make the map
`T -> B_T` a Markov-state quotient.  Strong lumpability would require the
multiset of successor `B`-states under the `M` individual arc flips to be
constant on every `B`-fibre.  An exact converse-merged audit finds zero split
fibres for `n=2,3,4`, but two at `n=5`.

The first explicit kernel pair has full-arc masks `8` and `10` (atlas states
`n5-a01` and `n5-a05`).  Both have

```text
A_T=(9,30,42,30,9),   B_T=(3/2,0),
H=9,   c3=3,   score sequence=(1,1,2,3,3),
E[Delta H]=-3/5.
```

Their individual `H`-gradient multisets are respectively

```text
(-8,-4,-4,0,0,0,0,2,2,6),
(-6,-6,-4,-4,0,0,2,4,4,4),
```

so their one-step diffusions are `14` and `78/5`; their target-`B`
histograms differ as well.  More directly, the same labelled flip `{1,2}`
sends their common stalk to `(-3/2,-1)` and `(-9/2,0)`, respectively.  Thus
even the complete forward polynomial fails the labelled transition square as
well as strong lumpability for the uniformly random arc-flip chain.

The operator distinction is exact.  If `P=I+S/M` is one uniform discrete
flip, then

```text
P B(z)=B(z)-(2/M)zB'(z),
P^t H_(2r)=(1-4r/M)^t H_(2r),
```

whereas only the Poissonized semigroup satisfies

```text
exp(tS)B(z)=B(exp(-2t)z).
```

No single discrete scale factor works once the degree-two and degree-four
modes both occur: the required fourth-power multipliers differ by
`(1-4/M)^2-(1-8/M)=16/M^2`.  Möbius radiality is therefore an exact Koopman
observable conjugacy, not a state operation-congruence.  It preserves every
future mean of `H` while losing successor incidence, diffusion, arc owners,
and wall chronology. ∎
