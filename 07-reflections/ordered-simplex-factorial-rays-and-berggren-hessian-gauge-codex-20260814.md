# Ordered-simplex factorial rays and the Berggren Hessian gauge

**Research synthesis, 2026-08-14.**  The proof sources are THM-3362,
THM-3367, and THM-3368; this reflection records the concept board, transfers,
hostiles, and next questions.  `FC(3)`, `JC(2)`, and `LRC(14)` all remain
open.

## Outcome first

One ternary branch alphabet produced two exact but different determinant
mechanisms.

1. The determinant of THM-3357's three-dimensional nodewise-transfer pencil
   is the split cubic

   ```text
   P=(x-y-z)(x+y-z)(x+y+z).
   ```

   Its ordered-simplex coordinates separate two arbitrary real odd profiles.
   THM-3362 proves an infinite, generally nonhomogeneous ambient `FC(3)`
   subclass: the first three factorial moments already detect every
   `A+B R+C R^2` in that family.

2. The determinant of the underlying two-dimensional spinor pencil is the
   Lorentz quadratic

   ```text
   q=a^2-b^2+c^2.
   ```

   THM-3367 proves that one fixed target gauge turns this pencil into the
   universal symmetric `2 x 2` pencil.  Row-integrability is exactly the
   Hessian condition.  Affine rank-two coefficient pullbacks cannot be
   Keller, every constant-determinant affine-line image is explicitly tame,
   and the highest Hessian layer of any remaining symmetric Keller map has
   only one isotropic coefficient direction.

The shared branch syntax is useful, but it is not an identification.  The
representation-level operation and the aggregation operation do not commute.
That failure is the sharp boundary preventing a spurious `FC -> JC` transfer.

## Inheritance pass

| role | inherited object | what was retained |
|---|---|---|
| closest proved mechanism | THM-3018 simplex moments; THM-3357 three-branch transfer | exact factorial integral and the labelled `L/M/R` matrices |
| canonical hostile | THM-3335's same compressed determinant data/different clocks; MISTAKE-237 | phase is not a scalar gate, and a symmetric binary JC subproblem is not general `JC(2)` |
| corrected near miss | MISTAKE-350 | `FC(3)` means three ambient variables and is not homogeneous-only |
| least-used sidecar | the ordered triangle `(U/T,V/T)` and the fixed target gauge `M^{-1}` | angular ordering for FC; row-integrability for JC |

The session kept four live concepts on the board:

- factorial moment orientation on a short polynomial subspace;
- the two-dimensional symmetric/Hessian Keller equation;
- the LRC determinant-gate Horn implication and its missing seed;
- representation before aggregation versus aggregation before
  representation.

Each pull was compared against all four rather than being treated as a local
identity.

## Pull 1: an ordered-simplex tensor-product compiler for `FC(3)`

Put

```text
T=x+y+z,       U=x-y-z,       V=x+y-z.
```

The exponential factorial integral sends the positive octant to the standard
simplex.  In the angular coordinates

```text
u=U/T,       v=V/T,
```

that simplex is the ordered triangle `-1<=u<=v<=1`.  For nonzero real odd
polynomials `h,k`, their homogenizations, and `e>=0`, define

```text
R=T^e H_h(U,T) H_k(V,T),       D=deg(h)+deg(k)+e.
```

The reflection `(u,v)->(-v,-u)` preserves the ordered triangle and exchanges
the two odd factors.  Consequently every factorial moment is compiled by

```text
calL(R^r)
  =(Dr+2)!/8 * integral(h(t)^r,-1,1) * integral(k(t)^r,-1,1).
```

This is stronger than the initially visible tensor square: `h` and `k` need
not agree.  The example `h=t,k=t^3,e=0` gives `R=UV^3` and

```text
[calL(R^r)]_(r=0)^4=[1,0,86400,0,49249028505600].
```

Odd powers vanish.  For the even moments `p,q,w` of orders `2,4,6`, strict
variance gives `q>p^2`, two Cauchy--Schwarz inequalities control the two
profiles, and the radial factorial gap gives

```text
p w/q^2 >= (7/5)^(2D) >= 2401/625 > 3.
```

The exact cubic eliminant for `A+B R+C R^2` is therefore positive.  The first
three factorial moments force `A=B=C=0`.  This is a genuine ambient `FC(3)`
result because a general member mixes degrees `0,D,2D`; it is not merely an
HFC statement.

Two boundaries matter.

- Two moments never suffice: `A=-p,C=1` and
  `B^2=-(q-p^2)/p` kill the first two, while the third is positive.
- Reality up to scalar phase cannot be deleted.  THM-3362 gives an exact
  complex odd profile for which moments one, two, and three vanish but the
  fourth does not.  This is a failure witness for the proposed short-window
  extension, not an `FC(3)` counterexample.

Common parity, not equality of the two profiles, is the factorization
sidecar.  The mixed pair `h=1,k=t` gives `R=V` and `calL(V)=1`, whereas the
false half-square formula predicts zero.  This hostile refutes the angular
factorization only; it does not assert that every mixed-parity three-slot
subspace fails first-three detection.

## Pull 2: the spinor pencil is exactly a Hessian pencil

For THM-3357's branch matrices `L,M,R`, write

```text
B(a,b,c)=aL+bM+cR.
```

With `D=diag(-1,1)` and `J` the coordinate swap,

```text
L=MD,       R=MJ,
B=M [[b-a,c],[c,a+b]].
```

Thus the fixed target gauge `M^{-1}` exposes the universal symmetric pencil,
and

```text
det B=a^2-b^2+c^2=-det Hess(H).
```

The coefficient pullback is the Jacobian matrix of a polynomial pair exactly
when the symmetric matrix is a Hessian.  Explicitly,

```text
a=(H_tt-H_ss)/2,  b=(H_tt+H_ss)/2,  c=H_st,
(P,Q)=(H_t,H_s+2H_t).
```

This supplies an exact source-to-target dictionary, not a reduction of
general `JC(2)` to the symmetric case.

Three useful consequences follow.

1. If `(a,b,c)` is affine in the two source variables and its affine direction
   space has rank two, then the quadratic part of `q` would be a totally
   isotropic plane in a nondegenerate ternary quadratic space.  Its Witt index
   is one, so this is impossible.  Any affine constant-`q` pullback has
   coefficient rank at most one.

2. If an integrable constant-`q` coefficient image lies in an affine line,
   then

   ```text
   H(r)=1/2 r^T C r+F(g^T r)+linear,
   g^T C^(-1)g=0.
   ```

   The inverse of the gradient map is explicit: after removing the fixed
   target gauge and the linear offset, set

   ```text
   sigma=g^T C^(-1)y,
   r=C^(-1)(y-F'(sigma)g).
   ```

   Hence the whole affine-line image class is tame.

3. If `det Hess(H)` is a nonzero constant and `deg H=d>=3`, then the top
   homogeneous piece satisfies `det Hess(H_d)=0`.  The binary homogeneous
   zero-Hessian lemma forces `H_d=A ell^d`.  Therefore the leading coefficient
   triple is one isotropic direction times `ell^(d-2)`.  A second coefficient
   direction can only enter below the top Hessian layer.

The last fact is a precise stopping result for the metallic-circuit idea: the
top layer has one channel, so a two-direction THM-3010 maximal-alternation
circuit cannot live there.  Lower-layer target-shear descent is the missing
sidecar.

## The operation-order hostile separating `FC` from `JC`

The cubic and quadratic determinants come from the same labelled branch
matrices but from different orders of operation:

```text
individual branches -> symmetric-square transfer -> weighted aggregation
individual branches -> weighted aggregation -> symmetric square.
```

These paths agree at a single branch and diverge on mixtures.  At weights
`(1,1,1)`,

```text
det Sym^2(L+M+R)=det(L+M+R)^3=1,
det(T_L+T_M+T_R)=P(1,1,1)=-3.
```

So the FC cubic is not the JC Hessian quadratic in disguise.  The preserved
data are branch labels and each individual representation; the destroyed
datum is compatibility with linear mixing.  The necessary sidecar for any
future transfer is the operation order itself.

## LRC anchor: a weighted Horn tariff and an exact no-seed boundary

THM-3368 extracts more quantitative content from THM-3357's Horn implication
than its original pass/pass corollary.  For child scores `G_X` and defects
`A_X=max(0,-G_X)`, its exact inequality immediately gives

```text
c G_M >= K-a A_L-b A_R.
```

Therefore

```text
a A_L+b A_R<=K  ==>  G_M>=0,
```

strictly when the bill is below `K`.  This tolerates two failing outer gates.
In the proved THM-2057 AP one-tail deck, parent `(353,356)` has

```text
(G_L,G_M,G_R)=(-169080,1066,-3893),
a A_L+b A_R=1338084208 < K=1606017978,
c G_M=K-a A_L-b A_R=267933770.
```

The tariff certifies the middle while THM-2057's independent clocks certify
both outers.  The full sibling packet closes even though the old pass/pass
Horn rule cannot fire.

The converse bridge is sharply false.  THM-3363 closes the physical row

```text
(1,2,3,4,6,12,168,36,60,108,132,156,180).
```

For the fixed child `d=(2,3)`, one can vary only the last column of a saturated
deck from `(90,0)` to `(90+3q,-2q)`.  The physical row `C_q^T d` remains
literal and saturation remains witnessed by a different determinant-one
pair, but

```text
Delta_(C_q)(d)=270+13q,
G_(91,C_q)(d)=-24557-1183q -> -infinity.
```

Thus even one fixed complement-clock-closed row supplies neither a Horn seed
nor a finite gate-defect bound after the saturated deck/gauge is forgotten.
The simpler logical hostile is the AP-deck root `(1,2)`: all three sibling
rows are clock-safe, yet all three gate scores are negative.

The lawful merge with proved THM-3366 is a typed decision pipeline.  Its
current exact support join independently deletes `298` refined `k=2` rows. Retain
the full labelled deck, child parameter, gate score and defect; separately
compile the physical row into its labelled body/tail support chart.  Apply the
weighted defect tariff first, then run complement-clock exits on every still
uncertified row.  Never replace a clock exit by `A_X=0`.  This improves finite
search pruning and sibling-packet closure; it does not prove `LRC(14)`.

## Connection ledger

| source -> target | map | preserved | lost / required sidecar |
|---|---|---|---|
| odd interval profiles -> `FC(3)` three-slot subspace | ordered-simplex homogenization | all target moments via an exact radial product; moment-zero predicates | separate factor moments, factor order, mixed-parity and nonseparable coupling |
| Berggren weights -> symmetric Keller subproblem | fixed target gauge `M^{-1}` | determinant and row-integrability | positivity chambers, branch words, ancestry; MISTAKE-237 blocks general-JC promotion |
| transfer cubic -> spinor quadratic | change representation and operation order | individual branch labels | weighted-sum compatibility; `(1,1,1)` is the hostile |
| complement clocks + Horn -> LRC search | defect tariff, then typed decision-DAG exits | middle gate certificate and pointwise row terminals | support quotient loses deck/gauge; retain labelled reconstruction and never set a clock-safe defect to zero |

## What genuinely moved

- **PROVED:** a new infinite ambient `FC(3)` subclass, including unequal odd
  profile pairs and every odd power of the THM-3357 determinant ray.
- **PROVED:** an exact Hessian gauge for the Berggren spinor pencil, a global
  affine rank-two coefficient no-go, explicit tameness of all affine-line
  constant-determinant images, and a top-layer one-direction theorem.
- **PROVED:** a weighted LRC determinant-defect tariff that can certify a
  middle sibling despite two failing outer gates, plus a fixed-row unbounded
  deformation proving complement-clock safety cannot seed that tariff.
- **REFUTED / narrowed:** equality of the FC angular profiles is not
  load-bearing; common odd parity is.  The FC cubic and JC quadratic are not
  interchangeable.  Complement-clock safety cannot supply a Horn gate seed.
- **UNCHANGED OPEN:** `FC(3)`, `JC(2)`, `DC(2)`, and `LRC(14)`.

## Cheapest next tests

1. **Factorial:** classify complex odd profile pairs by the three scalar
   inequalities behind the eliminant, rather than by coefficient reality.
   Test whether a Hermitian angular sidecar yields a larger phase-stable
   family without reviving the exact complex hostile.
2. **Jacobian:** after normalizing `H_d=A ell^d`, compute whether the next
   Hessian layer can always be removed by a target shear in any new bounded
   degree range.  Keep the symmetric-subproblem type explicit.
3. **LRC:** on an actual THM-2056 residual deck, record both the complement
   support word and the three sibling gate scores.  Search for a quantitative
   implication only inside that common labelled fibre; the global
   safe-to-gate implication is already false.
4. **Representation:** for other ternary trees, compare
   `sum rho(B_i)` with `rho(sum B_i)` before importing determinant identities.
   The one mixed-weight hostile should be mandatory.

The session's reusable procedural lesson is that two views built from the
same generators can share every vertex and still disagree after aggregation.
Always type the operation order before treating a representation-level
identity as a cross-frontier bridge.
