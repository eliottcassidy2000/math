# Normalized C3 vertical-equalizer audit

**Verdict: PASS for the new exact THM-3765-to-THM-3770 bridge; REPAIR one
proof-scope mismatch in the now-promoted THM-3770.**  Its unqualified
nonconstant-dressing no-go is true, but Section 3 currently states it only
after imposing birationality; a two-line application of Section 2 removes
that unnecessary hypothesis and should be made explicit.  The normalized
nonlinear rational-exact C3 family is exactly a nonconstant spectral dressing
of a birational equalized log-canonical pair.  Smoothness separates the old
component value from the appended dressing root, while polynomial regularity
would identify them.  The resulting nonzero vertical equalizer obstruction
can be moved from one component to the other but cannot be cancelled.  Two
canonical one-component rational primitives and a birational triangular
chart are additional exact consequences.

Audit surface: `origin/main` through `125aef69f`, including proved and
hostile-audited THM-3765 at `84bca3ba8`, the new THM-3770 proof at
`238ba8864`, blowdown sharpening at `9207339ec`, and promotion at
`b48410210`; reserved THM-3759 at `e0a1e09bc`, proved THM-3757 at
`783c40617`, and proved THM-3758 at `b778b12a9`.

## Exact C3 corollary

Work over an algebraically closed characteristic-zero field and use
`J(P,Q)=P_X Q_T-P_T Q_X`.  On the nonlinear smooth THM-3765 boundary take

```text
g,p,c != 0,
u=1+gT/2,
Q=h+pT+Xu^2,
q0=h-2p/g,
v=2p+gXu,
D=g(Q-q0)=uv.
```

Then the following statements are exact.

1. The exceptional fibre `Q=q0` is the disjoint reduced union `u=0` and
   `v=0`.  They are disjoint because `v-gXu=2p`.  The factor `u` is linear;
   `v` is primitive and degree one in `X` over `k[T]`, with
   `gcd(gu,2p)=1`, so both are irreducible.

2. For THM-3765's primitive `P0=cY/D`, where

   ```text
   Y=X-T[p+(g^2/4)XT],
   ```

   the two exact restriction identities are

   ```text
   Y-2p/g = u[X(2-u)-2p/g],
   Y+2p/g = (2-u)v/g.
   ```

   In the common target uniformizer `t=Q-q0`, the simple principal
   coefficients of `P0` are therefore

   ```text
   u=0: +2cp/g^2,                 v=0: -2cp/g^2.       (1)
   ```

3. A target correction with simple part `a/t` adds the diagonal vector
   `(a,a)` to `(1)`.  Simultaneous cancellation would require

   ```text
   a=-2cp/g^2  and  a=+2cp/g^2,
   ```

   which is impossible.  A higher-order target pole creates an uncancellable
   higher pole on both components, and a target function regular at `q0`
   changes neither coefficient.  Since THM-3765 proves that every rational
   mate is `P0+H(Q)`, the invariant equalizer defect of the complete mate
   torsor is

   ```text
   res_(u=0)-res_(v=0)=4cp/g^2 !=0.                   (2)
   ```

4. The obstruction has two canonical one-component representatives:

   ```text
   Pu=2c/(gu),                      J(Pu,Q)=c,
   Pv=2cX/v,                        J(Pv,Q)=c,          (3)
   Pv-Pu=-4cp/(gD)=-4cp/[g^2(Q-q0)] in k(Q),
   P0+c/g=(Pu+Pv)/2.
   ```

   `Pu` has residue vector `(4cp/g^2,0)` and a pole only on `u=0`; `Pv` has
   vector `(0,-4cp/g^2)` and a pole only on `v=0`.  Thus a target shear can
   clear either component, never both.  This is the minimal universal
   obstruction in this family.

5. The mate `(Pu,Q)` is explicitly birational.  If `R=Pu` and `L=Q`, its
   inverse is

   ```text
   u=2c/(gR),
   T=2(u-1)/g,
   X=(L-h-pT)/u^2.                                    (4)
   ```

   Hence the nonlinear boundary is rationally triangular after a reciprocal
   source coordinate, while its unavoidable finite pole prevents polynomial
   entry.

This also gives an exact family-wide dichotomy.  Inside the complete smooth
rational-exact locus of THM-3765, vertical cancellation is compatible with a
polynomial mate exactly on `g=0`, where `Q=X+h+pT` and the displayed linear
mates in THM-3765 apply.  If `g!=0`, smoothness forces `p!=0` and `(2)` is
nonzero.  The excluded `g!=0,p=0` boundary is singular and is not evidence
for either direction.

## Exact THM-3770 spectral-dressing realization

The two fibre factors themselves form the log-canonical pair required by
THM-3770.  Put

```text
U=u,                 W=v,                 m=-g^2/2.
```

Then

```text
J(U,W)=mU,
Spec_0(U,W)=(2p),
(W-2p)/U=gX,                      J(U,gX)=m.           (5)
```

Thus the underlying irreducible-`U` spectrum is equalized and the blowdown
quotient recovers the elementary Keller pair `(u,gX)`.  It is birational:

```text
T=2(U-1)/g,                       X=(W-2p)/(gU),
k(X,T)=k(U,W).                                            (6)
```

After the target translation, the nonlinear C3 member is exactly THM-3770's
nonconstant dressing

```text
Q-q0=U phi(W),                    phi(W)=W/g.          (7)
```

The appended root is `0`; the old component spectrum is `2p`.  THM-3765's
smoothness condition `p!=0` is therefore precisely THM-3770's root-avoidance
condition `phi(2p)!=0`, while regularity would require `2p=0`.  The two gates
are incompatible on every nonlinear smooth member.  THM-3770's canonical
dressed primitive, scaled to response `c`, is not merely analogous to `(3)`:

```text
-cW/[m(Q-q0)] = 2c/(gu)=Pu.                           (8)
```

This independently instantiates the theorem's principal coefficients:
`-2p/m=4p/g^2` on `U=0` and `-0/m=0` on `W=0` for response one.

For completeness, with `delta=J(-,Q)` the same factors obey

```text
delta(u)=-(g/2)u^2,             delta(v)=+(g/2)uv,
delta(log u)=-(g/2)u,           delta(log v)=+(g/2)u,
delta(uv)=0,
delta(u^-1)=g/2,                delta(X/v)=1/2.        (9)
```

So the successful rational channels are the reciprocal Darboux responses
`u^-1` and `X/v`; a constant-exponent formal log combination has either a
nonconstant response or zero and lies outside `k(X,T)`.  This last sentence is
only a side observation.  The typed THM-3770 dressing is `(7)`, and it passes
the exact C3 control.

## Precise bridges to adjacent canon

- **THM-3758 (PROVED).**  Both families pass generic rational exactness and
  then carry a nonzero anti-diagonal principal-part class modulo target
  diagonals on a two-component fibre.  The new C3 calculation strengthens
  the analogy by exhibiting the two one-component primitives and their exact
  target transition.  Both are now exact applications of THM-3770's general
  vertical recursion; only the C3 family is additionally identified here as
  its birational nonconstant spectral dressing.

- **THM-3759 (PROVED + INDEPENDENTLY HOSTILE-AUDITED).**  At the level of
  the displayed constant-flank ansatz alone, THM-3765's boundary has
  `psi=p+(g^2/4)z`.  It belongs to the constant-flank slice only if
  `psi'=g^2/4=0`, hence `g=0` in characteristic zero; then
  `chi=h+gz` is also constant and `Q` is linear.  Thus the nonlinear C3
  vertical wall and the proposed nonlinear constant-flank boundary are
  disjoint.  THM-3759's nonlinear one-flank cases fail even rational exactness
  by THM-3551, so they stop before THM-3770.

- **THM-3757 (PROVED).**  Every Pell-Chebyshev member stops at the earlier
  generic-fibre exactness gate: `dz/Y` is holomorphic of positive genus for
  `n>=2`, and has nonzero conic residues for `n=1`.  Also
  `deg chi_n=n+1>=2`, whereas the rational-exact THM-3765 boundary has affine
  `chi`; there is no overlap.  The split fibre of the `n=1` control must not
  be described as a vertical equalizer failure because there is no rational
  primitive torsor to regularize.

## Hostile audit of the new general THM-3770

The proof body passes its typed scope.

1. Smoothness makes every `Q-lambda` squarefree, so the common target
   parameter has valuation one on each vertical component.  The descending
   DVR recursion subtracts the unique diagonal scalar at each pole order and
   is necessary and sufficient.
2. A target rational function has no pole on a divisor where `Q` is
   nonconstant.  Independent target values assemble by partial fractions,
   and the normal-domain/height-one intersection then gives a polynomial.
3. For squarefree `U`, reduction of `J(U,W)=mU` on each irreducible factor
   makes `W` a derivation constant on a one-variable function field, hence a
   scalar over algebraically closed `k`.  Equal component scalars are
   equivalent to divisibility by `U` and therefore to a constant-Jacobian
   mate.  The blowdown equivalence and function-field degree are lossless.
4. Section 3 explicitly assumes `k(X,T)=k(U,W)`.  This gives
   `ker J(-,Q)=k(Q)` for the dressed target, the explicit rational primitive,
   and the complete vertical-torsor proof, exactly as the C3 realization
   `(5)--(8)` confirms.

There is one proof-scope mismatch on the current promoted truth surface.  The
YAML `status` states the nonconstant-dressing no-go without birationality,
while Section 3 introduces `(21)` before stating it.  The stronger status is
in fact correct and does not need to be weakened.  It follows immediately
from Section 2, using no rational-centralizer claim:

```text
Q=U phi(W)  gives  J(Q,W)=mQ.
If Q is smooth, Q is squarefree.
Hence a polynomial mate for Q would, by (13), equalize Spec_0(Q,W).
But Q=0 has old components U_i=0 with W=c_i and appended components
W-r=0 with W=r for every root r of phi.
Smoothness gives phi(c_i)!=0, so c_i!=r: the spectrum is not equalized.
```

This proves the no-go for every squarefree nonconstant `phi`, whether or not
`k(X,T)=k(U,W)`.  Birationality is needed only for the stronger conclusions
that `P0=-W/(mQ)` has complete torsor `P0+k(Q)`, that the mismatch is the
full vertical regularization obstruction, and that the constant-dressing
positive boundary is an automorphism.  The unresolved counterexample target
is therefore the **constant-dressing**, synchronized, nonbirational boundary,
not a nonconstant dressing.

Exact required addition before the current sentence “In addition to `(10)`,
assume `(21)`”:

```text
The polynomial no-go below does not require birationality.  Indeed
J(Q,W)=mQ, and smoothness makes Q squarefree.  Its zero-component spectrum
contains every old value c_i and every root r of phi; smoothness gives
c_i!=r.  Section 2 therefore forbids a polynomial mate.  We now add (21)
only to identify the complete rational-mate torsor and its vertical
principal coefficients.
```

Commit `b48410210` promoted the theorem but did not reconcile this
status/body scope mismatch.  No asserted mathematical conclusion is false,
and no change is needed to the blowdown equivalence.  After the proof
addition, the independently hostile-audited promotion is justified.  The
normalized C3 family should be recorded as its exact nonconstant-dressing
control.

## Exact replay

Tracked exact artifacts:

```text
04-computation/jc2_c3_vertical_equalizer_bridge_20260823.py
05-knowledge/results/jc2_c3_vertical_equalizer_bridge_20260823.out
```

Reproduce with

```powershell
python3 -B 04-computation/jc2_c3_vertical_equalizer_bridge_20260823.py
python3 -B -O 04-computation/jc2_c3_vertical_equalizer_bridge_20260823.py
```

Normal, optimized, and frozen stdout match byte-for-byte.  The assertion-free
probe performs 56 exact checks, including all identities above, both
birational inverses, the log-canonical spectrum/dressing identities, three
signed parameter controls, and exact rank gaps for polynomial mates through
total degree seven.  The bounded systems are hostile controls; the torsor and
valuation argument is all-degree.

```text
jc2_c3_vertical_equalizer_bridge_20260823.py   SHA256 0aeb85d302dc7e1c2b3735257c00e5987dfef61d5edf66124a20b4bd9ba4d152
jc2_c3_vertical_equalizer_bridge_20260823.out  SHA256 dfe887647424a07dc86b726b2878e0725cc0694b4a8cb8d3fd10d874cd134e2c
CHECKS=56
RESULT: PASS
```

The audit is recorded here as a reflection; no theorem namespace is promoted by this file.
