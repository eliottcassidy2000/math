# Independent audit: period-three global index

19 October 2026. **PASS.** Independent symbolic and proof audit of the final
staged period-three global-index derivation. All requested scope corrections are present in the audited
final prose.

Audited proof SHA256:
`00f8d1d53752c47617d090321f1fae0962043eb2c21783c0f30ddee0d1131207`.
Audited author-verifier SHA256:
`ea17e8ad2d44f1804078d9016615af4a3bbfe60bd95b24490d0aa85697e85023`.

## Findings

The main calculation is correct:

- For a polynomial return map, the sum of finite fixed-point residues is
  zero. Distinct orbit points, rather than root multiplicities, determine
  the cycle weight. At the distinguished parameter this gives
  `6/49 + 3 I = 0`, hence `I=-2/49` and the simple-parabolic iterative
  residue `1-I=51/49`.
- The direct independent local calculation gives
  `a=q(P')²`, `b=q'(P')²+q P'P''` at roots of `P`, with
  `gcd(a,P)=1` and `49b+2a²=0 mod P`. Thus the residue sign is independently
  confirmed by `I=b/a²`.
- `Res(q,Φ₃)=D(c)=16c²+4c+7`. At `c=-7/4`, `Res(q,P)=+7` and
  `D=Res(q,P)²=disc(P)=49`. The oriented cyclic gaps prove the squared
  identity. The supplied scaling counterexample correctly prevents
  treating the polynomial discriminant as a conjugacy invariant. The final
  clarity sentence correctly states that the underlying number field is
  unchanged, with field discriminant 49; only the scaled monic polynomial
  has discriminant 3136. Removing that exact added sentence reproduces
  the previously audited proof hash, confirming no other proof changes.
- In the trace coordinate, `A(v)=v²+v+7`,
  `Res(q,P_v)=A(v)`, `disc(P_v)=A(v)²`, and `1-μ(v)=v A(v)`.
  Consequently `D(c(v))=A(v)A(-v)` and the two per-cyclepoint indices sum
  to `-2/D(c(v))`.
- At `c=1/4`, the removable fixed-cluster value and the actual single-point
  index of the third iterate both equal `2/3`. At `D(c)=0`, fixed and
  primitive-period-three points collide; the separated rational formula
  must not be interpreted as the actual finite index of that one point.

Three wording corrections were requested and applied: specify that the general map
is polynomial; require exactly one multiple cycle before solving for its
individual iterative residue using only the simple-point budget; and
distinguish local-index conjugacy invariance from the finite-point
polynomial budget. The latter is not an arbitrary spherical-coordinate
invariant assertion.

## Local unfolding and first passage

Independent exact reduction gives
`α_r=∂_c f_c³(r)|_(c=-7/4)=2r²-11/2` and `a_r α_r=49`.
Thus in the fixed chart `z=a_r(x-r)`, with `ε=c+7/4`, the return map is
`g_ε(z)=z+z²+49ε+O(z³,εz,ε²)`.

For sufficiently small fixed positive `η`, and `ε→0+`, the first passage
from `z=-η` to `z≥η` satisfies `sqrt(ε) N_ε→π/7`. Here `N_ε` counts
third-return iterates. The corresponding base-iterate count is `3N_ε`.

To check the discretization and endpoints, take `0<δ<η` and then
`ε≪δ²`. On the local interval, the map is increasing and its displacement
is comparable to `z²+49ε`. Entry into the inner interval lies in
`[-δ,-δ+O(δ²+ε)]`; its exit lies in `[δ,δ+O(δ²+ε)]`.
The clock `T_ε(z)=atan(z/(7sqrt(ε)))/(7sqrt(ε))` changes by
`1+O(δ)` on every inner step. The endpoint clock difference multiplied
by `sqrt(ε)` tends to `π/7`, while each outer segment requires only
`O_δ(1)` steps. Sending `ε→0+` and then `δ→0` proves the limit.
This is a local passage theorem, not a critical-orbit or global-escape
claim.

## Provenance and executable evidence

The cyclotomic model, discriminant 49, and reversed-doubling/tripling
interpretation are inherited from the September 17 arithmetic-braids
geometry note and October 4 difference-visibility note. Their attribution
was checked directly. The explicit action `t→1-t-t²` agrees with
`ζ₇→ζ₇³`, including orientation.
The rational neutral-cycle exclusion is also correctly attributed to
Lemma R1 of the October 15 quadratic-cycle-phase note; that source and
its valuation argument were checked directly.

The linked [Milnor revised preprint](https://abel.math.harvard.edu/archive/118r_spring_05/docs/milnor.pdf#page=63)
was verified: Theorem 9.2 is on PDF page 63, and Problem 9-1 on PDF page 65.
Their section numbering belongs to the September 1991 revised preprint.

Run `python audit_period_three_index.py` and
`python -O audit_period_three_index.py`. Both pass all 32 independent
exact controls and produce byte-identical output. Failures use explicit
exceptions, so optimization does not disable the controls. The author's
final 62-control verifier was independently replayed in both modes with
byte-identical passing output. The final local first-passage proof,
including its monotonicity and endpoint overshoot treatment, was checked
independently. No unresolved mathematical or provenance issue remains
within this audit's stated scope.
