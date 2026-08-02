# Frontier portfolio synthesis -- 2026-08-02

**Status:** CURRENT SNAPSHOT / NOT A PROOF SOURCE.  Follow theorem links and
their literal statuses.  Corrections and promoted canon override this
synthesis.

## Executive map

| lane | inherited frontier | new result from this session | exact gain | still open |
|---|---|---|---|---|
| LRC(14) | direct `<=6`; projected caps `k=2:1579`, `k=3:229`; reflected sufficient-certificate residual `561` | [THM-3126](../01-canon/theorems/THM-3126-six-z229-prefix-divisor-universal-complete-cell-closure.md), proved and hostile-audited | six explicit `z1=229` prefix families close for every possible fourth tail; `792=744+48` divisor proof | THM-3113 screen exhaustion and `z1=228`; proved cap remains `229`; arbitrary `k<=1`; LRC(14) |
| planar Jacobian | `JC(2)`/`DC(2)` open; nonsplit response reduced to `F=VG^2`, `2VG'+V'G=2`; one heptic passport solved | [THM-3123](../01-canon/theorems/THM-3123-heptic-e3-remaining-accessory-classification-and-s7-monodromy.md), proved and hostile-audited | remaining two heptic accessory algebras are reduced length six; four new unmarked maps, all `S7` | Faber fluxes and Keller-chart entry; `JC(2)` |
| factorial conjectures | `FC(3)`/`SFC(3)` open; univariate slot atlases and product-Gamma width three proved in scope | [MISTAKE-350](../01-canon/MISTAKES.md), [THM-3124](../01-canon/theorems/THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census.md), [THM-3125](../01-canon/theorems/THM-3125-monomial-ray-first-window-factorial-closure-in-three-variables.md), proved and hostile-audited | corrects dimension/slot types; gives an `O(N)` quadratic moment recurrence and exact starts `1..200`; proves a genuine `FC(3)` monomial-ray subclass from product-Gamma | all resonant starts; translated rays; full `FC(3)`/`SFC(3)` |
| tournaments | scalar substitution product iff quotient transitive; exact cyclic kernel was open | [THM-3121](../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md), proved and hostile-audited | exact all-quotient path-cover transfer; `C3` reduces to three diagonal sums; `H(C3[C3,C3,C3])=3159` | computing arbitrary factor profiles remains hard; classify other low-state quotient kernels |
| integer sequences / closed forms | support and multiplicity must be separated; finite fitting is not proof | THM-3121, THM-3124, THM-3126, all proved | rational walk series, linear recurrence, and divisor-support bulk-plus-finite-correction laws arise from mechanisms | characterize when rational/multivariate transfers collapse to low-dimensional recurrences |

## Inheritance pass

### LRC(14)

- Closest proved mechanism: THM-2941's lossless projected implication plus
  THM-2984's sharp translated danger-band capacity.
- Canonical hostile: at denominator `d=2`, density alone ties the capacity;
  a same-size all-even carrier can fail although each actual carrier has full
  residue support.
- Corrected near miss: MISTAKE-347 invalidates the reflected split tails by
  representing one level ordering twice.  No later reflected cap is live.
- Least-used sidecar: the whole divisor profile
  `d -> |pi_d J|`, not one terminal denominator.

### Planar Jacobian

- Closest proved mechanism: THM-2796's balanced response normal form and
  THM-2840's remainder-ideal classification.
- Canonical hostile: an integer passport does not determine a map; the marked
  Nielsen orbit is load-bearing.
- Corrected near miss: degree-one inertia can be the identity; a
  distinguished-cycle shortcut is unsafe.
- Least-used sidecar: the marked cyclic chord/Nielsen class, which separates
  the three `(3,2,1,1)` unmarked maps.

### Factorial conjectures

- Closest proved mechanisms: THM-2824 for arbitrary univariate three-slot
  first windows and THM-3107 for anchored product-Gamma width three.
- Canonical hostile: `1+x` is a legal nonhomogeneous `FC(1)` input and has no
  single polar degree.
- Corrected near miss: by MISTAKE-350, `m` in `FC(m)`/`SFC(m)` is ambient
  dimension; `N(f)` is support size and window length.  Thus `1+x+x^2` is a
  three-slot `SFC(1)` input, while `1+x+y+z` is an `SFC(3)` input with four
  slots.
- Least-used sidecar: the boundary value `q(0)=a`, which forces the resonance
  in the quadratic recurrence.

### Tournaments and sequences

- Closest proved mechanisms: THM-1862's order-join product and THM-1975's
  path-cover profile as the required substitution sidecar.
- Canonical hostile: scalar `H(S_i)` cannot see how many maximal block runs a
  cyclic quotient revisits.
- Corrected near miss: a finite list of sequence terms does not prove a
  closed form (MISTAKE-342).
- Least-used sidecar: factorial-normalized component counts
  `F_T(r)=r! pc_T(r)`.

## New exact mathematics

### 1. Tournament path-cover substitution is a walk-content contraction

For quotient tournament `Q`, set `X=diag(x_i)` and

```text
W_Q(x)=1^T X(I-A_QX)^(-1)1.
```

If `pc_T(d)` counts spanning covers by `d` unordered directed paths, then

```text
pc_(Q[S_i])(d)
 = sum_c ((prod_i c_i!)/d!) [x^c]W_Q(x)^d prod_i pc_(S_i)(c_i).
```

The proof splits a global path cover at maximal block runs.  The factorials
assign distinct factor components to quotient-word occurrences; division by
`d!` is legal only after this decoration makes the permutation action free.

For `C3`, a quotient walk is periodic.  Its positive coefficient is `3` on
the diagonal, `1` when `max(c)-min(c)=1`, and `0` otherwise.  Consequently
the Hamiltonian count is a diagonal plus two nearest-shell convolutions and
is linear in the supplied factor-profile length.  This is an efficient
closed transfer, not a finite fit and not a general `#P` collapse.

### 2. The two missing heptic accessory strata are zero-dimensional

Put

```text
D=x^a(x-1)^b(x-lambda)(x-mu),
E=(1/7) T D'/D,              u=lambda+mu, v=lambda mu.
```

Dividing `D=S E^2+R` and using `E|D'` makes `R'/E` linear, so constant
remainder is exactly two equations.  The theorem obtains

```text
(4,1,1,1): v=(8u^2+9u+8)/7,
             100u^3+244u^2+237u+44=0;

(3,2,1,1): v=(24u^2-16u-16)/21,
             75u^3-89u^2-31u+61=0.
```

Both cubic discriminants and every collision/critical/accessory/resultant
gate are nonzero.  Reconstruction gives the balanced ODE exactly.  The first
passport has one unmarked map and the second three; explicit transpositions
with the seven-cycle prove all four monodromy groups are `S7`.  Together with
THM-2840, this is the proved five-map abstract heptic accessory atlas, not a
planar Keller classification.

### 3. Quadratic factorial moments have a boundary-forced recurrence

For `q(t)=a+bt+ct^2`, `M_n=L(q^n)`, and `Delta=b^2-4ac`, integration by parts
gives

```text
M_(n+1)=a^n(a+(n+1)b)
       +2(n+1)(2n+1)c M_n
       +n(n+1)Delta M_(n-1).
```

With `M_0=1`, `M_1=a+b+2c`, this computes `M_0,...,M_N` in `O(N)` ring
operations and constant recurrence state.  If three consecutive moments at
start `r` vanish and `abc!=0`, then necessarily

```text
b/a=-1/(r+2).
```

Thus a fixed exact quadratic has at most one possible bad window.  At the
resonance, badness is a common-root problem for two explicit polynomials in
`u=c/a`; two prime-field gcd certificates exclude every `1<=r<=200`.

### 4. A genuine sparse `FC(3)` family comes from Gauss multiplication

For nonzero `v in N_0^3`, `s=X^v`, and `d>=1`, THM-3125 considers

```text
F=A+B s^d+C s^(2d).
```

Along this ray,

```text
L_3((s^d)^n)=prod_i(d v_i n)!
            =G^n prod_(i:v_i>0) prod_(k=1)^(d v_i) (k/(d v_i))_n.
```

The geometric factor is absorbed by `t -> Gt`; the remaining weights are a
finite product of positive Gamma layers.  Proved THM-3107 then makes the
first three moments detect `A=B=C=0`.  The map preserves first-window
nullity but loses translated rays, off-ray supports, and shifted windows.

### 5. Six LRC prefix carriers have universal divisor escape

For each of six explicit `z1=229` prefixes, reconstruct its complete-cell
carrier `J subset Z/LZ`.  For every divisor `d>1` of its ruler,

```text
|pi_d J| > ceil(d/7).
```

Density proves `744` of the `792` carrier/divisor pairs.  The uniform tail
starts at `d=25`; the `48` tied small-modulus cases have full residue support,
with moduli in `{2,3,4,5,8,9,10,11,15,22}`.  THM-2984 then forces escape from
every translated strict danger band.  This closes the six prefix families
unconditionally once named; it does not prove that candidate THM-3113 has
found every `z1=229` residual.

## The common procedural object: quotient plus restoring sidecar

| source quotient | information lost | target predicate | minimal sidecar |
|---|---|---|---|
| LRC carrier `J -> |J|/L` | residue placement | escape every translated band | full divisor residue support at density ties |
| JC `(lambda,mu)->(u,v)` | order of simple poles | marked map/Nielsen class | one root of `z^2-uz+v` or marked chart |
| JC response data `->(F,V,G)` | source prefix and flux equations | Keller-chart entry | `A_src,B_src,d,s`, Faber bank, remaining fluxes |
| homogeneous `f -> g=f|Delta` | interactions among different radial degrees | full FC nullity | the complete nonhomogeneous radial-degree profile |
| tournament factor `S -> H(S)` | number of path components | cyclic substitution count | factorial path-cover profile `r!pc_S(r)` |
| sequence values `-> support` | multiplicity/collisions | indexed sums and operation response | collision tax or labelled fibres |

This is not an analogy claiming the objects are the same.  It is a reusable
audit question: identify the consumer, the quotient kernel, and the cheapest
sidecar that restores the consumer's predicate.

## Next decisive work

### Anchor -- LRC(14)

1. Give THM-3113 a full immutable independent screen replay.  Only that can
   connect the six-prefix closure to an exhaustive `z1=229` statement and
   validate the claimed `z1=228` zero residual.
2. Repair MISTAKE-347 at the assignment level: fix one ordered label pair and
   reverse only the original ratio inequality.  Re-audit all `561` reflected
   bodies before restoring any cone suffix.

### Niche -- factorial sequences

1. Seek a symbolic all-`r` proof that the resonant pair
   `P_(r,r),P_(r+1,r)` is coprime.  The recurrence says exactly where to look;
   more unconstrained coefficient search is lower value.
2. Identify whether the resonant polynomials form an orthogonal/Appell family
   with a Sturm, interlacing, or contiguous-relation certificate.
3. Extend the monomial-ray transfer first to translated Gamma layers with an
   honest tilt sidecar; do not infer translation invariance from the anchored
   theorem.

### Wildcard -- planar Jacobian

1. The independent hostile audit is complete; stop expanding the accessory
   atlas temporarily and attack
   the actual chart-entry map: retain `A_src,B_src,d,s`, form `M=SET`, and test
   `phi_Q=0`, `psi_Q in C`, `A_src K_Q=lambda M`.
2. Treat `S7` dessin monodromy and quartic `S4/C3` cofactor monodromy as
   differently typed covers unless an intertwining map is constructed.

### Tournament / closed-form lane

1. Classify quotient tournaments whose rational `W_Q` has coefficient
   support on finitely many bounded-width affine diagonals; these are the
   candidates for `C3`-style low-state transfers.
2. Compose THM-3121 recursively on modular decomposition trees and measure
   state width by the necessary path-cover profile, not vertex count alone.
3. For every proposed sequence formula, derive it from the rational kernel,
   recurrence, or divisor map and preregister the first hostile out-of-sample
   term.  Finite interpolation remains evidence only.

## Reproduction

```bash
python3 04-computation/tournament_path_cover_walk_kernel_thm3121.py
python3 -O 04-computation/tournament_path_cover_walk_kernel_thm3121.py

python3 04-computation/jc_heptic_e3_remaining_accessory_classification_thm3123.py
python3 -O 04-computation/jc_heptic_e3_remaining_accessory_classification_thm3123.py

python3 04-computation/factorial_quadratic_recurrence_census_thm3124.py
python3 -O 04-computation/factorial_quadratic_recurrence_census_thm3124.py

python3 04-computation/factorial_monomial_ray_fc3_thm3125.py
python3 -O 04-computation/factorial_monomial_ray_fc3_thm3125.py

python3 04-computation/lrc14_j7_k3_z229_six_prefix_divisor_universal_support.py \
  --output /tmp/thm3126.out
python3 -O 04-computation/lrc14_j7_k3_z229_six_prefix_divisor_universal_support.py \
  --output /tmp/thm3126.opt.out
cmp /tmp/thm3126.out /tmp/thm3126.opt.out
```

The exact outputs and hashes are recorded in the linked theorem files.  A
matching computation verifies the encoded finite claim; status promotion
also requires the theorem-level type, implication, and scope audit.
