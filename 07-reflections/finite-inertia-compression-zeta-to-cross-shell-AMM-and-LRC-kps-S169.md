# Finite inertia compression: from the 2026 zeta preprint to cross-shell AMM surgery

**Session:** kps-S169, 2026-08-12  
**Status:** synthesis with PROVED canon consumers THM-3337, THM-3338, and
THM-3340; the zeta paper remains an EXTERNAL PREPRINT CLAIM + FORMAL ARTIFACT
under partial local audit.

## 1. Inheritance pass

- **Closest proved mechanism:** THM-509's two-layer moment/integrality
  dichotomy, with THM-441's correlation-adjoint language.
- **Canonical hostile example:** THM-3285's 169 LRC middle origins have
  middle completions but empty outer co-support. Any compressed certificate
  that forces a physical owner there has forgotten the decisive sidecar.
- **Corrected near miss:** a positive Gram matrix can be full rank while
  completely missing integer Bernstein capacity. Selecting the columns
  `p^i q^(d_i)` gives a unitriangular minor for every AMM degree profile, so
  raw positive index is maximally blind.
- **Least-used sidecar:** indefinite block signature. A negative and positive
  defect may each be useful before the final quotient cancels them.

The live concept board was:

```text
finite compression | block inertia | cross-shell defect | cyclic rotation
owner-phase reflection | Macaulay rank | integer Bernstein capacity
```

## 2. What the external paper actually does

The linked 35-page paper is Claude, *More than two thirds of the zeros of the
Riemann zeta function lie on the critical line* (10 August 2026). It claims:

```text
simple on-line zeros / all zeros       >= 2/3,
distinct zeros / all zeros             >= 5/6,
optimized constants                    0.6725 and 0.83625,
primitive Dirichlet L-function analogue.
```

Its reusable engine is not generic spectral analysis:

1. compress Weil's Hermitian form to a finite Gabor family;
2. retain the functional-equation involution;
3. represent each on-line zero by a positive rank-one fixed block;
4. represent each off-line pair by a signature-`(1,1)` hyperbolic block;
5. compute prime-side trace and Frobenius-square asymptotics;
6. turn these into positive index via a Hermitian rank--trace inequality;
7. use pullback inertia to force fixed blocks.

The equality sentence in Lemma 3.2 has an upstream-reported scale typo: its
model should be `P=(c/2)Pi_1`. The inequality and applications are unchanged.
The formal companion is `anthropics/zeta-23-lean`. We inspected its source
tree and audit surface, but the local Lean toolchain download reset/timed out;
no local kernel-build claim is made.

## 3. The transfer contract

| Target | Map | Preserved predicate | Lost information | Verdict |
|---|---|---|---|---|
| AMM 12592 | finite Hamming-layer compression | total layer deficit and its sign | integer capacity under a raw Gram quotient | methodological transfer succeeds; raw spectral transfer fails |
| LRC(14) | phase reflection `b -> 2a-b` on `C_13` | one fixed phase plus six involution pairs | address, ancestry, outer co-support, current, exit margin | conditional exact test only |
| FC(3) | Macaulay map `Phi_D`, Hermitianization `Phi_D Phi_D*` | fixed-degree ideal rank | saturation, syzygy location, projective components | legitimate weak screen, not a proof engine |

The strict lesson is: compression helps only when the quotient retains the
coordinate that distinguishes the target from its hostile pair.

## 4. AMM: the paper's lesson became literal mathematics

### 4.1 Horizon eight: first cross-shell escape

THM-3032 proved `T(4)>=6` only for shell-balanced rules. THM-3337 abandoned
that quotient. Its one extractor has profile

```text
(2,3,5,5,7,7,8)
```

on critical values one through seven. The shells `{2,3}` and `{4,5,6,7}`
have nonzero exact opposite Hamming-layer deficits, so their global sum
bisects every length-eight layer. This proves `T_opt(4)=5`.

### 4.2 Horizon sixteen: one donor and two opposite shells

THM-3338 gives one extractor with profile

```text
(2,3,4,5,15,7,8,9,10,11,12,13,14,15,16).
```

Only `n=5` is delayed. All other `n<=15` stop at their first disagreement.
The shells `{1}` and `{2,3}` remain balanced, while the doubled deficits of
`{4,...,7}` and `{8,...,15}` are exact opposites. Exact enumeration of all
`65536` words and a separate coefficient audit pass. This one rule reaches
the floor simultaneously at every `n<=15` except five.

An independent audit found a complement-equivariant representative with the
same deadlines and a quadratic Boolean description of the floor rows. That
is useful normalization, but the canonical integer witness is already exact.

### 4.3 All pointwise values: cyclic rotation replaces MILP

The finite data exposed a much simpler theorem. Fix dyadic `M`. For
`2<=n<M`, color both initial-bit orientations heads when `n` is even and
tails when `n` is odd, stopping at flip `n+1`. Delay only `n=1` to flip `M`.

In a weight-`w` layer let

```text
A_w = C(M-2,w-1),
E_w = # even-critical words,
O_w = # odd-critical words with n>=3.
```

Left cyclic rotation bijects the odd words with the even words ending in
their initial bit. The unmatched even words end in the opposite bit. Rotating
such a word left by `n-1` places injects it into the `n=1` layer: the terminal
run of the image has length exactly `n-1`, so the map is reversible. Hence

```text
0 <= E_w-O_w <= 2A_w.
```

Lucas parity makes the defect even. Therefore

```text
r_w=A_w-(E_w-O_w)/2=C(M,w)/2-E_w
```

is an available integer count in one delayed donor layer. Selecting that
many `01` words bisects every Hamming layer. Splice the finite prefix into
THM-2225's unchanged constant-prefix tail.

THM-3340 consequently constructs, for every dyadic `M`, one fair extractor
with

```text
T_M(1)=M,       T_M(n)=n+1 for 2<=n<M.
```

Choosing `M>n` proves

```text
T_opt(n)=n+1 for every positive integer n.
```

This is pointwise: the method may depend on `n`.

### 4.4 Closed discrepancy factorization

The independent polynomial proof is stronger. With `m=M-2`,

```text
F_M(x)=sum_(n=2)^(M-1)(-1)^n(1+x^(n-1))(1+x)^(M-n-1)
      =A_m(x)+x^m A_m(1/x),
A_m(x)=((1+x)^m-1)/(x+2)
      =x sum_(i=0)^(m/2-1)(1+x)^(2i).
```

If `a_k=[x^k]A_m`, then the discrepancy coefficient is
`f_k=a_k+a_(m-k)`. Hockey-stick bounds give

```text
0<=f_k<=C(m,k),
```

and dyadic parity makes `f_k` even. Thus the donor repair actually satisfies

```text
C(m,k)/2 <= r_k <= C(m,k).
```

The recurrence

```text
2a_k+a_(k-1)=C(m,k),       1<=k<=m,
```

is a new finite-state handle on the uniform problem.

### 4.5 The donor can move, but its full-bit leak is exact

THM-3343 executes the first handoff.  For powers of two `d<M`, take the
annulus of length-`M` words whose first `d` bits are equal but which are not
constant.  Color `n>d` alternately by the parity of `n-d`.  Rotation by one
pairs every negative branch with a positive branch ending in its initial bit;
rotating each unmatched positive branch by `n-d` injects it into the
critical-`d` donor.  The donor repairs each Hamming layer exactly.

The annulus is composition-bisectable **if and only if both endpoints are
dyadic**:

```text
(1+x)^(M-d)(1+x^d)=1+x^M              in F_2[x].
```

Using consecutive endpoints first gives one uniform extractor

```text
T(n)=n+1 off the powers of two,       T(2^r)=2^(r+1).
```

THM-3344 then splits the two orientations at the upper boundary.  This
cancels `R(-1)` and lets the donor stop one bit earlier:

```text
T(1)=2;  T(n)=n+1 off powers of two;  T(2^r)=2^(r+1)-1 for r>=1.
```

The gain is exactly one bit in the whole shell-local floor-interior class.
Once `R(-1)=0`, the derivative receives magnitude `M-2` from the last row
and at most two from the penultimate row, so `R'(-1)!=0`.  No choice of
orientation signs can create a factor `(1+x)^2`.  Moving the donor therefore
concentrates the factor-two asymptotics on a zero-density spine but does not
reduce the worst slope.

## 5. What remains open in AMM

Pointwise difficulty is gone, and a simultaneous shifted donor now exists.
The remaining problem is to **split or smooth** its sparse factor-two cost.
THM-3344 proves that sign changes confined to floor-stopped interior rows can
save at most one additive bit; a slope improvement must let residual cross an
annulus boundary or spend slack on some interior rows.

The sharp next idea is a **donor handoff automaton**. Interpret
`A_m(x)` and its reciprocal as the two residual channels after a dyadic
prefix surgery. Instead of deciding the whole donor layer, defer one channel
to the next scale. The recurrence `2a_k+a_(k-1)=C(m,k)` is exactly the
two-state carry relation to compare with THM-3009's golden capacity wall.

Cheap probes:

1. split one donor across two or more future annuli and retain the signed
   residual channel instead of closing it immediately;
2. allow a sparse set of interior rows one extra bit and test whether their
   derivative contributions cancel the `M-2` top-row leak;
3. test whether those residual factors close under two states or require a
   growing state space;
4. derive the spectral radius of the exact carry matrix and compare it with
   `phi` and `C_arch=1+log_5(phi^2)`;
5. keep the pointwise/global quantifiers separate throughout.

Exploratory one-donor MILPs at `M=8,16,32` admit donor positions through
`4,5,7`, respectively, when every other row is forced to the floor. This is
only numerical discovery evidence; the `M=64` floating model became
unreliable because binomial coefficients exceed safe scaling. It is not a
canon claim.

## 6. LRC and FC transfers that survive

### LRC owner-phase reflection

For an owner candidate `a in C_13`, reflection `iota_a(b)=2a-b` has one fixed
phase and six pairs. A native Hermitian form can therefore have one positive
fixed block and six signature-`(1,1)` blocks. But a valid matrix must retain
literal common outer support at the same address and ancestry. THM-3285's
169 hostile middle origins are the mandatory first control: they have middle
phase completions but no outer co-support. Any positive-index certificate
that fires there is semantically invalid.

That hostile control has now been executed exactly.  On the three retained
addresses, the typed support matrix is

```text
0_169 direct-sum J_r direct-sum 0_169,
r(a,b)=(a,2a-b),
inertia=(91,78,338),       trace=13,       Frobenius^2=169.
```

The nonzero inertia is entirely the orbit count of the imposed involution:
13 fixed points and 78 two-cycles.  In both certified fields, exhaustive
shift/character searches find no finite-Gabor covariance of either middle
endpoint bank under reflection and no left/right exchange; only the 13 fixed
diagonal points match literally.  The outer semantic co-support is empty.
Therefore the paper's inertia mechanism does **not** currently transfer to
LRC(14).  It supplies a reproducible negative diagnostic, not an owner or row
certificate.  Reproduce with

```bash
python 04-computation/lrc14_c13_reflection_gabor_hostile_sidecar_20260812.py
```

### FC(3) Macaulay rank

For a moment ideal, `G_D=Phi_D Phi_D*` has

```text
n_+(G_D)=rank(Phi_D)
```

exactly. This makes trace/Frobenius stable rank a legal cheap screen at the
degree-29 cyclic-quartic frontier. Small Jacobian probes gave certificates
only in the range one to three because the rows are highly coherent. A full
rank screen would still need a saturation/Hilbert sidecar to imply projective
emptiness. Run it once; retire first-two-moment compression if the stable rank
is far below codimension one.

## 7. Outcome

The external paper did not import a zeta theorem into another problem. Its
successful transfer was a research move: preserve signed blocks until the
correct global compression. That move broke a shell-balanced AMM obstruction,
then exposed cyclic rotation and closed every pointwise deadline exactly.

The remaining high-value frontiers are crisp:

1. AMM: turn the donor recurrence into one uniform multiscale extractor;
2. LRC: run the `C_13` inertia form only with literal outer co-support;
3. FC: test the degree-29 Macaulay stable-rank screen with saturation kept
   explicit;
4. zeta paper: distinguish source/formal audit from independent local build
   and track the Lemma 3.2 equality correction.
