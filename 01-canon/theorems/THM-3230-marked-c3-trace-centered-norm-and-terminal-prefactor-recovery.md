---
id: THM-3230
title: "Marked C3 trace-centered norm and terminal-prefactor recovery"
status: >
  PROVED + VERIFIED-EXACT.  In a tame totally ramified cyclic cubic
  completion K(pi)/K with pi^3=epsilon*t, every primitive element x has a
  nonzero trace-centered part delta=pi*A+pi^2*B.  Its exact norm is
  epsilon*t*A^3+(epsilon*t)^2*B^3.  The two possible leading values are
  incongruent modulo three, so there is no leading cancellation: if h is
  the normalized value of delta and N_h the leading unit of its norm, then
  3 does not divide h and [N_h]=[epsilon]^h.  Thus a supplied local fixed
  sheet of a pure-C3 quartic, through the companion cubic and its centered
  norm, recovers the tame cubeclass.  In the THM-3081/3201 terminal scope,
  [L] alone is Bezout-gauge dependent; the intrinsic class is
  Lambda=[L*theta^(A-1)].  The centered norm recovers [K]=[epsilon]^-1 and
  hence Lambda.  It also gives the exact compatibility
  [q_0]^h[N_h]^m=1 with the graph-quartic leading coefficient.  This repairs
  the 3|m blindness of q_0 after a marked sheet is supplied, but supplies no
  global marking, cross-place gluing, polynomial realization, C3 exclusion,
  S4 exclusion, or JC(2) theorem.
source: jc-marked-c3-sidecar-2026-08-02
depends_on:
  - THM-3081-terminal-toric-residue-parameter-mobius-rigidity-and-autonomous-decoder
  - THM-3201-c3-local-resolvent-splitting-and-matching-newton-gate
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
  - THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan
  - THM-3077-pointed-norm-relative-line-lift-and-relation-carry-obstruction
script: 04-computation/jc_marked_c3_centered_norm_sidecar_thm3230.py
output: 05-knowledge/results/jc_marked_c3_centered_norm_sidecar_thm3230.out
script_sha256: 880ef2eada900a8cf5a9b33ad78ba3e9d7b8d06dbdb9777d41df9068d94d8054
output_sha256: 90b2eac43044968d8309c7cfdc5c65b7ee81ca1e0d70972f1b5bdef4c4eb5b93
hash_basis: LF-normalized bytes
---

# THM-3230 -- a marked cubic orbit recovers the missing terminal cubeclass

**PROVED + VERIFIED-EXACT.**

## 1. Result and information ledger

[THM-3201](THM-3201-c3-local-resolvent-splitting-and-matching-newton-gate.md)
shows that a pure `C3` place splits the local `V4` Kummer packet.  Its
fixed-plus-escaping graph calculation leaves

```text
[q_0]=[K]^m in C(theta)^*/C(theta)^(*)3,                 (1)
```

where `m` is the escaping pole order and `K` is the terminal target-leading
unit from
[THM-3081](THM-3081-terminal-toric-residue-parameter-mobius-rigidity-and-autonomous-decoder.md).
When `3|m`, equation `(1)` is blind to `[K]`.  Moreover, the phrase “the
missing class `[L]`” is not quite intrinsic: the terminal function `L`
changes under a different Bezout chart.

The missing operation is to mark the locally inertia-fixed sheet and then
**trace-center the moving cubic orbit before taking its norm**.  This gives a
second scalar whose first nonzero residue always lies in one of the two
nontrivial `C3` character degrees.  It recovers `[K]` even when `(1)` does
not.  The precise ledger is

```text
source:       a pure-C3 quartic completion plus its local fixed sheet;
operation:    divide off that sheet, trace-center the cubic root, take norm;
preserved:    the tame cubic unit cubeclass;
destroyed:    cyclic ordering and the individual character amplitude;
restored:     the intrinsic terminal prefactor class Lambda;
not supplied: a global/coherent fixed-sheet mark or polynomial realization. (2)
```

Thus this theorem closes a **local recovery problem**, not the global
`C3` branch.

## 2. The trace-zero Kummer normal form

Let

```text
F be a field of characteristic different from three containing zeta_3,
K=F((t)),                         v(t)=1,
L=K(pi),                          pi^3=epsilon*t,
epsilon in F[[t]]^*,              sigma(pi)=zeta_3*pi.  (3)
```

Assume that `epsilon*t` is not a cube, so `L/K` is a totally ramified cyclic
cubic extension.  Normalize the extension value by

```text
w(pi)=1,                          w restricted to K=3v. (4)
```

Every `x in L` has a unique expression

```text
x=C+pi*A+pi^2*B,                  A,B,C in K.            (5)
```

Since the two nontrivial character sums vanish,

```text
Tr_(L/K)(x)=3C.
```

Put

```text
delta=x-Tr(x)/3=pi*A+pi^2*B.                             (6)
```

The three conjugates of `(6)` give the exact identity

```text
N_(L/K)(delta)
 =product_(i=0)^2(zeta_3^i*pi*A+zeta_3^(2i)*pi^2*B)
 =(pi*A)^3+(pi^2*B)^3
 =epsilon*t*A^3+(epsilon*t)^2*B^3.                       (7)
```

This is the cubic analogue of retaining a centered vector before taking its
scalar shadow.  The mixed terms cancel by the character orthogonality, but
the two pure character cubes survive.

If `x notin K`, then `(A,B)!=(0,0)`.  Write

```text
a=v(A),                         b=v(B),
h=min(1+3a, 2+3b),                                      (8)
```

with the absent term assigned value `+infinity`.  The two finite candidates
in `(8)` cannot be equal: one is `1 mod 3`, the other `2 mod 3`.  Therefore
there is **no leading cancellation** in `(7)`, and

```text
w(delta)=h,                  v(N(delta))=h,
h not congruent to 0 mod 3.                              (9)
```

Write

```text
N(delta)=t^h(N_h+O(t)),                  N_h in F*.      (10)
```

If the first term of `(7)` leads, its residue is
`bar(epsilon)*bar(t^(-a)A)^3`; if the second leads, it is
`bar(epsilon)^2*bar(t^(-b)B)^3`.  Since the corresponding residues of `h`
are respectively one and two modulo three,

```text
[N_h]=[bar(epsilon)]^h in F*/(F*)^3.                     (11)
```

Because `h` is invertible modulo three, `(11)` recovers
`[bar(epsilon)]`.  This is stronger than a valuation gate: it recovers the
unit Kummer class which the valuation alone forgets.

The prime degree is useful here.  Every `x notin K` is primitive for `L/K`,
and conversely a primitive `x` makes `(6)` nonzero.  If `x in K`, then
`delta=0`; there is no unit `N_h`, and the theorem deliberately makes no
recovery claim.

## 3. A supplied fixed quartic sheet computes the norm

Let a monic separable quartic over `K` have, at the chosen pure-`C3` place,
the factorization

```text
f(T)=(T-alpha)g(T),
alpha in K,
g(T)=T^3+b_2T^2+b_1T+b_0 irreducible,                   (12)
```

and let `x` be a root of `g`.  The linear factor is the sheet fixed by local
inertia; the other three roots form the moving orbit.  Its mean and centered
norm are

```text
mu=Tr(x)/3=-b_2/3,
N_alpha=Norm(x-mu)=-g(mu)
       =-2b_2^3/27+b_1b_2/3-b_0.                        (13)
```

Thus `(10)--(11)` are computable from the quartic **after the local fixed
sheet has been supplied**.  If

```text
f=T^4+a_3T^3+a_2T^2+a_1T+a_0,
```

synthetic division gives

```text
b_2=a_3+alpha,
b_1=a_2+alpha*a_3+alpha^2,
b_0=a_1+alpha*a_2+alpha^2*a_3+alpha^3.                  (14)
```

Equations `(13)--(14)` use no cyclic ordering of the three moving roots.
They do use `alpha`.  The unmarked global `S4/V4` resolvent has forgotten
which quartic sheet supplies that linear factor.  A chosen valuation has a
unique locally fixed sheet, but carrying this choice coherently between
conjugate divisors or back to the affine source is additional global data.
Nothing in this theorem constructs that global mark.

## 4. The intrinsic terminal prefactor is not `[L]`

Now impose the coordinate-line terminal hypotheses of THM-3081 and the
pure-`C3` graph hypotheses of THM-3201.  Use `theta` for the Mobius residue
coordinate.  At the terminal toric stage write

```text
g,e>=1,                         gcd(g,e)=1,
A*g+B*e=1,
tau=rho^3 K(theta),
lead(U)=rho^(3-e)L(theta).                               (15)
```

Here `tau` is the leading coefficient of the target parameter relative to
the Kummer uniformizer.  In the normalization `(3)`,

```text
tau=bar(epsilon)^(-1),
[K]=[bar(epsilon)]^(-1) in C(theta)^*/C(theta)^(*)3.    (16)
```

All brackets from now on denote cubeclasses in `C(theta)`.

Every other Bezout solution is

```text
A'=A+k*e,                         B'=B-k*g.              (17)
```

THM-3081's value-one coordinate and leading scale transform as

```text
S'=S*Theta^(-k),                  rho'=rho*theta^(-k).   (18)
```

Keeping the physical functions `tau` and `lead(U)` fixed in `(15)` forces

```text
K'=theta^(3k)K,
L'=theta^(k(3-e))L.                                     (19)
```

Thus `[K]` is gauge invariant, but `[L]` need not be.  The exact invariant
replacement is

```text
Lambda=[L*theta^(A-1)].                                 (20)
```

Indeed, `(17)--(19)` change the exponent of `theta` in `(20)` by

```text
k(3-e)+k*e=3k,                                          (21)
```

which is a cube.  This proves rigorously that the missing object is
`Lambda`, not the chart-dependent class `[L]`.

If the Mobius coordinate is

```text
theta=(a_M*u+b_M)/(c_M*u+d_M),
Delta=a_M*d_M-b_M*c_M,
ell(theta)=a_M-c_M*theta,                               (22)
```

THM-3081's autonomous decoder is

```text
K/L=kappa/(3*Delta)*theta^(A-1)*ell(theta)^2.            (23)
```

Every nonzero complex constant is a cube.  Passing `(23)` to cubeclasses
therefore gives the gauge-invariant relation

```text
[K]=Lambda*[ell]^2,
Lambda=[K]*[ell]^(-2).                                  (24)
```

Equations `(11)`, `(16)`, and `(24)` prove the promised recovery.  Choose
`r in {1,2}` with

```text
r*h congruent to 1 mod 3.
```

Then

```text
[K]=[N_h]^(-r),
Lambda=[N_h]^(-r)*[ell]^(-2).                           (25)
```

The terminal prefactor has therefore ceased to be an unknown **local
cubeclass** once the fixed sheet and its first trace-centered orbit
coefficient are retained.  It remains an unconstrained class: `(25)` does
not prove that it is trivial or globally realizable.

## 5. Compatibility with the graph-quartic coefficient

For the graph-quartic pole order `m`, THM-3201 gives `(1)`.  Combining it
with

```text
[N_h]=[K]^(-h)                                           (26)
```

produces the chart-free compatibility law

```text
[q_0]^h*[N_h]^m=1.                                      (27)
```

Equivalently,

```text
q_0^h*N_h^m is a cube in C(theta).                       (28)
```

Since constants over `C` are cubes, `(28)` has the completely executable
divisor form

```text
h*div(q_0)+m*div(N_h) is coefficientwise divisible by 3. (29)
```

One may regard `div(f)` as a finitely supported signed atomic measure on
`P^1_theta`; `(29)` says that every atom, including infinity, has mass zero
modulo three.  This is a simultaneous horizontal multi-place test inside the
chosen vertical completion, not a new vertical-place globalization theorem.

The two lanes are now exact.

```text
3 does not divide m:
  q_0 already recovers [K]; N_h supplies an independent marked compatibility.

3 divides m:
  [q_0]=1, but h is prime to three and N_h still recovers [K] and Lambda. (30)
```

Thus trace-centering repairs precisely the divisibility lane in which the
leading graph coefficient loses the tame class.

## 6. Sharp `3|m` blindness hostile

The loss in `(30)` is attained by a one-parameter local Kummer family.  Let

```text
F=C(theta),                    K=F((t)),
pi^3=epsilon*t,
f_epsilon(T)=T*((T-t^(-1))^3-epsilon*t).                (31)
```

Its four roots are the fixed root zero and the moving orbit

```text
x_i=t^(-1)+zeta_3^i*pi.                                  (32)
```

With `w(pi)=1`, the escaping pole order is `m=3`.  Direct formal depression
of `(31)` gives the linear coefficient

```text
q=t^(-3)/8-epsilon*t.                                   (33)
```

Hence

```text
q_0=1/8                                                  (34)
```

for every `epsilon`: the leading graph coefficient is exactly independent
of the tame unit cubeclass, not merely congruent to it modulo cubes.

On the other hand, the moving-orbit mean is `t^(-1)`, so

```text
delta_i=zeta_3^i*pi,
Norm(delta)=epsilon*t,
h=1,                         N_h=epsilon.               (35)
```

Taking, for example, `epsilon=1` and `epsilon=theta` gives the same
`(m,q_0)=(3,1/8)` but distinct centered-norm cubeclasses, because `theta`
has a simple zero.  This is the minimal hostile proving that a `q_0`-only
argument cannot recover the missing class on the `3|m` lane.

The family `(31)` is a local cyclic-cubic quartic packet, not a global
`S4` polynomial Keller map.  It tests the information map and makes no
physical realization claim.

## 7. What has and has not changed

The local sidecar problem now has a precise answer:

```text
UNMARKED S4/V4 RESOLVENT:
  loses the V4 origin and the locally fixed quartic sheet.

LEADING GRAPH q_0:
  remembers [K]^m and is blind exactly when 3|m.

SUPPLIED FIXED SHEET + CENTERED CUBIC NORM:
  remembers [K]^(-h), with h automatically invertible modulo three.

TERMINAL DECODER:
  turns [K] into the intrinsic prefactor class Lambda.               (36)
```

This does **not** prove any of the following:

1. that an unmarked Keller map canonically supplies one fixed sheet across
   all conjugate divisors;
2. that the local fixed roots glue to a rational or polynomial global
   section;
3. that `Lambda` is trivial, has forbidden divisor, or satisfies the affine
   cofactor realization gate of THM-3064;
4. that an arbitrary Jelonek component is a target coordinate line; or
5. that pure `C3`, `A4`, `S4`, a degree-four Keller map, `JC(2)`, or `DC(2)`
   is excluded.

The next global obligation is sharper than “control `[L]`”: transport the
locally fixed sheet and the signed divisor of `N_h` between conjugate
components, then test the recovered `Lambda` against the true branchwise
inverse-different/cofactor class.  A second unmarked norm cannot perform that
transport.

## 8. Exact evidence

Run

```bash
python3 04-computation/jc_marked_c3_centered_norm_sidecar_thm3230.py
python3 -O 04-computation/jc_marked_c3_centered_norm_sidecar_thm3230.py
```

Both modes byte-match the stored transcript.  The companion checks:

1. the exact `C3` character-product identity in `(7)`;
2. the marked cubic formula `(13)--(14)`;
3. `361` hostile valuation pairs, including negative values, and the
   no-cancellation/congruence statement `(8)--(11)`;
4. `2,233` terminal Bezout-gauge cells for `(17)--(21)`;
5. `672` graph/centered-norm cubeclass compatibility cells; and
6. the exact depressed coefficient and centered norm in `(33)--(35)`.

Every truth-bearing gate uses an explicit exception and remains active under
optimized execution.

**QED.**
