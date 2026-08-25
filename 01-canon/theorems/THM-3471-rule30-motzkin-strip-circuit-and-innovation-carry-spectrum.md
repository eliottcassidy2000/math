---
id: THM-3471
title: "Rule 30 Motzkin source-strip circuit, innovation-carry spectrum, and charged macroblock compiler"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Whole radial source-depth
  strips admit one common quadratic Motzkin-field compiler.  The first three
  forced strips form an exact null circuit only after Green-slack
  scalarization, and a first transverse slack jet detects the circuit.
  Innovation coordinates turn phase translation into an odd-holonomy
  triangular odometer with maximal carry degree and full Walsh support; the
  physical center is a calibrated terminal-arc integral of that current.
  A separate charged Four-Russians-style macroblock construction gives an
  all-n O(n^2/log^2 n) word-RAM upper bound.  None of the three Rule 30
  prizes is claimed.
source: root-rule30-deeper-invariants-20260815
audit: >
  Three independent hostile audits rederived the strip/diagonal algebra,
  carry cohomology and marked calibration, and macroblock cost/cocycle laws.
  Ordinary and optimized replays of both companions match their stored
  outputs byte-for-byte, with no assert-dependent gates: ACCEPT.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
related:
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3459-rule30-ternary-intersection-factorial-truth-lift-and-keller-boundaries
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
external:
  - G. Christol, T. Kamae, M. Mendès France, and G. Rauzy, "Suites algébriques, automates et substitutions", Bulletin de la Société Mathématique de France 108 (1980), 401--419, DOI 10.24033/bsmf.1926 (PUBLISHED; CITED for the finite-field algebraic-series/automatic-sequence equivalence)
  - Stephen Wolfram, "Announcing the Rule 30 Prizes", https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/ (2019; CITED for the problem statements only)
  - Wolfram Rule 30 Prizes, https://rule30prize.org/ (CURRENT OFFICIAL LISTING checked 2026-08-15; active listing and submission status only)
script: 04-computation/rule30_motzkin_strip_carry_thm3471.py
output: 05-knowledge/results/rule30_motzkin_strip_carry_thm3471.out
script_sha256: 796cb70d53abdf351c29bae31f6a564e91258ca7c5308ea698b1dc090ac2af5a
output_sha256: cc1ba49a2b341d020bf897e9dbe8f79473882ccd29fdf06f69813d0032babd1b
secondary_script: 04-computation/rule30_macroblock_uniform_compiler_probe.py
secondary_output: 05-knowledge/results/rule30_macroblock_uniform_compiler_probe.out
secondary_script_sha256: 1ea8362b60c002bba2a5cfe4c54f9d77f75709a559ab1c226b5ea7675c012cc3
secondary_output_sha256: d10f1d06a3d7eaa31fb6e334b1bdf7a1e8035f224f8d1c56ce42f76a7cae2be4
hash_basis: raw bytes
---

# THM-3471 -- Rule 30 Motzkin source-strip circuit, innovation-carry spectrum, and charged macroblock compiler

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem continues the three Rule 30 fronts with one common lesson:
scalar state is often adequate only after the missing transverse or marked
coordinate has been proved harmless.  Here it is not harmless.  Summing a
Green transport grade creates an exact three-strip cancellation, independent
innovation bits hide spectrally full carries, and universal block rank
coexists with a charged subquadratic-in-the-cone simulation.

## 1. Inheritance, board, and conventions

The closest proved mechanism is THM-3468's radial decomposition

```text
c_t=B_t+sum_(u,v>=0) J_(u,v)(t),
B_t=H_1(t-2),
H_d(n)=[x^(n+d)](1+x+x^2)^n.                         (1)
```

Its canonical hostile is already internal: a nonperiodic summand need not
survive an infinite graded sum.  The corrected near miss is therefore not to
seek another nonperiodic scalar, but to retain the grade through the transfer.
The least-used sidecars are Green slack, the physical phase basepoint, and the
charged preprocessing model.

The live concept board is:

1. the ternary/Motzkin Green kernel;
2. radial source depth `u` and transport slack `v`;
3. the slack marker and its Hasse jets;
4. the innovation-cube odometer carry;
5. the marked terminal phase arc; and
6. rectangular space-time macroblocks.

All unlabelled algebraic sums are over `F_2`.  Put

```text
K(n,r)=[x^r](1+x+x^2)^n,
F_s(0)=e_s(0),
F_s(d)=e_s(d)+e_s(-d)                 (d>=1),         (2)
```

with the centered adjacent-`11` bond source `e_s` of THM-3468.  In its
two-slack coordinates,

```text
d=s-u,
n=d+v,
t=u+1+2d+v.                                           (3)
```

## 2. A common quadratic compiler for every whole source-depth strip

For fixed `u>=0`, define the source-edge word and its generating function

```text
alpha_u(d)=F_(u+d)(d),
alpha_u(d)=0 if u+d<2,
A_u(X)=sum_(d>=0) alpha_u(d)X^d.                      (4)
```

THM-3468 proves that `alpha_u` is eventually periodic, with a power-of-two
period dividing `2^(u+1)` after its finite head.  Hence `A_u(X)` is rational.

Retain transport slack with a formal marker `q`, and let

```text
R_u(z,q)
 =sum_(d,v>=0) alpha_u(d)K(d+v,v)
      z^(u+1+2d+v)q^v.                               (5)
```

The unmarked complete `u`-strip is `R_u(z):=R_u(z,1)`.

### 2.1 The diagonal lemma

Let `L(x)=1+x+x^2`, and let `xi` be the unique series with

```text
xi=qz L(xi).
```

The formal constant-term/diagonal identity gives, for every `d>=0`,

```text
sum_(v>=0) K(d+v,v)(qz)^v
 =L(xi)^d/(1+qz).                                    (6)
```

Indeed the general denominator is `1-qz L'(xi)`; in characteristic two,
`L'=1` and minus equals plus.  Set

```text
W_q=z^2 L(xi).
```

Then `W_q=z^2+O(z^3)` is the unique small root of

```text
q^2 W_q^2+(1+qz)W_q+z^2=0,                           (7)
```

and (5)--(6) give the exact whole-strip compiler

```text
boxed: R_u(z,q)=z^(u+1) A_u(W_q)/(1+qz).             (8)
```

At `q=1`, write `W=W_1`.  Every fixed strip, and every finite sum of fixed
strips, lies in the same field

```text
K=F_2(z,W),
W^2+(1+z)W+z^2=0.                                    (9)
```

The valuation word `G(z)=sum_(n>=0)H_1(n)z^n` also lies there:

```text
W=z(1+z)G,
B(z):=sum_(t>=2)B_tz^t=zW/(1+z).                     (10)
```

Since `G` is nonrational by THM-3468, (9) is a genuine quadratic extension.
Its nontrivial involution is

```text
sigma(W)=W+1+z.                                      (11)
```

Thus a finite-strip model is rational, equivalently eventually periodic,
if and only if it is fixed by `sigma`.  This is an exact finite decision
procedure after the finite source periods are compiled.  By the standard
finite-field algebraic-series correspondence, every such fixed finite-strip
model is also two-automatic.  The important point is stronger than the label:
all finite strips share one quadratic clock, so algebraic degree cannot
separate their cancellations.

## 3. The three forced boundary strips form one exact circuit

THM-3468's three forced source grades give, with the head values included,

```text
A_0(X)=X^2/(1+X),
A_1(X)=X^2/(1+X^2),
A_2(X)=X/(1+X^2).                                    (12)
```

Indeed `alpha_0(d)=1` for `d>=2`, `alpha_1(d)=1` for even `d>=2`, and
`alpha_2(d)=1` for odd `d>=1`, while `alpha_2(0)=F_2(0)=0`.

Substitution in (8) gives

```text
R_0=zW_q^2/((1+qz)(1+W_q)),
R_1=z^2W_q^2/((1+qz)(1+W_q^2)),
R_2=z^3W_q/((1+qz)(1+W_q^2)).                        (13)
```

On a common denominator, their sum is

```text
R_0+R_1+R_2
 = z(1+q)W_q^2((1+q)W_q+z)
   /((1+qz)(1+W_q^2)).                               (14)
```

After factoring `zW_q`, the remaining bracket is (7) subtracted from
`W_q^2+(1+z)W_q+z^2`.  Consequently

```text
boxed: R_0(z)+R_1(z)+R_2(z)=0.                       (15)
```

All three series are nonzero, so (15) is a minimally supported three-circuit.
The physically forced first three source-depth strips contribute nothing to
the unmarked center:

```text
boxed: c_t=B_t+sum_(u>=3) R_u(t).                    (16)
```

This upgrades “the three spines are not independent” to an exact boundary-
current cancellation.

### 3.1 Why all three branches cancel

There is a coefficientwise proof which exposes the mechanism.  Extend invalid
Green shells by zero.  At an even target `t=2m`, the even-source part of the
`u=0` strip equals the odd-source `u=1` part by

```text
H_(2k)(2r+1)=H_k(r),
H_(2k)(2r)=H_k(r),                                   (17)
```

while the remaining odd-shell/even-time terms vanish.  At an odd target
`t=2m+1`, put `r=m-k-1`.  The three odd-source terms are

```text
H_k(r)+H_(k+1)(r),
H_k(r),
H_(k-1)(r)+H_k(r).                                   (18)
```

Their XOR is

```text
H_(k-1)(r)+H_k(r)+H_(k+1)(r)=H_k(r+1),               (19)
```

which cancels the neighboring even-source `u=0` term.  Thus (15) is the
local ternary Green recurrence applied to all three boundary branches at
once, not a coincidental generating-function identity.  It is a circuit, not
a tournament and not a Berggren ancestry equivalence.

### 3.2 A single strip realizes the Mahler tail and the missing `1/6`

The inward strip has an especially sharp form:

```text
B(z)+R_2(z)=z^3/(1+z^2),
R_2(z)=z^3G(z^2).                                    (20)
```

The second equality is Frobenius, and the first is exactly the shifted Mahler
equation for `G`.  Coefficientwise,

```text
R_2(t)=H_1(t-2)+[t is odd, t>=3].                     (21)
```

Thus `B` and `R_2` are disjoint and partition the odd targets from time three:

```text
B_t R_2(t)=0,
dens(B)=1/3,
dens(R_2)=1/6.                                       (22)
```

The lone `u=2` strip would turn the valuation backbone into the balanced,
period-two odd-time word.  It satisfies THM-3468's density invoice exactly:
`rho_(R_2)-2rho_(B R_2)=1/6`.  But `R_0+R_1=R_2`, so the same carrier appears
twice in the physical boundary layer and cancels.  Any proof of balance must
recover the `1/6` from the genuine inward tail `u>=3`.

This is the sharp hostile to strengthening THM-3468's finite-rectangle tariff
to complete bounded-source strips: unbounded Green slack at the single depth
`u=2` already cancels the entire nonrational component of `B`.

## 4. The slack marker repairs the invisible circuit at first order

Equation (14) shows more than unmarked cancellation.  Its order at `q=1` is
exactly one, because after division by `1+q` the specialization is

```text
z^2W^2/((1+z)(1+W^2)) !=0.                           (23)
```

Therefore the first Hasse/ordinary derivative in characteristic two is

```text
boxed:
  (d/dq)(R_0+R_1+R_2)|_(q=1)=R_1(z).                 (24)
```

Coefficientwise, (24) says that weighting each event by `v mod 2` changes the
zero circuit into the middle strip.  Since

```text
v=t-1-2s+u = t-1+u mod 2,                            (25)
```

this first jet is a cheap checkerboard source-depth color at a fixed target.
Finer transport information requires higher Hasse jets, or the full marker
`q`.  The exact lesson is that setting `q=1` and then taking the associated
grade do not commute.

## 5. Innovation is odd holonomy, not a free Bernoulli time process

Retain THM-3468's phase readouts and transition current, for `k>=2`,

```text
T_k(h)=b_k(h+k),
T_k(h+1)+T_k(h)=Q_k(h),
xor_(h mod P_k) Q_k(h)=epsilon_k.                    (26)
```

Let `X_k=F_2^(I_(k-1))` be the innovation coordinates below depth `k`, and
let

```text
Gamma:Z/P_k Z -> X_k,
tau=Gamma o (h |-> h+1) o Gamma^(-1).                (27)
```

Pull the current back to the cube:

```text
g_k(Gamma(h))=Q_k(h).                                 (28)
```

Then `g_k` is a one-cochain on the single `P_k`-cycle and

```text
xor_(x in X_k) g_k(x)=epsilon_k.                      (29)
```

If `epsilon_k=0`, there are exactly two primitives

```text
f(tau x)+f(x)=g_k(x),                                 (30)
```

differing by the constant one; the physical owner `T_k(0)=c_k` selects one.
If `epsilon_k=1`, no primitive exists on `X_k`.  Adjoin `z=T_k` and use

```text
tau'(x,z)=(tau x,z+g_k(x)).                           (31)
```

After one base cycle, `z` flips, so (31) is one cycle of length `2P_k`.
This is exactly the next innovation double cover.  Innovation is therefore
the nonzero holonomy class of the transition current; noninnovation is a
coboundary plus one calibrated owner bit.

### 5.1 Every innovation carry has maximal degree and full spectrum

Enumerate innovation depths `kappa_1<kappa_2<...`.  Put `q_1=1` on the
zero-variable cube.  For `r>=2`, at `k=kappa_r`, write

```text
q_r:F_2^(r-1)->F_2
```

for the pullback (28).  Equation (29) says that `q_r` has odd support.  The
top ANF coefficient of a Boolean function is the parity of its truth table,
so

```text
deg(q_r)=r-1,                                         (32)
```

every earlier innovation bit is essential, and its exact deterministic
decision-tree depth is `r-1`.

For `r>=3`, every unnormalized real Walsh coefficient

```text
What_qr(S)=sum_x (-1)^(q_r(x)+S dot x)                (33)
```

is congruent to `2 mod 4`, hence is nonzero.  For nonempty `S`, the pure
character sum is zero and (33) is minus twice an odd signed support sum; for
`S=0`, `2^(r-1)` is divisible by four and the same congruence holds.  The
one-variable `r=2` case is the sharp exception.

Thus the Haar coordinates are independent as **spatial coordinates**, while
phase addition has carries of maximal algebraic degree and full Walsh support.
If the physical phase origin is forgotten, the same dynamical object
`(X,tau)` can be marked at any `h_0 in Z_2`; `Gamma(h_0)` ranges over every
infinite innovation word, while `tau`, its carry functions, Haar law, and all
their invariants stay fixed.  Finite marks are already arbitrary because
`Gamma_K` is bijective.  Those unpointed invariants alone therefore cannot
determine the marked innovation word.  Physical calibration is a load-bearing
sidecar.

## 6. The center is a calibrated terminal-arc integral

The owner is not arbitrary once the light-cone boundary is retained.  For
`k>0`,

```text
T_k(-k)=b_k(0)=0,
T_k(0)=c_k.                                           (34)
```

Telescoping (26) gives the exact marked formula, for every `k>=2`,

```text
boxed: c_k=xor_(h=-k)^(-1) Q_k(h).                    (35)
```

If `r=k mod P_k`, periodicity and (29) reduce it to

```text
c_k=((floor(k/P_k) mod 2) epsilon_k)
    +xor_(h=-r)^(-1)Q_k(h).                           (36)
```

Thus total support parity, balance, or full Walsh support does not close the
problem: one needs the parity of a moving terminal arc.

There is a fully physical left-front version.  Put

```text
ell_j(t)=a_t(-t+j),       ell_j=0 for j<0.             (37)
```

Changing variables in (35) yields

```text
c_k=
 xor_(0<=s<=k, s congruent k mod 2)
   (ell_(s-1)((k+s-2)/2) or ell_s((k+s-2)/2)).        (38)
```

The prize-density target is therefore a marked diagonal parity of aligned
left-front adjacent-parent OR events.  This is stronger and more actionable than
“typical phase is not the origin,” but it is not yet a discrepancy bound for
that arc.

The light-cone calibration actually selects the basepoint uniquely.  If
`x=Gamma(h_0)` and `f_k` is the depth-`k` readout, then

```text
f_k(tau^(-k)x)=T_k(h_0-k)=b_k(h_0).                  (38a)
```

The boundary conditions `f_k(tau^(-k)x)=0` for every `k>=1` say that the
entire packed edge state at phase `h_0` is the seed `1`.  Through a finite
cutoff `K`, the exact seed period forces

```text
h_0=0 mod P_(K+1).                                   (38b)
```

The periods are unbounded, so the inverse limit forces `h_0=0` in `Z_2`.
The owner is therefore not free; the remaining debt is the ordered current
on the moving arc `[-k,0)`.

Downstream THM-4050 uses the physical light-cone zero to shorten this arc to
`[-floor(k/2)-1,0)`, then stops it at the nearest actual zero. The resulting
address is exactly a marked `1 0^(2r-2)` cylinder. Its Haar law is explicit,
but deterministic temporal discrepancy remains open.

## 7. A charged uniform macroblock compiler

This section fixes a model.  Use a uniform unit-cost random-access word-RAM
with word size

```text
w=ceil(log_2(n+2)),                                   (39)
```

standard binary `n` in one word, unit-cost Boolean operations, shifts,
addition, and addressed word load/store, a constant program, no advice or
oracle, and all runtime preprocessing charged.

For `s,h>=1`, define the raw-cone block map

```text
M_(s,h):{0,1}^(s+2h) -> {0,1}^s                      (40)
```

as the middle `s` cells after `h` Rule 30 steps.  If `s+2h<=w`, enumerate its
table at runtime.  Each entry costs `O(h)` packed word operations, so table
construction costs

```text
O(h 2^(s+2h)).                                        (41)
```

Index the raw input as positions `0,...,s+2h-1`; the table stores positions
`h,...,h+s-1` after `h` steps.  Exterior zeros cannot enter those outputs.
At macro time `rh`, partition the needed future interval into globally aligned
`s`-cell blocks.  Locality makes each output block one lookup on its expanded
`s+2h` predecessor window; a window crossing a machine-word boundary still
uses `O(1)` reads and shifts.  An explicitly zero-padded packed interval
containing the full time-`n` cone prevents artificial boundary contact.
Summing the `O((rh+h)/s+1)` calls over macrosteps and handling the final
fewer-than-`h` rows by ordinary packed updates gives

```text
T(n;s,h)=O(h2^(s+2h)+n^2/(sh)+n/h+n),
S(n;s,h)=O(2^(s+2h)+n/w) words.                       (42)
```

For `w>=16`, take

```text
h=floor(w/8),       s=2h.                             (43)
```

Then `s+2h=4h<=w/2`, the table has `O(sqrt(n))` entries, and

```text
boxed:
T(n)=O(n^2/log^2 n) word operations,
S(n)=O(n/log n) words.                                (44)
```

Small `n` use direct simulation.  Charging `O(w)` bit operations per word
operation gives the weaker accounting `O(n^2/log n)`; no equivalence with a
Turing-machine cost model is asserted.

To emit the entire center prefix, store with each entry the `h`-bit vertical
trace of the first site in its output block.  Align the origin there and read
one trace per macrostep.  Since `s+h=3h<=w`, this has the same asymptotic cost
as computing `c_n` alone.

### 7.1 The time-varying cocycle and rank hostile

The tables obey the exact width-changing cocycle

```text
M_(s,h+k)=M_(s,k) o M_(s+2k,h).                       (45)
```

Time compression invoices a wider predecessor state; it is not iteration of
one fixed finite map.

Left permutivity gives, by induction,

```text
(F^h x)_j=x_(j-h)+G_h(x_(j-h+1),...,x_(j+h)).        (46)
```

After reindexing the block input,

```text
M_(s,h)(x)_r=x_r+H_(r,h)(x_(r+1),...,x_(r+2h)).      (47)
```

For every fixed final `2h` input bits, the prefix-to-output map is triangular
with diagonal one and hence a permutation.  Its graph-incidence matrix is a
`2^s` permutation matrix over every field, and a uniform prefix gives a full
Walsh cube.  Maximal universal block rank therefore coexists with (44): this
is an internal hostile to promoting universal-input rank into a fixed-seed
time lower bound.

Within the explicitly restricted raw-cone rectangular lookup model, an
`L`-bit address satisfies `s+2h<=L`, hence `sh<=L^2/8`; the choice `s=2h`
is geometrically optimal for the number of rectangular calls.  This is not a
general lower bound.  The companion finite audit also finds that for `s=2h`
the cumulative one-seed aligned windows through macro times `0,h,...,t`
cover every table address (including the exterior-zero address) at

```text
(h, first cumulative full-coverage macro time)
 =(1,11),(2,152),(3,1161),(4,6144).                   (48)
```

That is a finite-exact hostile to naive on-demand memoization, not an
asymptotic statement or simultaneous coverage in one row.  At the square
scale `s=h`, an address is literally three adjacent `h`-bit blocks.  The
analogous finite-exact cumulative first full-coverage times are

```text
(h,t)=(1,3),(2,52),(3,210),(4,980),(5,3740),(6,14718). (49)
```

This is a simultaneous three-block saturation control, not an all-`h`
ternary-tree law.

### 7.2 Exact block-CA renormalization

The square scale has a lossless dynamical meaning.  Let `pi_h` group a binary
configuration into consecutive `h`-bit symbols, and set

```text
C_h=pi_h o F^h o pi_h^(-1)                           (50)
```

on the alphabet `A_h=F_2^h`.  Then `C_h` is a radius-one cellular automaton
whose local rule is exactly

```text
M_(h,h):A_h^3 -> A_h.                                (51)
```

The three blocks act simultaneously; they are neighboring parents of one
supercell, not three descendants in a rooted tree.  Equation (47) says that
this alphabet rule is left-permutive in its first parent.

If `G_(m,h)` groups `m` consecutive `h`-symbols into one `mh`-symbol, then

```text
boxed: C_(mh)=G_(m,h) o C_h^m o G_(m,h)^(-1).        (52)
```

This follows immediately from `pi_(mh)=G_(m,h)pi_h`.  Time rescaling is exact,
but the alphabet grows from `2^h` to `2^(mh)`; (52) is a precise
time-versus-state invoice, not a fixed finite renormalization.

Left permutivity also proves that every finite output cylinder of `C_h` has
uniform Bernoulli measure: choose the two right boundary input symbols freely
and solve the remaining inputs right-to-left.  Perfect universal block
Bernoulli/Walsh behavior at every scale therefore coexists with one exceptional
marked seed orbit.  It gives no center-density theorem.

## 8. Preservation, loss, and cross-frontier meaning

| map | preserves | destroys | required sidecar / boundary |
|---|---|---|---|
| `(d,v)->R_u(z,q)` | every fixed-`u` target coefficient and slack weight | source time after summation | `u`, oriented bond source |
| `q->1` | scalar center contribution | transport slack and its Hasse moments | first jet detects the boundary circuit |
| three forced strips -> zero | exact unmarked Green response | which branch carries the Mahler tail | labelled `u` or `q` |
| finite strips -> quadratic field | exact GF and rationality test | unbounded source-depth completion | full `u>=3` tail |
| phase -> innovation cube | cylinders and Haar measure | native cyclic group law | nonlinear odometer `tau` |
| current -> holonomy bit | innovation/noninnovation | primitive/owner when holonomy vanishes | `c_k` or terminal arc |
| carry -> ANF/Walsh spectrum | maximal degree and full spectral support | marked address | physical cone calibration |
| raw cone -> macro table | exact `h`-step local evolution | no fixed homogeneous state map | block width, halo, preprocessing cost |
| binary CA -> `C_h` | exact sampled spacetime and Bernoulli law | fixed alphabet size | alphabet `F_2^h` grows with scale |

The three-strip relation is the Rule 30 analogue of a boundary/Stokes circuit:
all three branches must be retained before scalarization.  It resembles the
simultaneous three-branch bookkeeping in THM-3357, but it is generated by the
Motzkin characteristic equation and preserves no Berggren node ancestry.
The slack jet is analogous to the radial/face sidecars which repair factorial
quotients; it supplies no factorial-conjecture transfer.  The carry cochain is
a lawful owner/current analogy with LRC, but has no circle-covering semantics.
No intrinsic pairwise orientation appears, so tournament language would add
nothing.

## 9. Verification and exact open frontier

The intended exact audits are

```bash
python3 04-computation/rule30_motzkin_strip_carry_thm3471.py
python3 -O 04-computation/rule30_motzkin_strip_carry_thm3471.py
python3 04-computation/rule30_macroblock_uniform_compiler_probe.py
python3 -O 04-computation/rule30_macroblock_uniform_compiler_probe.py
```

The strip/carry companion independently compares centered Rule 30 rows with
ternary Green transport, checks (15), (20), and (24), and audits the finite
carry, cochain, Walsh, and marked-arc controls.  The macroblock companion uses
independent list, packed, centered, and right-edge implementations; its
bounded address-coverage claims are explicitly finite exact.  Universal
claims are the proofs above, not extrapolations from those ranges.

The three sharpened prize-facing problems are now:

1. control the `u>=3` completion while retaining at least the slack jet;
2. bound the marked terminal-arc discrepancy of the spectrally full carry;
3. beat or lower-bound the charged macroblock cocycle in a fixed binary-index
   model, rather than inferring time from universal block rank.

The official listing was last checked on `2026-08-15`; on that dated evidence
the repository treats all three prize questions as open.  This theorem claims
no prize solution and no literature novelty.
