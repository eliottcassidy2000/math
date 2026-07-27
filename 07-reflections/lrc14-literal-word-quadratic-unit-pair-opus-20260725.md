---
source: codex-2026-07-25-literal-word-quadratic-pair
status: >
  PROOF CANDIDATE + VERIFIED-EXACT CONSTANTS + INDEPENDENTLY
  HOSTILE-AUDITED. On the strictly shallower selected-owner return branch
  of THM-2302/2305, one exact terminal word has a literal source handoff
  subset with two genuine same-root-character Fourier atoms whose
  difference is m times the deepest blocker, gcd(m,91)=1. The first
  multiplier bounds are 2,013,431 (strict) and 16,656,658
  (repeated-first, conditional on a named deeper excluded blocker). The
  selected restricted atoms need not survive in the unrestricted owner
  spectrum, so this is not yet THM-2302's marked source edge or LRC(14).
depends_on:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2305-canonical-blocker-word-handoff-hypergraph
related:
  - THM-2286-endpoint-prony-lift-and-phase-no-go
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2312-sparse-root-bispectrum-positive-word-current
---

# An exact handoff word carries a literal unit-coloured pair

## 1. Why the quadratic object should be restricted first

THM-2302's quadratic no-go concerns the invalid inference

```text
word-restricted covariance
  -> an edge in the unrestricted owner spectrum.                    (1)
```

That inference is false because multiplication by a word mask is
convolution in frequency. But the restricted coefficient itself is not
lost. If the mask is transported back through the actual clock before
deconvolution, it is the Fourier transform of one literal source subset.

This gives a useful result strictly between THM-2305 and the stronger
ancestry target:

```text
unit-coloured pair in one literal handoff subset,

but not necessarily in the full exclusive-owner source.              (2)
```

## 2. Setup

Fix the selected owner `j` at depth `lambda`, its prescribed clock

```text
K=lambda+1,
```

and a positive canonical THM-2305 terminal word `Q`. Put

```text
E=E_j,

G=P^lambda 1_E,

E_Q=E intersection T^(-K)Q,

rho_Q=measure(E_Q)=integral_Q P G>0.                 (3)
```

For the thirteen predecessor values define

```text
v_r(y)=G((y+r)/13),

S_1(y)=sum_r v_r(y),

M_a(y)=sum_(r=0)^12 v_r(y)zeta^(-ar),     1<=a<=12,

N_a(y)=exp(-2*pi*i*a*y/13)M_a(y),

W_a=1_Q N_a.                                        (4)
```

Every active fibre has at most two entries in `[0,1]`.

Assume a named strictly deeper blocker

```text
c_*=13^K s,                   s in Z_(>0),

E intersection D_(c_*)=empty.                       (5)
```

On the strict THM-2302 shallow-owner branch one takes `c_*=c_3`.
Repeated-first conclusions below are conditional on the same named
deeper-exclusion hypothesis; THM-2296/2305 alone can select the deepest
owner there.

## 3. A sharper positive quadratic density

Define

```text
X_Q
 =1_Q sum_(a=1)^12 |M_a|^2
 =sum_(a=1)^12 |W_a|^2.                             (6)
```

If the two possible sheet masses are `x,y in [0,1]`, cyclic Parseval gives

```text
sum_(a=1)^12 |M_a|^2
 =13(x^2+y^2)-(x+y)^2
 =12(x^2+y^2)-2xy.                                  (7)
```

The exact lower remainder is

```text
12(x^2+y^2)-2xy-(11/2)(x+y)^2
 =(13/2)(x-y)^2>=0.                                 (8)
```

The quadratic in (7) is convex on the unit square, so its maximum is the
largest corner value, namely `22` at `(1,1)`. Therefore, a.e.,

```text
(11/2)1_Q S_1^2<=X_Q<=22.                           (9)
```

The word lies in a blocker danger set, so `measure(Q)<=1/7`. Since

```text
integral_Q S_1=13rho_Q,
```

Cauchy--Schwarz and (9) give

```text
mu_Q:=integral X_Q
 >=(11/2)(13rho_Q)^2/measure(Q)
 >=(13013/2)rho_Q^2.                               (10)
```

If `X_Q(y)>0`, some `z in E` has `T^Kz=y`. Equation (5) gives

```text
||s y||=||c_* z||>=1/14,
```

so

```text
support(X_Q) subset D_s^c.                          (11)
```

## 4. Shell extraction and exact bounds

Apply THM-2293's conditional-covariance shell lemma to the nonnegative
density `X_Q`, using the cap `22`. It supplies

```text
d=ms,

m!=0,                   gcd(m,91)=1,

Fourier(X_Q,d)!=0,                                  (12)
```

with

```text
0<|m|<=floor[33800*22/(961mu_0)]+1.                (13)
```

For the strict THM-2305 word floor,

```text
rho_Q>=39002430583/160481782761300,

mu_0
 =1521189591381733719889
   /3958257526818579244260000,

0<|m|<=2013431.                                     (14)
```

For the repeated-first floor, conditional on (5),

```text
rho_Q>=13560199813/160481782761300,

mu_0
 =183879018968485234969
   /3958257526818579244260000,

0<|m|<=16656658.                                    (15)
```

These are the first integers accepted by this exact tail ledger:
`mu_0/22-33800/(961m)` is positive at the displayed endpoint and
nonpositive one integer earlier.

## 5. Autocorrelation becomes a literal source pair

From (6) and (12), some character `a in {1,...,12}` satisfies

```text
Fourier(|W_a|^2,d)!=0.                              (16)
```

Because `W_a in L^2`, the autocorrelation series is absolutely convergent
for each fixed `d`:

```text
Fourier(|W_a|^2,d)
 =sum_(h in Z)
    Fourier(W_a,h) conjugate(Fourier(W_a,h-d)).     (17)
```

Thus some integer `h` makes both factors in (17) nonzero.

The mask must now be transported, rather than removed. The exact identity

```text
1_(T^(-1)Q) P^lambda 1_E
 =P^lambda 1_(E_Q)                                 (18)
```

gives

```text
Fourier(W_a,h)
 =13 Fourier(1_(E_Q),13^lambda(a+13h)).             (19)
```

Consequently `1_(E_Q)` has two genuine nonzero Fourier atoms

```text
A =13^lambda(a+13h),

A'=13^lambda(a+13(h-d)).                            (20)
```

Both have exact thirteen-adic grade `lambda`, the same nonzero root
character `a`, and

```text
A-A'
 =13^(lambda+1)d
 =m c_*,

gcd(m,91)=1.                                        (21)
```

Their coefficient product transforms under a source translation `tau`
by the unit deepest phase

```text
exp(-2*pi*i*m c_* tau).                              (22)
```

This is a literal, affine-covariant, unit-coloured pair on the exact
pure/fork handoff word.

## 6. The exact stopping boundary

The theorem proves more arity-wise than THM-2312's literal triad: it has
two legs, the same root character, and much smaller multiplier bounds.
It proves less about ancestry.

Multiplication by `1_Q` can create frequencies. The two selected atoms in
(20) need not both occur in the unrestricted rooted gauge `N_a`, or
equivalently in the spectrum of the full exclusive-owner indicator
`1_E`. Positive covariance

```text
<N_a,W_a>=integral_Q |N_a|^2>0
```

does guarantee **some** common `N_a,W_a` frequency, but does not force
that common frequency to be an endpoint of the pair (20).

The seven-periodic hostile from the companion bispectrum audit is exact:
there are `N_a` and `W_a=1_QN_a` for which `Fourier(|W_a|^2,1)!=0`, while
the unrestricted spectrum of `N_a` occupies one residue class modulo
seven and has no unit-coloured pair. It witnesses the first failed
implication in (1).

Therefore (21) is not yet:

```text
a unit edge incident to THM-2302's marked full-owner atom;

a THM-2303 relative terminal-component phase edge;

or a contradiction to the scalar cover.                            (23)
```

The next decisive target is an incidence theorem, not another existence
theorem:

```text
force one endpoint of (21) into
  support Fourier(N_a) intersection support Fourier(W_a),

or bound the spectral leakage from the other exact words strongly enough
  to prevent cancellation there.                                  (24)
```

This distinction is the same source/target sidecar that repeatedly appears
in THM-2286, THM-2299, and THM-2302. LRC(14) remains open.

## 7. Information ledger

```text
source:
  one positive exact THM-2305 word and a named deeper excluded blocker;

map:
  restrict the quadratic root energy first, apply the whole shell
  extractor, then transport the mask back through the actual clock;

preserved:
  source owner, prescribed clock, exact target word, root character,
  literal source subset, deepest unit colour, and affine phase;

destroyed:
  survival of either selected atom in the unrestricted owner spectrum;

needed sidecar:
  restricted-to-unrestricted spectral incidence or a quantitative leakage
  inequality;

hostile control:
  a seven-periodic unrestricted gauge whose word restriction creates a
  unit pair.                                                        (25)
```

## 8. Exact reproduction

The standard-library companion and transcript are

```text
04-computation/lrc14_literal_word_quadratic_pair_probe.py

05-knowledge/results/lrc14_literal_word_quadratic_pair_probe.out.
```

Their LF byte hashes are

```text
script:
  baac18ee0205c77ad84c2a3d6c41e34fd8c0bd7aa8252944e53ca295702ae8bf;

output:
  75bb3390118b0d44678a569d7bf7ebf3e79bdf0efc9a4e7d1188e4f91d85d01a.
```

Normal and optimized runs are byte-identical to the stored output. The
companion checks the exact two-sheet remainder, sharp cap, both mass
floors, both first accepted multiplier bounds and their preceding failures,
and the grade-transport scale. Equations (11)--(21) are the analytic proof,
not computational assumptions.
