---
id: THM-2326
title: "Vertexwise septimally primitive c3 degree"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Let f be a
  step function supported outside the danger comb D_c, with J jumps, and
  let f-hat(A) be nonzero with c not dividing A. Modulating the exact
  disjointness f*1_(D_c)=0 by A forces an incident coefficient
  f-hat(A+mc) with 1<=m<=7J-1 and 7 not dividing m. Applied to either
  shallow canonical exclusive source E_j and c=c_3, every marked shell
  vertex has a same-grade, same-root-character c_3-neighbour with
  septimally primitive multiplier m<=14S-1. Thus THM-2323's common
  bare/word vertex has a canonical incident edge whose only remaining
  forbidden colour can be 13. A sharp thirteen-periodic real indicator
  shows that endpoint recurrence alone cannot remove that factor.
  No 91-unit incident edge, scalar-row exclusion, or proof of LRC(14) is
  claimed.
source: codex-2026-07-25-vertexwise-septimal-degree
depends_on:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
related:
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
  - THM-2322-local-hostile-c3-toothpick-ladder
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
script: 04-computation/lrc14_vertexwise_septimal_c3_degree_thm2326.py
output: 05-knowledge/results/lrc14_vertexwise_septimal_c3_degree_thm2326.out
script_sha256: e0a8fda969ad141822b837b555a6b97fd37b348c70d50d05a0fea15cb0dbc375
output_sha256: 7334acaadfd1f20c179917418c773e65aca7571b6b26ce97c416bac7b856ae77
hash_basis: working-tree bytes (LF)
---

# THM-2326 -- every marked canonical vertex has septimally primitive c3 degree

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2302 gives every marked shell vertex an incident `c_3`-multiple edge,
but its endpoint recurrence does not control the multiplier colour.
THM-2293 separately gives a `91`-unit edge somewhere in the same
root-character graph. The quantifiers do not compose.

The canonical source has one extra property not used by the generic
endpoint recurrence: either shallow exclusive owner is disjoint from the
deepest danger comb. Modulate that exact disjointness by the marked
frequency. The zeroth Fourier term is then the marked atom itself, while
every other surviving comb term is a `c_3`-neighbour whose multiplier is
not divisible by seven.

The result is vertexwise and canonical:

```text
marked canonical source vertex
  -> modulated source/deepest-comb disjointness
  -> incident c_3 edge with 7 not dividing its multiplier
  -> only the possible factor 13 remains.                          (1)
```

## 1. The danger-comb Fourier needle

On the circle put

```text
D_c={x:||cx||<1/14}.                                (2)
```

It is the union of `c` centered intervals of half-width `1/(14c)`.
For the convention

```text
h_hat(n)=integral_T h(x)exp(-2*pi*i*n*x)dx,         (3)
```

the root-of-unity sum gives

```text
(1_(D_c))_hat(n)=0                     if c does not divide n,

(1_(D_c))_hat(cm)
 =sin(pi*m/7)/(pi*m)                  if m!=0,

(1_(D_c))_hat(0)=1/7.                              (4)
```

In particular, among nonzero multiples of `c`, the coefficient in (4)
vanishes exactly when `7` divides `m`. The deepest danger comb is therefore
a Fourier needle supported on the septimally primitive directions together
with the zero mode.

## 2. The modulated disjointness identity

Let `f` be a complex-valued step function with `J` nonzero circle jumps.
Assume

```text
f 1_(D_c)=0 almost everywhere,                     (5)

f_hat(A)!=0,
c does not divide A.                               (6)
```

The second condition in (6) guarantees that no frequency `A+cm` is zero.
Define

```text
W_A(x)=exp(-2*pi*i*A*x)f(x).                        (7)
```

Then

```text
(W_A)_hat(n)=f_hat(A+n),
W_A 1_(D_c)=0.                                     (8)
```

Parseval, followed by (4), gives the exact identity

```text
0
 =sum_(n in Z)(W_A)_hat(n)
        conjugate((1_(D_c))_hat(n))

 =(1/7)f_hat(A)
   +sum_(m!=0, 7 does not divide m)
      f_hat(A+cm) sin(pi*m/7)/(pi*m).               (9)
```

There is no conditional rearrangement in (9). Both functions are in
`L^2`, so Cauchy--Schwarz makes the Parseval series absolutely convergent:

```text
sum_n |(W_A)_hat(n)(1_(D_c))_hat(n)|
 <=||W_A||_2 ||1_(D_c)||_2<infinity.                (10)
```

The first term on the right of (9) is nonzero. Therefore some signed
integer

```text
m_0!=0,
7 does not divide m_0                              (11)
```

satisfies

```text
f_hat(A+c m_0)!=0.                                 (12)
```

This is already an incident septimally primitive edge. The next step lands
one at a positive bounded multiplier.

## 3. Residuewise endpoint landing

Let `r in {1,...,6}` be the residue of `m_0` modulo seven. For `h in Z`
define

```text
X_h
 =(A+c(r+7h)) f_hat(A+c(r+7h)).                    (13)
```

Distributional differentiation of the step function gives

```text
2*pi*i X_h
 =sum_(x in Jump(f))
    Delta_f(x) exp(-2*pi*i(A+cr)x)
               [exp(-2*pi*i*7cx)]^h.               (14)
```

After equal endpoint nodes are combined, (14) is an exponential sequence
on at most `L<=J` distinct nodes. It is not the zero sequence: write
`m_0=r+7h_0` and use (6) and (12).

A nonzero `L`-node exponential sequence cannot vanish at `L` consecutive
integers. Applying the Vandermonde lemma to

```text
h=0,1,...,L-1
```

therefore gives an index with `X_h!=0`. Put

```text
m=r+7h.
```

Then

```text
1<=m<=7J-1,
7 does not divide m,
f_hat(A+mc)!=0.                                    (15)
```

This proves the general **septimally primitive marked-degree lemma**.
Unlike the generic marked-degree lemma in THM-2302, it uses the actual
source/deepest-comb disjointness before applying endpoint recurrence.

## 4. Canonical LRC specialization

On a strict LRC(14) scalar profile, let `E_j` be either shallow canonical
exclusive source used by THM-2293:

```text
j in {1,2},
f_j=1_(E_j).
```

By exclusivity,

```text
E_j intersection D_(c_3)=empty.                    (16)
```

The canonical source has at most

```text
J_(f_j)<=2S                                        (17)
```

jumps. Let `A` be any nonzero canonical source atom of exact shallow
grade `b` and root character `kappa`:

```text
(f_j)_hat(A)!=0,

A=13^b a,                 13 does not divide a,

A/13^b=kappa mod 13.                               (18)
```

The deepest blocker has

```text
c_3=13^c u_3,             c>b,
13 does not divide u_3.                             (19)
```

Thus `c_3` cannot divide `A`, so Sections 2--3 apply with `c=c_3`.
There is an integer

```text
1<=m<=14S-1,
7 does not divide m                                (20)
```

such that

```text
(f_j)_hat(A+m c_3)!=0.                              (21)
```

The new atom remains in exactly the same shell:

```text
A+m c_3
 =13^b[a+13^(c-b)u_3m],

nu_13(A+m c_3)=b,

(A+m c_3)/13^b=kappa mod 13.                       (22)
```

Equations (20)--(22) prove that **every** marked vertex in either canonical
shallow root-character graph has a bounded incident local arithmetic edge
whose multiplier is a seven-unit. No averaging over vertices and no
choice of a different root character occurs.

For THM-2302's marked service vertex, the original endpoint retains its
named current-service coefficient. For the middle-owner vertex supplied
by THM-2323, the original endpoint is simultaneously a bare-source atom
and a literal positive-word atom at the same `91`-colour and gauge.
THM-2326 attaches to that endpoint a canonical bare-source neighbour with
`7`-unit multiplier. It does not assert that the neighbour retains the
word coefficient.

## 5. The exact remaining colour and a sharp recurrence hostile

The multiplier in (20) can still satisfy

```text
13 divides m.                                      (23)
```

The modulated comb identity (9) cannot see this colour: the Fourier zeros
of `1_(D_(c_3))` detect multiples of seven, not multiples of thirteen.

Endpoint recurrence alone genuinely cannot repair (23). Let

```text
F=union_(r=0)^12 [r/13,(r+1/4)/13).                (24)
```

This is a nonconstant real indicator with `26` jumps and period `1/13`.
Writing `I=[0,1/52)`, one has

```text
1_F=sum_(r=0)^12 1_(I+r/13).
```

Therefore

```text
(1_F)_hat(n)=0 unless 13 divides n,                (25)

(1_F)_hat(13)=(1-i)/(2*pi)!=0.                     (26)
```

If `s` is coprime to thirteen, then

```text
(1_F)_hat(13+ms)!=0
   only if 13 divides m.                            (27)
```

For example, at `s=1` the next positive surviving multiplier is `m=13`,
because `(1_F)_hat(26)!=0`. Thus a marked endpoint sequence can
concentrate every visible incident multiplier in the forbidden
thirteen-colour.

The hostile (24) is a boundary for endpoint recurrence, not a canonical
LRC row and not a counterexample to a future theorem that also uses the
full canonical endpoint field, the THM-2323 primitive-colour diagonal, or
another canonical sidecar.

## 6. Frontier effect

The exact coloured-degree ledger is now

```text
THM-2302:
  every marked canonical vertex has a nonzero-multiplier neighbour;

THM-2326:
  every marked canonical vertex has a 7-unit-multiplier neighbour;

THM-2293:
  every selected character graph has a 91-unit edge somewhere;

THM-2323:
  every primitive colour has a common bare/word marked vertex;

THM-2323 conditional incidence branch:
  under its explicit low-cofactor criterion, two common bare/word marks
  themselves form a 91-unit c_3 edge;

remaining outside THM-2323's conditional branch:
  remove a possible factor of 13 from the vertexwise THM-2326 edge,
  while retaining the word mark and target-plane address.             (28)
```

The faithful carrier is an undirected edge-coloured graph. Reversing the
edge replaces `m` by `-m`; no target predicate or time arrow orients it.
A tournament encoding would erase the load-bearing colour in (20).

No scalar profile is excluded. The scalar ledger remains `165`, and
LRC(14) remains open.

## 7. Exact verification

The companion checks the danger-comb support and zero predicates, every
bound in the residuewise landing, `100000` grade/root-preservation rows,
the canonical `14S-1` ledger, and the exact thirteen-periodic hostile
support through `200000` signed multipliers. All load-bearing tests raise
explicitly under ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_vertexwise_septimal_c3_degree_thm2326.py
python3 -O 04-computation/lrc14_vertexwise_septimal_c3_degree_thm2326.py
```

The transcripts must match

```text
05-knowledge/results/lrc14_vertexwise_septimal_c3_degree_thm2326.out
```

byte-for-byte after LF normalization.
