---
id: THM-2377
title: "Septimal valuation-collision and Bockstein carry gate"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. A carrier
  made from finitely many fixed radius-1/14 danger/safe or radius-1/7
  guard-safe factors, one higher septimal graft, and one still higher
  target tooth is exactly circulant whenever the fixed-factor speeds
  have pairwise distinct 7-adic valuations. Thus every nonzero
  high-target mode needs a repeated lower valuation layer. At a
  two-term minimal layer the active Fourier indices satisfy an exact
  chi_7 colour-transfer law, and after leading cancellation their
  quotient is the Bockstein carry which must climb to the graft and
  target layers. Distinct-level factors are scalar spectators until the
  first repeated layer enters. This is a necessary carrier gate, not a
  sufficient target certificate or scalar-row exclusion; the ledger
  remains 165 and LRC(14) remains open.
source: codex-2026-07-25-septimal-valuation-collision
depends_on:
  - THM-1156-tooth-seam-chi7-bipartition
  - THM-2372-hard-septimal-signed-stalk-and-toothpick-divisibility
related:
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2369-complete-line-target-dirichlet-and-balanced-observable-no-go
script: 04-computation/lrc14_septimal_collision_carry_thm2377.py
output: 05-knowledge/results/lrc14_septimal_collision_carry_thm2377.out
script_sha256: 92ab72ce9fa0c3e5cfb135d5dfa7b30924339423bfe3b3a4e11301dab6c23ea6
output_sha256: 75ddc06b25f773e00e898fc433004da64f9c5e92f386fc30b2c13437b3a7b3f1
hash_basis: working-tree bytes (LF)
---

# THM-2377 -- septimal valuation-collision and Bockstein carry gate

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2372 shows that the scalar cover forces a `98`-nested
low/high-blocker pair, yet the associated minimal lawful carrier is
exactly circulant. The next question is not merely whether an extra mask
is present. It is whether two retained masks occupy the same septimal
layer.

The answer is exact:

```text
all retained lower valuations distinct
  -> every lower factor is a scalar spectator
  -> the high-target tensor is circulant;

nonzero high-target mode
  -> at least two active lower factors share the minimum valuation
  -> their leading F_7 cancellation creates a Bockstein carry.       (1)
```

This identifies the faithful vertices for the next tournament/Fano
probe: same-layer collision pairs with their leading indices and carry,
not raw runner labels.

## 1. Interval factors and their septimal spectrum

Put

```text
d_0(x)=1_(||x||<1/14),

g(x)=1-d_0(x),

u_2(x)=1_(||x||>=1/7).                              (2)
```

For a factor `phi in {d_0,g,u_2}`, its mean is respectively

```text
mu(phi) in {1/7,6/7,5/7}.                           (3)
```

For every nonzero integer `a`,

```text
d_0^hat(a)= sin(pi a/7)/(pi a),

g^hat(a)=  -sin(pi a/7)/(pi a),

u_2^hat(a)=-sin(2pi a/7)/(pi a).                    (4)
```

Consequently

```text
phi^hat(a)!=0       implies       7 does not divide a. (5)
```

The same implication holds after an arbitrary fixed translate of
`phi`; translation changes only the coefficient phase.

For `L in {1,2}`, retain the notation

```text
u_L(x)=1_(||x||>=L/14),

mu(u_L)=(7-L)/7.                                    (6)
```

It has the same nonzero-index support law (5).

## 2. All-distinct lower layers are exactly circulant

Let

```text
v_1,...,v_m,q,d
```

be positive integers, let `theta_i in T`, and choose
`phi_i in {d_0,g,u_2}`. Assume

```text
nu_7(v_i)<nu_7(q)<nu_7(d)               for every i,

nu_7(v_1),...,nu_7(v_m)                 pairwise distinct. (7)
```

For `r,t in F_13`, define

```text
T_L(r,t)
 =integral_T
   product_(i=1)^m phi_i(v_i x+theta_i)
   d_0(dx-r/13)g(dx-t/13)
   u_L(qx+t/13) dx.                                  (8)
```

Then

```text
T_L(r,t)
 =mu(u_L) product_i mu(phi_i) J(r-t),                (9)

J(a)
 =integral_T d_0(y-a/13)g(y)dy,

J(0)=0,

J(+1)=J(-1)=1/13,

J(a)=1/7                         otherwise.          (10)
```

In particular `T_L` is exactly circulant.

### Proof

Use Fejer approximants for every bounded interval factor. They are
uniformly bounded, converge in `L^1`, and make the following spectral
argument finite. A contributing monomial has integer indices

```text
a_1,...,a_m,k,n
```

obeying

```text
sum_i a_i v_i+kq+nd=0.                              (11)
```

Whenever `a_i` or `k` is nonzero and its interval coefficient survives,
(5) gives

```text
nu_7(a_i v_i)=nu_7(v_i),

nu_7(kq)=nu_7(q).                                   (12)
```

These valuations are pairwise distinct and all strictly below
`nu_7(nd)` when `n!=0`. If any lower index were nonzero, the left side
of (11) would have a unique term of minimum `7`-adic valuation and could
not vanish. Hence

```text
a_1=...=a_m=k=0,
```

and then (11) forces `n=0`. Only the constant coefficient of every lower
factor survives. Passing to the `L^1` limit gives (9); Haar invariance
under `x -> dx` gives (10). No infinite Fourier product is rearranged.
QED.

The contrapositive is the valuation-collision gate:

> If a tensor of the form (8) is noncirculant, then at least two fixed
> lower speeds `v_i,v_j` have
>
> ```text
> nu_7(v_i)=nu_7(v_j).                               (13)
> ```

The theorem is phase-uniform and permits any finite number of fixed
factors.

## 3. The exact one-extra-factor hostile

Take

```text
(c,s,q,d)=(13,7,49,8918),

d=686c,

(nu_7(c),nu_7(s),nu_7(q),nu_7(d))=(0,1,2,3).        (14)
```

For `J,L in {1,2}`, specialize (8) to

```text
T_(J,L)(r,t)
 =integral_T
   d_0(cx)u_J(sx)
   d_0(dx-r/13)g(dx-t/13)
   u_L(qx+t/13)dx.                                  (15)
```

Equation (9) gives

```text
T_(J,L)(r,t)
 =((7-J)(7-L)/343)J(r-t).                           (16)
```

On the common exact endpoint grid

```text
D=182 lcm(c,s,q,d)=1623076,
```

the four difference profiles, in order
`0,1,...,12`, are

```text
(J,L)=(1,1):
  (0,13104,24336,24336,24336,24336,24336,
     24336,24336,24336,24336,24336,13104),

(J,L)=(1,2) or (2,1):
  (0,10920,20280,20280,20280,20280,20280,
     20280,20280,20280,20280,20280,10920),

(J,L)=(2,2):
  (0,9100,16900,16900,16900,16900,16900,
     16900,16900,16900,16900,16900,9100).           (17)
```

Thus merely retaining one additional lower mask does not escape the
zero-drift branch. What matters is collision of valuation layers.

## 4. The first repeated layer has a `chi_7` law

Consider any contributing tuple before the distinct-layer hypothesis is
imposed. Let `e` be the minimum valuation among its active lower terms.
The ultrametric inequality forces that minimum to occur at least twice.

Suppose exactly two terms attain it, at speeds

```text
v_i=7^e V_i,              v_j=7^e V_j,

7 does not divide V_i V_j,                            (18)
```

with nonzero Fourier indices `a,b`. Reduction of (11) at the leading
digit gives

```text
aV_i+bV_j=0                    in F_7.               (19)
```

Let

```text
epsilon(v)=chi_7(v/7^nu_7(v))
```

be THM-1156's canonical speed colour. Since
`chi_7(-1)=-1`, equation (19) gives the exact transfer law

```text
chi_7(b/a)
 =-epsilon(v_i)epsilon(v_j).                         (20)
```

Therefore:

- opposite speed colours force the two Fourier indices to have the same
  `chi_7` colour;
- equal speed colours force opposite index colours.

This is a genuine map from a runner/tooth seam colour to the spectral
colour which can drive a target. The sidecar `(a,b)` is essential; a
runner-only bipartition loses it.

### A sharp positive same-layer control

The collision gate is attained by a physical tensor. Take

```text
(c,s,q,d)=(13,1,7,1274),

(nu_7(c),nu_7(s),nu_7(q),nu_7(d))=(0,0,1,2),
```

and use `d_0(cx)g(sx)`, `L=1` in (8). Exact interval counting gives the
rank-one difference/co-shift formula

```text
T_1(r,t)=J(r-t)N_t/637,

N_t
 =78                         if t=0,
 =74                         if t=+1 or -1,
 =71                         otherwise.             (20a)
```

In particular, on the line `r=t+1`,

```text
8281 T_1(1+t,t)
 =78                         if t=0,
 =74                         if t=+1 or -1,
 =71                         otherwise.             (20b)
```

One short proof averages the thirteen roots above `z=13x`. The factor
`g(x)` removes root `0`; the moving `g(7x+t/13)` removes root `-2t` on
the central `z` interval and one additional translated root on each
outer interval. Under the deep scale `98z`, those three regions contain
exactly `2,6,6` full periods, with weights

```text
1/637, 3/637, 3/637.
```

The three root counts give (20a)--(20b). For every nonzero target colour `b`,
the `t`-Fourier coefficient is proportional to

```text
7+3(zeta_13^b+zeta_13^(-b))
 =7+6cos(2 pi b/13)
>0.                                                   (20c)
```

Thus all twelve nonzero target colours survive. Repeated valuation is
not sufficient in general, but it is the sharp first place where genuine
drift becomes possible.

## 5. The Bockstein carry is the next digit

For the two-low-factor carrier, write

```text
c=7^e C,                 s=7^e S,

q=7^f Q,                 d=7^g D,

e<f<g,                   7 does not divide CSQD.    (21)
```

A contributing nonzero target monomial satisfies

```text
ac+bs+kq=nd.                                          (22)
```

After division by `7^e`, the first necessary congruence is

```text
aC+bS=0                         mod 7^(f-e).         (23)
```

Define its exact carry

```text
h=(aC+bS)/7^(f-e).                                  (24)
```

The remaining equation is

```text
h+kQ=0                         mod 7^(g-f).          (25)
```

For adjacent layers `f=e+1`, the residue

```text
beta(a,b)=(aC+bS)/7               mod 7             (26)
```

is the Bockstein sidecar invisible to the leading `chi_7` colour. The
next digit is

```text
beta(a,b)+kQ=0                    mod 7.             (27)
```

Thus the faithful collision vertex is

```text
(same-layer speed pair, index ratio b/a,
 leading chi_7 transfer, Bockstein carry beta),     (28)
```

not a bare edge in a tournament.

There is also a canonical recursive descendant. Choose the unique
balanced

```text
rho in {-3,-2,-1,1,2,3}
```

with

```text
S congruent rho C                  (mod 7).
```

The allowed Fourier pair `(-rho,1)` has nonzero interval coefficients
and produces

```text
-rho c+s=7^(e+1)h_1,

h_1=(S-rho C)/7,

|h_1|<=(S+3C)/7.                                  (29)
```

Iteration of (29), with the successive carries retained, is a signed
`7`-adic Euclidean/toothpick ladder. It is a much smaller state space
than the original integer speeds.

## 6. Consequence for the hard scalar lane and deletion order

In THM-2372's hard lane,

```text
nu_7(c_*)<M=nu_7(q_*)<nu_7(d)
```

for the forced nested high blocker `d`. Apply (13) to whatever fixed
guard/unit/blocker factors are retained in the carrier:

```text
high-target drift
 -> two retained lower factors share a septimal valuation.          (30)
```

Accordingly one may order target-neutral fixed-factor insertions by
valuation class. Every prefix containing at most one representative of
each class is exactly circulant by Section 2. Within such an order, the
first THM-2370 deletion layer capable of carrying drift is the insertion
of a second fixed factor into an already occupied class. Unique-layer
fixed factors before that point are provably drift-neutral scalar
spectators. Lawfully co-shifted factors outside the model (8) require
their own gauge audit.

This does not say that every repeated layer drifts. Even the canonical
pair (29) must climb through (25), survive its interval coefficients and
thirteen co-shift phases, and avoid cancellation after the remaining
owner masks are restored. The theorem identifies the first possible
carrier and the missing coordinates; it does not supply the final
charged target certificate.

No scalar profile is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 7. Exact companion

The dependency-free exact companion:

- checks the all-distinct spectral exclusion on finite hostile banks;
- reproduces all four rational profiles in (17);
- exhausts the `F_7^*` solutions of (19) and verifies the `chi_7`
  transfer (20);
- checks the two-stage carry congruences (23)--(27);
- enumerates the balanced descendants (29) and their contraction bound;
- exactly recounts the same-layer tensor (20a)--(20b) and all twelve
  nonzero target colours in (20c); and
- records explicit controls showing that a repeated layer is necessary
  but not asserted sufficient.

Run

```bash
python3 04-computation/lrc14_septimal_collision_carry_thm2377.py
python3 -O 04-computation/lrc14_septimal_collision_carry_thm2377.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_septimal_collision_carry_thm2377.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
