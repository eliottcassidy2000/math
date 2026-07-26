---
id: THM-2421
title: "All-clock septimal ancestry endpoint-event detector"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED; INDEPENDENT AUDIT
  PENDING. For every integer clock R and rational finite-interval source
  E and target Q, the seven inverse-branch ancestry counts form an exact
  step vector. A signed source endpoint x creates the labelled interior
  jump sign(x)e_(floor(Rx) mod 7) at terminal phase {Rx}. Hence one base
  profile plus O(#boundary(E)) signed events reconstructs the complete
  seven-colour ancestry word without enumerating R branches. Its exact
  target-restricted cyclic Dirichlet energy Gamma is positive exactly
  when a positive-measure part of Q has nonflat ancestry, and Gamma is at
  least twice that nonflat mass. The theorem applies abstractly to every
  THM-2305 E_j,Q_(j,sigma) packet, but no explicit genuine global
  scalar-cover packet is presently available for the decisive
  computation. No scalar row is excluded and LRC(14) remains open.
source: codex-2026-07-26-septimal-ancestry-events
depends_on: []
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2399-physical-one-clean-forty-nine-orbit-sharpness
  - THM-2409-unfiltered-septimal-source-completion-and-word-phase-boundary
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
  - THM-2419-valuation-normalized-homogenization-of-affine-sideband-shells
script: 04-computation/lrc14_septimal_ancestry_endpoint_events_thm2421.py
output: 05-knowledge/results/lrc14_septimal_ancestry_endpoint_events_thm2421.out
script_sha256: 39347ab51b9f056962cf31e8b992ee89b27af53a67c5f09a7b2abc309ac91d5c
output_sha256: d3dd8addcf5b5ad25b72a38893e152ea0fd7416253814a98a9f9b214375a5a25
hash_basis: working-tree bytes (LF)
---

# THM-2421 -- all-clock ancestry is a signed endpoint-event word

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED; independent audit
pending.**

THM-2418 identifies the exact base-thirteen carry but also gives a fixed
real/even rank-one source hostile. The missing question is therefore not
whether the unweighted carry matrix has charged rank. It is whether the
actual exclusive source, on the actual terminal word, uses the seven
carry residues unequally on a set of positive measure.

There is a much cheaper exact carrier than the `R=13^k` prefix list:

```text
source endpoints
  -> signed terminal events labelled by prefix residue modulo seven
  -> complete ancestry step vector
  -> cyclic Dirichlet/translation-Gram surplus.                (1)
```

The construction works for every integer `R`, including clocks far too
large for prefix enumeration.

## 1. The ancestry vector

Work on the half-open coordinate `[0,1)`. Let `E` and `Q` be finite
disjoint unions of half-open intervals with rational endpoints, and let
`R>=1` be an integer. For `c in F_7` and `0<y<1`, define

```text
N_c(y)
 =#{0<=n<R:
       n=c mod 7 and (y+n)/R belongs to E}.                    (2)
```

Here `y` is the terminal coordinate, `n` is the inverse-branch prefix,
and `c=n mod 7` is its septimal carry colour. Boundary values are
irrelevant to every measure statement below.

If

```text
E=disjoint_union_i [a_i,b_i),
```

then an exact branch-free evaluation at any non-event point is

```text
L_i(y)=max(0,ceil(Ra_i-y)),
U_i(y)=min(R-1,ceil(Rb_i-y)-1),

N_c(y)
 =sum_i #{L_i(y)<=n<=U_i(y):n=c mod 7}.                        (3)
```

The last cardinality is one quotient-and-remainder calculation. Thus a
single chamber value costs `O(7 #intervals(E))`, not `O(R)`.

## 2. Signed endpoint-event law

Give every left endpoint of `E` sign `+1` and every right endpoint sign
`-1`. For an endpoint `x` with

```text
{Rx}!=0,
```

put

```text
y_x={Rx},
c_x=floor(Rx) mod 7.                                         (4)
```

Then the jump of the complete vector at `y_x` is

```text
N(y_x+)-N(y_x-)
 =sum_(x':{Rx'}=y_x) sign(x') e_(c_x').                       (5)
```

In particular, coincident endpoint images must be **added with their
signs and labels**; replacing them by an unlabelled breakpoint loses the
ancestry current.

To prove (5), fix one endpoint `x` and write

```text
Rx=floor(Rx)+{Rx}.
```

Exactly the inverse branch

```text
n=floor(Rx)
```

crosses `x` as `y` increases through `{Rx}`. A left endpoint enters `E`
and a right endpoint exits it, giving the claimed signed basis vector.
Different endpoints with the same terminal phase contribute additively.
No other branch changes there.

Endpoints with `{Rx}=0` lie at the chosen terminal cut `0=1`. They can
change null boundary values but no open chamber or integral. Formula
(3) initializes the first open chamber exactly, after which (5)
reconstructs every remaining chamber. Therefore:

```text
# interior event phases <= 2 #intervals(E),                    (6)
```

and the complete step word is computable with bit complexity depending
logarithmically on `R`, rather than by visiting all `R` prefixes.

## 3. The exact nonflatness detector

Define the target-restricted ancestry energy

```text
Gamma_Q(E;R)
 =integral_Q sum_(c in F_7)(N_(c+1)(y)-N_c(y))^2 dy.           (7)
```

All event and target endpoints are rational, so `Gamma_Q(E;R)` is an
exact nonnegative rational number. Put

```text
eta_Q(E;R)
 =measure{y in Q:N_0(y),...,N_6(y) are not all equal}.         (8)
```

Then

```text
Gamma_Q(E;R)>0
  iff eta_Q(E;R)>0,                                           (9)

Gamma_Q(E;R)>=2 eta_Q(E;R).                                  (10)
```

Indeed, the integrand in (7) vanishes exactly on constant vectors
because the residue cycle is connected. On a nonconstant integer
vector, the seven integer differences sum to zero and include at least
one positive and one negative value. Their squared sum is therefore at
least two. This proves both claims and identifies the equality boundary:
the pointwise lower bound is attained exactly when the cyclic difference
word has one entry `+1`, one entry `-1`, and five zeros.

The functional form is a nearest-neighbour translation-Gram distance:

```text
sum_c(N_(c+1)-N_c)^2
 =2(sum_c N_c^2-sum_c N_c N_(c+1)).                           (11)
```

Thus `Gamma` is the seven-colour ancestry analogue of THM-2365's
translation-Gram drift. It measures exact distance from the rank-one
constant line, not merely total source mass.

## 4. All six charged colours, and the alternating clock

Let `zeta=zeta_7` and

```text
Nhat_e(y)=sum_(c=0)^6 N_c(y) zeta^(-ec).
```

Parseval gives the pointwise identity

```text
sum_c(N_(c+1)-N_c)^2
 =1/7 sum_(e=1)^6
      |zeta^e-1|^2 |Nhat_e|^2.                               (12)
```

Moreover, because the coefficients `N_c(y)` are rational, vanishing of
even one charged value `Nhat_e`, `e!=0`, forces the polynomial

```text
sum_(c=0)^6 N_c X^c
```

to be a rational multiple of `Phi_7=1+X+...+X^6`. Hence

```text
N(y) nonflat
  iff Nhat_e(y)!=0 for every e in F_7^*.                       (13)
```

At a base-thirteen clock `R=13^k`, THM-2418's physical root law is

```text
l=n+(-1)^k r mod 7.                                          (14)
```

Changing `k` parity reflects the root coordinate because
`13=-1 mod 7`; changing affine origin cyclically permutes it. The energy
(7) is invariant under both operations. This is the exact sense in
which the limiting ancestry object is alternating while its scalar
detector is all-clock: the parity is retained in the labelled event
word and harmlessly quotiented by a dihedral-invariant energy.

## 5. Sharp hostile and positive controls

THM-2418's single fixed source

```text
G=[3/13,10/13)
```

has, for every `R=13^k` and almost every terminal `y`,

```text
N(G;R)(y)
 =(13^(k-1),...,13^(k-1)).                                   (15)
```

Both endpoints map to the terminal cut, so there are no interior events,
and

```text
Gamma_Q(G;13^k)=0
```

for every target `Q`. This proves that the detector does not silently
assume away the canonical rank-one hostile.

At the opposite boundary, one prefix cylinder

```text
E=[0,1/R)
```

has `N=(1,0,0,0,0,0,0)` almost everywhere and therefore

```text
Gamma_Q(E;R)=2 measure(Q).                                   (16)
```

The cheap scalar test

```text
sum_c N_c not congruent to 0 mod 7
```

is sufficient for nonflatness but very far from necessary. The exact
seven-cylinder source

```text
E=[0,6/R) union [7/R,8/R)
```

has

```text
N=(2,1,1,1,1,1,0),
sum_c N_c=7,
Gamma_[0,1)(E;R)=6.                                         (17)
```

Thus the full vector, not total ancestry modulo seven, is the correct
test.

There is a second, stronger boundary. At `R=13`, take one half of each
of the first seven prefix cylinders: the first half for prefix residues
`0,1,2` and the second half for residues `3,4,5,6`. Explicitly,

```text
E_bal
 =union_(c=0)^2 [c/13,(2c+1)/26)
   union union_(c=3)^6 [(2c+1)/26,(c+1)/13).                  (18)
```

On the first terminal half its ancestry profile is
`(1,1,1,0,0,0,0)`; on the second it is `(0,0,0,1,1,1,1)`.
Consequently

```text
integral_0^1 N_c(y)dy=1/2             for every c,
Gamma_[0,1)(E_bal;13)=2.                               (19)
```

Thus even the **integrated** seven carry histogram can be perfectly
flat while positive-measure pointwise ancestry is nonflat. The squared
event detector retains information which a linear integrated carry
coefficient can cancel.

Finally, for the exact terminal word `Q=D_7`, of mass `1/7`, the flat
source (15) still gives zero while the one-cylinder control gives

```text
Gamma_(D_7)([0,1/R);R)=2/7.                                 (20)
```

Terminal restriction therefore neither manufactures nor automatically
destroys the ancestry signal.

## 6. Relation to the fixed source--terminal classifier

For a fixed source `G=1_E`, the weighted carry masses in THM-2418 are

```text
a_(R,c)
 =integral_T 1_E(x) Q({Rx})
    1_(floor(Rx)=c mod 7)dx.
```

The inverse-branch disintegration gives the exact identity

```text
R a_(R,c)=integral_Q N_c(y)dy.                               (21)
```

Hence `Gamma_Q(E;R)=0` forces the centred weighted carry histogram to
vanish. The converse is false by (18)--(19): the integrated histogram
for `E_bal,Q=1,R=13` is flat although `Gamma=2`. The THM-2418
scale-periodic classifier and the present detector therefore sit at two
different levels:

```text
linear integrated carry histogram
  <- can cancel across terminal chambers;

pointwise ancestry event word / quadratic Gamma
  <- retains chamberwise nonflatness.                         (22)
```

This distinction is load-bearing when a source owner changes its carry
imbalance across the target word.

## 7. The genuine scalar-cover application remains open

For a THM-2305 scalar-cover row, the exclusive owner `E_j` and every
terminal word `Q_(j,sigma)` are rational finite-interval sets. At its
base-thirteen clock `R`, define

```text
N_(j,c)(y)
 =#{0<=n<R:
      n=c mod 7 and (y+n)/R belongs to E_j}.                   (23)
```

Equations (3)--(7) give the cheapest genuine scalar-cover computation:

```text
compute the exact endpoints of E_j and Q_(j,sigma);
push each E_j endpoint to ({Rx},floor(Rx) mod 7,sign);
sweep the joint E_j-event/Q_(j,sigma)-endpoint chambers;
test Gamma_(Q_(j,sigma))(E_j;R)>0.                            (24)
```

A positive value is exactly a positive-measure nonflat ancestry vector
on the actual terminal word. By (12)--(13), it retains every charged
carry colour in `L^2`. It does **not** by itself prove a relation-address
landing, prevent cancellation of a chosen linear Fourier coefficient,
or provide THM-2419's same-affine-shell residue-zero reference.

No explicit genuine global scalar-cover `E_j,Q_(j,sigma)` packet is
currently present in canon. THM-2399 supplies a strict physical packet
whose scalar cover is only local on one forty-nine-orbit and which has a
global safe point; feeding it into (20) as though it were global would
erase the theorem's main boundary. Consequently (20), not a numerical
LRC row verdict, is the present advance.

## 8. Connection and information-loss ledger

```text
source:
  the exclusive owner E_j, exact terminal word Q_(j,sigma), and clock R;

target:
  the seven-vector ancestry step word N(y), then its nonnegative
  translation-Gram energy Gamma;

map:
  inverse branch n -> carry colour n mod 7, and signed endpoint
  x -> ({Rx},floor(Rx) mod 7,sign(x));

preserved by the event word:
  source/target sets, clock, carry residue, endpoint sign, coincident
  event multiplicity, chamber order, base profile, and 13-adic parity;

destroyed by N:
  the individual prefix n inside one residue class and the corresponding
  exact physical preimage;

destroyed by Gamma:
  event order, affine phase, which residue is large, linear character
  phase, and individual endpoint ancestry;

needed continuation sidecar:
  retain the labelled signed event word and base profile, not Gamma alone;
  for physical word landing also retain the quotient/physical word pair,
  and for affine-sideband composition retain a same-shell reference.     (25)
```

The endpoint word is therefore the minimal exact transport stalk for this
question. `Gamma` is its cheapest decisive scalar **only for
nonflatness**, not a substitute for the stalk in later composition.

## 9. Exact hostile audit

The dependency-free companion has two independent in-file evaluation
paths:

1. definition-level enumeration of every endpoint/branch breakpoint and
   all `R` inverse branches on every resulting chamber; and
2. integer-range initialization followed by the signed event sweep.

It compares them on 192 exact random interval-set fixtures, 2,062
chambers, 1,113 nonzero aggregate jumps, circle-cut controls, and a
coincident mixed-sign event. It then checks (15) through `k=8`, where
`R=815730721`, without an `R`-loop; checks (16)--(20); confirms the
divisible-total hostile (17); and verifies the flat-integrated but
pointwise-charged control (18)--(19). Reproduce with

```bash
python3 04-computation/lrc14_septimal_ancestry_endpoint_events_thm2421.py
python3 -O 04-computation/lrc14_septimal_ancestry_endpoint_events_thm2421.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_septimal_ancestry_endpoint_events_thm2421.out
```

after LF normalization. Every load-bearing check remains active under
optimized Python. Independent hostile audit remains pending. QED.
