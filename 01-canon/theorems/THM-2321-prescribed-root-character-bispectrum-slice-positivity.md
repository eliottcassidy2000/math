---
id: THM-2321
title: "Prescribed-root-character bispectrum slice positivity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For a
  nonzero nonnegative vector on F_p supported on at most two sites, every
  fixed generic bispectrum slope lambda has complete scalar-fibre sum at
  least (p-4)/4 times mass cubed. Dually, for every prescribed nonzero
  first character k, the real sum over its second legs is at least
  (p-4)(1-cos(pi/p))/4 times mass cubed. At p=13, on every positive
  canonical blocker word and for every prescribed k, one cubic current has
  first character k and real part greater than (5733/22)rho^3. Whenever
  THM-2318's one-shot corner exists, its collision atom and the cubic word
  current can share one exact normalized root colour. Projectivized
  character pairs canonically label all fourteen directions in THM-2315's
  target plane, and every direction has positive aggregate current.
  This selects every abstract gain class, but an individual prescribed
  root-colour/gain cell may have negative real part and no actual relation
  address is proved incident to the current. LRC(14) remains open.
source: codex-2026-07-25-prescribed-character-bispectrum
depends_on:
  - THM-2312-sparse-root-bispectrum-positive-word-current
  - THM-2318-one-shot-three-prime-mobius-amplifier
related:
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
script: 04-computation/lrc14_prescribed_character_bispectrum_thm2321.py
output: 05-knowledge/results/lrc14_prescribed_character_bispectrum_thm2321.out
script_sha256: f9d1f2f55e575975877abcca9d90d242ca68da2c3360a36d5427bd6c2fe0e568
output_sha256: 5e6695c8c0c691983052e2830627986a54fe39994309c54c00ed567fc07d386c
hash_basis: working-tree bytes (LF)
---

# THM-2321 -- every bispectrum row and column survives

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2312 proved positivity only after summing all `132` nontrivial
bispectrum cells. It then selected one ordered pair `(k,l)`. That statement
appears to lose both:

```text
the projective shape lambda=l/k,
the first root character k.                         (1)
```

For a root fibre supported on at most two sheets, neither loss is real.
The `132` cells form a `12 x 11` rectangle

```text
F_13^* x (F_13 minus {0,-1}),
(k,lambda) -> (k,lambda k).                         (2)
```

Every column has a sharp positive sum. Less obviously, the real sum of
every row is also strictly positive. Consequently the first character can
be prescribed before selecting the second leg.

This aligns THM-2318's exact collision colour with a cubic current on the
same canonical word. Projectivizing the ordered pair also gives the exact
coordinate shape of THM-2315's target plane. Neither representation-level
alignment makes the ordinary collision atom a cubic factor or realizes an
actual relation address.

## 1. The bispectrum rectangle

Let `p` be an odd prime, let `zeta` be a primitive `p`th root of unity, and
let

```text
v=(v_r)_(r in F_p),          v_r>=0,

M_k=sum_r v_r zeta^(-kr),    k in F_p.              (3)
```

Assume `v` is nonzero and has at most two occupied sites. Write their
masses as `a,b>=0`,

```text
V=a+b>0,
x=ab/V^2<=1/4.                                    (4)
```

For nonzero `k,l,k+l`, put

```text
C_(k,l)=M_k M_l conjugate(M_(k+l)).                 (5)
```

Every allowed pair has the unique form

```text
(k,l)=(k,lambda k),

k in F_p^*,
lambda in Lambda_p:=F_p minus {0,-1}.               (6)
```

Thus the complete face is a disjoint union of `p-2` free scalar fibres,
or dually `p-1` first-character rows.

## 2. Every fixed projective shape is positive

For `lambda in Lambda_p`, define

```text
S_lambda(v)
 =sum_(k in F_p^*) C_(k,lambda k).                  (7)
```

> **Fixed-shape identity.**
>
> ```text
> S_lambda(v)=p sum_r v_r^3-(sum_r v_r)^3.          (8)
> ```

### Proof

Expand one cell using physical sites `r,s,t`. Summing over nonzero `k`
gives

```text
sum_(k!=0) zeta^(-k(r+lambda s-(1+lambda)t))
 =p[r+lambda s=(1+lambda)t]-1.                     (9)
```

Translate and rescale the two occupied sites to `0,1`. For binary
`r,s,t`, the equation in (9) has only

```text
(r,s,t)=(0,0,0),(1,1,1).                           (10)
```

Indeed, the six mixed triples would force one of

```text
1=0,       lambda=0,       1+lambda=0.
```

All are excluded. The `p` part of (9) therefore contributes
`p(a^3+b^3)`, while the `-1` part contributes `-V^3`. This is (8). QED.

In terms of (4),

```text
S_lambda(v)/V^3
 =p(1-3x)-1
 =(p-1)-3px.                                       (11)
```

Hence, for every `p>=5` and **every prescribed shape**,

```text
S_lambda(v)>=(p-4)V^3/4>0.                         (12)
```

The constant is sharp. Equality in the non-strict floor occurs exactly
when both sheets have equal positive mass.

This decomposes THM-2312's whole-face identity exactly:

```text
B_p(v)=sum_(lambda in Lambda_p)S_lambda(v).         (13)
```

The old sharp whole-face floor is simply `p-2` copies of (12). No shape
pigeonhole is needed.

## 3. Every prescribed first character is positive

Fix any `k in F_p^*` and define its row sum

```text
R_k(v)
 =sum_(l!=0,-k) C_(k,l).                            (14)
```

Let

```text
N_k=sum_r v_r^2 zeta^(-kr).                         (15)
```

Complete correlation in the second character gives

```text
sum_(l in F_p) M_l conjugate(M_(k+l))
 =p conjugate(N_k).                                 (16)
```

The omitted terms `l=0,-k` are both
`V conjugate(M_k)`. Therefore

```text
R_k
 =p M_k conjugate(N_k)-2V|M_k|^2.                  (17)
```

Translate the first occupied site to zero, let `d` be the nonzero
difference to the second, and put

```text
theta=2*pi*k*d/p.                                   (18)
```

An exact expansion of (17) gives

```text
Re R_k/V^3
 =(p-2)+x[-3p+4+(p-4)cos(theta)].                   (19)
```

The coefficient of `x` is at most `-2p<0`. Since `x<=1/4`,

```text
Re R_k
 >=(p-4)(1+cos(theta))V^3/4.                        (20)
```

For odd `p`, no nonzero `p`th root is `-1`. Among nonzero residues, the
smallest cosine is `-cos(pi/p)`. Thus, for every prescribed `k`,

```text
Re R_k
 >=(p-4)(1-cos(pi/p))V^3/4
 >0.                                                (21)
```

The bound is sharp for each `p`: take equal masses and

```text
k d congruent +(p-1)/2 or -(p-1)/2 mod p.          (22)
```

Equations (14) and (21) prove that for every `k!=0`, some admissible `l`
satisfies

```text
Re C_(k,l)>0.                                       (23)
```

This is not a consequence of the fixed-shape formula: it is the
orthogonal row positivity of the same `12 x 11` matrix.

## 4. Quantitative canonical-word current at every character

Use THM-2312's positive THM-2305 word `Q`, prescribed owner and clock.
For its root fibre write

```text
v_r(y)=G_j((y+r)/13),
V(y)=sum_r v_r(y),                                  (24)

C_(k,l)(Q)
 =integral_Q M_k(y)M_l(y)conjugate(M_(k+l)(y))dy.  (25)
```

Every active fibre has at most two occupied roots,

```text
integral_Q V=13rho_Q,
measure(Q)<=1/7.                                    (26)
```

Hölder gives the exact common input

```text
integral_Q V^3
 >=49(13rho_Q)^3
 =107653rho_Q^3.                                   (27)
```

### 4.1 Prescribe the slope

For every `lambda in F_13 minus {0,-1}`, equations (12) and (27) give

```text
sum_(k!=0) C_(k,lambda k)(Q)
 >=(968877/4)rho_Q^3.                              (28)
```

There are twelve terms, so some `k` satisfies the same sharp per-cell
floor as THM-2312:

```text
Re C_(k,lambda k)(Q)
 >=(322959/16)rho_Q^3.                             (29)
```

This now holds for **every prescribed slope**, not merely one unknown
shape.

### 4.2 Prescribe the first character

For every `k in F_13^*`, equations (21) and (27) give

```text
sum_(l!=0,-k) Re C_(k,l)(Q)

 >=(968877/4)(1-cos(pi/13))rho_Q^3.                (30)
```

There are eleven terms. Hence some `l` obeys

```text
Re C_(k,l)(Q)
 >=(968877/44)(1-cos(pi/13))rho_Q^3.               (31)
```

There is also a convenient rational floor. Concavity of sine on
`[0,pi/2]` gives

```text
sin(pi/(2p))>1/p,
1-cos(pi/p)>2/p^2.                                 (32)
```

At `p=13`, equations (30)--(32) yield

```text
sum_(l!=0,-k) Re C_(k,l)(Q)
 >(5733/2)rho_Q^3,                                 (33)

some l:
 Re C_(k,l)(Q)>(5733/22)rho_Q^3>0.                 (34)
```

The prescribed `k` may depend on any other lawful sidecar; only `l` is
selected after positivity.

## 5. Exact alignment with the collision colour

Suppose the same positive canonical word satisfies THM-2318's one-shot
corner hypothesis. That corner supplies a common word-current/source atom

```text
N=c13^(s-1)n,
gcd(n,91ell)=1.                                    (35)
```

Write

```text
c=13^lambda u,             13 does not divide u,

kappa_res=n mod 13,
kappa_N=u n mod 13
       =N/13^(nu_13(N)) mod 13.                    (36)
```

The first coordinate is the residual core colour after the owner is
factored; the second is the actual normalized root character of `N`.
Both are nonzero and both may be prescribed in Section 4. To align the
cubic with the Fourier atom itself, prescribe `k=kappa_N` in (34). The
resulting decorated pair

```text
(ordinary atom N, cubic current C_(kappa_N,l)(Q))   (37)
```

now shares:

- the selected source owner;
- the prescribed clock;
- the complete pure or fork word `Q`;
- the exact nonzero root character `kappa_N`.

On THM-2318's multiplier-four hostile carrier,

```text
N=13^3 n,
gcd(n,91)=1,
kappa_res=kappa_N=n mod 13.                        (38)
```

Thus the grade-three, `91`-unit collision sidecar and the cubic current can
be marked by the same exact root colour. Root-colour selection is no longer
the stopping reason.

Equation (37) is deliberately a **decorated pair**, not a factorization.
Nothing proves that `f_hat(N)` is one of the three Fourier factors inside
`C_(kappa_N,l)`, that the cubic has ordinary spatial frequency `N`, or that
either object is incident to THM-2293's `c_3` edge.

## 6. Comparison with the 91-root unit needle

The candidate THM-2319 supplies a complementary CRT object. Reduction

```text
(Z/91Z)^* -> F_13^*                                (39)
```

has six primitive lifts above every `kappa_N`. THM-2319's rational needle
saturation makes every such colour nonzero separately for the bare source
and the word-restricted source, at bounded but potentially different gauge
indices.

The present theorem selects only the mod-thirteen first character. It does
not select one of the six septimal lifts, equalize the bare/word gauge
indices, or turn THM-2319's separately nonzero coefficients into a mixed
current. THM-2319's exact negative mixed-polarization witness shows why
this last step cannot follow from separate unit-face positivity.

Thus the two gains compose only at the level of retained coordinates:

```text
THM-2321: prescribe the root colour modulo 13;
THM-2319: expose every primitive lift separately;

still open:
  one common lift and gauge incident to a bare-source shell edge.       (40)
```

This section uses THM-2319 only as a status-marked comparison, not as a
proved dependency of the slice identities.

## 7. The complete projective pair line labels every abstract target gain

The generic rectangle has a canonical projective completion. For

```text
d=[r:s] in P^1(F_p),
```

define its scalar-fibre current

```text
S_d(v)
 =sum_(t in F_p^*)
   M_(tr)M_(ts)conjugate(M_(t(r+s))).               (41)
```

If `r,s,r+s` are all nonzero, normalize `r=1`. Then (41) is the fixed-shape
sum (7), so

```text
S_d(v)>=(p-4)V^3/4.                                 (42)
```

Exactly one of `r,s,r+s` vanishes at the three remaining directions

```text
[1:0],       [0:1],       [1:-1].
```

On each, Parseval gives the same positive boundary current

```text
S_d(v)
 =M_0 sum_(t!=0)|M_t|^2
 =V(p sum_r v_r^2-V^2)
 >=(p-2)V^3/2.                                     (43)
```

Thus every direction of `P^1(F_p)` has positive aggregate current. At
`p=13`, equations (27), (42), and (43) give:

```text
generic direction:
 aggregate >=(968877/4)rho_Q^3,
 one scalar term >=(322959/16)rho_Q^3;

boundary direction:
 aggregate >=(1184183/2)rho_Q^3,
 one scalar term >=(1184183/24)rho_Q^3.             (44)
```

Now use THM-2309's canonical target quotient

```text
K/L=F_13 e_a direct_sum F_13 e_b.
```

The ordered character legs give the literal coordinate map

```text
Phi:P^1(F_13)->P(K/L),

Phi([r:s])=[r e_a+s e_b].                           (45)
```

It sends the two character axes to the pure target axes and sends
`[1:g]` to the mixed target gain `g`. Common rescaling of the root
character generator multiplies `r,s` together and disappears
projectively. Leg swap and target swap agree:

```text
[r:s]->[s:r],
g->g^(-1).                                         (46)
```

The third analytic mark `[1:-1]`, where the output character is trivial,
lands on the honest mixed gain `-1`; equation (43) proves that it is
positive rather than missing.

Consequently the exact word types have matching projective currents:

```text
pure word {a}:     choose [1:0];
pure word {b}:     choose [0:1];
fork word {a,b}:  for every g in F_13^*, choose [1:g].               (47)
```

For every algebraic target-plane direction compatible with the word,
there is therefore a nonzero cubic current carrying the same **abstract
projective label**. Abstract target-gain selection is no longer an
unstructured `132 x 12` search.

### The remaining address-incidence gap

Equation (45) is a canonical identification of two coordinate
representations. It is not a theorem that a relation address in the coset

```text
L+r e_a+s e_b                                      (48)
```

contributes to the Fourier current `S_[r:s]`. Both nontrivial cubic input
legs in (41) are formed from the same root-word function; neither is
polarized by target `a` or target `b`.

THM-2309's target-graft freedom shows why its pivot minor cannot supply
incidence by itself. Replacing one exact graft row

```text
r_(k_b)+t_b
```

by

```text
r_(k_b)+alpha t_b,
alpha mod 13 in F_13^*,                             (49)
```

leaves the owner pivot, its determinant, every non-`b` coordinate modulo
thirteen, and complete brightness unchanged, while the relative `b`-column
scale sweeps all twelve nonzero values. The packet rows remain in the zero
class `L`; this family does not construct twelve mixed addresses. It proves
that pivot/minor data alone do not make a chosen coset visible to the
current.

A lawful incidence theorem must target-polarize the two cubic legs or
produce one visible relation address in (48) whose coefficient participates
in the word current. THM-2319's exact negative mixed-polarization control
warns that such a theorem will need more than separate nonvanishing.

## 8. Sharp hostile boundaries

Both slice theorems have exact stopping examples.

1. At the threshold `p=4`, equal two-sheet masses give zero in the
   fixed-character formula for every allowed phase; an antipodal pair with
   `k=1` already has `M_k=0`. At `p=3`, two equal adjacent masses give real
   row sum `-1` after unit normalization.

2. Positive row and column sums do not make every cell positive. At
   `p=13`, take equal masses at sites `0,1` and prescribe

   ```text
   (k,l,k+l)=(1,6,7).
   ```

   The cell is exactly

   ```text
   8 cos(pi/13)cos(6pi/13)cos(7pi/13)<0.            (50)
   ```

   Therefore this real-positivity argument may prescribe `k` **or**
   `lambda` before selection, but it cannot force a positive real part in
   an arbitrary prescribed intersection cell.

3. THM-2318's hostile carrier still kills the ordinary base pair
   coefficient. The new cubic at the same character does not resurrect
   that linear coefficient or terminal component phase.

These examples prevent the row/column theorem from being misread as
entrywise total positivity or pair-edge incidence.

## 9. Connection and loss ledger

```text
source:
  THM-2312's at-most-two-sheet nonnegative root fibres;

new representations:
  the 12 x 11 (first character x projective shape) bispectrum matrix,
  and its three-marked projective completion;

operations:
  sum one fixed shape over scalar fibres, or sum one fixed first character
  over second legs before selecting a positive cell;

preserved:
  owner, prescribed clock, exact canonical word, cubic complex phase
  closure, any prescribed root colour, every compatible abstract
  target-plane direction, and quantitative margins;

new composition:
  choose THM-2318's normalized atom colour as the prescribed first
  character, or choose a THM-2315 target direction as the prescribed
  projective shape;

destroyed or unselected:
  the second character, an individual prescribed matrix cell, a common
  ordinary Fourier frequency, a common septimal lift/gauge, target-polarized
  input legs, a realized relation address in the labelled gain class,
  pair-edge incidence, and component base phase;

exact remaining object:
  a target-polarized mixed tensor making one visible target-plane relation
  address participate in the correspondingly labelled cubic current. (51)
```

The last object would turn the label map (45) into an incidence theorem.
Neither a tournament head selector nor a bare support hypergraph contains
that witness.

No scalar profile is excluded, and LRC(14) remains open.

## 10. Exact verification

The companion uses exact integer, `Fraction`, and cyclotomic group-algebra
arithmetic. It checks the `132=12*11` factorization, `47664` fixed-shape
rows, `4704` fixed-character rows, all constants in (28)--(34), the
three projective boundary shapes, the `p=3/4` controls, the negative cell
(50), and the full twelve-value target-graft gauge while the owner pivot
and every non-target coordinate stay fixed.

Reproduce with

```bash
python3 04-computation/lrc14_prescribed_character_bispectrum_thm2321.py
python3 -O 04-computation/lrc14_prescribed_character_bispectrum_thm2321.py
```

Every load-bearing check raises explicitly, so optimized mode executes the
same audit. QED.
