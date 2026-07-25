---
id: THM-2327
title: "Two-colour triangle for marked unit c3 incidence"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. On every
  canonical middle-owner profile and for every nonzero root character,
  there is a c_3-multiple edge with multiplier coprime to 91 incident
  to a bare-source vertex carrying the same positive-word coefficient.
  The first edge is obtained at modulus 13D from the gcd-normalized
  bare/word carrier and has multiplier 1 modulo 13. THM-2326 gives an
  incident edge whose multiplier is nonzero modulo 7. If neither edge
  is already a 91-unit, their third side is. The resulting exact
  word/deepest-comb/bare cross-bispectral monomial is nonzero. The
  theorem does not select a THM-2315 target-plane gain or terminal
  component phase, excludes no scalar profile, and does not prove
  LRC(14).
source: codex-2026-07-25-two-colour-marked-unit-triangle
depends_on:
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2326-vertexwise-septimally-primitive-c3-degree
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2315-target-gain-corolla-and-pair-shadow-kernel
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
script: 04-computation/lrc14_two_colour_marked_unit_triangle_thm2327.py
output: 05-knowledge/results/lrc14_two_colour_marked_unit_triangle_thm2327.out
script_sha256: accf57d1160adcff66750228c6f821ff2609d93919ca12f1e45554188bf85820
output_sha256: c0dc2b32e541180fdc06d7c13cb933338b27ffeb1154e0e94ed51de820aa6475
hash_basis: working-tree bytes (LF)
---

# THM-2327 -- two-colour marked unit `c_3` triangle

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

The middle-owner obstruction had become a quantifier mismatch. THM-2323
marks a bare Fourier vertex by the full positive word. THM-2326 gives
every such vertex an incident `c_3` edge whose multiplier avoids seven.
THM-2293 gives a `91`-unit edge somewhere, but not necessarily at the
marked vertex.

The missing move is not another average over the graph. It is a
three-vertex colour completion. A gcd-normalized root needle produces a
marked edge whose multiplier avoids thirteen. The modulated deepest-comb
identity produces an incident edge whose multiplier avoids seven. If the
first edge is still divisible by seven and the second is still divisible
by thirteen, their third side avoids both primes.

```text
                 t: 13 does not divide t
       marked A ------------------------- marked B
                \                         /
                 \ s: 7 does not divide s
                  \
                   bare C

one of t, s, or s-t is a unit modulo 91.             (1)
```

At least one endpoint of the selected unit edge is one of the marked
vertices `A,B`. This closes the unit-colour/word-mark incidence gap on
the canonical middle-owner carrier. It does not yet connect that edge to
a prescribed target-plane gain or terminal-component phase.

## 1. The abstract two-colour completion

Let `t,s` be nonzero integers such that

```text
13 does not divide t,
7 does not divide s.                                (2)
```

Then at least one of

```text
t, s, s-t                                           (3)
```

is coprime to `91`. Indeed:

- if `7` does not divide `t`, use `t`;
- otherwise, if `13` does not divide `s`, use `s`;
- in the remaining case, `7|t` and `13|s`, while

  ```text
  s-t congruent s !=0 mod 7,
  s-t congruent -t!=0 mod 13.                       (4)
  ```

This is the Chinese-remainder triangle lemma. More generally, the same
argument works for any two coprime primes. Its carrier is an undirected
edge-coloured triangle, not a tournament: orienting the edges would add
no information, while forgetting either divisibility colour destroys the
proof.

## 2. The gcd-normalized middle-owner carrier

Work in the canonical strict LRC(14) setup with middle owner `E=E_2` and
deepest blocker speed `c_3`. Write

```text
c_2=13^b u_2,
c_3=13^c u_3,                  c>b,

g=gcd(c_2,c_3),
a=c_2/g,
D=c_3/g.                                           (5)
```

Then

```text
gcd(a,D)=1,
13 does not divide a,
13 divides D.                                      (6)
```

Let `Q` be any positive canonical word supplied by THM-2305 and put

```text
E_Q=E intersection T^(-(b+1))Q.                    (7)
```

Thus `E_Q` is a positive-measure literal subset of `E`; a nonzero
coefficient of `1_(E_Q)` retains the complete word ancestry. Normalize
by the common speed factor rather than by all of `c_2`:

```text
F_Q=P_g 1_(E_Q),
F_E=P_g 1_E.                                       (8)
```

Positivity of the Perron operator gives

```text
0<=F_Q<=F_E,
integral F_Q=measure(E_Q)>0.                        (9)
```

The support is not generally one interval. It lies in the `a`-tooth
pullback

```text
D_a={x:||a x||<1/14}.                              (10)
```

Indeed, if the inverse branch `(x+r)/g` contributes to (8), then it lies
in `E subset D_(c_2)=D_(ga)`, and hence

```text
||a(x+r)||=||a x||<1/14.                           (11)
```

This normalization is the decisive quotient repair. Normalizing by
`c_2` puts the spectrum on the sublattice of `c_2`-multiples and forces
every `c_3`-edge multiplier to contain the quotient `a`. Normalizing by
`g` retains the missing coset coordinate.

The jump ledger is unchanged. Write

```text
g=13^b v,                  13 does not divide v.
```

THM-2319 proves

```text
G=P^b 1_E                       has at most 2S jumps,
P^b 1_(E_Q)=G(1_Q after T)      has at most 6S jumps. (12)
```

Since `P_g=P_v P^b` and Perron transport cannot increase the number of
nonzero image jumps,

```text
J_E:=#Jump(F_E)<=2S,
J_Q:=#Jump(F_Q)<=6S,
L:=J_E J_Q<=12S^2.                              (13)
```

## 3. Automorphic straightening of the toothpick needle

The fixed-colour argument of THM-2323 extends from one danger arc to
`D_a`. Let `N` be any positive integer with

```text
gcd(a,N)=1,                                        (14)
```

and let `f,h` be nonzero rational step functions satisfying

```text
0<=f<=h,
support(f),support(h) subset D_a.                   (15)
```

Expose the `N` root branches and freeze a primitive character `k`. The
physical cross-correlation is

```text
C_k=sum_(d mod N)c_d zeta_N^(kd),                  (16)
```

where every `c_d` is nonnegative rational and

```text
c_0=N integral f h>0.                              (17)
```

If `c_d` is nonzero, two contributing points lie in `D_a`. After applying
the map `x -> ax`, both lie in the single open danger arc `D_1`.
Consequently the least signed residue `e_d` satisfying

```text
e_d congruent a d mod N
```

obeys

```text
|e_d|<N/7.                                         (18)
```

Suppose (16) vanished. Because both `a` and `k` are units modulo `N`, the
Galois automorphism

```text
zeta_N -> zeta_N^(a k^(-1))                        (19)
```

would give

```text
0=sum_d c_d zeta_N^(a d).
```

Every term on the right has argument of absolute value less than
`2*pi/7<pi/2`, and (17) supplies a positive term. Its real part is
therefore strictly positive, a contradiction.

Parseval and the endpoint-product Vandermonde argument of THM-2323 now
give, for every primitive `k` and every integer block origin `H`, some

```text
H<=h<=H+J_f J_h-1
```

such that

```text
f_hat(k+Nh)h_hat(k+Nh)!=0.                         (20)
```

Neither `7|N` nor a totient lower bound is needed. Multiplication by `a`
is the automorphism that straightens every tooth of (10) into the same
acute phase sector.

## 4. A marked edge avoiding thirteen

Apply Section 3 to (8) with

```text
N=13D.                                             (21)
```

By (6), `gcd(a,N)=1`. Fix any desired nonzero middle-shell root character
`kappa modulo 13`. Since

```text
g=13^b v_g,                 13 does not divide v_g,
```

choose a unit class `K_0 modulo D` satisfying

```text
v_g K_0 congruent kappa mod 13,                    (22)
```

and take its representative `1<=K_0<D`. Such a class exists by CRT:
use the nonzero class in (22) at thirteen and any unit class at every
other prime power dividing `D`.

Put

```text
K_1=K_0+D.                                         (23)
```

The prime divisors of `N=13D` are exactly those of `D`. Thus both `K_0`
and `K_1` are primitive modulo `N`; they are distinct, lie in
`{1,...,N-1}`, and have the same nonzero residue modulo thirteen.

Apply (20), with block origin zero, to each `K_i`. There are

```text
0<=h_i<=L-1,
q_i=K_i+N h_i                                     (24)
```

such that both `F_E` and `F_Q` have nonzero coefficient at `q_i`. By the
Perron Fourier identity, the physical atoms

```text
A=g q_0,
B=g q_1                                            (25)
```

satisfy

```text
(1_E)_hat(A)(1_(E_Q))_hat(A)!=0,
(1_E)_hat(B)(1_(E_Q))_hat(B)!=0.                   (26)
```

Both are exact grade-`b` atoms of root character `kappa`. Moreover,

```text
B-A
 =g[D+13D(h_1-h_0)]
 =t c_3,

t=1+13(h_1-h_0).                                   (27)
```

In particular,

```text
t!=0,
13 does not divide t,
|t|<=13L-12<=156S^2-12.                            (28)
```

Thus `AB` is a same-word marked `c_3` edge carrying the first complementary
colour. It may still have `7|t`.

## 5. Complete the unit edge at the marked vertex

Apply THM-2326 to the marked bare atom `A`. Since

```text
nu_13(A)=b<c=nu_13(c_3),
```

there is an integer

```text
1<=s<=14S-1,
7 does not divide s                                (29)
```

and a bare-source atom

```text
C=A+s c_3,
(1_E)_hat(C)!=0.                                   (30)
```

The atoms `A,B,C` have the same grade and root character. Apply the
two-colour lemma:

```text
if 7 does not divide t:
    choose X=A, Y=B, m=t;

else if 13 does not divide s:
    choose X=A, Y=C, m=s;

else:
    choose X=B, Y=C, m=s-t.                         (31)
```

In every case,

```text
Y=X+m c_3,
gcd(m,91)=1,                                       (32)

(1_(E_Q))_hat(X)!=0,
(1_E)_hat(X)(1_E)_hat(Y)!=0.                       (33)
```

Thus `X` is a common bare/word vertex and `XY` is a canonical
unit-coloured `c_3` edge. In the first case both endpoints retain the
word coefficient. In the other cases only the displayed endpoint `X` is
proved word-marked.

The uniform bound is

```text
|m|
 <=|t|+s
 <=13L+14S-13
 <=156S^2+14S-13.                                  (34)
```

No condition on `a`, no bound on the unit cofactor of `D`, and no scalar
row restriction remains.

## 6. The exact mixed Fourier triangle

Let

```text
D_(c_3)={x:||c_3 x||<1/14}.
```

Its exact Fourier coefficient is

```text
(1_(D_(c_3)))_hat(m c_3)
 =sin(pi*m/7)/(pi*m)                               (35)
```

for nonzero `m`. Equation (32) makes (35) nonzero. Therefore (32)--(33)
give the heterogeneous cross-bispectral monomial

```text
B_(E_Q,D_(c_3),E)(X,m c_3)

 :=(1_(E_Q))_hat(X)
   (1_(D_(c_3)))_hat(m c_3)
   conjugate((1_E)_hat(Y))

 !=0,                     X+m c_3=Y.               (36)
```

This is genuine ordinary-frequency incidence, not merely residue-label
availability or an abstract graph edge. It is one nonzero monomial, not
a claim that a sum of mixed bispectral terms has positive real part.
The cancellations in THM-2319's mixed-polarization hostile example remain
valid.

## 7. What this closes and what remains

The marked middle-owner ledger is now:

```text
THM-2323:
  a complete positive word and the bare source coexist at one Fourier
  atom of every prescribed primitive colour;

gcd-normalized 13D needle:
  two such marked atoms form a c_3 edge avoiding thirteen;

THM-2326:
  every marked atom has an incident c_3 edge avoiding seven;

THM-2327:
  a three-vertex colour circuit forces a 91-unit c_3 edge incident to a
  word-marked bare atom, and hence the nonzero mixed monomial (36).     (37)
```

This removes the explicit `gcd(a,91)=1` boundary from THM-2323's marked
incidence problem. The old conditional branch remains stronger where it
applies because it makes both endpoints of the selected unit edge carry
the same word.

The following are not consequences:

- (36) does not choose a projective target gain in the THM-2315 corolla;
- it does not retain the terminal-component phase/current of THM-2303;
- it does not identify a full-lattice relation address from THM-2325
  with this analytic frequency triangle;
- it does not exclude any of the `165` live scalar profiles.

The remaining object is therefore narrower than “find a unit edge.” It
is an address-decorated, target-polarized mixed Fourier circuit whose
source leg is the marked endpoint in (36). A tournament shadow is too
coarse: it forgets the two prime colours, the undirected third-side law,
and the word mark that make (31) work.

LRC(14) remains open.

## 8. Exact verification

The companion exhausts `419328` signed two-colour residue pairs,
constructs primitive `K_0,K_0+D` lifts across representative deep
cofactors (including `7|a` controls), checks the Galois multiplier and
physical `g`-normalized identity, verifies `3175200` quantitative bound
rows and `53280` grade/root rows, and checks the exact deepest-comb
nonvanishing predicate. All load-bearing tests raise explicitly under
ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_two_colour_marked_unit_triangle_thm2327.py
python3 -O 04-computation/lrc14_two_colour_marked_unit_triangle_thm2327.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_two_colour_marked_unit_triangle_thm2327.out
```

byte-for-byte after LF normalization.
