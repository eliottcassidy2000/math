# Fixed Berggren drift doubles the Fibonacci clock and harmonic averaging smooths parity

**Status: PROVED in THM-3497 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**
The companion is
`04-computation/berggren_fixed_drift_word_language_parity_probe_20260816.py`.
This is a finite branch-language theorem, not a physical LRC
current or a Jacobian statement.

## 1. The quantifier ladder

The variable-translation language asked whether, for each word `w`, *some*
translation could make its projective action affine with the true composite
linear part.  THM-3339 supplies more structure: one fixed calibrated branch
cocycle

```text
Atilde(u)=J u+p,       Btilde(u)=J u+p,       Ctilde(u)=u+p. (1)
```

Write `rho_fix(w)` for the affine permutation obtained by composing (1) in
root-to-child order.  Call `w` **fixed-drift calibratable** when there is a
point gauge `f_w:P1(F_3)->V4` satisfying

```text
f_w rho_3(w) f_w^(-1)=rho_fix(w).                         (2)
```

The point gauge may still depend on the whole word, but the target drift no
longer does.  The four relevant quantifier levels are therefore

```text
same 24-point atlas;
wordwise gauge and wordwise translation;
wordwise gauge with the fixed THM-3339 drift;              (this note)
one simultaneous branch-monoid gauge.                     (impossible) (3)
```

Keeping these levels separate resolves the apparent conflict between the
positive word census and the generatorwise no-go.

## 2. The exact 192-state automaton

In the pinned order `(0,p,q,r)`, put

```text
G(u)=Ju+p=(0 1 3 2),              T(u)=u+p=(0 1)(2 3).    (4)
```

Then

```text
G^4=T^2=1,                        TGT=G^(-1),              (5)
```

so `<G,T>` is THM-3339's order-eight affine `D4`.  A word has the pair state

```text
(rho_3(w),rho_fix(w)) in S4 x D4,                         (6)
```

and (2) holds exactly when the two components have the same cycle type.

The three paired generators generate the full direct product:

```text
image=S4 x D4,                       |image|=192.          (7)
```

The exact search finds all `24` pure-source states and all `8` pure-target
states, proving (7), not merely surjectivity of the two projections.
Partition refinement proceeds through

```text
9, 85, 172, 192, 192                                      (8)
```

blocks.  Thus every state is distinguishable: the DFA is minimal and its
syntactic monoid is the regular `S4 x D4` action.

The `B` spine remains a hostile positive control.  Since `B_3=G` in the
natural four-point gauge,

```text
B^n passes (2) for every n>=0.                            (9)
```

The fixed drift narrows the rest of the language; already at levels one and
two only `B` and `BB` pass.

## 3. The bipartite character and parity-split densities

The abelianization of `D4` has a character

```text
chi(G)=chi(T)=1.                                          (10)
```

Every branch letter therefore flips `chi`, so the 192 states split into two
classes of 96 and level `n` lies entirely in class `n mod 2`.  The accepting
state counts are

```text
even class: 16 of 96;
odd class:  18 of 96.                                    (11)
```

The nine two-letter increments generate the full 96-state kernel of `chi`.
Moreover, `B^4` and `C^6` are identity returns of two-step lengths two and
three.  Thus the two-step walks are irreducible and aperiodic, hence mixing,
on their respective classes.  Consequently,
if `a_n` counts passing depth-`n` words, then

```text
a_(2m)/3^(2m)       -> 1/6,
a_(2m+1)/3^(2m+1)   -> 3/16.                              (12)
```

The exact first terms are

```text
1,1,1,4,13,46,113,421,1121,3667,9801,33166,...           (13)
```

They satisfy the order-fourteen recurrence whose nonzero lags are

```text
a_n=7a_(n-2)+7a_(n-4)+71a_(n-6)+213a_(n-8)
       +189a_(n-10)+1701a_(n-12)-2187a_(n-14).            (14)
```

The companion verifies more than 192 consecutive residuals.  Because the
residual itself is a scalar observation of the same 192-dimensional transfer
operator, Cayley--Hamilton promotes that finite gate to all `n`.

The generating-function denominator factors as

```text
(x-1)(x+1)(3x-1)(3x+1)(3x^2+1)
 (3x^2-2x+1)(3x^2+2x+1)(9x^4-2x^2+1).                   (15)
```

The normalized polar coefficients

```text
lim_(x->1/3)(1-3x)F(x)=17/96,
lim_(x->-1/3)(1+3x)F(x)=-1/96
```

are not analytic residues.  Every other characteristic root has modulus at
most `sqrt(3)`.  Hence

```text
a_n=(17/96)3^n-(1/96)(-3)^n+O(3^(n/2)),                 (16)
```

which proves (12) with an exact error scale.

## 4. Natural density fails, harmonic density survives

Use the same shortlex enumeration of the ternary tree as in the
variable-translation result.  Complete levels dominate the counting
function.  Summing the alternating densities geometrically gives two
different endpoint limits:

```text
lim_(m->infty) |S_fix intersect [1,N_(2m)]|/N_(2m)=11/64,
lim_(m->infty) |S_fix intersect [1,N_(2m+1)]|/N_(2m+1)=35/192. (17)
```

Thus `S_fix` has no ordinary natural density.

Harmonic weighting behaves differently.  The two-step gate above gives
uniform equidistribution from every prefix state in its 96-state class.  A
fixed prefix is a triadic rank cylinder, so cylinder step-functions
approximate the continuous level weight `t -> 1/(1/2+t)`.  The harmonic mass
of one complete shortlex level therefore tends to

```text
(1/6)log 3       on even levels,
(3/16)log 3      on odd levels.                           (18)
```

Two-level Cesaro averaging and the fact that one partial level lies in an
integer interval of endpoint ratio at most three, hence contributes only
`O(1)`, then give

```text
sum_(m<=N, m in S_fix) 1/m
       =(17/96)log N+o(log N).                            (19)
```

This is an explicit instance where a subset of the harmonic series has a
stable logarithmic density although its underlying subset of natural numbers
has no natural density.  Harmonic averaging preserves the mean of the two
parity chambers and destroys which chamber occurred at a particular depth.

## 5. Restriction to the Fibonacci three rays

On THM-3339's addresses

```text
(BA)^r,                 A(BA)^r,                 C(BC)^r, (20)
```

the fixed-drift gate is

```text
(BA)^r passes       iff r=0 mod 4;
A(BA)^r passes      iff r=3 mod 4;
C(BC)^r passes      iff r=3 mod 4.                       (21)
```

Translated to the consecutive Fibonacci index,

```text
k mod 12 in {0,1,2}.                                     (22)
```

The variable-translation language gave the same three-residue window modulo
six.  Fixing the drift doubles that clock to twelve.  This is an exact finite
echo of THM-3487's six-versus-twelve cycle repair, but it is not an
identification of the two bundles: here gauges vary by composite word, while
THM-3487 studies a base-preserving connection on one periodic path.

On Fibonacci indices, (22) gives

```text
sum_(k<=N, k mod 12 in {0,1,2}) 1/k
       =(1/4)log N+C+o(1).                               (23)
```

As before, reciprocal Fibonacci values and reciprocal triple coordinates
converge exponentially; (23) concerns the index subset only.

## 6. What the drift preserves and destroys

- source: the projective composite action `rho_3(w)`;
- target: THM-3339's fixed affine composite `rho_fix(w)`;
- map: a word-dependent point conjugacy `f_w`;
- preserved predicate: full four-point cycle type;
- destroyed information: the chosen conjugacy and its compatibility with
  neighboring words;
- retained sidecar: the full `D4` drift, not just the linear bit;
- cheapest hostile: `A` and `C` fail, while every `B^n` passes;
- global boundary: no common `f` exists for all three generators.

The fixed language is therefore a precise transport of the drift instance,
not a flattening of it.  Its mod-twelve and harmonic shadows are genuine; a
physical LRC word-current, D5 cohomology map, or Jacobian invariant would
still require an independently typed source/target map.
