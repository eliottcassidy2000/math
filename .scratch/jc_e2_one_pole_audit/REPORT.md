# Independent audit: balanced `e=2,h=1` one-pole classification

Status: **PASS, with one nonfatal control-strength correction recommended.**

I independently rederived the draft from THM-2796, ran its companion in
ordinary and optimized modes, built a separate assertion-free exact audit,
and checked the new Chebyshev interpretation.  I found no mathematical
counterexample, missing class, hidden affine equivalence, sign error, or
factor-typing failure.  The claimed all-degree classification and the
degree-26 six-class conclusion are sound.

## 1. Re-derivation and signs

For `e=2,h=1`, balance gives `s=N-4`, the unique pole has multiplicity
`N`, and after translating it to zero

```text
T=x,  D=x^N,  E=x^2+ux+v,  P=SE^2.
```

Multiplying THM-2796 equation (2) by `E` gives, with no sign change,

```text
xP'-NP=CE.
```

Since `P` is monic of degree `N`, coefficient comparison gives exactly

```text
P=x^N-C[x^2/(N-2)+ux/(N-1)+v/N].
```

Putting `E=(x-1)(x-rho)` and imposing `P(1)=P(rho)=0` yields

```text
C=N(N-1)(N-2)/(N-(N-2)rho),

H_N(rho)=rho^(N-1)(N-(N-2)rho)-N rho+N-2=0.
```

The audit recovers `H_N` independently by eliminating `C` between the two
evaluation equations.  The apparent extra factor in that raw elimination
is `rho`, which is a unit on the admissible scheme.

Both possible normalized denominators are disjoint:

```text
N-(N-2)rho,              N rho-(N-2).
```

The first is checked directly by substituting `rho=N/(N-2)` in `H_N`; the
second follows by anti-reciprocity, or by direct substitution.  Neither can
vanish for `N>=4`.  Also `rho`, `rho-1`, and `C` are units.

At either root `z` of `E`, the first integral and `P(z)=0` give `P'(z)=0`.
Differentiation gives

```text
zP''(z)=CE'(z),
```

so the multiplicity is exactly two because `z,C,E'(z)` are all nonzero.
In normalized coordinates the two exact certificates are

```text
P''(1)=C(1-rho),
P''(rho)=C(rho-1)/rho.
```

Any repeated root of `P` outside `E` would make `P=P'=0`, and the first
integral would force `E=0`, a contradiction.  Thus `S=P/E^2` is monic of
degree `N-4`, squarefree, and disjoint from `E`.  Finally
`P(0)=-C rho/N` is a unit, so the pole is disjoint.  This proves the full
converse typing, not only divisibility.

## 2. Ratio scheme and simplicity

The root `rho=1` has exact multiplicity three:

```text
H_N(1)=H_N'(1)=H_N''(1)=0,
H_N'''(1)=-N(N-1)(N-2).
```

If another root were multiple, `H_N'=0` would give

```text
rho^(N-2)[(N-1)-(N-2)rho]=1.
```

Combining this with `H_N=0`, using `rho!=0`, reduces exactly to

```text
-(N-2)(N-1)(rho-1)^2=0.
```

Hence

```text
Q_N=H_N/(rho-1)^3
```

is squarefree of degree `N-3`, and every one of its roots is admissible.

The draft companion's function `reduce_coefficient_mod_q` has one generic
validity-gate weakness: checking only that a denominator is *not identically
zero* modulo a reducible squarefree `q` does not prove that it is a unit.
For example,

```text
q=(rho-2)(rho-3),       denominator=rho-2
```

passes that weaker test while being undefined on one component.  This does
**not** invalidate the present result.  Every denominator actually arising
in the draft's `P(rho)` and `E^2` remainder checks is coprime to `Q_N`; the
independent audit verifies this coefficient by coefficient using
`gcd(Q_N,denominator)=1`.  I recommend strengthening the helper from
“nonzero remainder” to a gcd/unit check before canonization.

## 3. Exact Chebyshev bridge

The new bridge is correct, including its sign.  With
`rho=exp(2 i theta)`,

```text
exp(-iN theta) H_N(exp(2 i theta))
 =2i[N sin((N-2)theta)-(N-2)sin(N theta)].
```

Meanwhile

```text
d/dtheta U_(N-2)(cos theta)
 =- [N sin((N-2)theta)-(N-2)sin(N theta)]
    /(2 sin(theta)^2).
```

Thus the equations have exactly the same interior zeros.  The polynomial
`U_(N-2)` has `N-2` simple roots in `(-1,1)`, so its derivative has exactly
`N-3` simple interlacing roots there.  They produce `N-3` distinct ratios
on the unit circle, exhausting `Q_N`.  This is a genuinely independent
all-degree proof of the root count, simplicity, and unit-circle locus.

For maximum precision the draft should display the last scalar derivative
identity, rather than saying only that differentiation “shows” the
equivalence.  Both `theta=0` and `theta=pi` represent the same discarded
ratio `rho=1`; this is merely a parameter-endpoint duplication.

## 4. Equivalence relation and class count

After the unique pole is translated to zero, every affine source
equivalence is a scaling.  Once one double zero is normalized to `1`, the
unordered set of double zeros leaves only

```text
rho -> rho               or               rho -> 1/rho.
```

The independent audit verifies coefficientwise on the full reduced ratio
scheme that the substitution `x=rho*y`, followed by monic normalization,
takes the map at `rho` to the map at `1/rho`.  Conversely, a scaling taking
one normalized unordered pair to another has one of these two forms.
There are no further affine identifications.

The three branch partitions are distinct, including at `N=4`, so allowing
the usual branch-compatible target equivalence introduces no hidden
permutation of branch values.  As a separate combinatorial check, I fixed
the pole `N`-cycle, enumerated every product of two disjoint transpositions
whose third cycle type is `(N-2,1,1)`, and quotiented by the full
centralizer of the `N`-cycle.  For every `4<=N<=14`, this independent dessin
count agrees with the ratio-orbit count.

Anti-reciprocity,

```text
rho^N H_N(1/rho)=-H_N(rho),
```

makes inversion act on the `N-3` admissible ratios.  Its only possible fixed
points are `+1` and `-1`; `+1` is discarded, while `-1` occurs exactly for
even `N`.  Burnside therefore gives

```text
# affine classes = floor((N-2)/2).
```

## 5. Split boundary and degree 26

For `N=4`, `S=1` and

```text
rho=-1,  C=4,  P=(x^2-1)^2,  V=v x^6,
```

so the unique class is split.  For every `N>=5`, either `N` is odd (the
unique pole multiplicity is odd) or `N` is even and `s=N-4>0`; THM-2796
criterion (4) therefore makes every class genuinely nonsplit.

Since

```text
deg V=s+N+2=2N-2,
```

degree `26` is exactly `N=14`.  Here `Q_14` has degree `11`, is squarefree,
has all roots on the unit circle, and factors as the claimed `(rho+1)`
times the reciprocal degree-ten factor.  There is one inversion-fixed root
and five two-element inversion orbits, hence exactly **six** affine
classes.  This is a response-layer classification only; it does not prove
Faber-flux intersection, chart entry, or JC(2).

## 6. Exact controls

Ordinary and optimized transcripts agree for both companions.

```text
draft script:
  6bbe1f4b76fbc4067ad7dce454660d7b41fa65f032fabb77456bd16230bbf06f
draft transcript:
  613eb7c15c2edaadeda5767080c233cfe1ecc00b1019bf42a0679dea5d1a346f

independent audit:
  7e08cd9d62d64ffeef440ae0f85e15aa92e2861a278e1e343234dbfc9ce472d9
independent transcript:
  cbeb0cfb92a9d7016849eb7179900346d3c69fff5e22274d9da2c63a0e750d2d
```

Reproduce with:

```bash
python3 .scratch/jc_e2_one_pole/e2_one_pole.py
python3 -O .scratch/jc_e2_one_pole/e2_one_pole.py
python3 .scratch/jc_e2_one_pole_audit/audit.py
python3 -O .scratch/jc_e2_one_pole_audit/audit.py
```

Both scripts contain zero Python `assert` nodes.
