# Uniform stability and the unique exponential closure of signed duplication

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
Separate continuation of the frozen
[sharp all-degree theorem](signed_degree5_empty_core_next_sep06_uniform.md).
No theorem ID or external-priority claim.

## 1. A degree-independent sequential equivalence

Let `r^(j)` be finite real lists, of possibly varying length, satisfying

```text
p_1=sum_i r_i=1,       p_2=sum_i r_i^2=1,
G_j(s)=prod_i(1+r_i s),
E_j=e_2(r_1^2,...,r_n^2)=(1-p_4)/2>0,
D_j=-[s^4]G_j(s)^2=(5-8p_3+3p_4)/6,
c_*=(13-8sqrt(2))/3.                                      (1)
```

Zero entries and both signs are allowed. The following are equivalent:

1. `D_j/E_j -> c_*`.
2. After a permutation of each list and padding with zeros,
   `r_1,r_2 -> 1/sqrt(2)` and `sum_(i>=3) r_i^2 ->0`.
3. Locally uniformly for complex `s`,

```text
G_j(s) -> (1+s/sqrt(2))^2 exp[(1-sqrt(2))s].                (2)
```

In this case `E_j->1/4`, `D_j->c_*/4`, and the dust has the prescribed
signed first moment

```text
sum_(i>=3) r_i -> 1-sqrt(2).                               (3)
```

In particular the number of nonzero dust entries tends to infinity.
No conclusion that the dust has only one sign or equal entries is valid.

Equivalently, (1) has uniform qualitative stability: for every `epsilon>0`
there is a `delta>0`, independent of the length, such that
`D/E<c_*+delta` forces, after permutation,

```text
(r_1-1/sqrt(2))^2+(r_2-1/sqrt(2))^2+sum_(i>=3)r_i^2
  < epsilon^2.                                            (4)
```

This statement supplies an existence modulus, not an explicit rate.

## 2. Inheritance, target, and the lost coordinate

The inherited mechanism is the exact integer-multiplicity theorem in the
linked all-degree note: each finite-degree minimum has two equal positive
reciprocal roots and equal negative dust. Its sharp uniform constant is an
unattained infimum. The closest hostile boundary is `(1,0,...)`, where the
energy vanishes. The new target is rigidity of arbitrary near minimizers,
not only the explicitly minimizing rows. The least-used sidecar is the
first moment of roots that disappear in square norm.

The map is from the normalized cancellation sphere `p_1=p_2=1` to its
closure under coordinatewise convergence after sorting by absolute value.
Third and fourth moments survive this map; the first moment need not.
The nonzero lost first moment is kept separately and becomes the linear
exponential factor. The live board is: actual square coefficient; root
energy; compactness with variable degree; the one-root zero-energy orbit;
mixed-sign dust; exponential closure. This map is confined to real-rooted
ordinary cores after normalization. It supplies no generic positivity or
transport statement for LRC masks or arbitrary Laurent rows.

The apparent shortcut "the limiting moment gap vanishes, hence there are
two atoms" is false. That gap also vanishes at one positive atom. Section
4 retains the ratio to exclude this boundary. Section 6 gives a separate
exact hostile family with dust of both signs.

## 3. A sharp sphere inequality with two equality types

Put `A=2-sqrt(2)` and `B=sqrt(2)-1=1-A`. For every finite or countable real
sequence with `sum r_i^2<=1`,

```text
H(r):=sum_i(r_i^3-A r_i^4) <= B.                           (5)
```

Equality holds precisely at permutations of

```text
(1,0,...),       (1/sqrt(2),1/sqrt(2),0,...).              (6)
```

All sums in (5) converge absolutely. Replacing a negative entry by its
absolute value strictly increases `H`, so equality requires nonnegative
entries. Write `x_i=r_i^2`, and `h(x)=x^(3/2)-A x^2`. For `0<x<=1`,

```text
h'(x)=sqrt(x)(3/2-2A sqrt(x))>0.
```

If the total square mass is less than one, append its positive square
root; this strictly increases `H`. It therefore suffices first to maximize
`sum h(x_i)` on each finite probability simplex.

A maximizer cannot have at least three positive coordinates. Its two
smallest positive coordinates `x,y` satisfy `x+y<=2/3`. Along the feasible
redistribution `(x,y)->(x+t,y-t)`, the second derivative is

```text
h''(x)+h''(y)
  = (3/4)(1/sqrt(x)+1/sqrt(y))-4A
 >= (3/2)sqrt(3)-4A >0.                                  (7)
```

For the last strict sign, `A<3/5` and `12/5<(3/2)sqrt(3)` suffice.
A twice differentiable function cannot have this positive second
derivative at an interior local maximum.

One-coordinate support gives `h(1)=B`. For two-coordinate support, set
`t=sqrt(x)+sqrt(1-x)`, so `1<=t<=sqrt(2)`. Exact elimination gives

```text
B-h(x)-h(1-x)
  = (t-1)^2 (sqrt(2)-t)(A t+1)/2 >=0.                     (8)
```

Its only equality cases are the single-coordinate endpoint or `x=1/2`.
This proves the finite inequality and classification.

For countably many entries, apply the finite inequality to truncations
and pass to the absolutely convergent limit. Equality cannot have square
mass less than one, since appending the missing mass would contradict the
inequality. Nor can equality have infinitely many positive coordinates:
choose two sufficiently small positive square coordinates with sum at
most `2/3` and use the same finite redistribution (7). The remaining
finite-support equality cases have already been classified. This proves
(5)--(6) without a first-moment constraint.

## 4. Compactness keeps the moments; the quotient selects two atoms

Define `g=B-p_3+A p_4>=0`. The exact identity is

```text
6(D-c_* E)=8g,      D/E-c_* =4g/(3E).                     (9)
```

This already gives an independent elementary proof of the strict uniform
bound `D/E>c_*`. Equality in (9) would require one of (6); the one-atom
case has `E=0`, and the two-atom case has `p_1=sqrt(2)`, contrary to the
current normalization. The mixed-sign family in Section 6 proves
sharpness of the uniform constant. The finite-degree theorem is therefore
inheritance and context, not a required dependency for the present
uniform bound or closure result.

Consequently `D/E->c_*` implies `g->0`, since `E<=1/2`.
Sort each list by nonincreasing absolute value and pad with zeros. Then
`|r_k|<=1/sqrt(k)`. A diagonal subsequence converges coordinatewise to a
sequence `y` with `sum y_i^2<=1`. Uniform tail estimates give

```text
sum_(i>K)|r_i|^3 <=1/sqrt(K+1),
sum_(i>K)r_i^4   <=1/(K+1).                               (10)
```

Thus `p_3,p_4` converge to the corresponding moments of `y`; (9) forces
`H(y)=B`. Section 3 leaves just the two alternatives (6). Each has full
square norm one, so coordinatewise convergence is also convergence in
square norm after these permutations.

The one-atom alternative is excluded by the quotient. Suppose its first
coordinate `a->1`, set `sigma^2=sum_(i>=2)r_i^2=1-a^2`, and write
`u_3=sum_(i>=2)r_i^3`, `u_4=sum_(i>=2)r_i^4`. Positive energy ensures
`sigma>0`; the tail bounds are `|u_3|<=sigma^3`, `0<=u_4<=sigma^4`.
Then

```text
D/E = [5-8a^3+3a^4-8u_3+3u_4]
      /[3(1-a^4-u_4)].                                   (11)
```

The scalar numerator factors as

```text
5-8a^3+3a^4=(1-a)(-3a^3+5a^2+5a+5).
```

Divide numerator and denominator of (11) by
`sigma^2=(1-a)(1+a)`. The scalar numerator tends to `6`, its tail terms
vanish, and the denominator `3(1+a^2-u_4/sigma^2)` tends to `6`.
Hence `D/E->1`, contradicting `c_*<1`.

Every subsequence of a near-minimizing sequence therefore has a further
subsequence converging in square norm to the two-atom configuration.
This proves condition 2 along the whole sequence (use the sorted order,
or choose minimizing permutations). Conversely condition 2 gives
`p_3->1/sqrt(2)`, `p_4->1/2`, and (1) gives condition 1. The uniform
formulation (4) follows by contradiction from this sequential equivalence
and the strict lower bound `D/E>c_*` just proved.

The signed dust sum (3) follows from `p_1=1`. If it has `N` nonzero
entries and square mass `eta^2`, Cauchy--Schwarz gives
`N eta^2 >= (sum dust)^2 ->(sqrt(2)-1)^2>0`. Since `eta->0`, necessarily
`N->infinity`. This retains the signed first moment without imposing
either a sign partition or a multiplicity pattern on the dust.

## 5. The polynomial closure is locally uniform

Under condition 2, the largest absolute dust entry tends to zero. For a
fixed complex disk `|s|<=R`, eventually every dust factor satisfies
`|r_i s|<=1/2`. The analytic logarithm near one then obeys

```text
|log(1+r_i s)-r_i s| <= R^2 r_i^2.
```

Sum over the dust, use its vanishing square mass and (3), and exponentiate.
The dust product converges locally uniformly to
`exp[(1-sqrt(2))s]`; the two remaining factors give (2). Conversely local
uniform convergence gives convergence of the first four coefficients by
Cauchy's formula. Newton's identities then give the same limits of
`p_3,p_4`, so `E->1/4>0` and `D/E->c_*`. Thus all three statements are
equivalent, including their actual polynomial consequence.

For an unnormalized finite real list with `e_2(r)=0` and positive energy,
`p_1^2=p_2>0`; divide all entries by `p_1` to reach (1). A negative sum
therefore flips all signs before interpreting the two positive atoms.
For a real ordinary core `a_0 prod(1+r_i s)`, divide by `a_0` and scale
the variable by `1/p_1`. Monomial factors must first be removed with the
coefficient index shifted accordingly. A zero reciprocal entry is a
degree drop (an ordinary root at infinity), not an ordinary root at zero.
Complex gauges require this real normalization before a signed ratio is
assigned, as in the inherited degree-five and Laurent-transport notes.

## 6. Exact hostile controls: mixed signs, zeros, and the energy boundary

For an integer `L>=6`, put

```text
kappa=5/L,
beta=2/[2+sqrt(2+2kappa)],
u_L=(1,1, [-2beta/L repeated L], [beta/L repeated L]).      (12)
```

The dust has sum `-beta` and square sum `kappa beta^2`. The defining
quadratic `(1-kappa)beta^2-4beta+2=0` gives `e_2(u_L)=0` exactly.
Normalize by the positive sum `2-beta`. Since `beta->2-sqrt(2)`, this
family satisfies condition 2 and hence approaches `c_*`. Yet its positive
dust sum tends to `sqrt(2)-1`, while its negative dust sum tends to
`-2(sqrt(2)-1)`. Thus even positive dust with a nonvanishing total first
moment is compatible with sharp near equality. Any number of zero
entries may be appended without changing the statement or the ratio.

Separately, the lists `(1,t,-t/(1+t))` for positive `t->0` have exact
`e_2=0`, positive energy, and `D/E=1`. After normalization they approach
the one-atom zero-energy boundary and their gap `g` tends to zero. This
is an exact hostile control against dropping the denominator in (9).

## 7. Reproduction and audit manifest

[Source](../../04-computation/signed_uniform_stability_empty_core_next_sep06.py)
and [output](signed_uniform_stability_empty_core_next_sep06.out).

```bash
python3 -B 04-computation/signed_uniform_stability_empty_core_next_sep06.py
python3 -B -O 04-computation/signed_uniform_stability_empty_core_next_sep06.py
```

The declared universe contains six sphere controls, the five exact
mixed-sign rows `L=6,10,25,100,400`, and the four one-atom boundary rows
`t=1/2,1/10,1/100,1/1000`. All **99 explicit gates** pass. The source
checks the entire two-support polynomial factorization in `Q(sqrt(2))`,
the curvature signs, both sphere equality types and their hostile signs,
the exact mixed-sign cancellation, complete energy and ordinary square
coefficient, zero-entry padding, and the normalized limiting square
coefficient. Compressed binomial products are cross-checked against
factor-by-factor multiplication for the first three mixed rows. Both
normal and optimized execution reproduce exactly; no assertion removal
can suppress the gates. These are bounded algebraic controls, not a
census or a substitute for the analytic compactness proof.

The root agent independently audited the simplex curvature, the exact
two-support factorization, countable equality classification, third- and
fourth-moment tail continuity, ratio limit at the one-atom boundary,
whole-sequence argument, Cauchy--Schwarz dust count, local-uniform
logarithm argument and converse, and the mixed-sign tuning quadratic:
**PASS**. Independent referee `certificate_audit` separately checked the
full analytic argument, reconstructed both decisive factorizations, the
degree-three ratio-one hostile, and both mixed-dust first-moment limits.
They also audited the exact quadratic arithmetic and compressed/literal
coefficient paths, independently replayed all 99 gates in normal and
optimized modes, and verified byte-for-byte output agreement and both
raw hashes: **PASS**, with no mathematical correction. The proof does not claim an
explicit stability rate or a unique dust sign pattern.

```text
source SHA256 08ff85688ca267cbb19829be1946c8c663d4823f22342757f970196e41632b0b
output SHA256 1e9299f556fa60bc06a4459b143a7667ecfd4fca14050d7e7a7ebffd68a63c2c
trace  SHA256 edbb2b086b5981ffcd7f1fb0befd197550e2ee628f20862b30aa034441667167
```
