---
id: THM-1091
title: Deep-sheet resonance-ramification identity
status: PROVED — exact normalized finite-cyclic Fourier identity, THM-769 containment corollary, and finite-cover character-energy law; no nonzero-mode cancellation bound or deep-branch emptiness is claimed
source: codex-2026-07-18-S67 global-bridge audit and proof
depends_on: [THM-769, THM-1006]
related: [THM-761, THM-774, THM-994, THM-1072, THM-1090, THM-1075, THM-1092, HYP-6820]
---

# THM-1091 — Deep-sheet resonance-ramification identity

The exact information missing from sheet capacity is the Fourier phase of the
individual sheet masks.  This theorem makes that assertion literal.  A
ramified runner has a short sheet period, hence Fourier support in the
annihilator of that period.  The simultaneous uncovered-sheet count is the
zero-sum convolution of those supported transforms.

Everything below is finite.  There is no convergence or summation-order issue.

## 1. Statement and normalization

Fix an integer `s>=2` and write

```text
G_s = Z/sZ,                 e_s(x) = exp(2 pi i x/s).
```

Let `I` be a nonempty finite index set, let `w_i` be an integer not divisible
by `s` for every `i in I`, and fix `tau in R`.  Define the closed danger mask
and its Boolean safe complement by

```text
M_i(tau) = {j in G_s : ||w_i(tau+j)/s|| <= 1/13},
f_i^tau(j) = 1_(G_s minus M_i(tau))(j).                       (1)
```

Put

```text
g_i = gcd(w_i,s),           D_i = s/g_i >= 2.                 (2)
```

For a function `f:G_s -> C`, use normalized Haar Fourier transform

```text
fhat(r) = (1/s) sum_(j in G_s) f(j)e_s(-rj),
f(j)    =       sum_(r in G_s) fhat(r)e_s(rj).                (3)
```

Divisibility `g_i | r_i` for a residue `r_i in G_s` means that `r_i` belongs
to the well-defined subgroup `g_i Z/sZ`.  Define the ramified resonance module

```text
R_s(w) = {(r_i) in product_i G_s :
          sum_i r_i = 0 in G_s and g_i | r_i for every i}.    (4)
```

Then:

> **Deep-sheet resonance-ramification identity.**  For every `tau in R`,
>
> ```text
> (1/s)|G_s minus union_i M_i(tau)|
>   = sum_(r in R_s(w)) product_i fhat_i^tau(r_i).             (5)
> ```
>
> In particular, if the masks cover every sheet, then
>
> ```text
> sum_(r in R_s(w), r != 0) product_i fhat_i^tau(r_i)
>   = - product_i (1-|M_i(tau)|/s).                           (6)
> ```

Every factor on the right of (6) is strictly positive.  Thus an exact sheet
cover cannot be certified by the zero modes alone: supported nonzero
resonances must cancel their positive product exactly.

There is an immediate arithmetic obstruction.  Let `H_i=g_i Z/sZ`, and put
`d=gcd_i g_i`.  The sum homomorphism `product_i H_i -> G_s` has image
`d Z/sZ`, so

```text
|R_s(w)| = d product_i D_i/s.                                (6a)
```

If the masks cover, (6) is nonzero and at least one nonzero tuple in
`R_s(w)` has `product_i fhat_i^tau(r_i) != 0`.  In particular

```text
d product_i D_i > s                                          (6b)
```

is necessary.  This is only a support obstruction; (6b) does not ensure that
any allowed coefficient is nonzero or that their sum has the required sign.
For one mask, (6a) equals `1`, so the obstruction already rules out a
one-tightener sheet cover.  For two masks it reads

```text
|R_s(w)|=s/lcm(g_1,g_2),
```

so a cover requires `lcm(g_1,g_2)<s` before any metric information is used.

For comparison with an unnormalized discrete Fourier transform, put

```text
F_i(r) = sum_j f_i^tau(j)e_s(-rj) = s fhat_i^tau(r).
```

If `m=|I|`, (5) is equivalently

```text
|G_s minus union_i M_i(tau)|
   = s^(1-m) sum_(r in R_s(w)) product_i F_i(r_i).             (7)
```

The powers of `s` in (5) and (7) are therefore completely fixed by the
Fourier convention.

## 2. Periodicity and ramified Fourier support

For every `j in G_s`,

```text
w_i(tau+j+D_i)/s
  = w_i(tau+j)/s + w_i/g_i.
```

The last summand is an integer.  Distance to the nearest integer is unchanged,
so

```text
M_i(tau)+D_i = M_i(tau),       f_i^tau(j+D_i)=f_i^tau(j).      (8)
```

This is the ramification statement: `f_i^tau` is pulled back from
`Z/D_i Z` along reduction `G_s -> Z/D_i Z`.

Decompose `j=a+kD_i`, where `0<=a<D_i` and `0<=k<g_i`.  Using (8),

```text
fhat_i^tau(r)
 = (1/s) sum_(a=0)^(D_i-1) f_i^tau(a)e_s(-ra)
           sum_(k=0)^(g_i-1) exp(-2 pi i rk/g_i).             (9)
```

The last geometric sum is zero unless `g_i|r`, and equals `g_i` when
`g_i|r`.  Hence

```text
fhat_i^tau(r)=0 unless g_i|r.                                 (10)
```

No converse about nonvanishing is intended: a coefficient allowed by (10)
may still vanish because of the exact location of the mask.

The masks are proper.  Indeed,
`gcd(w_i/g_i,D_i)=1`, so as `j` varies the phases in (1) form a translate of
a `D_i`-grid, with every grid point repeated `g_i` times.  A
closed arc of length `2/13` contains at most

```text
floor(2D_i/13)+1
```

points of such a grid.  Therefore

```text
|M_i(tau)| <= g_i(floor(2D_i/13)+1) < g_iD_i=s,               (11)
```

where the strict inequality follows from
`floor(2D/13)+1 <= 2D/13+1 < D` for `D>=2`.  Consequently

```text
fhat_i^tau(0)=(1/s)sum_j f_i^tau(j)=1-|M_i(tau)|/s>0.         (12)
```

This proves both the support restriction and the positivity assertion used
in (6).

## 3. Fourier inversion, orthogonality, and proof of the identity

The character orthogonality relation is

```text
(1/s) sum_(j in G_s) e_s(rj) = 1 if r=0 in G_s, and 0 otherwise.  (13)
```

It gives inversion in the normalization (3):

```text
sum_r fhat(r)e_s(rj)
 = (1/s)sum_k f(k)sum_r e_s(r(j-k))
 = f(j).                                                       (14)
```

Because every `f_i^tau` is Boolean, the product

```text
product_i f_i^tau(j)
```

is exactly the indicator that sheet `j` is in none of the danger masks.
Apply (14) to each factor, multiply, and average over `j`:

```text
(1/s) sum_j product_i f_i^tau(j)
 = sum_((r_i) in product_i G_s)
     product_i fhat_i^tau(r_i)
     [(1/s)sum_j e_s(j sum_i r_i)]
 = sum_(sum_i r_i=0) product_i fhat_i^tau(r_i).               (15)
```

By (10), all summands outside the extra conditions `g_i|r_i` vanish.  The
last expression is therefore the right side of (5).  Its all-zero term is,
by (12),

```text
product_i fhat_i^tau(0)=product_i(1-|M_i(tau)|/s).            (16)
```

If the masks cover `G_s`, the left side of (5) is zero.  Separating (16) from
the remaining terms proves (6).  Complex summands occur in conjugate pairs,
so their total is the negative real number in (6).

### Independence of the representative of `tau`

Replacing `tau` by `tau+a`, `a in Z`, merely relabels the sheets:

```text
f_i^(tau+a)(j)=f_i^tau(j+a),
fhat_i^(tau+a)(r)=e_s(ar)fhat_i^tau(r).                       (17)
```

On a resonance tuple the product of the phase factors in (17) is
`e_s(a sum_i r_i)=1`.  Thus both sides of (5) are intrinsic to
`tau in R/Z`, even though a real representative labels the individual sheets.

## 4. The exact THM-769 deep-packet corollary

Now specialize to THM-769.  Let `A` be a primitive twelve-speed set with a
global maximum of height `1/13` at a reduced point

```text
t* = p/(13s),        s>1,
```

and use its on-sheet/off-sheet split

```text
E={v in A:s|v}=sU,        F=A minus E.
```

For

```text
G_U={tau in R/Z: phi_U(tau)>1/13},                            (18)
```

THM-769's exact containment says that, for every `tau in G_U`, all `s` lifts

```text
t_j=(tau+j)/s,       j in G_s,
```

must be blocked by the off-sheet speeds.  Directly, an on-sheet speed `su`
has

```text
||su t_j||=||u(tau+j)||=||u tau||>1/13,
```

while tightness says that some member of `A` has clearance at most `1/13` at
every `t_j`.  That member must lie in `F`.  Hence, with the masks (1),

```text
union_(w in F) M_w(tau)=G_s       for every tau in G_U.       (19)
```

Applying (6) pointwise gives the exact necessary condition

```text
sum_(r in R_s(F), r != 0) product_(w in F) fhat_w^tau(r_w)
 = - product_(w in F)(1-|M_w(tau)|/s) < 0                    (20)
```

for every `tau in G_U`.

Equation (20) is the promised bridge.  THM-769's capacity inequality replaces
each sheet mask by a cardinality upper bound, hence retains only
zero-frequency information and no position.  The exact metric statement
remembers which lifts are covered; in (20), that information is precisely the
phase and amplitude of the nonzero coefficients.  The arithmetic gcd data do
not disappear: they cut the possible cancellation down to the ramified
support `g_w|r_w`.

In particular, a proposed deep packet is impossible at a given `tau` if its
supported nonzero resonance sum cannot equal the strictly negative right side
of (20).  The theorem does not supply a uniform estimate proving that
impossibility.

## 5. Character energy of a finite cover

There is a complementary exact identity for cover multiplicity.  It is useful
at the capacity-equality edges but should not be confused with the product
identity (5).

Let `G` be a finite abelian group of order `c`, let `A_1,...,A_m` be subsets
covering `G`, and suppose

```text
sum_i |A_i| = c+Delta.
```

For a character `chi` of `G`, use the normalized transform

```text
hhat(chi)=(1/c)sum_(x in G)h(x)conj(chi(x)).                  (21)
```

Here `1` denotes the trivial character.

Let `m(x)=sum_i 1_(A_i)(x)` and `e(x)=m(x)-1`.  Covering makes `e` a
nonnegative integer-valued function, and `sum_x e(x)=Delta`.  Parseval gives

```text
sum_(chi != 1) |sum_i widehat(1_(A_i))(chi)|^2
 = [c sum_x e(x)^2-Delta^2]/c^2
 >= Delta(c-Delta)/c^2.                                     (22)
```

Indeed, at every nontrivial character the transform of the constant function
`1` vanishes, so the sum inside the absolute value is `ehat(chi)`.  Normalized
Parseval says

```text
sum_chi |ehat(chi)|^2=(1/c)sum_x e(x)^2,
```

and the trivial coefficient is `ehat(1)=Delta/c`.  This proves the equality
in (22).  Since `e(x)` is a nonnegative integer,
`sum_x e(x)^2>=sum_x e(x)=Delta`, proving the inequality.  The lower bound is
informative for `0<=Delta<=c`; the identity itself has no such restriction.

At exact capacity, `Delta=0`, a cover is a partition and every nontrivial
combined character coefficient vanishes.  Conversely, even without assuming
a cover, if `sum_i|A_i|=c` and all the nontrivial combined coefficients vanish,
then Fourier inversion makes `m` the constant function `1`; hence the masks
partition `G`.  With unnormalized transforms, (22) is the equivalent identity

```text
sum_(chi != 1) |sum_i sum_(x in A_i)conj(chi(x))|^2
 = c sum_x e(x)^2-Delta^2.                                  (23)
```

Thus sheet capacity is the trivial character, character energy measures the
multiplicity forced into the nontrivial characters, and (5) retains the still
higher joint product needed to decide whether any sheet is uncovered.

## 6. The two- and three-sheet equality edges

### `s=2`: opposite parity is one nonzero Fourier cancellation

In THM-769's two-tightener edge, both exceptional speeds are odd.  Whenever
an odd speed `w` is eligible, meaning `||w tau||<=2/13`, its danger mask is a
singleton.  If `N_w(tau)` is the unique nearest integer to `w tau`, then

```text
M_w(tau)={N_w(tau) mod 2}.                                   (24)
```

Indeed, `j` is dangerous exactly when `w tau+wj` is within `2/13` of an even
integer; because `w` is odd, this is equivalent to
`j=N_w(tau) (mod 2)`.  This is precisely THM-769's parity colour.

For a singleton danger mask `{a}` in `Z/2Z`, its safe complement has

```text
fhat(0)=1/2,              fhat(1)=-(1/2)(-1)^a.               (25)
```

Two such masks cover both sheets exactly when their colours are opposite.
Then the two zero modes contribute `1/4`, while the only nonzero resonance
contributes `-1/4`.  Equation (6) is therefore exactly the persistent
opposite-parity condition of THM-769, not merely an analogy with it.

### `s=3`: three-colour ownership is character cancellation

At THM-769's three-tightener equality edge, every exceptional speed is a unit
modulo `3`.  Eligibility `||w tau||<=3/13` again makes the danger mask a
singleton,

```text
M_w(tau)={a_w(tau)},
a_w(tau)=-w^(-1)N_w(tau) (mod 3),                            (26)
```

which is THM-769's sheet colour.  For a singleton `{a}` in `Z/3Z`,

```text
fhat(0)=2/3,
fhat(r)=-(1/3)e_3(-ra),       r=1,2.                         (27)
```

Three eligible singleton masks cover at total capacity `3` exactly when their
colours are all distinct.  Equivalently, the two nontrivial character sums of
their danger indicators vanish, by (22) with `Delta=0`.  Equivalently again,
the nonzero resonances in (5) cancel the complement zero mode `(2/3)^3=8/27`.
This recovers the persistent three-colour condition of THM-769.

The word *persistent* remains essential in both cases: these equalities must
hold for every `tau in G_U`, not only at the original binding point.

## 7. Common-scale flags and the continuum resonance dictionary

The support proof is not special to danger intervals.  If `D|c` and a
function `f:Z/cZ -> C` is pulled back through

```text
Z/cZ -> Z/DZ,
```

then `f(j+D)=f(j)`.  With `g=c/D`, its Fourier transform is supported on

```text
g Z/cZ,                                                        (28)
```

the annihilator of the kernel of that quotient.  Conversely, inversion shows
that support in (28) implies `D`-periodicity.

This is the harmonic form of the quotient-fibre flags used in THM-994,
THM-1072, and THM-1090.  Retaining an order-`D` pullback mask means retaining coefficients
on the corresponding annihilator subgroup; discarding transverse fibre
position discards the remaining supported phases.  Equation (28) explains
why gcd/prime-power ramification and Fourier sparsity are the same finite
object.

Equation (5) is also the exact finite-cyclic counterpart of THM-1075's lattice,
whose continuum identity is now an absolutely convergent raw series, as well
as a Fejer limit, by THM-1092.  Circle integration there imposes an integer equality
`sum_i n_i a_i=0`; averaging over `G_s` here imposes a congruence
`sum_i r_i=0 (mod s)`.  The finite identity has no convergence issue and does
not by itself bound either its own nonzero resonance sum or the continuum
higher-order terms.

## 8. Carrier and tournament audit

The proof object is most naturally either of the following equivalent
incidence carriers:

1. the bipartite sheet-mask hypergraph with sheet vertices `j in G_s` and
   hyperedges `M_w(tau)`; or
2. the weighted zero-sum frequency hypergraph whose vertices are supported
   Fourier modes, whose hyperedges are tuples in `R_s(w)`, and whose weights
   are the products of Fourier coefficients in (5).

The first preserves literal cover, closed endpoints, and mask multiplicity.
The second preserves the same predicate through (5), while making gcd
ramification visible.  Runners alone, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, and proof obligations were all considered as possible vertices.
Only sheet incidence or the full resonance hypergraph preserves the joint
uncovered-sheet predicate without adding data back in.  In particular, a
pairwise runner relation loses resonances supported on three or more masks.

There is one exact tournament sidecar at the `s=3` equality edge.  Its vertices
are the three eligible exceptional runners, and its pairwise observable is
the colour difference

```text
a_v(tau)-a_u(tau) in Z/3Z.
```

The sheet-origin gauge translates every colour together and leaves this
difference unchanged.  Direct `u -> v` when the difference is `+1`; equal
colours are ties, oriented by one fixed ordering of runner labels, whose
restriction supplies the tie Hamiltonian path.  Three distinct colours give
the directed triangle: score histogram `{1,1,1}`, one strongly connected
component, one directed 3-cycle, and three directed Hamiltonian paths.  A
colour collision gives a transitive tournament: score histogram `{0,1,2}`, no
directed cycle, three singleton strongly connected components, and one
Hamiltonian path.  Thus on this particular edge the tournament cycle is
equivalent to the three-colour cover.

At `s=2`, the unlabelled fingerprints of every two-vertex tournament are the
same single arc.  One can encode equal versus opposite parity into its labelled
direction only by stipulating an extra switch rule; that merely renames the
binary colour predicate and adds no structural invariant.  For general `s`,
masks need not be singletons and there is no canonical binary orientation.
One may orient pairs by a chosen correlation key and break ties along a
Hamiltonian path, but its score sequence, cycles, and SCCs do not determine
(5).  The challenged assumption is therefore that runners should be
tournament vertices.  The exact general carrier is a ramified hypergraph; the
three-colour tournament is a sharp low-sheet quotient, not the global proof
object.

## 9. Honest scope

THM-1091 proves an identity and two exact low-sheet interpretations.  It does
not prove any of the following:

- that the nonzero resonance sum in (20) is too small to cancel the zero mode;
- that a putative cancellation cannot persist over the continuum `G_U`;
- the content inequality `val(A)<=gcd(A)` isolated in THM-1006;
- emptiness of the deep primitive twelve-speed branch; or
- the full LRC14 residual obligation.

The new reduction is nevertheless strict: capacity plus primitivity, which
THM-1006 shows cannot close the deep branch, is replaced by the exact target
(20).  Any successful continuation must exploit either phase variation in
`tau`, incompatibility among the ramified support subgroups, or a genuinely
higher-order bound on the supported zero-sum convolution.  Those are now the
only places where the metric sheet positions can enter.
