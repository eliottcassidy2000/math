---
id: THM-4285
title: "Universal two-Eisenstein-norm shells and the exact contact trichotomy"
status: >
  PROVED RELATIVE TO THM-4280 + FINITE-EXACT AUDIT PASS. The A2 direct-sum
  theta series has coefficient 12(sigma(n)-3 sigma(n/3))=12 sigma(m)>0
  for n=3^a m with 3 not dividing m. Consequently every positive multiple-of-
  four degree shell has a V-slice, and THM-4280's minimal uniform noncollapse
  contact is exactly 5Q on Eisenstein-norm quotients and 3Q on nonnorm
  quotients; its existing 2Q branch covers nonmultiples of four. This is a
  shell zero-test, not pairwise recovery, raw A23 descent, or JC(2).
source: codex-continuation-frontiers-20260829
depends_on:
  - THM-4280-integral-three-channel-fat-contact-observer-and-sharp-five-jet-bound
related:
  - THM-4279-four-channel-formal-log-hasse-observer-for-e0-hom-at-fat-contact
  - THM-4284-a23-conductor-defect-and-degree-shell-first-character-nondescent
primary_script: 04-computation/jc23_two_eisenstein_norm_universal_shell_thm4285.py
primary_output: 05-knowledge/results/jc23_two_eisenstein_norm_universal_shell_thm4285.out
primary_script_sha256: 95132bdd059d8199f5b57323c4ae7b67368a6c8142e66f583c87da5ac689fe22
primary_output_sha256: d45dfe9b247be56d0d80acd4ac7f27a7855c2b5dbc1e7c9cf21bb1d1e0729fc2
hash_basis: raw LF bytes
audit: >
  PASS. Direct complete A2 enumeration through n=4096, followed by exact
  convolution, agrees entrywise with an independent divisor sieve. The audit
  freezes the n=0 exception, y-to-minus-y convention, powers of three,
  inert-prime nonnorms, norm-shell strata, and stable vector digests. The
  finite computation controls but does not replace the all-n modular proof.
---

# THM-4285 -- universal two-Eisenstein-norm shells and the exact contact trichotomy

**PROVED RELATIVE TO THM-4280 + FINITE-EXACT AUDIT PASS. THE PLANAR
JACOBIAN CONJECTURE REMAINS OPEN.**

## 1. Statement

Let

```text
O=Z[omega],                 omega^2+omega+1=0,
N(x+y omega)=x^2-xy+y^2,                              (1)
```

and count ordered, signed coordinate quadruples by

```text
r(n)=#{(x,y,z,w) in Z^4:N(x+y omega)+N(z+w omega)=n}.  (2)
```

Then `r(0)=1`, while for every `n>=1`,

```text
r(n)=12(sigma(n)-3 1_(3|n) sigma(n/3)).                (3)
```

If `n=3^a m` with `3` not dividing `m`, the parenthesis in `(3)` is exactly
`sigma(m)`. Hence

```text
r(n)=12 sigma(m)>0.                                   (4)
```

In particular, every nonnegative integer is a sum of two Eisenstein norms.
The value at zero is stated separately: substituting `n=0` into a divisor-sum
formula is not meaningful.

Now retain THM-4280's exact characteristic-zero `W=Lambda=0` fibre, based
actual global Hom classes, marked point `Q`, target-value sidecar, and
ramification convention. For a positive degree `d` whose actual shell is
nonempty, let `J(d)` denote the least truncation length such that every
nonconstant class in that shell stays nonconstant on `J(d)Q`. Then

```text
J(d)=2,                          if 4 does not divide d;
J(4n)=5,                         if n is an Eisenstein norm;
J(4n)=3,                         if n is not an Eisenstein norm.    (5)
```

This is the exact **minimal uniform noncollapse contact**. It is not pairwise
faithfulness or recovery of a whole degree shell.

There is also an exact multiplicity sidecar on THM-4280's consumer
`V=O u direct-sum O v`. Put

```text
rho(n)=#{e in O:N(e)=n}.                                  (6)
```

The degree-`4n` part of `V` has `r(n)` classes. For `n>=1`, exactly `rho(n)`
have local ramification four, and the remaining `r(n)-rho(n)` have local
ramification two. Thus on a nonnorm quotient every one of the `r(n)` classes
in this `V`-slice has exact ramification two.

## 2. The A2 direct-sum theta identity

Put `q=exp(2 pi i tau)` and

```text
Theta(tau)=sum_(x,y in Z) q^(x^2-xy+y^2).                (7)
```

The Gram matrix of the quadratic form in `(7)` is

```text
B=[[2,-1],[-1,2]],       det(B)=3,
3 B^(-1)=[[2,1],[1,2]].                                 (8)
```

Thus `B` is an even positive-definite rank-two lattice of level three. The
even-lattice theta transformation, obtained by Poisson summation on the
lattice and its dual, gives

```text
Theta in M_1(Gamma_0(3),chi_(-3)).                       (9)
```

Squaring removes the quadratic character, so

```text
Theta^2 in M_2(Gamma_0(3)).                              (10)
```

For completeness, define the quasimodular Eisenstein series

```text
E_2(tau)=1-24 sum_(n>=1) sigma(n)q^n.                    (11)
```

The anomalous term in the transformation law of `E_2` cancels in

```text
G(tau)=(3E_2(3tau)-E_2(tau))/2,                          (12)
```

so `G` is also in `M_2(Gamma_0(3))`.

Here `[SL_2(Z):Gamma_0(3)]=4`; the modular curve has two cusps, no elliptic
point of order two, and one elliptic point of order three. Its genus is

```text
1+4/12-0/4-1/3-2/2=0.                                   (13)
```

Consequently `S_2(Gamma_0(3))=0`, while the weight-two Eisenstein subspace
has dimension `2-1=1`. Therefore

```text
dim M_2(Gamma_0(3))=1.                                   (14)
```

Both `Theta^2` and `G` have constant coefficient one, so `(10)--(14)` force

```text
Theta(tau)^2=(3E_2(3tau)-E_2(tau))/2.                    (15)
```

The coefficient of `q^n` on the left is exactly `(2)`. Expanding the right
side using `(11)` proves `(3)` for every `n>=1` and also gives the separate
constant coefficient `r(0)=1`.

To see positivity without any asymptotic argument, write `n=3^a m`, with
`3` not dividing `m`. If `a=0`, the correction term in `(3)` is absent. If
`a>=1`, multiplicativity of `sigma` gives

```text
sigma(n)-3sigma(n/3)
 =sigma(m)(sigma(3^a)-3sigma(3^(a-1)))
 =sigma(m).                                             (16)
```

Equations `(3)` and `(16)` prove `(4)` and universality.

## 3. Exact degree-shell contact

THM-4280 proves the actual local ramification spectrum `{1,2,4}` and the
degree-shell criterion

```text
e_Q(h)=4 iff h=e v with e in O-{0},
deg(e v)=4N(e).                                         (17)
```

It also proves that every map on a nonempty shell with `4` not dividing its
degree has ramification one. Such a map is invisible on `1Q`, which retains
only the target value, and visible on `2Q`; this gives the first line of `(5)`.

Suppose next that `d=4n`. If `n=N(e)`, the class `ev` in `(17)` has exact
ramification four. THM-4280's uniform `5Q` bound is therefore both sufficient
and necessary, proving the second line of `(5)`.

If `n` is not an Eisenstein norm, `(17)` excludes ramification four, so the
spectrum makes `3Q` sufficient. By `(4)`, choose

```text
n=N(a)+N(e).                                             (18)
```

Neither `a` nor `e` can vanish: otherwise `n` itself would be a norm. The
actual global class

```text
h=a u+e v                                                (19)
```

has degree `4n`, and THM-4280's integral channel matrix gives

```text
c_1(h)=0,                    c_2(h)=a!=0.                (20)
```

Thus `(19)` has exact ramification two. It is invisible on `2Q`, so `3Q` is
necessary as well as sufficient. This proves the final line of `(5)` and
upgrades THM-4280's earlier infinite family and four examples to every
nonnorm quotient.

Finally, the `O`-basis `(u,v)` makes `(a,e)` in `(19)` unique. Equations
`(17)` and `(20)` say that its ramification is four exactly when `a=0` and
two exactly when `a!=0`. Counting the two strata gives `(6)` and the
multiplicity assertion.

## 4. Exact audit

The optimization-safe script

```text
04-computation/jc23_two_eisenstein_norm_universal_shell_thm4285.py
```

checks the first `4097` coefficients (`0<=n<=4096`) by two independent exact
paths:

1. Since
   `2N(x+y omega)=(x-y)^2+x^2+y^2>=x^2+y^2`, it completely enumerates the
   A2 norm histogram in the square `|x|,|y|<=91` and convolves that histogram
   with itself.
2. Independently, it builds `sigma(n)` by a divisor sieve and evaluates the
   right side of `(3)`.

The arrays agree entrywise. Controls also freeze `r(0)=1`, the equivalent
`x^2+xy+y^2` sign convention, the maximally cancelling powers
`1,3,...,2187`, the nonnorm rows `2,5,8,11`, and the order-two/order-four
stratum digests. The stored transcript is

```text
05-knowledge/results/jc23_two_eisenstein_norm_universal_shell_thm4285.out
```

Reproduce it from the repository root with

```bash
python3 -B 04-computation/jc23_two_eisenstein_norm_universal_shell_thm4285.py
python3 -B -O 04-computation/jc23_two_eisenstein_norm_universal_shell_thm4285.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_two_eisenstein_norm_universal_shell_thm4285.py
```

The finite audit is a hostile control of `(15)`; it is not the all-integer
proof, which is the modular identity in Section 2.

## 5. Scope and connection ledger

The exact connection is

```text
source:              the A2 direct-sum lattice O^2 and its theta series;
target:              THM-4280's integral V-consumer degree shells;
map:                 (a,e) maps to the actual Hom class a u+e v;
preserved predicate: degree 4(N(a)+N(e)) and first formal-log exponent;
destroyed data:      H-components, target value, and pairwise shell identity;
restoring sidecars:  the marked target value, actual-global-Hom membership,
                     the fixed W=Lambda=0 fibre, Q, b, and omega_E;
sharp hostiles:      r(0), pure ev, powers of three, and inert nonnorms.
                                                               (21)
```

No claim here supplies a transverse `W`-jet, saturated raw A23 descent,
pairwise recovery on a degree shell, new LRC incidence deletion, physical
exact-`M=12` entry, `JC(2)`, or `DC(2)`.
