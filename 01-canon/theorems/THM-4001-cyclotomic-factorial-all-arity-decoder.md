---
id: THM-4001
title: "Coordinatewise cyclotomic factorial moments decode every fixed arity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over a
  characteristic-zero field containing the dth roots of unity, the
  coordinatewise (mu_d)^k Reynolds projector followed by the k-variable
  factorial functional sends the dmth power of a linear form to the complete
  homogeneous polynomial h_m(a_1^d,...,a_k^d).  For known k, the first k
  projected moments recover the unordered multiset of dth powers, including
  repetitions and zeros, by the elementary/complete-homogeneous Newton
  recurrence.  Recovering the coefficients themselves requires an injective
  dth-root branch.  A diagonal cyclic average is insufficient, the projector
  is not multiplicative, and the result has no FC, HFC, GMC, JC, or LRC
  consequence.
source: root + audit_factorial_all_arity_decoder / all-frontiers continuation, 2026-08-24
audit: >
  INDEPENDENT MATHEMATICAL AND HOSTILE AUDIT PASS.  The coefficient proof,
  coordinatewise selector, Newton inversion, repeated/zero boundary, root
  torsor, known-arity requirement, and nonmultiplicativity were independently
  reconstructed.  The exact companion checks 64 formal selector cases,
  24,704 exact response rows on 1,018 signed coefficient vectors, 4,072
  Newton/recurrence rows, and 3,136 nonnegative unordered finite controls.
  Diagonal and one-coordinate selector hostiles, the 1729 taxicab collision,
  and normal/optimized replay are explicit.  Normal and optimized LF streams
  equal the frozen output at 108,399 active gates.
depends_on: []
related:
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3825-prime-colour-valuation-two-cube-decoder
script: 04-computation/cyclotomic_factorial_all_arity_decoder_thm4001.py
output: 05-knowledge/results/cyclotomic_factorial_all_arity_decoder_thm4001.out
script_sha256: 45088e0169befc91ccb8097534e4108fe152718c7471283ea29ea6afd18425bb
output_sha256: 810ee433370f82988f92acdb185279c3cf49c40551889345d0117e99c5ad4a61
semantic_sha256: 6d53d78ad129ae4bb403ac01e1dfcd65a121782815f98d899e49bd661d350f57
hash_basis: raw LF bytes
---

# THM-4001 -- coordinatewise cyclotomic factorial moments decode fixed arity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is an
elementary response-decoder theorem.  Its selector and root-branch sidecars
are part of the statement.

## 1. Exact setting and projector

Fix integers

```text
d>=2,                         k>=2,
```

and let `K` be a characteristic-zero field containing a primitive `d`th root
of unity `zeta`.  On `K[x_1,...,x_k]`, define the factorial functional

```text
L(product_i x_i^e_i)=product_i e_i!                         (1)
```

and extend it linearly.  Define the **coordinatewise** Reynolds projector

```text
Pi_d(G)(x_1,...,x_k)
 =d^(-k) sum_(r_1,...,r_k in Z/dZ)
          G(zeta^r_1 x_1,...,zeta^r_k x_k).                 (2)
```

The adjective coordinatewise is load-bearing.  Character orthogonality gives

```text
Pi_d(product_i x_i^e_i)
 =product_i x_i^e_i     if d divides every e_i,
 =0                     otherwise.                          (3)
```

Because the inputs below have total degree divisible by `d`, averaging any
fixed `k-1` coordinates independently is equivalent to `(2)`: divisibility of
the last exponent then follows from the total degree.  A diagonal action or a
single-coordinate action for `k>2` is not equivalent.

Let

```text
F=sum_(i=1)^k a_i x_i,              a_i in K,
u_i=a_i^d,                                                (4)
```

and, for every `m>=0`, define

```text
M_m=L(Pi_d(F^(dm)))/(dm)!.                               (5)
```

Characteristic zero makes the denominator in `(5)` invertible.  The literal
formula is not asserted in finite characteristic, where factorials and the
Reynolds denominator may vanish.

## 2. Projected factorial response

For every `m>=0`,

```text
M_m
 =sum_(j_1+...+j_k=m) product_i a_i^(d j_i)
 =h_m(u_1,...,u_k).                                      (6)
```

Here `h_m` is the complete homogeneous symmetric polynomial and `h_0=1`.

### Proof

Expand `F^(dm)`.  By `(3)`, the only surviving exponent vectors have the form

```text
(e_1,...,e_k)=(d j_1,...,d j_k),       sum_i j_i=m.       (7)
```

For such a vector, the multinomial coefficient and the factorial readout
cancel exactly:

```text
binom(dm;d j_1,...,d j_k) product_i (d j_i)!=(dm)!.       (8)
```

After division by `(dm)!`, the channel has coefficient one and monomial
`product_i u_i^j_i`.  Summing `(8)` over the weak compositions in `(7)` proves
`(6)`.  This is a coefficient identity; no positivity or finite computation
is used.  QED

Equivalently, the whole response sequence has generating function

```text
H(t):=sum_(m>=0) M_m t^m
     =product_(i=1)^k (1-u_i t)^(-1).                    (9)
```

Common scaling is visible at every arity:

```text
(a_1,...,a_k)->(c a_1,...,c a_k)
       implies M_m->c^(dm) M_m.                          (10)
```

## 3. First k moments recover the dth-power multiset

Let `e_r=e_r(u_1,...,u_k)` be the elementary symmetric polynomials, with
`e_0=1`.  Equation `(9)` is exactly

```text
H(t) (sum_(r=0)^k (-1)^r e_r t^r)=1.                    (11)
```

Comparing coefficients recursively gives, for `1<=r<=k`,

```text
e_r=sum_(j=1)^r (-1)^(j-1) e_(r-j) M_j.                 (12)
```

Thus `M_1,...,M_k` recover `e_1,...,e_k`, and hence the monic polynomial

```text
P(T)=T^k-e_1 T^(k-1)+e_2 T^(k-2)-...+(-1)^k e_k
    =product_(i=1)^k (T-u_i).                            (13)
```

The root multiset of `(13)` is therefore the unordered multiset

```text
{a_1^d,...,a_k^d},                                      (14)
```

with multiplicity.  Repeated roots and zero roots require no genericity
hypothesis.  Factoring may be performed over an algebraic closure, although
the displayed roots already lie in `K`.

Once the `e_r` are known, the response continues by the order-`k` recurrence

```text
M_m=e_1 M_(m-1)-e_2 M_(m-2)+...+(-1)^(k-1)e_k M_(m-k),  (15)
```

for `m>=1`, with `M_0=1` and `M_n=0` for `n<0`.  Equations `(12)--(15)` are
the complete-homogeneous/elementary Newton inversion; no division by a
Vandermonde discriminant occurs.

### Repetitions, zeros, and known arity

For example, at `d=3`,

```text
(a_1,a_2,a_3,a_4)=(0,-1,2,2),
(u_1,u_2,u_3,u_4)=(0,-1,8,8),
(M_1,M_2,M_3,M_4)=(15,177,1871,18609),                  (16)
```

and `(12)` reconstructs

```text
P(T)=T^4-15T^3+48T^2+64T=T(T+1)(T-8)^2.                (17)
```

The arity `k` must be supplied.  Adding a zero root changes neither `(9)` nor
any positive moment:

```text
h_m(1,8)=h_m(0,1,8)             for every m>=0.         (18)
```

Therefore the response cannot infer an unknown number of zero coefficients.
No stronger unknown-arity claim is made.

### Root-branch boundary

The intrinsic output is `(14)`, not the original coefficient multiset.  For
arbitrary `xi_i in mu_d`,

```text
(a_1,...,a_k)->(xi_1 a_1,...,xi_k a_k)                 (19)
```

preserves every `M_m`.  Consequently:

- nonnegative real coefficients are recovered uniquely for every `d`;
- real coefficients are recovered uniquely when `d` is odd;
- positive integer coefficients, including the two-cube application, are
  recovered by their positive integer `d`th roots;
- for even `d`, real signs are lost; and
- over `C`, the independent `mu_d^k` root torsor is lost.

Positivity is unnecessary for recovery of `(14)`; it is only one lawful
sidecar for lifting `(14)` back to the `a_i`.

## 4. Selector and multiplicativity hostiles

The coordinatewise definition cannot be weakened cosmetically.  At

```text
k=3,       d=2,       F=x+y+z,       m=1,              (20)
```

the normalized readout is

```text
coordinatewise (mu_2)^3 projector:   3,
diagonal mu_2 projector:              6,
one-coordinate mu_2 projector:        4.                (21)
```

The diagonal action fixes all of the even-degree polynomial `F^2`; averaging
only `z` also leaves the mixed term `2xy`.  Only the coordinatewise action
retains precisely the exponent vectors divisible by two.

The projector is a linear idempotent, not an algebra homomorphism.  For
distinct `i,j`,

```text
[x_i^d x_j^d] Pi_d(F^d)^2
   =2 a_i^d a_j^d,
[x_i^d x_j^d] Pi_d(F^(2d))
   =binom(2d,d) a_i^d a_j^d.                            (22)
```

These coefficients differ for every `d>=2`.  Thus

```text
Pi_d((F^d)^2) != Pi_d(F^d)^2                            (23)
```

whenever the displayed cross coefficient is nonzero.

## 5. Taxicab control and all-arity interpretation

At `d=3,k=2`, the first projected moment is the ordinary two-cube value.  The
canonical collision is split by the second moment:

```text
(a,b)=(1,12):   (M_1,M_2)=(1729,2987713),
(a,b)=(9,10):   (M_1,M_2)=(1729,2260441).              (24)
```

More generally, for every known arity `k`, the first `k` projected moments
separate unordered positive `k`-tuples through their `d`th powers.  This is a
response decoder, not a statement about how often the first power sum
collides or how its support is distributed.

THM-3825 gives a different tradeoff at arity two: one arithmetic scalar
decodes a restricted inert/cube-free carrier.  THM-4001 uses `k` factorial
response scalars and imposes no inertness, coprimality, cube-free shell, or
factorization hypothesis.

## 6. Strict nonconsequences and nearby canon

The sequence can be written

```text
M_m=(L o Pi_d)((F^d)^m)/(dm)!.                          (25)
```

It therefore uses the modified functional `L o Pi_d`.  By `(23)`, it cannot
be replaced by original factorial moments of the fixed invariant polynomial
`Pi_d(F^d)`.  The nonzero normalization does not affect vanishing, but the
change of functional does.  In particular, when all `a_i` are positive every
term in `(6)` is positive, so this packet supplies no null sequence.

The relation to the nearest factorial results is sharply typed:

- THM-2810 excludes a fixed finite multiplicative carrier for the entire
  characteristic-zero factorial tower.  The per-input order-`k` response
  recurrence `(15)` belongs to one projected family and is not such a carrier.
- THM-3018 supplies the factorial functional and its homogeneous simplex
  interpretation.  It applies the original functional and does not contain
  the coordinatewise projector or decoder `(12)`.
- THM-3300 identifies the torus-invariant Gaussian/factorial class and finite
  cyclic selection boundary.  It does not state this all-arity response or
  recover a root multiset.
- THM-3825 is the restricted one-scalar arithmetic decoder described above;
  it neither subsumes nor is subsumed by the unrestricted `k`-scalar tradeoff.

Accordingly, THM-4001 proves or refutes no case of `FC`, `HFC`, `SFC`, `GMC`,
`NC2`, `JC(2)`, `DC(2)`, or `LRC(14)`.  It does not recover coordinate
placement, orientation, global sign, other runner speeds, phase, owner,
arrival, or loneliness.

## 7. Exact verification

Run

```text
python -B 04-computation/cyclotomic_factorial_all_arity_decoder_thm4001.py
python -B -O 04-computation/cyclotomic_factorial_all_arity_decoder_thm4001.py
```

and compare LF-normalized output with

```text
05-knowledge/results/cyclotomic_factorial_all_arity_decoder_thm4001.out.
```

The assertion-free companion checks:

- `64` formal coefficient dictionaries for the coordinatewise selector,
  containing `484` retained and `11,053` rejected monomials;
- `24,704` exact response rows on `1,018` signed coefficient vectors;
- `4,072` Newton inversions and order-`k` recurrences;
- `3,136` nonnegative unordered tuples, including repeats and zeros;
- the diagonal, one-coordinate, root-branch, unknown-arity,
  nonmultiplicativity, and 1729 hostiles.

Normal and optimized runs equal the frozen LF stream at `108,399` active
gates.  The universal quantifiers are proved by `(7)--(13)`, not inferred from
the finite controls.

**QED.**
