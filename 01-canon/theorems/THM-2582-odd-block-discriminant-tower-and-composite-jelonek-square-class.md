---
id: THM-2582
title: "Odd-block discriminant tower and composite Jelonek square class"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  For a separable finite index algebra, the discriminant square class of a
  norm-product of degree-m conjugate polynomial blocks is the norm of the
  block discriminant times the index discriminant to the power m.  The
  mechanism is the parity of the cross-block resultant: symmetric for even
  m and alternating for odd m.  For the fixed sporadic cubic Keller map, a
  corrected inverse-chart root product gives Norm(L)=H/(64L).  Hence the
  degree-nine x-eliminant of F o F has global discriminant square class H;
  its sole odd irreducible divisor is H, while L has even valuation.  This is
  a field/divisor-parity theorem, not an exact multiplicity, higher-iterate,
  general Keller, JC(2), or GMC(2) result.
source: root-frontier-final-2026-07-28
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
related:
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
  - MISTAKE-290
script: 04-computation/keller_odd_block_discriminant_square_class_thm2582.py
output: 05-knowledge/results/keller_odd_block_discriminant_square_class_thm2582.out
script_sha256: afcf930e050e2cae525e07fef6127c8e57b522898d0365f78832110166260a9c
output_sha256: d9a2b77f76c02d8072793decf12a96e8a99bd6277ba25f7f5333a6d7d91c1ff4
hash_basis: working-tree bytes (LF)
---

# THM-2582 -- odd blocks carry the outer discriminant class

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

The three exact slices in HYP-9033 showed that the old Jelonek factor `L`
occurs to even order in the degree-nine composite discriminant, while the new
factor `H` occurs to odd order.  This theorem proves that pattern globally.
The missing mechanism is a standard but easily hidden sign representation:
the product of all cross-block resultants is alternating when the blocks have
odd degree.  Its square contributes one copy of the discriminant of the
indexing polynomial.

The fixed-map calculation is then only the square-class multiplication

```text
outer cubic:                 [-L],

norm of three inner cubics:  [-H/L],

composite:                   [-L][-H/L]=[H].                 (1)
```

Thus the old factor cancels in parity rather than disappearing from the full
discriminant.

## 1. Discriminant of a norm-product

Let `K` be a field of characteristic different from two.  Let `g(Y)` be a
separable degree-`n` polynomial, with roots `y_1,...,y_n` in a splitting
field.  For each root let

```text
f_i(X)=f(y_i;X)                                             (2)
```

be a degree-`m` polynomial.  Assume that the `f_i` are separable and pairwise
coprime, and put

```text
P(X)=product_(i=1)^n f_i(X).                                (3)
```

The coefficients of `P` lie in `K`; equivalently, it is the polynomial norm
of the conjugate block `f` from the finite separable index algebra.  Then in
`K^*/K^{*2}`,

```text
[Disc_X(P)]
 =[Norm(Disc_X(f))] [Disc_Y(g)]^m.                          (4)
```

The statement is unchanged if `P` is multiplied by any nonzero element of
`K`, because its degree is `mn` and

```text
Disc(cP)=c^(2mn-2) Disc(P).                                 (5)
```

### Proof

The product-discriminant identity gives

```text
Disc(P)
 =product_i Disc(f_i)
  *product_(i<j) Res(f_i,f_j)^2.                            (6)
```

The first product is `Norm(Disc(f))`.  Write

```text
R=product_(i<j) Res(f_i,f_j).                               (7)
```

Interchanging two index roots interchanges two degree-`m` blocks.  Since

```text
Res(f_j,f_i)=(-1)^(m^2) Res(f_i,f_j)=(-1)^m Res(f_i,f_j),   (8)
```

`R` is symmetric when `m` is even and alternating when `m` is odd.

If `m` is even, `R` is a symmetric rational function of the `y_i`, hence is
in `K`; its square contributes no square class.  If `m` is odd, divide by the
Vandermonde

```text
V=product_(i<j)(y_i-y_j).                                  (9)
```

The quotient `R/V` is symmetric and belongs to `K`.  If `a_g` is the leading
coefficient of `g`, then

```text
Disc(g)=a_g^(2n-2)V^2,                                     (10)
```

and the power of `a_g` is itself a square.  Therefore `[R^2]=[Disc(g)]` for
odd `m`, and `[R^2]=1` for even `m`.  These are exactly the two parities of
`[Disc(g)]^m`, proving (4). QED.

This proof needs the pairwise-coprime hypothesis only to keep the displayed
classes nonzero.  A single squarefree specialization proves generic
coprimality in the fixed application below.

## 2. The corrected inverse-chart norm of the Jelonek factor

Work over

```text
K=Q(a,b,c).                                                (11)
```

Use THM-2576's notation

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
T=4-3bc,
S=27ac^2-9bc+8,

E(X)=LX^3+TX-2c.                                           (12)
```

For a root `w` of `E`, the rational inverse section produces the middle
fibre point

```text
q(w)=(w,Y(w),Z(w)).                                        (13)
```

THM-2576 defines polynomials `D,N,P,Q` satisfying

```text
s=wY(w)=-N(w)/D(w),

w^6 L(q(w))=P(w,s),

Q(w)=D(w)^4 P(w,-N(w)/D(w)).                               (14)
```

MISTAKE-290 repairs only the convention sign in the large resultant.  With
the standard Sylvester/root-product convention, the exact identity is

```text
Res_X(E,Q)=a^8 c^18 S^8 H.                                 (15)
```

The small companion resultant is

```text
Res_X(E,D)=a^2 c^3 S^2.                                    (16)
```

Let `w_1,w_2,w_3` be the roots of `E`.  Since `deg Q=15` and `deg D=2`,
the root-product convention gives

```text
product_i Q(w_i)=a^8 c^18 S^8 H/L^15,

product_i D(w_i)=a^2 c^3 S^2/L^2,

product_i w_i=2c/L.                                        (17)
```

Using (14),

```text
product_i L(q(w_i))
 = [product_i Q(w_i)]
   /[(product_i D(w_i))^4 (product_i w_i)^6]

 = [a^8c^18S^8H/L^15]
   /[64a^8c^18S^8/L^14]

 = H/(64L).                                                (18)
```

Thus the function-field norm law previously isolated by the independent
inverse-section calculation follows directly from THM-2576's other chart
after the resultant convention is corrected.  Its sign is positive.

At `(a,b,c)=(1,1,1)`, the exact control is

```text
L=25,                    H=951326441195,

H/(64L)=190265288239/320.                                  (19)
```

The same specialization gives a degree-nine squarefree composite `x`-core,
with leading coefficient `-H`; hence the conjugate cubic blocks are pairwise
coprime generically and the nonzero hypothesis of Section 1 holds.

## 3. Global square class of the degree-nine eliminant

For each middle point `q_i=q(w_i)`, let

```text
E_i(xi)=L(q_i)xi^3+T(q_i)xi-2q_(i,3).                      (20)
```

THM-2473 gives

```text
Disc_xi(E_i)=-4 S(q_i)^2 L(q_i).                           (21)
```

The polynomial norm

```text
E_2(xi)=product_(i=1)^3 E_i(xi)                            (22)
```

is the generic degree-nine `x`-eliminant of `F o F`, up to multiplication
by a nonzero element of `K`.  Apply (4) with index degree `n=3` and odd block
degree `m=3`.

First, (12) gives

```text
[Disc_X(E)]=[-4S^2L]=[-L].                                 (23)
```

Second, (18) and (21) give

```text
[Norm(Disc(E_i))]
 =[(-4)^3 Norm(S(q_i))^2 H/(64L)]
 =[-H/L].                                                  (24)
```

The square-class law now gives

```text
[Disc_xi(E_2)]
 =[-H/L][-L]
 =[H]                in K^*/K^{*2}.                        (25)
```

This is the promised global identity.  Since
`Q[a,b,c]` is a unique-factorization domain and THM-2576 proves that `H` is
irreducible, (25) says that in any denominator-cleared primitive polynomial
representative of the degree-nine discriminant:

```text
valuation_H(Disc(E_2)) is odd,

valuation_R(Disc(E_2)) is even for every irreducible R!=H. (26)
```

In particular `valuation_L` is even because `gcd(L,H)=1`.  Equivalently, the
odd square-free divisor is `H` alone.  This proves globally what the three
HYP-9033 slices had shown finitely.

## 4. Mechanism and higher-level boundary

The cancellation in (25) has two distinct sources:

```text
Norm(L)=H/(64L)        contributes the inverse old component;

odd block alternation contributes one outer discriminant L. (27)
```

They cancel exactly in parity.  Omitting the alternating cross-block factor
would predict `LH`; omitting the denominator in the norm would also predict
the wrong tower.  The conductor/image denominator and the sign
representation are both load-bearing.

Equation (4) is all-degree, but (25) is only a level-two fixed-map
application.  At level three one would still need:

1. the next image divisor and its norm law;
2. generic separability/coprimality of the next block family; and
3. control of whether successive image divisors are distinct.

Therefore this theorem does not prove the HYP-9033 prediction that every
level has only the newest odd factor, nor its component-count prediction.
It gives the exact induction mechanism and closes the first nontrivial rung.

## 5. Exact companion and scope

Run

```bash
python3 04-computation/keller_odd_block_discriminant_square_class_thm2582.py
python3 -O 04-computation/keller_odd_block_discriminant_square_class_thm2582.py
```

The companion checks:

- odd cubic and even quadratic block controls for (4), with nontrivial square
  quotient `(t+2)^2`;
- the block-swap parity `(-1)^m`;
- the exact small resultant (16), the exponent ledger behind (18), and the
  transported coefficient hash of `H`;
- MISTAKE-290's exact PRS/Sylvester sign hostile at `(1,1,1)`;
- a squarefree degree-nine composite core at that target; and
- the final square-class cancellation (25).

The result concerns the discriminant divisor parity of one coordinate
eliminant for one fixed three-dimensional Keller map.  It gives no exact
positive multiplicities, classification of the Keller monoid, higher tower
component count, planar Jacobian conclusion, or GMC(2) consequence.

**QED.**
