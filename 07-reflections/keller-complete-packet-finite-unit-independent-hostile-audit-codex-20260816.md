# The finite branch is the beta-homogeneous pullback prime

**Status: INDEPENDENT HOSTILE AUDIT ACCEPTS THE PROPOSED THM-3529
FINITE-SHEET UNIT MECHANISM AT FIXED-MAP COMPLETE-PACKET SCOPE.  THE
THM-3529 FILE REMAINS RESERVED PENDING OWNER CANONIZATION.**

For the fixed sporadic Keller map of
THM-2473 -- `sporadic-keller-branch-tower-depressed-trisection-anatomy`, put

```text
B=F^*L=L(F_1,F_2,F_3).
```

The proposed conclusion is correct:

```text
s_L(P)=0
```

for every nonzero polynomial `P` carrying a complete packet `A(e,m)`.  By
THM-3528 --
`fixed-keller-all-level-cleared-norm-polynomiality-and-finite-sheet-defect`,
every positive-level raw cleared-norm rung is consequently coprime to the old
`L` (the seed itself is `P_0=L`).  This
does not prove that any later rung is irreducible or a new image prime.

The audit used the following inheritance board.

1. THM-3528 supplies the proved divisor identity
   `s_L(P)=ord_(C_fin)(P)`.
2. THM-2473 supplies the unique generic finite affine inverse point and the
   everywhere-etale fixed map.
3. THM-3506 -- `fixed-keller-five-face-norm-transform-and-271-99-boundary`
   supplies the complete minimum-`beta` face.
4. MISTAKE-415 is the corrected near miss: a visible face alone is not a
   closed norm state.  Here the complete face and the finite divisor are
   both retained.
5. The underused sidecar is the literal source pullback divisor `F^*L`, not
   another evaluation at a chosen finite point.

## 1. Exact pullback and beta grading

Independent expansion gives

```text
B = -9x^4y^2z^2 -54x^3y^3z -18x^3yz^2 -81x^2y^4
    -72x^2y^2z -9x^2z^2 -54xy^3 +6xyz +63y^2 +16z.
```

It has ten terms, content one, `deg_x(B)=4`, and

```text
B(0,y,z)=63y^2+16z != 0.                              (1)
```

For

```text
beta(x^i y^j z^k)=i-j-2k,
```

all ten monomials have weight `-2`.  Thus `B` is literally
`beta`-homogeneous; there is no sign change between the candidate and
THM-3506 conventions.

## 2. The Laurent quadratic proves that B is prime and reduced

Localize at `x` and introduce

```text
p=xy,              q=x^2z.
```

This is an isomorphism

```text
Q[x^{+-1},y,z] = Q[x^{+-1},p,q].
```

Direct substitution gives `x^2B=C(p,q)`, where, as a quadratic in `q`,

```text
C = -9(p+1)^2 q^2
    -2(27p^3+36p^2-3p-8)q
    -9p^2(9p^2+6p-7).                                 (2)
```

The three coefficients in (2) have gcd one in `Q[p]`, and

```text
disc_q(C)=64(3p+4).                                    (3)
```

The valuation of (3) at the linear prime `3p+4` is one, so it is not a
square in `Q(p)`.  The primitive quadratic is therefore irreducible in
`Q[p,q]` by Gauss's lemma.  It remains prime after adjoining and inverting
`x`, so `B` is irreducible in `Q[x^{+-1},y,z]` up to a Laurent unit.

If `B=UV` factored in `Q[x,y,z]`, one factor would become a Laurent unit,
hence would be `c x^r`.  Equation (1) rules out `r>0`; polynomiality rules
out `r<0`.  The factor is therefore constant.  Thus

```text
(B) is a height-one prime and V(B) is irreducible and reduced.            (4)
```

This direct argument is independent of the exact factorization routine,
which returns the same one-factor answer.

## 3. No hidden affine divisor component exists

There is also a geometric cross-check.  The target polynomial `L`, viewed
as a quadratic in `a`, has discriminant

```text
disc_a(L)=4(4-3bc)^3,
```

which is nonsquare in `Q(b,c)`.  Hence `D=V(L)` is integral.  Since
`det J_F=-2`, `F` is etale and therefore quasi-finite everywhere.  Base
change along `D` identifies

```text
F^{-1}(D)=V(F^*L)=V(B),                                (5)
```

and makes `V(B)->D` etale.  Thus the pullback is reduced.  Every divisor
component has dimension two; quasi-finiteness forces its image closure to
have dimension two, so it must dominate the two-dimensional integral
divisor `D`.  No divisor can hide over the exceptional empty-fibre curve or
another lower-dimensional stratum.

At the generic point of `D`, THM-2473's core cubic becomes

```text
(4-3bc)w-2c,
```

with the single rational root `w=2c/(4-3bc)`.  The inverse formulas of
THM-3495 --
`level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness`
recover `y,z` rationally there.  Hence the generic finite fibre is one
reduced residue-degree-one point.  The geometric component argument and
the Laurent irreducibility proof agree.

## 4. The finite branch closure is exactly V(B)

THM-3528 already proves that the closure `C_fin` of this generic finite
branch is a source divisor and that

```text
s_L(P)=ord_(C_fin)(P).                                 (6)
```

Because `L(F(q_fin))=0`, one has `C_fin subseteq V(B)`.  Both are
irreducible closed subsets of dimension two, while (4) says `V(B)` has only
one component.  Therefore

```text
C_fin=V(B),              s_L(P)=ord_B(P).              (7)
```

So yes: once the proved THM-3528 divisor statement is retained, literal
containment plus irreducibility of `B` immediately gives the desired
identification.  Global quasi-finiteness supplies an independent exclusion
of vertical divisor components but is no longer load-bearing after (4).

## 5. The complete minimum-beta face excludes B

The possibly negative weight vector `(1,-1,-2)` still gives a genuine
`Z`-grading of `Q[x,y,z]`.  Every individual polynomial has a least occupied
weight.  Since the associated graded ring is a domain, nonzero initial forms
multiply without cancellation.

Suppose `B|P` and write `P=BQ`.  Homogeneity of `B` gives

```text
in_min-beta(P)=B in_min-beta(Q).                       (8)
```

For a complete packet `A(e,m)`, THM-3506 prescribes the nonzero face

```text
in_min-beta(P) = c y^(3e-2m) z^(e-m)
  (y^2+27z)^(2m/3)(y^2+108z)^(m/3),      c!=0.         (9)
```

It lies in `Q[y,z]` and has `x`-degree zero.  But over the coefficient
domain `Q[y,z]`, equation (8) has

```text
deg_x(B in_min-beta(Q))=4+deg_x(in_min-beta(Q))>=4,    (10)
```

contradicting (9).  Hence `B` does not divide `P`.  Equations (7)--(10)
prove

```text
s_L(P)=0                                                (11)
```

for every complete packet.

Neither negative weights nor cancellation create a loophole: only finitely
many grades occur in each polynomial, and the product of the two lowest
nonzero homogeneous pieces is nonzero in a domain.

## 6. Consequence and exact boundary

THM-3528 proves for `P` with packet `A(e,m)` that

```text
ord_L(L^eN(P))=s_L(P).
```

Together with (11), this gives `gcd(L,L^eN(P))=1`.  Since THM-3528 also
proves that every raw rung is a complete packet, every raw rung after the
seed is old-`L`-coprime.

This proves only:

- the finite branch is a unit for every complete packet;
- all raw cleared norms in this fixed tower are coprime to the old `L`.

It does not prove later-rung irreducibility, squarefreeness, primitive
integral normalization, a new image equation or image prime, separability,
full eliminant degree, distinct Jelonek components, discriminant recursion,
an arbitrary Keller-map statement, or any general Jacobian-conjecture claim.

## 7. Certificate repair found by the audit

The incoming exact companion correctly proves every symbolic identity but
its first bounded census used only `0<=2m<=e`, the stronger raw-orbit cone.
Complete packets allow the full cone

```text
0<=m<=e,             3|m.
```

For example `A(3,3)` has nonnegative exponents and was omitted.  This was not
a theorem counterexample because the symbolic argument above is uniform in
`e,m`.  The companion was repaired to test all `15,250` positive-`e`
packets through `e=300`, with `A(3,3)` as an explicit boundary hostile.
Normal, optimized, and stored-output replays agree after the repair.
