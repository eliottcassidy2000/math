# Independent audit of the valuation-twelve response classification

**Status: ACCEPTED: analytic scope and independent FINITE-EXACT reconstruction.**
The accepted theorem is the complete finite response with source valuation
at least twelve, low source weight at most 23 plus the designated `p³y⁶`
column, perturbations in rows eleven through fifteen, and the five literal
background rows in the
[report](planar_jc_long_20260906_memory_valuation12.md).
It includes the sharp weight-22 threshold, affine five-dimensional
replacement family, affine ten-dimensional terminal fibres, and the
consumer over the actual inherited `G_m`. It does not classify earlier
source valuations, later rows, full polynomial depth, or JC(2).

## 1. Independent mathematical route

The primary engine uses 155 free coefficients for the two perturbation
polynomials and eliminates the literal source, bracket, and depth equations.
This audit instead solves every bracket row first. The row operator is

```text
n[(x/2)delta C+(3/8)(x²+2)delta A].
```

Its complete kernel is `theta(A0',C0')`, with `deg theta<=n`: `x` and
`x²+2` are coprime, and the degree caps impose precisely that bound.
A different particular solution from the inherited pivot routine is
obtained, for target `nF`, by
`delta A=4F(0)/3` and
`delta C=2[F-(3/8)(x²+2)delta A]/x`. The numerator vanishes at zero and
the degree caps hold. Thus the audit's 70 tangent parameters represent
*every* admissible capped bracket response, with no chosen terminal
subfamily hidden in the parameterization.

It then computes the actual source coefficients in rows twelve through
fifteen, giving 58 equations. For depth it enumerates every literal
generator `x^a u^b p^c y^e` with `a+b<=d` and first row at most 15,
separates the genuine diagonal summands, and computes their left-nullspaces.
This independently supplies all 91 depth constraints. It does not use
the primary's binomial annihilator formula or import any repository
mathematical implementation.

Exact elimination of this 149-equation system gives tangent rank 60,
source rank five, joint rank 65, and kernel dimension ten. The five
source rows agree entry for entry with the primary's result; agreement
of just the ranks would not have been sufficient. All entries and pivots
are rational constants, so the classification holds over every
characteristic-zero field without an exceptional boundary locus.

## 2. Complete source universe, family, and minimality

The monomial enumeration independently recovers exactly the fourteen low
columns in the report and the fifteenth high column. At fixed valuation,
the leading powers `x^b` are distinct, so lower valuations cannot cancel
to produce omitted sources. The weight bound also excludes invisible
terms of valuation above fifteen. Therefore the tested source universe
is the entire stated filtered polynomial space.

For weight caps 18, 19, 20, 21, 22, 23 the independently computed pairs
`(low rank,rank after adjoining the high column)` are respectively

```text
(1,2),(1,2),(3,4),(3,4),(5,5),(5,5).
```

This proves the sharp threshold for a *nonzero* high coefficient. The
explicit signed convention is `R=L-p³y⁶`, hence `h=-1`; it agrees with
every coefficient in the displayed `L0`. At weight at most 22 there are
ten low source coordinates and five independent equations, hence exactly
five free replacement parameters. Four are the stated odd monomials;
the fifth is the even `K12` direction. At weight at most 23 there are
four additional odd directions. These are complete affine families.

The audit constructs six source columns: the affine base replacement
and the five free directions. It independently solves for their complete
tangent lifts and substitutes the whole six-column bank into **every**
original equation. This verifies every parameter value simultaneously,
not just a few rational points. It separately constructs and verifies
all ten homogeneous terminal directions.
The exact slice `v=52578/117601` is separately checked to recover the
earlier valuation-thirteen replacement, including its signs and all
source coordinates.

## 3. Fifth-row inheritance and the actual boundary consumer

The newly needed row is not inferred from the primary's numerical matrix.
I read its construction in the proved
[THM-4426 / source-normal-row-fourteen-weight-eighteen-memory-repair](../../01-canon/theorems/THM-4426-source-normal-row-fourteen-weight-eighteen-memory-repair.md)
companion and checked its dependence on the earlier source. With
`Phi=eta=0`, `Delta=896/15`, `Theta=512/75`, and `K=-32/5`,
THM-4308 gives the row-five target

```text
G5=-731648/2025+(6144/25)x²-(1952/45)x⁴.
```

The audit independently substitutes the report's `A4,C4` into the entire
row-four bracket and the predicted row-five source and obtains zero
residuals. If two capped pairs solved both equations, their difference
would be `theta_4(A0',C0')`, annihilated by the next-row `m=5` Student
operator. That operator is injective by
[THM-4308 / source-normal-bracket-hasse-truncation-through-row-eight](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md).
Thus the fifth row is unique. Later source faces begin after row five,
so they cannot alter it anywhere on the actual boundary. No later
boundary-producer import is needed for this check.

For a perturbation beginning at row eleven, source quadratic terms start
at row 22 and bracket quadratic terms at row 21. The entire tested finite
system is therefore exactly linear, and only background rows zero through
four can enter. Adding the classified response to any actual THM-4438
boundary point preserves all tested equations and projected depth.
Fixing its prefix through row ten gives the stated terminal fibre.

The source coordinates `z,h` of valuations nine and ten remain unchanged
and recover the boundary parameter; the four new odd coefficients and
the `p²y⁶` coefficient then recover the five replacement coordinates.
This proves the actual source parameter space `G_m x A^5`. A zero of
the original high coefficient changes neither the rank nor the family;
only the nonzero-response minimum becomes vacuous there. The inherited
irreducible-quartic hostile makes that coefficient nonzero at every
rational nonzero boundary parameter, exactly as stated in the report.

## 4. Reproduction and manifest

Standalone audit
[source](../../04-computation/planar_jc_long_20260906_memory_valuation12_audit.py)
and [output](planar_jc_long_20260906_memory_valuation12_audit.out):

```bash
python3 -B 04-computation/planar_jc_long_20260906_memory_valuation12_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_valuation12_audit.py
```

Both modes pass **5,146 live gates** and agree byte for byte. No producer,
repository mathematical module, random sample, floating-point computation,
or optimization-disabled assertion is used.

```text
source SHA256:
c54ba58ad5311cdcfdf354fbeaf8a290570f2d1d33561bdd92bae4728fb747ad
output SHA256:
851865545a2fea41ae7e2c05555ccd5a60ed135db467e65909c86a34a9f19867
semantic SHA256 (source columns, full relation bank, ranks, row counts):
85cb53d0490576e2d6508ffd6aa1357f5e2c5f0763b2d58809e73b3ce72162e6
```

The primary source now has1,646 explicit live controls, including all six
affine family columns in all271 original equations. Its normal/optimized
streams agree. Final byte hashes are in the session manifest; its semantic
matrix remains the one independently reproduced above.
