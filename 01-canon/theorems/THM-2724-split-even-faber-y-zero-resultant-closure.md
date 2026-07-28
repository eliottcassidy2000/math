---
id: THM-2724
title: "Split even-Faber y-zero resultant closure"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT, AWAITING
  INDEPENDENT HOSTILE AUDIT.  On the y identically zero, lambda!=0 chosen-
  sheet split polynomial exact-square-prefix degree-twenty-two even-Faber
  boundary, the first two fluxes have a degree-23 q-resultant with
  parameter-independent leading coefficient 96059601 and constant term a
  nonzero multiple of lambda^3.  Thus q and u lie in the constant field, so
  the third Faber observable is constant, contradicting R_Q'=kappa/U.  The
  y-nonzero chart, every odd Faber seed, lambda=0, the full split branch, and
  JC(2) remain outside this candidate.
source: root/split-even-y-zero-resultant-2026-07-28
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
related:
  - THM-2718-split-prime23-five-pole-rational-primitive-closure
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2725-split-even-faber-nonzero-first-flux-unified-closure
script: 04-computation/jc2_split_even_faber_y_zero_resultant_thm2724.py
output: 05-knowledge/results/jc2_split_even_faber_y_zero_resultant_thm2724.out
script_sha256: faee30ff7afb307447fb6082dc16e6ec056dc27ef123ab66a1171f0ce30502ae
output_sha256: e597f26ab2cdd99cb6ec0d971f4a087f47b8d0ad1d75040a6af9fab7a6074c8f
hash_basis: working-tree bytes (LF)
---

# THM-2724 -- the split even-Faber `y=0` boundary is empty

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT, AWAITING
INDEPENDENT HOSTILE AUDIT.**  Nothing may depend on this candidate until its
status is promoted.

## 1. Exact statement

Let `C` be an algebraically closed field of characteristic zero and put
`K=C(x)`.  Consider a polynomial Keller pair in the split polynomial
exact-square-prefix terminal chart whose reduced, target-translated mate is

```text
Q=E_22+B E_14+C E_10+D E_6+E E_2,                  (1)
```

with `B,C,D,E in C`; write `W in C` for its second flux constant.  Thus all
eleven odd Faber coefficients are zero.  On one chosen split sheet assume

```text
Phi_Q=lambda in C*,              Psi_Q=W.            (2)
```

Use THM-2411's centered coordinates

```text
y=11s_ctr,             u=d_ctr T,             Z=T^2,
q^2=T.                                                   (3)
```

Then the boundary

```text
y=0 identically                                            (4)
```

contains no physical Keller trajectory.

The constants in (1) may vanish arbitrarily; no `B!=0` hypothesis is used.
The theorem does not treat `y!=0`, `lambda=0`, or any odd Faber seed.

## 2. The two exact fluxes at `y=0`

THM-2411 gives the primitive integer flux numerators `N1,N2`.  Setting
`y=0` yields

```text
N1=1331(616B-1089u)Z+4(2342560Cu-3748096E),          (5)

N2=15944049Z^2-206145280CZ
   +1443016960Bu^2-1978994688Du-1319329792W
   -1190488992u^3.                                   (6)
```

The chosen-sheet first-flux equation and the second constant flux are

```text
q N1=-7496192 lambda,                   N2=0.          (7)
```

Because `Z=T^2=q^4`, define

```text
G1(q,u)=q N1(q^4,u)+7496192 lambda.                   (8)
```

It is linear in `u`; explicitly,

```text
G1=
-14641q(-640C+99q^4)u
+117128(7Bq^5-128Eq+64lambda).                        (9)
```

The polynomial `N2(q^4,u)` is cubic in `u` with fixed leading coefficient

```text
[u^3]N2=-1190488992!=0.                              (10)
```

Equation (9) also records immediately that `q=0` is impossible when
`lambda!=0`.

## 3. The degree-23 eliminant never vanishes identically

Direct Sylvester elimination in `u` gives

```text
Res_u(G1,N2)=505447028499293771 P23(q),               (11)
```

where the primitive factor, grouped by `q`-degree, is

```text
P23=
 96059601 q^23
-3104956800 C q^19
+(1483603968 B^3-6744342528 BD+36130406400 C^2-7948689408 W)q^15
+(-17983078400 B^3C-30519853056 B^2E+87199580160 BCD
   -181665792000 C^3+154156400640 CW+123325120512 DE)q^11
+(15259926528 B^2-61662560256 D)lambda q^10
+(657666867200 B^2CE-281857228800 BC^2D-372051542016 BE^2
   +335544320000 C^4-996566630400 C^2W-1594506608640 CDE)q^7
+(-328833433600 B^2C+372051542016 BE+797253304320 CD)lambda q^6
-93012885504 B lambda^2 q^5
+(-6012954214400 BCE^2+2147483648000 C^3W
   +5153960755200 C^2DE+7937099563008 E^3)q^3
+(6012954214400 BCE-2576980377600 C^2D
   -11905649344512 E^2)lambda q^2
+(-1503238553600 BC+5952824672256 E)lambda^2 q
-992137445376 lambda^3.                               (12)
```

The proof needs only the two pivots

```text
[q^23]P23=96059601!=0,
[q^0]P23=-992137445376 lambda^3!=0.                   (13)
```

The first makes `P23` a nonzero polynomial for every specialization of
`B,C,D,E,W`; the second makes its nonzero-`q` boundary explicit.  The exact
support and size are

```text
supp_q(P23)={0,1,2,3,5,6,7,10,11,15,19,23},
terms(P23)=34.                                        (14)
```

The companion computes (11) directly and independently reconstructs the
same polynomial, up to the standard resultant sign convention, by
substituting the unique linear root of (9) into the cubic (6).  Its canonical
polynomial digest is

```text
e50e733ec85df8c1aff2bf805a6442424b82b3c3dbb4cb5563cea363c9f7b6e0.
```

Thus every physical solution of (7) satisfies

```text
P23(q)=0.                                             (15)
```

## 4. Constant-field descent

All coefficients of `P23` lie in `C`.  Equation (15) makes
`q in C(x)` algebraic over `C`.  Because `C` is algebraically closed, the
constant field of `C(x)` is exactly `C`, so

```text
q in C*.                                              (16)
```

The nonzero assertion follows either from (13) or directly from (7).
Consequently

```text
T=q^2 in C*.                                          (17)
```

Now (6) is a nonzero cubic polynomial over `C` in the single rational
function `u`, with the fixed leading coefficient (10).  Hence `u` is
algebraic over `C`, and again

```text
u in C,                  d_ctr=u/T in C,
s_ctr=y/11=0.                                           (18)
```

Every centered quartic coordinate and every coefficient entering the third
Faber observable is therefore constant.  It follows that

```text
R_Q in C,                       R_Q'=0.                (19)
```

## 5. The third flux gives the contradiction

On the chosen split sheet write the exact-square-prefix leading coefficient
as `V_src=U^2` and use the depressed coordinate

```text
w=Uz+B_src/(2U).                                      (20)
```

Here `U` is a nonzero element of `C(x)` (indeed it may be chosen in `C[x]`).
THM-2129's exact Hamiltonian identity and the Keller equation give

```text
J_(x,w)(P,Q)
 =(w^2+p/4)Phi_Q'+w Psi_Q'+R_Q'
 =kappa/U,                         kappa in C*.       (21)
```

The first two fluxes in (2) are constants, so (21) reduces to

```text
R_Q'=kappa/U!=0.                                      (22)
```

This contradicts (19).  Therefore the physical locus (4) is empty.

No nonsplit parity enters the proof: the first flux is the nonzero split
constant `lambda`, and only its derivative vanishes in (21).

## 6. Exact boundary and reproduction

The conclusion is precisely

```text
{split degree-22 even-Faber Keller trajectories with
 y identically zero and lambda!=0}=empty.             (23)
```

The proof does **not**:

1. treat the `y!=0` scale chart;
2. treat `lambda=0`;
3. retain any of the eleven odd Faber seeds;
4. close the full split degree-twenty-two branch; or
5. prove `JC(2)` or `DC(2)`.

Run

```bash
python3 04-computation/jc2_split_even_faber_y_zero_resultant_thm2724.py
python3 -O 04-computation/jc2_split_even_faber_y_zero_resultant_thm2724.py
```

Both modes byte-match the declared output.  The script is self-contained,
uses exact rational polynomial arithmetic, contains no Python `assert`, and
checks (5)--(14) by the two independent eliminant paths described above.
