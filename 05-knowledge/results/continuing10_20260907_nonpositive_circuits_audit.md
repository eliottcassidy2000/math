# Independent referee: joint native D moments and the closed nonpositive octant

**Status: PASS after one prose-only transcription repair.** The complete continuous-region theorem, its weak boundaries, and every stated finite control are accepted. The general two-negative full-model classification remains open.

Primary proof: [continuing10_20260907_nonpositive_circuits.md](continuing10_20260907_nonpositive_circuits.md). The independent verifier first reconstructed the native algebra and full certificate without importing or executing the producer. The producer source was read only after that path passed.

## 1. Correction and native normalization

The first referee read found `-887y/7` in the displayed sixth D moment. Literal division gives `-915y/7`. The producer source, JSON, H4 determinant, and both sign polynomials already used the correct coefficient. The author corrected only the prose; no mathematical certificate or output changed.

With the precise monic B and D from the primary, the correct last displayed moment is

```
m6 = 7x^2/12 + 8969x/12 - 58463 - 915y/7 + 15z.
```

The referee reconstructs all native moments and computes H3 and H4 by literal permutation expansions of their determinants. This reproduces the complete polynomial identities, including

```
J3=-(343x^2-67788x+2592y+3157056)/1008,
partial_z^2 J4=-6.
```

The circuit formulas use the actual ordinary coefficients `(1,13,55,x,y,z)`, normalized by their binomial coefficients. Thus fixing the root sum and square sum at 13 and 59 fixes `R1=338/275`, and the three quotients are exactly those in the primary. The polynomial interlacer C and the circuit symbols C2,C3,C4 are not interchangeable.

For a weakly interlacing real-rooted monic denominator B, each repeated root of multiplicity m is a numerator root of multiplicity at least m-1. After cancellation, the proper rational function has only simple poles with nonnegative residues. Its entire ordinary Hankel sequence is positive semidefinite. Hence J3,J4 are necessary nonnegative quantities even with repeated beta roots or common-factor cancellation. The proof only uses this necessary direction; it does not infer full geometry from these two determinants.

## 2. Domain reduction and all boundary cases

The cleared circuit conditions force

```
x>=x0=831875/8788>0,
y>=ell(x)=13x^3/166375>0,
z>=k(x,y)=44y^3/x^3>0.
```

Consequently every denominator in the circuit interpretation and subsequent normalization is strictly positive. J3>=0 gives `y<=h(x)` exactly. Combining these bounds yields the stated cubic g(x)<=0. The independently reconstructed g(99) is positive, and all coefficients of g'(99+u) are strictly positive. This proves x<99, including exclusion of equality, on the entire unbounded x>=99 range; no sampled cutoff is used.

Every feasible x,y pair is covered by the primary closed-square substitution. At a collapsed fibre `h=ell`, v=0 suffices. At values with `h<ell` there is no feasible pair, but the polynomial certificate is valid even on that extra part of the image. Therefore inclusion of irrelevant reversed intervals does not remove a possible feasible point or change an inequality. Enlarging the x chart to include 99 is harmless.

The referee reconstructs the two positive normalizations exactly:

```
P0=-7112448*x^6*J4(x,y,k),
P1=-1008*x^3*partial_z J4(x,y,k).
```

Both are actual primitive integer polynomials. In particular the signs are not inferred from an unspecified rational denominator or an absolute-value normalization. P0 has all 15 stated terms and P1 has all six terms.

## 3. Independent full Bernstein identity

The independent path does not use the producer's power-to-Bernstein triangular formula. It expands every supplied tensor Bernstein basis polynomial back into the monomial basis, using exact fraction arithmetic, and sums the complete array. The result is then compared coefficient-for-coefficient with the native determinant polynomial after the independently reconstructed curved chart substitution.

Both full identities pass. The actual bidegrees are (18,6) and (9,3), with precisely 133 and 40 index pairs and no omissions or duplications. Every one of the 173 coefficients is strictly positive. The exact minima are

```
min P0 coefficients = 38259412467502741144816725/262144,
min P1 coefficients = 4049307478755/256.
```

The saved power packets, degree metadata and minimum values also match the independent reconstruction. Positivity of the Bernstein basis and its partition of unity on the closed square then proves strict positivity everywhere, including every edge, every corner, and every collapsed feasible fibre. This is a complete polynomial certificate; neither a numerical grid nor a check of a few vertices substitutes for the identity.

## 4. Unbounded z and exact consequence

The two positive normalizations give `J4(x,y,k)<0` and `partial_z J4(x,y,k)<0`. Literal expansion independently verifies

```
J4(x,y,k+w)=J4(x,y,k)+w*partial_z J4(x,y,k)-3w^2.
```

It is strictly negative for every w>=0. This contradicts J4>=0 on the entire proposed region, including the equality z=k. No upper bound on z, product, coefficient height or root separation is needed.

Thus the theorem excludes the **closed all-three-nonpositive octant**: every positive-product full B/C/D model point has at least one circuit quotient strictly above one. It does not assert that two quotients must be above one. At z=0 the finite fourth quotient is undefined; the primary neither defines it artificially nor needs to, because the cleared inequalities already force z>0. This is the correct boundary typing.

## 5. Hostiles and positive controls

All displayed rational controls are independently regenerated from their native x,y,z coordinates. Their four Newton ratios, three circuit quotients and all five C and D leading minors agree exactly with the certificate.

- The fourth-determinant-alone hostile has three negative circuits, all strict Newton inequalities, J4>0 and J3<0. It refutes only the discarded one-determinant shortcut.
- The inherited strict C-only all-negative point has J3>0 and J4<0. Independent Sturm isolation verifies its positive simple B roots and strict C interlacing. The all-one center has the same determinant sign split, confirming that ties do not escape the proof.
- `(95,68,1)` and `(86,50,9)` have two negative circuits and a positive-definite ordinary D matrix through moment degree six; both fail its degree-eight determinant. No B/C geometry is inferred for the first point. For the second, separate rational Sturm intervals certify all five positive simple B roots and strict C interlacing, as claimed.
- All four strict full-model positive controls retain independent B/C/D root ordering and the exact words `+++`, `++-`, `+-+`, and `-++`.

These controls preserve the remaining gap. In particular actual B/C geometry plus the full D matrix through degree six still permits `+--`. A proof excluding that word in the full model needs more information than this relaxation. No phase, inverse carry, integer support or original-response sign follows from the current theorem.

## 6. Reproduction, lineage, and connection

The independent source contains its own rational Sturm and interval routines, literal determinant expansion, ratio reconstruction, inverse Bernstein transformation and complete index-universe checks. Frozen source/output/JSON pins prevent a changed producer from silently entering the audit. Inputs may be adjacent to the referee or in the filed results directory. No producer is imported or executed.

```
python3 -B 04-computation/continuing10_20260907_nonpositive_circuits_audit.py
python3 -B -O 04-computation/continuing10_20260907_nonpositive_circuits_audit.py
```

The operative transfer is native rational-function moments -> their joint determinant region -> a curved coefficient chart -> exact positive polynomial certificates. Higher D moments and B/C geometry are discarded only as a relaxation, so an infeasibility proof still applies to the full model. The two-negative survivors identify precisely why that relaxation cannot be used as a complete word decoder.

The source, target, denominators, fibre width, tensor degree and z-tail are all retained. The strongest new consequence is the closed-octant exclusion, while the Anchor/full-model two-negative question and the original Laurent response remain open. The parent owns promotion and filing.

Both frozen referee runs pass **129 always-active exact gates**, with byte-identical actual LF stdout.

- Referee source SHA256: `ea3908f02331f64838659a50221d85d64f41f90acfe423d7a75666eb079941a6`.
- Normal and optimized stdout SHA256: `3ee32a5c264c8f8afc9fe467ed0764a2885f395aa971df14a5db903b7067fd92`.
- Audited producer source SHA256: `d392980046171db2ea1df6cb759c04882fd57305ae811fe3d7f5746df6b6a973`.
- Audited producer stdout SHA256: `36e5c21f2113d4a3dfbc825af03ddb0d3570a3e359c58ac5a3a78f164db9908b`.
- Certificate SHA256: `be4762cf527f1add2887ee2b29a98440a83f0c68454571b54d182b7395b58092`.
- Corrected producer prose SHA256 before promotion: `192a02eebe3c75193a1662b7d7238d37b26c4136bae017d10fa82744be696df0`.
