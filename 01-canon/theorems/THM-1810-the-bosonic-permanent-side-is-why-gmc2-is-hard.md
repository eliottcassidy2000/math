---
id: THM-1810
title: "Wick permanents and tournament score products are distinct unsigned lifts"
status: >
  CORRECTED 2026-10-23. PROVED algebraic identities: normalized complex Gaussian
  moments equal a! = per(J_a); det(J_a)=0 for a>=2; the unsigned all-tournament
  score polynomial is product_(i<j)(x_i+x_j), whereas the Vandermonde permanent
  counts only transitive tournaments. The former equality of these two unsigned
  polynomials is REFUTED minimally at n=3. Positive combinatorial Wick signs do
  not prohibit coefficient cancellation. Degree-uniform detection lower bounds
  are inherited only from the repaired THM-1790. Complexity and bosonic/fermionic
  comparisons are heuristic analogies, not reductions or constraints on all proofs.
source: klein-2026-07-20-S387; corrected by session37 independent inheritance audit, 2026-10-23
depends_on:
  - THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels
  - THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one
related:
  - THM-1800
  - THM-1475
  - THM-1770-gmc2-detection-depth-grows-with-radial-degree
attribution: >
  The original fermionic/transitivity framing is credited to THM-1805 and boxeph's
  THM-1800; the Gaussian permanent complement was proposed by klein-S387.
  This correction distinguishes two different unsigned lifts and withdraws the
  invalid no-cancellation and complexity implications. The repaired projective
  detection-depth theorem belongs to the existing THM-1790/MISTAKE-246 lineage.
---

# THM-1810 — distinguish the two unsigned lifts

> **Correction, 2026-10-23.** The previous version identified the Vandermonde
> permanent with the all-tournament score sum and inferred absence of moment
> cancellation from positive Wick signs. Both assertions are false. The exact
> survivors, minimal witnesses and corrected scope are recorded below.
> [Independent proof and exact receipts](../../05-knowledge/results/session37_inheritance_audit_20261023/REPORT.md).

## 1. The Gaussian repeated-state identities survive

For a normalized circular complex Gaussian Z,

```
E[Z^a conjugate(Z)^b] = 0                 if a!=b,
E[|Z|^(2a)] = a! = per(J_a),
det(J_a) = 0                            for a>=2.
```

Here J_a is the a-by-a all-ones matrix. Wick expansions have positive
combinatorial matching signs; this is the precise bosonic/permanent analogy.
It does not say that arbitrary polynomial coefficients or covariance products
are positive, or that polynomial moments cannot cancel.

## 2. All tournaments versus transitive tournaments

Fix the column convention `V_ij=x_i^(n-1-j)` for 0<=i,j<n. Then

```
det V = product_(i<j)(x_i-x_j)
      = sum_(all tournaments T) (-1)^(#upsets(T)) x^score(T),

U_n(x) = product_(i<j)(x_i+x_j)
       = sum_(all tournaments T) x^score(T),

per V = sum_(transitive tournaments T) x^score(T).
```

An upset means the larger-index vertex beats the smaller-index vertex. If
ascending exponent columns are used, the determinant acquires the global
factor `(-1)^(n choose 2)`; its permanent is unchanged.

**Proof.** The two edge products choose a winner variable for every edge.
The permanent instead sums assignments of the distinct exponents
0,1,...,n-1 to vertices. These are precisely transitive tournament score
vectors: the maximum-score vertex beats all others, and removal gives
induction. Conversely a total order supplies one such tournament and score
vector. Therefore `U_n-per V` is exactly the generating polynomial of
intransitive tournaments, with nonnegative coefficients.

The unsigned equality holds only for n<=2. Its smallest failure is

```
n=3: U_3(x) = per V + 2*x_0*x_1*x_2.
```

At all ones, `U_n=2^(n choose 2)` whereas `per V=n!`; already the values are 8 and 6 at n=3.
Independent exact coefficient dictionaries verify all three identities for
n=1,...,5. Thus the signed cancellation identity survives, but stripping the
edge-product signs is not the same operation as replacing det(V) by per(V).

## 3. Positive Wick signs do not prevent cancellation

Let `L(s^k)=k!`. The degree-one polynomial `B=s-1` gives

```
L(B)=1!-0!=0,              L(B^2)=1.
```

So the previous sentence "no sign, so nothing cancels" is false for the
polynomial moment functional. The signs or phases can come from coefficients.
The already-repaired THM-1790 supplies a sharper exact control:

```
H=s(s-2),  P=W+Z(ZW-2),
E[P]=E[P^2]=E[P^3]=0,      E[P^4]=48.
```

Its general detection-depth lower bounds use homogeneous projective
existence in dimension d+1 and eventual moment nonvanishing. They do not
follow from signlessness. Keep the correction in MISTAKE-246 and
`THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one.md` as the
load-bearing source, including its distinction between lower bounds and
unproved general exact upper depths.

## 4. Scope and historical boundary

The original note was a structural synthesis, not a complexity reduction.
The general determinant/permanent complexity contrast does not prove a
complexity bound for GMC(2), and does not establish that every possible proof
must be bosonic or exclude all fermionic methods. Nor is constant-term
projection automatically a determinant functional. Those identifications
are withdrawn as mathematical deductions; they may be retained only as
explicitly heuristic comparisons with no claimed preservation theorem.

This corrected note makes no new assertion about the current status of
GMC(2). Its old "remains OPEN" wording described its historical date and is
not a current-status authority. Consult current canon for that separate topic.

The original `bosonic_fermionic_klein_S387.py` and output were unavailable
in the audited snapshot. The earlier claim that the false equality had been
"verified n=3,4" is superseded by the exact minimal witness and replacement
script in the [Session 37 packet](../../05-knowledge/results/session37_inheritance_audit_20261023/REPORT.md). Attribution and the valid Gaussian,
signed-tournament and projective-depth results are preserved.
