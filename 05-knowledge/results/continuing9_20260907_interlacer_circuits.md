# C-only interlacing still permits every circuit word at the anchored moments

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
At the fixed degree-five moments \(S=13,Q=59\), every ternary Newton-circuit
word occurs with five distinct positive beta roots and strict C interlacing.
There is an open neighborhood of the all-one circuit vector with those
properties that is entirely excluded by D. Four exact points certify the
strict full-model words \(+++,++-,+-+,-++\). No other full-model word is
excluded by this report.

[Independent audit](continuing9_20260907_circuit_interlacing_audit.md).

## 1. Native coefficient pencil and the inherited map

Retain exactly
\[
 B(v)=v^5-13v^4+55v^3-xv^2+yv-z,
\]
\[
 C(v)=v^4-12v^3+45v^2-\frac23xv+\frac37y,\qquad
 D(v)=v^4-11v^3+36v^2-\frac5{12}xv+\frac17y.
\]
The source is the fixed-first-two-moment chart in
[the main fixed-moment theorem](continuing9_20260907_fixed_moment_circuits.md).
Its \(d=5,q=275/338,\mu=13/5\) coefficients give the positive beta parameters,
whose ordinary symmetric coefficients are \(13,55,x,y,z\).
The exact circuit quotients are
\[
 C_2=\frac{831875}{8788x},\qquad
 C_3=\frac{13x^3}{166375y},\qquad
 C_4=\frac{44y^3}{x^3z}.
 \tag{1}
\]
The symbol C for an interlacer is distinguished from its subscripted circuit
quotients. No original Laurent phase or response is imposed in this sidecar.

The second source is the native degree-eight moment decoder,
[continuing4 moment packet](continuing4_20260906_moments_packet.md).
Write \(A/B=\sum_{j\ge0}m_jv^{-j-1}\), for \(A=C,D\), and
\(H_A=(m_{i+j})_{i,j=0}^4\). At simple real B roots \(\beta_i\),
\[
 H_A=V\,\operatorname{diag}
 \left(\frac{A(\beta_i)}{B'(\beta_i)}\right)V^T,
 \quad V_{ji}=\beta_i^j.
 \tag{2}
\]
The Vandermonde matrix is invertible. Thus \(H_A\) is positive definite
exactly when all five residues are positive, equivalently when A strictly
interlaces B. The source independently counts five positive simple B roots
for every declared point, then uses exact leading principal minors and
Sylvester's criterion. It never replaces the native moments by free ones.

## 2. Exact center and an open C-only neighborhood

The all-one circuit center has
\[
 x_0=\frac{831875}{8788},\quad
 y_0=\frac{3460080078125}{52206766144},\quad
 z_0=\frac{791547381622314453125}{52414354446428935168}.
 \tag{3}
\]
All five leading principal minors of \(H_C\) are strictly positive.
Their exact values, including the two large final fractions, are frozen
in the certificate. The first three are
\[
 1,\quad2,\quad\frac{10176837020375}{822256566768}.
 \tag{4}
\]
The D leading minors have signs \(+,+,+,-,-\). In particular its ordinary
moment packet through degree four is positive definite, with third minor
\(189400546772609/13156105068288\), but its ordinary degree-six determinant is
\[
 \det(m_{i+j}^{D})_{i,j=0}^3
 =
 -\frac{22803919898205048846211421605604713663080078125}
        {76327251049212984252697527408027867930427392}<0.
 \tag{5}
\]
This is a sufficient exact rejection of D; no approximation of an individual
negative residue is needed. It is not claimed to be the first failed
condition among every possible shifted or localizing packet.

The main chart is continuous in the circuit coordinates near \((1,1,1)\),
fixes \(S,Q\), and preserves positive simple beta roots on an explicit
neighborhood. The five strict C minor inequalities and the strict D
inequality (5) persist on a possibly smaller open neighborhood. Equation (2)
therefore proves:

**Theorem.** Some open neighborhood of \((C_2,C_3,C_4)=(1,1,1)\) is feasible
with five positive simple beta roots and strict C interlacing, while every
point in that neighborhood fails the D interlacing predicate.

This supplies all ternary signs by choosing small independent positive,
negative, or zero perturbations of the three circuit coordinates. The
additional explicit certificate below avoids leaving the witnesses merely
existential. It also pinpoints the lost predicate in an attempted transfer
from the fixed-moment theorem to the full Laurent model.

## 3. Explicit witnesses for all 27 C-only words

For each \(\sigma\in\{-1,0,1\}^3\), prescribe
\[
 C_j=1+\sigma_j/2048
\]
and use the exact coefficient chart of the main theorem. All 27 resulting
rational triples \((x,y,z)\) are saved in the certificate. They are the same
27 main-theorem degree-five targets, with the same scaling, rather than
independently fitted roots.

Exact computation proves, separately for each of these 27 points:

- B has five distinct positive roots;
- all five leading C-Hankel minors are positive;
- the ordinary degree-six D-Hankel determinant is negative;
- the circuit quotients are exactly the prescribed values.

Thus every ternary word has an explicit strict C-only witness, and each
displayed witness fails D. This finite certificate does **not** assert that
the entire closed radius-\(1/2048\) cube is C-admissible. The smaller open
neighborhood theorem in Section 2 follows from continuity, not from checking
these finitely many points.

## 4. Four exact full-model positive controls

The following rational points have five positive simple beta roots and
strict simultaneous C,D interlacing. The terminating decimals in the table
denote exact rational numbers.

| Strict circuit word | x | y | z |
|---|---:|---:|---:|
| +++ | 77.454 | 8.902 | 0.02558 |
| ++- | 77.3613 | 8.6001 | 0.17694 |
| +-+ | 86.2333 | 51.3919 | 8.6469 |
| -++ | 97.7028 | 70.6021 | 14.5020 |

At each point all five leading minors of each native Hankel are positive.
The complete moments, minors, and quotient signs are retained in the
certificate. The parent's numerical full-domain scout supplied these
locators; their acceptance here is wholly exact.

The statement that no full-model point has two negative circuits remains
**OPEN**. The four positive controls do not establish a complete sign
classification. No inference from the number or distribution of numerical
samples enters either theorem above.

## 5. Reproduction and connection boundary

The source-to-target map preserves \(13,55,x,y,z\), root positivity, native
C residues, and the exact circuit values. It loses D unless checked
separately. Even simultaneous interlacing still does not impose an original
phase, inverse carry, or response sign. The strongest proved survivor is:
two fixed moments plus one of these interlacers cannot exclude a circuit
sign word; the second interlacer removes an entire open neighborhood where
all 27 words otherwise occur.

Run the [standalone source](../../04-computation/continuing9_20260907_interlacer_circuits.py):

~~~text
python continuing9_20260907_interlacer_circuits.py
python -O continuing9_20260907_interlacer_circuits.py
~~~

The [certificate](continuing9_20260907_interlacer_circuits_certificate.json)
contains the center, 27 C-only witnesses, four full-model controls, native
moments through degree eight and all relevant leading minors. Its input
pin is the frozen main coefficient certificate; no producer is imported.
Both [normal output](continuing9_20260907_interlacer_circuits.out) and
[optimized output](continuing9_20260907_interlacer_circuits_optimized.out)
pass 214 always-active exact gates and have identical raw LF bytes.
The source sets stdout newline handling explicitly.

The parent owns audit promotion and filing. All claims about the absence
of additional full-model circuit words remain open.
