---
id: THM-4445
title: "LRC14 signed (1,1,1) sharp obstruction classification"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT + INDEPENDENTLY
  AUDITED. The signed (1,1,1) family is exactly a+b=c. The row (1,4,5)
  has min E=physical mass=1/28; every other row lies sharply between
  31/392 and 6/55, with lower equality (1,7,8) and upper equality
  (1,10,11). Thus (1,4,5) is the only row at or below 6/77, and there is no
  equality row. This classifies a persistent local obstruction; entry,
  synchronization, and LRC(14) remain open.
source: root additive-circuit continuation + independent clean-room referee, 2026-09-06
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
related:
  - THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits
  - THM-4442-lrc14-bounded-ten-body-parity-free-scale-three-completion
  - THM-4444-lrc14-signed-112-sharp-one-ray-classification
primary_script: 04-computation/lrc14_signed_111_sharp_obstruction_thm4445.py
primary_output: 05-knowledge/results/lrc14_signed_111_sharp_obstruction_thm4445.out
primary_script_sha256: b6580784528167045fe5704180b93b738cbdf4ccee6bfd2496a12adf877c4ee1
primary_output_sha256: 188eb03f9b9c31fc9330282df26cd07e56ab5e69c7cd9ccdcd4934f11f141d75
independent_script: 04-computation/lrc14_signed_111_sharp_obstruction_thm4445_independent.py
independent_output: 05-knowledge/results/lrc14_signed_111_sharp_obstruction_thm4445_independent.out
independent_script_sha256: d564e61e92c9195fd098b86077efe144b6852d3110e07c846cd7e9be9a7966cd
independent_output_sha256: 53a47cade027fae563be6f75b040cb5d05c144c517b76c95a81f05012c67c248
report: 05-knowledge/results/lrc14_signed_111_sharp_obstruction_thm4445.md
audit: 05-knowledge/results/lrc14_signed_111_sharp_obstruction_thm4445_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4445 -- LRC14 signed (1,1,1) sharp obstruction classification

**PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT + INDEPENDENTLY
AUDITED.** This completes the sharp local circuit classification, but proves
that the additive circuit is a persistent obstruction rather than removing
it. `LRC(14)` remains **OPEN**.

## 1. Sharp statement

Let \(w=(a,b,c)\) be primitive, sorted, distinct, positive, and ternary-unit,
and suppose it has a signed relation with coefficient magnitudes
\((1,1,1)\). Then necessarily and sufficiently
\[
                              a+b=c.                               \tag{1}
\]

For \(w=(1,4,5)\),
\[
                       \mu(F_w)=\min_iE_i(w)={1\over28}.            \tag{2}
\]
Every other eligible row satisfies the sharp band
\[
 {31\over392}\le \mu(F_w)\le\min_iE_i(w)\le {6\over55}.             \tag{3}
\]
Equality at the lower endpoint in either observable occurs only for
\(w=(1,7,8)\); equality at the upper endpoint in either occurs only for
\(w=(1,10,11)\).

Since
\[
 {31\over392}-{6\over77}={5\over4312}>0,
\]
the row \((1,4,5)\) is the only signed \((1,1,1)\) row at or below \(6/77\)
for either observable. There is no equality row at \(6/77\).

## 2. Complete ray

Modulo global reversal, positivity and strict sorting leave only the sign
vector \((1,1,-1)\), proving (1). Put \(u=(1,1,-1)\). If a complete carrier
is \(C=(x,y,z)\), then
\[
                 a(x+z)+b(y+z)=0.
\]
Primitivity gives \(\gcd(a,b)=1\), so for an integer \(h\),
\[
 x+z=bh,\qquad y+z=-ah,\qquad x-y=ch.
\]
The strict \(x,y\) roofs imply
\[
 |ch|\le |x|+|y|
 <{3[(b+c)+(a+c)]\over14}={9c\over14}<c.
\]
Hence \(h=0\), and the entire carrier set is exactly
\[
 \Lambda(w)=\{k(1,1,-1):0<|k|<3c/14,\ 3\nmid k\}.                  \tag{4}
\]
This preserves the deleted-third address and every strict endpoint.

## 3. Profiles and two-sided quadrature

Put \(r=3/14\), \(t=a/c\), and extend the profiles by zero for \(s\ge r\):
\[
\begin{aligned}
 f_a(s)&=\min\!\left(2r,{r(2-t)-s\over1-t}\right),\\
 f_b(s)&=\min\!\left(2r,{r(1+t)-s\over t}\right),\\
 f_c(s)&=\min\!\left(2r,{r-s\over t(1-t)}\right).
\end{aligned}                                                     \tag{5}
\]
The \(b\)-profile is pointwise at least the \(a\)-profile. Exact integration
gives
\[
\begin{aligned}
 {4\over3}\int f_a&={9+3t\over98},\\
 {4\over3}\int f_c&={12(1-t+t^2)\over98},\\
 {4\over3}\int\min_i f_i&={9\over98}.                              \tag{6}
\end{aligned}
\]
The minimum of the first two continuum projections ranges from \(9/98\) to
\(39/392\), with upper crossing at \(t=1/4\). The physical bulk is the
ratio-independent \(9/98\).

For \(R_<(T)=\#\{1\le k<T:3\nmid k\}\), three-term blocks give
\[
 {2T\over3}-{2\over3}\le R_<(T)<{2T\over3}+{2\over3}.
\]
Layer cake across profile height \(2r\) yields, for every profile in (5),
\[
 L_f-{4\over7c}\le {2\over c}\sum_{\substack{k\ge1\\3\nmid k}}f(k/c)
 <L_f+{4\over7c},\qquad L_f={4\over3}\int f.                       \tag{7}
\]
The lower half of (7) is the new decisive sidecar. It proves that both
\(\mu(F_w)\) and \(\min E_i(w)\) exceed \(6/77\) from the first admissible
height \(c=43\), and exceed \(31/392\) beyond the finite sharp head.

## 4. Exact head and persistent obstruction

The complete raw head \(c\le44\) has 77 rows. It proves (2) and the lower
edge of (3), including the unique equality \((1,7,8)\). The upper half of
(7), together with the continuum projection ceiling \(39/392\), reduces the
global \(6/55\) upper bound to a finite head; its unique maximizer is
\((1,10,11)\).

For an orthogonal parity check, every primitive additive row has one even
coordinate. A 1,901-row head through \(c=223\), followed by the exact upper
tail, gives these sharp maxima:

| even coordinate | sharp \(\max\min E\) | witness | sharp \(\max\mu\) | witness |
|---|---:|---|---:|---|
| \(c\) | \(2946/28861\) | \((7,31,38)\) | \(222/2275\) | \((1,25,26)\) |
| \(b\) | \(6/55\) | \((1,10,11)\) | \(6/55\) | \((1,10,11)\) |
| \(a\) | \(223/2156\) | \((4,7,11)\) | \(102/1001\) | \((2,11,13)\) |

Each parity sheet has an infinite primitive ternary-unit family with
\(a/c\to1/4\). By (7), their physical masses tend to \(9/98\), while their
minimum projections tend to \(39/392\). Thus the obstruction is cofinal and
cannot be repaired by a larger finite head.

## 5. Independent audit

The clean-room referee imports neither the primary script nor a repository
geometry engine. It uses a modular one-coordinate kernel solver, a separate
ray compiler, and a generic rational breakpoint integrator. Across all 1,901
rows it reproduces the complete raw ray, profiles, two-sided bound, threshold
locus, sharp band, parity leaders, and first strict tail cutoffs in 417,669
checks.

The arithmetic lower cutoff against \(6/77\) is \(c=42\), which is excluded
because \(3\mid42\); hence \(43\) is exactly the first admissible height.
Weak lower and strict upper endpoints in (7) are both audited.

## 6. Consequence and scope

Combining THM-4437, THM-4441, THM-4444, and this theorem gives the exact
local \(6/77\) classification:

- all additive \((1,1,1)\) rows except \((1,4,5)\) are strict hostiles;
- \((2,11,20)\) is the sole nonadditive strict hostile;
- \((1,5,11)\) is the sole nonadditive boundary packet;
- every other primitive ternary-unit triple is strictly below \(6/77\).

This does not settle whether a body good set intersects the physical
failure set. THM-4442 handles every tail only after entry into its bounded
ten-body chart. Arbitrary entry, synchronization, and LRC(14) remain open.

## 7. Reproduction

```powershell
python -B 04-computation/lrc14_signed_111_sharp_obstruction_thm4445.py
python -B -O 04-computation/lrc14_signed_111_sharp_obstruction_thm4445.py
python -B 04-computation/lrc14_signed_111_sharp_obstruction_thm4445_independent.py
python -B -O 04-computation/lrc14_signed_111_sharp_obstruction_thm4445_independent.py
```

All theorem decisions are exact integer or rational comparisons.
