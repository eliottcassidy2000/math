# Projective receipt rays and a source-owned nonzero-clock family

Started 2026-10-09; completed 2026-10-11. **PROVED:** the elementary all-height join, source-owned normalization, projective/lattice distinction, and repeated-family identities below. **FINITE-EXACT:** the retained 1,949-label scan and specified successor guards. **INHERITED FINITE-EXACT:** the common-future class of THM-4605, not replayed here. **OPEN:** ROOT for its giant sources and for the new terminal ports.

Artifacts: [checker](../../04-computation/experiments/collatz_projective_receipts_20261009.py) and [output](collatz_projective_receipts_20261009.out).

## 1. Inheritance and the retained coordinates

The [child-cut affine-family compiler](child_cut_obligation_partition_20261004.md), §2, already proves that compatible ordered valuation words define an exact arithmetic ray of common endpoints and sources. This note uses that mechanism; affine parameterization itself is not new. The inherited [223/233 packet](collatz_prime_partition_223_20261007.md), §5, supplies the actual endpoint 425 and the word from 223 used below.

The [native-cell observer](collatz_debt_observer_integration_20261009.md) requires source evaluation to distinguish a singleton meeting from equality on a whole cell. The [fibre ledger](collatz_fibre_integration_20261009.md) retains both actual clocks. Its normalized Mersenne defect is not reused for a child that is no longer Mersenne. The [three-generator normal form](collatz_child_normal_forms_20261007.md) supplies the exact $G_1,G_5,G_{17}$ successor guards.

The canonical hostile here is an integral convex midpoint satisfying both endpoint equations and payment, but violating the supplied words. The corrected near miss is treating the real cone as the complete integer certificate. The least-used sidecar is the marked arithmetic lattice. The live board is: incidence plane; native odd lattice; immutable source; two clock differences; terminal alternatives; common-future normalization.

| Source → target | Map / preserved predicate | Loss and required sidecar | Decisive test |
|---|---|---|---|
| Receipt ray → simplex segment | Normalize homogeneous coordinates; real incidence/payment persist | Integer parameter and first-hit domain | Midpoint $(445,351,668)$ |
| Supplied source → run child | Exact ternary division; native words and original source retained | ROOT still needs a child proof | $283\to425\leftarrow223$ |
| Several run choices → maximal run | Actual $G_1$ chain preserves the terminal obligation | Retain connecting words and clock costs | Equation (7) |
| Giant exponent → symbolic child | Full-denominator modular reader | Not a ROOT finder | 36 exact exponent labels |
| Two joins → composed join | Match the same middle prefix | Independent endpoint choices are invalid | Two-return ray (12) |

## 2. A nongenealogical paid family

Use the accelerated odd map $U(n)=(3n+1)/2^{v_2(3n+1)}$. Let $r\ge1$, let $t>0$ be odd, and impose
\[
t\equiv5\,3^{-(r+1)}\pmod8.
\]
Thus $t\equiv7\pmod8$ for even $r$, and $t\equiv5\pmod8$ for odd $r$. Set
\[
n=\frac{3^rt-1}{2},\qquad m=2^{r+1}t-1,\qquad
J=\frac{3^{r+1}t-1}{4}.                                      \tag{1}
\]
These are positive odd integers with the actual diagram
\[
n\xrightarrow{(1)}J\xleftarrow{(1^r,3)}m.                     \tag{2}
\]

After $j=0,\ldots,r$ one-valuations on the child its state is
\[
m_j=3^j2^{r+1-j}t-1.
\]
For $j<r$, its next numerator has exact 2-adic valuation one. At $j=r$, that numerator is $2(3^{r+1}t-1)$, whose valuation is exactly three by the modulo-8 condition. Its endpoint is $J$. The same condition makes $n$ odd with next valuation one.

The exact payment difference is
\[
n-m=\frac{(3^r-2^{r+2})t+1}{2}.
\]
For $r=1,2,3$ it is negative; for every $r\ge4$ it is positive. Thus the sharp paid regime is $r\ge4$, and there
\[
1<m<n<J.                                                     \tag{3}
\]
The source rises on its supplied edge. The child prefix does not pass through the parent: $n\equiv1\pmod3$, whereas every child endpoint after at least one of its one-valuations, and its final $J$, is $2\pmod3$. Initially $m<n$. All child prefixes exceed 1, so both prefixes are strict pre-ROOT.

The independent carry identity is
\[
F_{(1^r,3)}(m)=
\frac{3^{r+1}m+3^{r+1}-2^{r+1}}{2^{r+3}}.
\]
The raw source-minus-child clock differences are
\[
(\text{Terras cost difference},\text{odd-edge difference})
=(-r-2,-r).                                                   \tag{4}
\]
If either source is rooted, the same differences hold for the complete first-hit clocks because the common suffix cancels. Without a ROOT proof, (4) still measures the actual finite prefixes; it is not an operation on infinite stopping times.

The least $r=4$ example is $283\to425\leftarrow223$. A supplied, authenticated ROOT word for 223 has odd rank 24 and cost 46; stripping its five-letter prefix and prepending $(1)$ yields ROOT(283), rank 20 and cost 40. No all-height ROOT claim follows.

## 3. Source-owned fuel and exact normalization

For a supplied positive odd $n$, put $R=v_3(2n+1)$. The rule applies exactly for
\[
n\equiv3\pmod4,\qquad4\le r\le R.
\]
Its source domain is precisely
\[
n\equiv283\pmod{324}.                                       \tag{5}
\]
No new source is chosen. The offered child is
\[
m_r=\frac{2^{r+1}(2n+1)}{3^r}-1,\qquad
m_{r+1}+1=\frac23(m_r+1).
\]
Hence $r=R$ gives the smallest child.

There is **no new source-guard mass over the inherited $G_5$ bank**: (5) is contained in $n\equiv4\pmod9$, where
\[
G_5(n)=\frac{8n-5}{9}=4\,3^{r-2}t-1.
\]
For $r\ge4$, $2^{r-1}<3^{r-2}$, and therefore
\[
m_r<G_5(n)<n.                                               \tag{6}
\]
Moreover $v_3(n+5)=2$ and $G_5(n)+5=8(n+5)/9$ is a 3-unit, so repeated $G_5$ cannot reproduce these ports. At source 283 the old rule chooses 251 instead of 223.

Maximal-run selection is stronger than a size heuristic. For every $r<R$, the cofactor $(2n+1)/3^r$ is divisible by 3, so $m_r\equiv2\pmod3$ and
\[
\boxed{G_1(m_r)=\frac{2m_r-1}{3}=m_{r+1}.}                   \tag{7}
\]
Equivalently $m_R\xrightarrow{1^{R-r}}m_r$ is actual. Thus every offered run choice has the same ROOT obligation, with an explicit connecting prefix. Selecting $R$ loses no terminal alternative **within this family**. That conclusion depends on (7); mere guard containment or smaller size would not prove it.

## 4. Projective compactification and its arithmetic boundary

For fixed $r\ge4$, let $t=t_r+8s$, $s\ge0$, with $t_r=7$ or 5 as above. Then
\[
(n,m,J)=(n_0,m_0,J_0)+s(4\,3^r,\ 2^{r+4},\ 2\,3^{r+1}).       \tag{8}
\]
Homogenize as $z=(n,m,J,h)$. Incidence is the pair of linear equations
\[
3n+h-2J=0,\qquad
3^{r+1}m+(3^{r+1}-2^{r+1})h-2^{r+3}J=0,                     \tag{9}
\]
and payment is $n-m>0$. Write $z_0=(n_0,m_0,J_0,1)$ and let $d$ be the period vector in (8) with last coordinate zero. There is an exact simplex decomposition
\[
\frac{z_0+sd}{\mathbf1^T(z_0+sd)}
=\lambda_s\frac{z_0}{\mathbf1^Tz_0}
 +(1-\lambda_s)\frac d{\mathbf1^Td},\qquad
\lambda_s=\frac{\mathbf1^Tz_0}{\mathbf1^Tz_0+s\mathbf1^Td}.    \tag{10}
\]
It preserves real incidence and payment on a compact segment. The endpoint $d$ has $h=0$: it is a limiting direction, not an ordinary finite source.

The segment does not retain the arithmetic lattice $s\in\mathbb Z_{\ge0}$, odd endpoint, or actual words. The midpoint of the first two $r=4$ rows is $(445,351,668)$. It consists of positive integers, satisfies (9) with $h=1$, and obeys $351<445<668$, but source 445 has actual valuation three rather than one. Convex mixing preserves these real equalities while destroying the certificate. Choosing source and child coordinates independently from different rows—a box operation—can destroy incidence itself.

The positive repair is exact: retain $h=1$, the compatible odd-integer lattice of the inherited affine-family compiler, and strict first-hit cuts. Positive odd endpoint integrality then forces all valuations in the supplied inverse words. Equation (10) is a chart for those marked receipts, not an arithmetic existence theorem.

## 5. A genuine return to the same family

For maximal run $R$, write $t=(2n+1)/3^R$, a 3-unit. The child already satisfies $m_R\equiv3\pmod4$. It returns to this paid family exactly when
\[
3^4\mid2m_R+1=2^{R+2}t-1
\quad\Longleftrightarrow\quad
t\equiv2^{-(R+2)}\pmod{81}.                                 \tag{11}
\]
This is compatible with the modulo-8 condition by CRT. Among the 54 unit cofactor classes modulo 81 exactly one has the return. This is a parameter-frequency statement, not a probability about one supplied source.

A concrete all-height two-return ray is
\[
n=13891+26244s,\quad
m=10975+20736s,\quad
z=8671+16384s,\qquad s\ge0.                                  \tag{12}
\]
The first run is exactly maximal $R=4$; the second has at least four units of fuel, so choosing $r=4$ is legal but need not be maximal. Both edges are paid dependencies. Matching the actual middle prefix produces
\[
n\xrightarrow{(1)}(3n+1)/2
\xleftarrow{(1^4,3,1^3,3)}z,                                \tag{13}
\]
with raw clock differences $(-12,-8)$.

Every chosen legal run has $m+1<(64/81)(n+1)$. Since $2\cdot64^3<81^3$, three reductions more than halve the shifted source. Any ordinary iteration stops in fewer than $3\,\operatorname{bitlength}(n+1)$ reductions. The terminal can simply be outside (5); termination of this router is not ROOT. The chain constructor retains the original source word $(1)$ and aligns each middle prefix without new orbit discovery.

## 6. Exact symbolic ports from the inherited giant class

For $M_E=2^E-1$, $E\ge2$, the binary source condition is automatic. Elementary LTE gives
\[
v_3(2^{E+1}-1)=
\begin{cases}0,&E\text{ even},\\1+v_3(E+1),&E\text{ odd}.\end{cases}
\]
Thus the rule applies exactly for $E\equiv53\pmod{54}$, and emits
\[
m_{E,r}=\frac{2^{r+1}(2^{E+1}-1)}{3^r}-1,\qquad4\le r\le R.   \tag{14}
\]
This child is not Mersenne: its shifted odd cofactor is $5$ or $7\pmod8$, hence exceeds 1. Its leading size factor is
\[
\alpha_r=\frac{2^{r+2}}{3^r},\qquad
m_{E,r}=\alpha_rM_E+\alpha_r/2-1<\alpha_rM_E.
\]

The checker consumes the exact [fibre packet](collatz_fibre_integration_20261009.md) label set, not its interval hull. Of 1,949 labels, **36 are eligible**, with maximal-fuel histogram
\[
\{4:24,\quad5:8,\quad6:3,\quad7:1\}.
\]
Retaining all legal $r$ gives **53 pairs**; (7) normalizes them exactly to those 36 maximal ports. The least eligible exponent is 99,708,991,793. The unique $R=7$ source is 99,708,992,495, giving the best factor $512/2187<1/4$.

All 36 children nevertheless exceed the old minimum $M_{99,708,991,741}$. The least eligible exponent is 52 bits higher; exact comparisons use only the small exponent gaps. These are different non-Mersenne terminal ports, not a new minimum-height record.

To compute (14) modulo $q$, the reader computes its numerator modulo $3^r q$, divides by $3^r$, and subtracts 1. The full denominator is retained even when $3\mid q$. No giant source or ROOT word is materialized.

### Further independent paid guards

The inherited $G_1,G_5,G_{17}$ guards are disjoint. Five maximal ports have a successor in that bank:

| Original exponent | Selected inverse word on its maximal child |
|---:|---|
| 99,708,992,225 | $G_5$ |
| 99,708,992,711 | $G_5,G_5$ |
| 99,708,992,819 | $G_5$ |
| 99,708,993,197 | $G_5,G_1$ |
| 99,708,993,683 | $G_5$ |

The exact source residues authenticate every inverse word, and each resulting child is still smaller than the original Mersenne source. Prefixing its actual forward word to the child side of (2) preserves that original obligation.

All five final children are divisible by 3: they have no ordinary odd predecessor, no generator in this bank, and no run return. The other 31 maximal ports have no first generator. None of the 36 ports or their intermediate selected states satisfies (11). For all 53 initial run choices, the recorded route is exactly $G_1^{R-r}$ followed by the maximal port's route.

This is a finite exact stopping statement for the declared run/$G_1,G_5,G_{17}$ bank at these specified symbolic sources. It does not exclude other inverse words, sibling lifts, forward prefixes, or ROOT proofs. All giant terminal obligations remain **OPEN**.

## 7. API and verification scope

The ordinary join authenticates source, requested run, actual first-hit prefixes, and strict payment. Discharge checks a supplied child ROOT word and its common-future prefix before substitution. Symbolic ports and extensions retain exact exponents, native residues, full denominators, and ordered words. Boolean/float aliases and forged packets are rejected.

The declared controls are: 845 family instances $r=4,\ldots,16$, parameters 0 through 64; every positive odd source through 20,000 for (5); supplied ROOT(223); midpoint and independent-coordinate hostiles; 100 lifts of (12); every exponent 2 through 1,000 and literal comparisons for its 18 eligible Mersennes; all 1,949 inherited labels, 53 run choices and their successor guards; malformed and ROOT-boundary cases.

    python -B -X utf8 04-computation/experiments/collatz_projective_receipts_20261009.py
    python -B -O -X utf8 04-computation/experiments/collatz_projective_receipts_20261009.py

The all-height claims have the algebraic proofs above. Label/guard counts are the specified finite computation. Geometry, payment, or finite successor exhaustion does not discharge an unspecified ROOT obligation.

Normal and optimized runs agree exactly with the saved output: **25,595 checks**.
LF-normalized SHA256:
`74f2babe50523668532786d60a0b6c20390b9ab2bb92ec9b8e09d8f51ca0172b`.
The complete elementary proof and production APIs also passed two independent peer audits.
