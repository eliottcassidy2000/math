# The equality sextic: exact rational braid corroboration

**Status: PROVED analytic transport + FINITE-EXACT, INDEPENDENTLY AUDITED.**
The [independent audit](planar_jc48_sep06_boundary_braid_audit.md) passes
the complete proof, exact source, and normal/optimized frozen replays. This is a second, separate route
to the concrete exclusion already established by the
[marked infinity-plumbing proof](planar_jc48_sep06_boundary_plumbing.md).
It does not compute the whole affine complement group or claim a general
Jacobian-conjecture result.

## 1. Actual polynomial, inherited target, and limited consumer

Let C be the image of
\[
 U(t)=t^4+t^3+t^2,\qquad V(t)=2t^6+3t^5+2t^3+2t^2.
\]
The [resolution supplier](planar_jc48_sep06_resolution_budget.md) proves
that C is an irreducible sextic with its actual birational normalization,
one finite (2,5) cusp, four ordinary nodes and an infinity (2,9) branch.
Its resolved Nori budget is equality, \(D^2=2N=8\), so the strict
Nori criterion does not decide this curve.

The exact monic equation, of degree four in the vertical coordinate v, is
\[
\begin{aligned}
F(u,v)={}&16u^6-47u^5-60u^4v-218u^4-8u^3v^2-53u^3v-112u^3\\
&+68u^2v^2+16u^2v+84u^2+15uv^3+26uv^2-84uv\\
&+v^4-3v^3+21v^2.
\end{aligned}\tag{1}
\]
This equals \(\operatorname{Res}_t(U(t)-u,V(t)-v)\) exactly. Its vertical
discriminant is
\[
-u^5(256u^2+11u+12)
 (1225u^4+4639u^3+14945u^2+26313u+15435)^2.
\tag{2}
\]
Both identities can be reconstructed by ordinary symbolic resultant and
discriminant operations, independently of the path certificate. In
particular, no abstract quartic model is declared to be the actual curve.

The inherited [odd-cusp passport](planar_jc48_sep06_odd_cusp.md) would
force a hypothetical Keller map with whole nonproperness support C to have
six-sheet transitive monodromy, with every positive smooth-curve meridian
mapped to a single three-cycle. Each vertical-fibre meridian has this
same cycle type. The present finite certificate proves that no such
transitive representation satisfies four necessary actual braid
relations. It does not need the stronger inherited assertion that the
image is all of \(A_6\).

The small concept board retains five different coordinates:

| Concept | Preserved predicate | Missing datum supplied here |
|---|---|---|
| Semilocal cusp/node passport | Six sheets and marked meridian cycle type | Relations along actual global paths |
| Nori equality boundary | Actual resolved curve and numerical budget | A different global obstruction |
| Numerical root tracking | A candidate labelled path | Exact disjoint-root tubes |
| Exact rational braid | Actual homotopy in configuration space | Marked fibre-to-affine transport |
| Finite permutation bank | All raw three-cycle assignments | Necessity, not sufficiency, for a Keller map |

The closest hostile is a nonabelian marked \(A_5\) representation of the
infinity plumbing group: even a valid group obstruction must not be
silently upgraded to abelianity. The corrected near miss is an uncertified
numerical braid. The least-used sidecar is the monic-polynomial section:
it turns the braid relation into exact equality of marked tuples, rather
than equality only up to a free simultaneous conjugation.

## 2. A rational Rouché tube certifies actual paths

All verification arithmetic is in \(\mathbb Q(i)\). For a Gaussian
rational \(z\), write
\[
 L(z)=\max(|\Re z|,|\Im z|),\qquad
 U(z)=|\Re z|+|\Im z|.
\]
Thus \(L(z)\le |z|\le U(z)\).

Consider a straight base segment from \(u_0\) to \(u_1\), and four
candidate root centres \(z_{0i}\) at its start. Define
\[
 r_i={1\over32}\min_{j\ne i}L(z_{0i}-z_{0j}).
\tag{3}
\]
The open disks \(D(z_{0i},r_i)\) are pairwise disjoint. Expand exactly
\[
 F(u_0+h,z_{0i}+z)=\sum_{j=0}^6\sum_{k=0}^4B_{jk}h^jz^k.
\]
Put \(d=U(u_1-u_0)\). The verifier requires
\[
\begin{split}
L(B_{01})r_i >{}& U(B_{00})+
 \sum_{k=2}^4U(B_{0k})r_i^k\\
&+\sum_{j=1}^6\sum_{k=0}^4U(B_{jk})d^jr_i^k.
\end{split}\tag{4}
\]
For every point of the base segment, Rouché's theorem compares the
polynomial on \(|z|=r_i\) with \(B_{01}z\). It gives exactly one
root in each of the four disjoint disks. Since (1) is monic of degree
four, these are all roots, each simple. No unaccounted vertical root or
singularity can occur along an accepted base segment.

Endpoint matching is also paid. Let \(z_{1i}\) be the proposed endpoint
centres and let \(r'_{i}\) be (3) at that endpoint. The verifier requires
\[
 U(z_{1i}-z_{0i})<r_i/2
\]
and checks (4), with no base variation, in the smaller disk
\[
 D(z_{1i},\epsilon_i),\qquad
 \epsilon_i=\min(r_i,r'_i)/64.
\tag{5}
\]
The unique endpoint root therefore belongs to both the preceding large
disk and the next step's large disk. This is a labelled transport
certificate, not four unrelated counts of roots at successive samples.
Initial root isolation is checked separately in radius \(r_i/64\).

On each accepted segment the actual roots and the straight segments of
candidate centres remain in the same four convex disjoint disks. The
small endpoint displacements in (5) agree across adjacent segments, so
linear interpolation within these disks gives a consistent homotopy
between the actual configuration path and the rational polygonal one.
The initial configuration is common to all loops. Each final unordered
centre set is exactly the initial centre set. Thus the endpoint
displacement homotopies close with the actual endpoint permutation, and
all four loops use the same identification of their base fibre.

This argument preserves the full actual path in configuration space.
It does not infer any braid from approximate discriminant locations.

## 3. Four verified loops and their exact words

All loops are based at \(u_*=2+i\). The frozen witness follows rational
subdivisions of four diamonds and their common type of access path. Their
centres and radii are
\[
\begin{array}{c|c|c}
j&c_j&r_j\\\hline
0&0&1/16\\
1&-22528/2^{20}+i\,225902/2^{20}&1/32\\
2&-22528/2^{20}-i\,225902/2^{20}&1/32\\
3&-612673/2^{20}+i\,2796037/2^{20}&1/16.
\end{array}
\]
For each row the polygon vertices are
\[
 u_*,\ c_j+r_j,\ c_j+ir_j,\ c_j-r_j,
 c_j-ir_j,\ c_j+r_j,\ u_*.
\]
The numbers of accepted straight segments are respectively
\(15987,1330,1824,2362\), totalling **21,503**.
The centres were originally proposed by floating root calculations and
then rounded to denominator \(2^{38}\). Their authority comes entirely
from (3)--(5); the default verifier neither imports NumPy nor uses floating
arithmetic. The generating mode is an optional witness finder.

Multiply vertical coordinates by the oriented rational linear factor
\(1+i/4\). For each polygon segment, real-coordinate crossings occur at
rational times. The verifier checks that endpoint real coordinates are
distinct, crossing times are distinct, the crossing strands are adjacent
in current real order, and their imaginary separation is nonzero. The
positive letter denotes a counterclockwise half-turn: the strand starting
to the left is below the other at the crossing. Only adjacent inverse
letters are cancelled. The resulting chronological words are

```text
W0 = [2,3,2,-1,-2,1,1,1,1,2,1,-3,-2]
W1 = [2,3,2,-1,2,1,-2,-3,-2]
W2 = [2,3,-2,-1,-3,-2,-3,2,3,-1,-3,1,3,2,3,1,2,-3,-2]
W3 = [2,3,2,1,2,2,-1,-2,-3,-2].
```

It is unnecessary to prove that a diamond encloses exactly one critical
value, or that these four loops generate the entire regular-base group.
Each separately certified closed loop gives a necessary relation. This
removes an otherwise unnecessary algebraic-root-location obligation.

## 4. Why actual marked tuples satisfy these relations

Write \(B=\mathbb C\setminus\{\operatorname{disc}_vF=0\}\). Over B,
(1) gives the usual configuration-space pullback bundle with fibre a
plane punctured at four points. A fibre has free fundamental group
\(F_4\), generated by four positive meridians.

There is an actual section outside all vertical roots. For example, the
coefficients of the monic quartic give a continuous Cauchy root bound
\(R(u)>\max\{|v|:F(u,v)=0\}\); the real point \(v=R(u)\) defines a
continuous section over the entire u-plane, including critical values.
Over the compact region filling any of our four loops, one can instead
use a single sufficiently large constant \(v_*\). Its base loop is
null-homotopic in the actual affine complement. Consequently the moving
fibre meridians return to the same elements of that complement group,
without an undetermined simultaneous conjugation.

The vertical \(F_4\) surjects onto the affine complement group. Indeed,
over B the fundamental group is generated by fibre generators and lifts
of the punctured-base generators; the section supplies such lifts.
Restoring the finitely many critical fibres kills those section lifts.
Every loop in the full complement can be perturbed to avoid the finitely
many critical vertical fibres, so no additional generator is missed.
This argument uses the monic equation and a finite discriminant, and does
not require a generic transverse projection at projective infinity.

In a geometric meridian basis the elementary Hurwitz coordinate change
is, in the convention used by the source,
\[
 H_i^+(\ldots,a,b,\ldots)
   =(\ldots,aba^{-1},a,\ldots),\qquad
 H_i^-=(H_i^+)^{-1}:
 (a,b)\longmapsto(b,b^{-1}ab).
\tag{6}
\]
Here is an explicit convention check. Take meridian access stems from
below the ordered punctures. A counterclockwise half-twist transports the
old left loop to the final right loop, and the old right loop to
\(x_{\rm right}^{-1}x_{\rm left}x_{\rm right}\). Thus its map on loops
is the inverse-Artin map. The representation tuple is transported
contravariantly, by composing with the inverse loop map, which is exactly
\(H_i^+\) in (6). Successive crossings therefore apply (6) in chronological
order, as the source does. This separates geometric loop transport from
representation-coordinate transport.

The finite control below additionally checks global word reversal,
global sign reversal, and both operations. These are supplemental
convention controls; the declared below-stem convention already gives
the actual marked action used in the proof.

The actual representation must therefore supply four single three-cycles
which are fixed by the four marked actions, using one consistent
convention. This is an exact tuple condition. Merely asking whether each
word preserves a simultaneous-conjugacy class would lose the section
information and would be a weaker test.

## 5. Complete finite consumer and its scope

There are forty single three-cycles on six labels. Simultaneous relabeling
can send the first one to \((012)\), so it suffices to retain **all
\(40^3=64,000\)** assignments to the remaining three meridians. This
normalization loses no representation or transitivity information. The
source computes (6) with literal permutations on all six labels.
Successive survivor counts after imposing \(W_0,W_1,W_2,W_3\) are
\[
 30,400,\quad760,\quad19,\quad1.
\tag{7}
\]
The remaining tuple consists of four copies of the first three-cycle.
It is a positive consistency control, but its generated action fixes
three labels. There is **no transitive tuple**.

A separate cheap replay of the same literal consumer with every word
reversed, with every sign reversed, and with both operations performed
returns the identical counts (7) and the identical constant survivor.
This tests the convention boundary without modifying the frozen rational
path source or calling a numerical braid computation. A reproduction is

```python
import gzip, importlib.util, json
p = '04-computation/planar_jc48_sep06_boundary_braid.py'
s = importlib.util.spec_from_file_location('braid', p)
m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
w = [x['word'] for x in json.loads(gzip.decompress(m.CERT.read_bytes()))['loops']]
for ww in (w, [x[::-1] for x in w], [[-i for i in x] for x in w],
           [[-i for i in x[::-1]] for x in w]):
    counts, rows, transitive = m.group_census(ww)
    if counts != [30400,760,19,1] or rows != [(0,0,0,0)] or transitive:
        raise RuntimeError('marked convention control failed')
```

Since the four fibre meridians generate the affine group, a transitive
six-sheet representation cannot have this sole surviving tuple.
Together with the actual passport, the certificate excludes C as the
whole nonproperness support of a Keller map. This is a **finite-exact
corroborating route**, with an analytic transport proof. It neither
presents the full affine group nor proves that it is cyclic: other
meridian cycle types and other finite targets were not enumerated.

## 6. Frozen witness, reproduction and independent audit boundary

The frozen files are

* [exact verifier](../../04-computation/planar_jc48_sep06_boundary_braid.py);
* [compressed rational path witness](planar_jc48_sep06_boundary_braid_certificate.json.gz);
* [normal replay output](planar_jc48_sep06_boundary_braid.out).

The compressed witness has 1,393,532 bytes, decompressing to 5,621,433
bytes of JSON. Its size records actual accepted paths, not an implicit
claim about unverified intermediate samples. The verifier accepts every
segment only by the explicit strict inequalities above.

```text
python3 -B 04-computation/planar_jc48_sep06_boundary_braid.py
python3 -B -O 04-computation/planar_jc48_sep06_boundary_braid.py
source SHA256 459fc3a3920f2335331ee69704accb170a3b3c5ca31f6903c009672faf6d6046
compressed witness SHA256 a06d7434da196de9a0f721fd913362dac7fd678358673c8d929d7f2c2a2c6951
raw JSON SHA256 5e4c61c9a64e1a30ab841de524be966583ee2a29929c1544aea6f4e2f4675868
output SHA256 295b812bdfb38a6b0b2a80dddf3fab26de80d3ed00a6d30414f1cc05951a3b3d
```

The completed normal replay has **172,168 always-active gates**, with a
635-byte output. The independent referee has reproduced both normal and
optimized outputs byte-for-byte and passed the complete final text. The
audit hash is `598fdccb4c893959e45d39f9e801465e3d5da92389eeed750a8453ac01bc38a5`. The frozen source deliberately retains its candidate header and
prints that the analytic transport proof is required separately; neither
code execution nor the original numerical witness finder replaces that
proof.
