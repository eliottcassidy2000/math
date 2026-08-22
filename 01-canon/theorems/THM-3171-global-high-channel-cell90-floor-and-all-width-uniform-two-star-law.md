---
id: THM-3171
title: "Global high-channel cell-90 floor and all-width uniform two-star law"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On every high primitive channel of the reflected cell-90 two-star for
  H=(1,2,3,4,6,12), every ordered label lane has overlap at least 1/105,
  except the unique label-{1,6}, primitive-3:5 lane, whose sharp floor is
  2030/280393.  The exceptional lane forces at least three other regular
  high edges.  Consequently the fixed uniform law on the nine two-star
  edges closes cell 90 for every six-distinct cap-two level assignment with
  minimum level at least six.  This is a fixed-cell sufficient certificate,
  not LRC(14).
audit: >
  The proof partitions the infinite channel universe into an exact
  level-10000 floor-moment run, a symbolic finite-state tail for Q<=26,
  two monotone centered-grid regions, and a Dirichlet d-block proof of the
  four residual primitive label-12 lanes.  Separate audits reconstruct
  literal arcs, the Euclidean moment recursion, the 11608 affine resonance
  channels, and the exceptional ratio graph.  Every Python gate uses
  require/RuntimeError rather than assert and has byte-identical normal and
  optimized output.
source: root/frontier-synthesis/lrc-uniform-tail/2026-08-02
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3150-cap-two-two-star-skeleton-and-common-dilation-tail
  - THM-3156-uniform-two-star-law-cap-two-boundary-through-width-100
related:
  - THM-3135-directed-cycle-weak-order-lane-cover-and-reflected-h-boundary
finite_source: 04-computation/lrc_high_channel_floor_referee_thm3171.cpp
finite_output: 05-knowledge/results/lrc_high_channel_floor_referee_thm3171.out
finite_source_sha256: 48b635746237b15c3a0ab2874b6d163c42b511019aaa8ae5aa3ba036ecaed1a6
finite_output_sha256: e530eb3ae6f825108463950f58ab8d09968f554be439d5c75af40925a802095a
small_tail_script: 04-computation/lrc_small_channel_symbolic_tail_thm3171.py
small_tail_output: 05-knowledge/results/lrc_small_channel_symbolic_tail_thm3171.out
small_tail_script_sha256: d73273a4cf4b88bea2890e001166d96cb07dd9b61f3a248ff1538ec44579796a
small_tail_output_sha256: 1e582a60215f97009b371faa2db8426d696a6feca45e2e40ce3594b60651db8f
label12_core_dependency: 04-computation/lrc_label12_multiblock_core_thm3171.py
label12_core_dependency_sha256: 7cff3a20e4d874f1143d096d3c9f4fa37448eb76c0ae0408513dd81780d8382b
label12_script: 04-computation/lrc_label12_dirichlet_referee_thm3171.py
label12_output: 05-knowledge/results/lrc_label12_dirichlet_referee_thm3171.out
label12_script_sha256: ae4669304e0d5a81ea6eaf75f385e2294fba6899ad72938c3ded8736e3ff5007
label12_output_sha256: 3e78cdbe97be0975532c6119425ad7e427eeed13138b2947a54b0eda780df43a
label12_independent_script: 04-computation/lrc_label12_dirichlet_independent_audit_thm3171.py
label12_independent_output: 05-knowledge/results/lrc_label12_dirichlet_independent_audit_thm3171.out
label12_independent_script_sha256: 89e7745b2e82f254421a4b9077a9c9a81812e58c0e14290f4ec1c506e660f234
label12_independent_output_sha256: ad64555eb75291129d0634c58624b863785f5bd942184e1748a32ec9b09ea7ad
floor_independent_script: 04-computation/lrc_high_channel_floor_independent_audit_thm3171.py
floor_independent_output: 05-knowledge/results/lrc_high_channel_floor_independent_audit_thm3171.out
floor_independent_script_sha256: 44095fdf32838930477332074a448c1321ddb6208b347ef7f7311b33bfb59b17
floor_independent_output_sha256: cb66d4c2bb444d628268928d62c1ef9a7c657879793b840731848e4b9ffb8bd4
graph_script: 04-computation/lrc_exceptional_edge_high_count_audit_thm3171.py
graph_output: 05-knowledge/results/lrc_exceptional_edge_high_count_audit_thm3171.out
graph_script_sha256: d5695add5d8b5a71210aafa826f00d4bbc0605a6735ed981724f66b168366835
graph_output_sha256: 10ff465c818fcf53898092677c1e1bce5f6ba0a9481a54519859c3fa4f75293a
hash_basis: LF-normalized bytes
---

# THM-3171 -- global high-channel cell-90 floor and all-width uniform two-star law

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**

## 1. Statement

Work in the sufficient reflected \(k=1\) family of THM-2941 and fix

\[
 H=(1,2,3,4,6,12),\qquad L=168,\qquad j=90.                 \tag{1}
\]

The pivots carry labels \(1,2\).  Their nine-edge two-star is

\[
 E_*=\{12,13,14,16,1\,12,23,24,26,2\,12\}.                 \tag{2}
\]

For distinct physical levels \(p<q\le2p\), put
\(g=\gcd(p,q)\), \(p=gP\), \(q=gQ\).  A channel is **high** when
\(P+Q\ge8\).  If label \(e\) occupies \(p\), label \(f\) occupies \(q\),
write \(I_{e,f}(p,q)\) for the exact cell-90 reflected overlap.  For every
ordered lane belonging to (2), every \(p\ge6\), and every high channel,

\[
 I_{e,f}(p,q)\ge
 \begin{cases}
  2030/280393,&\{e,f\}=\{1,6\},\ (P,Q)=(3,5),\\[2mm]
  1/105,&\text{otherwise}.
 \end{cases}                                                \tag{3}
\]

The first constant is sharp.  Equality occurs when label \(6\) occupies
level \(6\) and label \(1\) occupies level \(10\).

Now let \(q_1,\ldots,q_6\) be six distinct integers, let
\(D=\min_i q_i\ge6\), and assume \(\max_iq_i\le2D\).  Put the six labels of
\(H\) on these levels in any order.  Then the one fixed probability law
which gives weight \(1/9\) to every edge of \(E_*\) has expected overlap
strictly larger than the full singleton debt.  Hence some two-star edge
closes reflected cell 90 by Bonferroni.

This is an all-width theorem only for this fixed cell and this inherited
sufficient family.  It neither classifies physical survivors nor proves
\(LRC(14)\).

## 2. Exact located-overlap formula

Put

\[
 z=168p-e,\quad w=168q-f,\quad
 r=90e\bmod168,\quad s=90f\bmod168,\quad
 B=\frac{rw-sz}{168}.                                       \tag{4}
\]

The last quantity is an integer: reduce its numerator modulo \(168\).
Directly solving the reflected modular inequality gives the disjoint arcs

\[
 A_e(p)=\bigcup_{k=0}^{p-1}
 \left[\frac{r+168k-12}{z},\frac{r+168k+12}{z}\right].       \tag{5}
\]

For two intervals with cleared radii \(12w,12z\) and cleared centre
distance \(168d\), their overlap numerator is

\[
 [12(z+w)-168d]_+-[12(w-z)-168d]_+.                         \tag{6}
\]

For a fixed \(k\), only
\(\ell=\lfloor(B+kw)/z\rfloor\) or \(\ell+1\) can meet (5), and the
cap-two inequalities make at most one of them contribute.  Therefore

\[
 I_{e,f}(p,q)=\frac1{zw}\sum_{k=0}^{p-1}
 \left(
 [12(z+w)-168d_z(B+kw)]_+
 -
 [12(w-z)-168d_z(B+kw)]_+
 \right),                                                   \tag{7}
\]

where \(d_z(x)=\min(x\bmod z,z-(x\bmod z))\).  This proves the integer
digital-tent formula used by every referee.

Counts and first moments of residues below a threshold reduce to the three
floor moments

\[
 F_\nu(n;m,a,b)=\sum_{k=0}^{n-1}k^\nu
 \left\lfloor\frac{ak+b}{m}\right\rfloor,\qquad \nu=0,1,2.  \tag{8}
\]

Splitting \(a,b\) into quotient and remainder parts and transposing the
lattice region under the floor gives the Euclidean recursion
\((n,m,a,b)\mapsto(y,a\bmod m,m,m-b\bmod m+a\bmod m-1)\).
The modulus strictly decreases.  Thus (8), and hence (7), is evaluated
exactly in \(O(\log p)\) arithmetic steps.  The C++ referee contains this
recursion and uses signed 128-bit integers; its declared universe stays
well inside that range.

## 3. A centered-grid lower bound

The tent in (7), as a periodic function of its normalized centre
displacement, has height \(24/w\), mean \(24/(7z)\), and Lipschitz constant
\(168/w\).  Partition the \(gP\) orbit samples into \(g\) blocks of \(P\)
samples and compare each block with the exact \(P\)-grid of step \(Q/P\).

For any phase, layer-cake counting on a \(P\)-grid differs from interval
length times \(P\) by at most two points.  Integrating the layers gives

\[
 \sum_{r\bmod P} f(x+rQ/P)\ge \frac{24P}{7z}-\frac{48}{w}.  \tag{9}
\]

The actual and ideal steps differ by

\[
 \frac wz-\frac QP=\frac{Qe-Pf}{Pz}.                        \tag{10}
\]

Choose centered representatives modulo \(P\).  Their absolute values sum
to \(\lfloor P^2/4\rfloor\).  Applying the Lipschitz constant to (9), block
by block, proves

\[
 I_{e,f}(gP,gQ)\ge
 \frac{24gP}{7z}-\frac{48g}{w}
 -\frac{168g\,|Qe-Pf|\,\lfloor P^2/4\rfloor}{Pwz}.          \tag{11}
\]

For distinct \(e,f\in H\),

\[
 |Qe-Pf|\le\frac{23Q}{2}.                                  \tag{12}
\]

Since \(z\ge168gP-12\), \(w\ge168gQ-12\), and the first term of (11) is
strictly larger than \(1/49\), it follows that

\[
 I_{e,f}(gP,gQ)>
 \frac1{49}-E(P,Q,g),                                      \tag{13}
\]
\[
 E(P,Q,g)=
 \frac{48g}{168gQ-12}
 +\frac{483gPQ}{(168gP-12)(168gQ-12)}.                     \tag{14}
\]

Each factor \(x/(168gx-12)\) decreases in \(x\), and (14) also decreases
in \(g\).  Since \(P>Q/2\),

\[
\begin{aligned}
 Q\ge27,\ g\ge57:\quad
 E(P,Q,g)&\le E(14,27,57)=\frac{4190811}{385086712},\\
 \frac1{49}-E(14,27,57)-\frac1{105}
 &=\frac{447611}{283038733320}>0,                           \tag{15}\\
 Q\ge124,\ g\ge2:\quad
 E(P,Q,g)&\le E(63,124,2)=\frac{133019}{12238746},\\
 \frac1{49}-E(63,124,2)-\frac1{105}
 &=\frac{15667}{999497590}>0.                               \tag{16}
\end{aligned}
\]

Thus every lane in either region has mass \(>1/105\).

When neither label is \(12\), the sharper bound
\(|Qe-Pf|\le11Q/2\) replaces \(483\) in (14) by \(231\).
For \(g=1,Q>10000\), hence \(P\ge5001\), monotonicity gives

\[
 I_{e,f}(P,Q)-\frac1{105}>
 \frac1{49}
 -\frac{48}{168(10001)-12}
 -\frac{231(5001)(10001)}
 {(168(5001)-12)(168(10001)-12)}
 -\frac1{105}
 =\frac{102641018749}{38426702262480}>0.                    \tag{17}
\]

This closes every residual primitive lane not involving label \(12\).

## 4. The four primitive label-12 lanes

It remains to treat

\[
 (e,f)\in\{(1,12),(12,1),(2,12),(12,2)\},\quad
 g=1,\quad Q>10000.                                        \tag{18}
\]

Write \(h=Q-P\).  Dirichlet's pigeonhole argument applied to
\(h/P\) gives integers

\[
 1\le d\le8,\quad0\le a\le d,\quad
 c=dh-aP,\quad |c|\le P/9.                                 \tag{19}
\]

Here \(c\ne0\): otherwise \(\gcd(P,h)=1\) would imply \(P\mid d\), while
\(P\ge5001>d\).  Put \(D=d+a\).  Then

\[
 dQ-DP=c,\qquad
 d\frac{w-z}{z}-a=\frac{A}{z},\qquad
 A=168c+eD-df.                                              \tag{20}
\]

Group the first \(d\lfloor P/d\rfloor\) nonnegative terms of (7) into
\(d\)-blocks, discarding fewer than \(d\) remaining terms.  If \(G_d\) is
the sum of the \(d\) shifted tents, the block origins form a rotation of
step

\[
 \rho=|A|/z,\qquad n=\lfloor P/d\rfloor,\qquad T=n\rho.      \tag{21}
\]

The label perturbation in (20) satisfies
\(|eD-df|\le184\).  From (19), \(P\ge5001\), and \(z\ge168P-12\),

\[
 \rho\le
 \frac{(56/3)P+184}{168P-12}
 \le\frac{23384}{210039}<\frac9{80}.                        \tag{22}
\]

For a periodic piecewise-linear \(G\), the exact composite-trapezoid
identity on a lifted path of length \(T\) is its path integral divided by
\(\rho\), plus the endpoint half-difference, plus one Peano term for each
derivative jump.  Every jump at fractional cell position \(u\) contributes
\(\Delta G'\rho u(1-u)/2\), whose absolute value is at most
\(|\Delta G'|\rho/8\).

If \(T\ge5\), discard the nonnegative fractional-path integral and use
\(\lfloor T\rfloor/T\ge4/5\).  The complete blocks contain at least
\(P-7\) samples; the endpoint loss is at most \(96/w\); and the four kink
families in each of at most eight tents give the displayed Peano loss.
With \(R=9/80\), the resulting uniform bound is

\[
 \frac{24(P-7)}{1176P}\frac45
 -\frac{96}{168(P+1)-12}
 -\frac{84}{168(P+1)-12}(PR^2+8R).                         \tag{23}
\]

Its positive term increases with \(P\), its two loss terms decrease, and at
\(P=5001\) it equals

\[
 \frac{1073260228297}{109824296467200}
 =\frac1{105}+
 \frac{9104849219}{36608098822400}>\frac1{105}.             \tag{24}
\]

If \(T<5\), then \(n\ge(P-7)/d\) and (21) imply

\[
 |A|<\frac{840dP}{P-7}.
\]

For \(P\ge5001,d\le8\), integrality gives \(|A|\le6729\), and (20) gives

\[
 1\le|c|\le
 \left\lfloor\frac{6729+184}{168}\right\rfloor=41.           \tag{25}
\]

There are exactly 11,608 compatible tuples
\((d,a,c,P\bmod d,e,f)\).  Along each resulting affine ray
\(P=p_0+dn,\ Q=q_0+(d+a)n\), the canonical piecewise-linear core integrates
the limiting \(d\)-block tent exactly and proves

\[
 I_{e,f}(P,Q)\ge J_0-\frac KP.                              \tag{26}
\]

All constants are rational and \(K\ge0\), so it suffices to check the first
admissible \(P\) on each ray.  The least certificate is attained by

\[
 (d,a,c,p_0,q_0,e,f,A,n_0)=(1,0,41,1,42,1,12,6877,5000)
\]

and equals

\[
 \frac{1008404155769629793}{71368830483641671344}
 =\frac1{105}+
 \frac{1643505041531878901}{356844152418208356720}
 >\frac1{105}.                                               \tag{27}
\]

Equations (24) and (27) prove all four lanes in (18).

## 5. Exhaustion of the high-channel universe

The exact C++ floor-moment referee checks every
\(6\le p<q\le10000\), \(q\le2p\), with \(P+Q\ge8\): 24,989,160 level
pairs and 449,804,880 ordered edge checks.  It applies (3)'s weaker
threshold only to the 3:5 label-\(\{1,6\}\) lane and reports `PASS`.

The symbolic referee separately proves every \(Q\le26\) ray for all
admissible dilations.  On each residue class of \(g\) modulo
\(|Qe-Pf|\), all active floor endpoints eventually become affine; the
cleared margin is then an exact quadratic.  It checks every finite head,
all 126,892 stabilized residue tails, and the quadratic integer minimum.
Its maximum stabilization level is 594, and its unique sharp class is the
exception in (3).

It remains only to verify that no channel lies between these gates.  Suppose
\(Q\ge27\) and \(q=gQ>10000\).

- If \(g\ge57\), (15) applies.
- If \(2\le g\le56\) and \(Q\ge124\), (16) applies.
- If \(2\le g\le56\) and \(Q\le123\), then
  \(q\le56\cdot123=6888\), contradicting \(q>10000\).
- If \(g=1\), (17) handles lanes without label \(12\), while Section 4
  handles the four label-12 lanes.

Together with the finite and small-channel gates, these cases prove (3).

For the sharp orientation
\((e,p;f,q)=(6,3g;1,5g)\), direct summation of (7) gives

\[
 I_g=\frac{4284g^2-2520g+C_{g\bmod3}}
 {(504g-6)(840g-1)},\qquad(C_0,C_1,C_2)=(0,-336,84).         \tag{28}
\]

The three cleared consecutive differences are positive multiples of

\[
\begin{aligned}
1787436g^2+2073474g+17,\quad
1211238g^2+1314819g+139285,\quad
999558g^2+964791g-34808.
\end{aligned}
\]

They are positive for \(g\ge1\); hence (28) increases, and its admissible
minimum \(g=2\) is \(2030/280393\).  The independent floor audit checks the
literal arcs against (7), the Euclidean recursion against 50,000 brute
moment controls, the full high-channel universe through level 300, and
(28) through \(g=1000\).

## 6. Exceptional edge and the uniform law

A distinct cap-two channel is low exactly when its oriented ratio lies in

\[
 S=\{1/2,2/3,3/4,4/3,3/2,2\}.                              \tag{29}
\]

Without the exceptional edge, THM-3150's exact ratio graph gives at least
three regular high edges, so

\[
 \sum_{uv\in E_*}I_{uv}\ge3/105=1/35.                       \tag{30}
\]

Suppose the exceptional label-\(\{1,6\}\) edge occurs.  Normalize the
label-1 pivot level to one.  Its label-6 leaf is
\(x=3/5\) or \(5/3\).  Normalize the label-2 pivot level to \(t\).
A leaf \(y\) is low to both pivots only if \(y\in S\) and \(y/t\in S\),
so \(t\in S/S\).  Inside the cap,

\[
 S/S=\{1/2,9/16,2/3,3/4,8/9,1,9/8,4/3,3/2,16/9,2\}.        \tag{31}
\]

If \(t\notin S/S\), the pivot edge is high and each of four leaves has low
degree at most one; at least five of the nine edges are high.  If
\(t\in S/S\), distinctness removes \(t=1,x\).  The exact graph referee
enumerates the 12 feasible \((x,t)\) rows and every subset of at most three
additional positive-degree leaves; the minimum high-edge count is four.
Thus an exceptional edge always has at least three other regular high
edges.  By (3),

\[
 \sum_{uv\in E_*}I_{uv}\ge
 \frac3{105}+\frac{2030}{280393}>\frac1{35}.                \tag{32}
\]

Therefore (30) holds for every six-distinct cap-two assignment, and the
uniform nine-edge law has expected overlap at least \(1/315\).

For a label \(e\) at level \(q_e\), THM-2941 gives singleton excess

\[
 \delta_e(q_e)=\frac{e}{7(168q_e-e)}.                        \tag{33}
\]

For \(D=\min_eq_e\), put
\(B(D)=\sum_{e\in H}e/[7(168D-e)]\).  This decreases in \(D\), and

\[
 \frac1{315}-B(101)=
 \frac{94550632925919195101}{32173637491076453316945}>0.    \tag{34}
\]

Hence the uniform expected overlap exceeds the full debt
\(\sum_e\delta_e(q_e)\) for every \(D\ge101\).  THM-3156 proves the same
strict inequality for every \(6\le D\le100\).  Thus one edge has overlap
larger than the debt for every \(D\ge6\).  Pairwise Bonferroni gives

\[
 \mu\!\left(\bigcup_{e\in H}A_e(q_e)\right)
 \le\frac67+\sum_e\delta_e(q_e)-I_{uv}<\frac67.             \tag{35}
\]

Cell 90 is body-safe, so THM-2941's reflected completion criterion closes
this fixed cell.  This proves the theorem.

## 7. Reproduction and validity boundary

No binary artifact is canonicalized.  Rebuild and replay the finite gate in
a temporary directory:

~~~bash
LRC3171_BUILD="$(mktemp -d)"
clang++ -std=c++17 -O3 -pthread \
  04-computation/lrc_high_channel_floor_referee_thm3171.cpp \
  -o "$LRC3171_BUILD/referee"
"$LRC3171_BUILD/referee" 10000 12
~~~

Compare that transcript with
`05-knowledge/results/lrc_high_channel_floor_referee_thm3171.out`.

The five Python gates are:

~~~text
04-computation/lrc_small_channel_symbolic_tail_thm3171.py
04-computation/lrc_label12_dirichlet_referee_thm3171.py
04-computation/lrc_label12_dirichlet_independent_audit_thm3171.py
04-computation/lrc_high_channel_floor_independent_audit_thm3171.py
04-computation/lrc_exceptional_edge_high_count_audit_thm3171.py
~~~

Run each with both `python3` and `python3 -O`; each pair is
byte-identical and equals its declared output.  The Dirichlet referee
hash-pins the multiblock core.  Its independent audit hash-pins both that
core and the referee.  The floor audit hash-pins the C++ source.  The YAML
header pins every canonical source, dependency, and output under
LF-normalized SHA-256.

The theorem proves an all-width fixed-cell/cap-two sufficient certificate.
It does not show that every LRC(14) survivor enters this reflected family,
does not close other cells, and is not a proof of \(LRC(14)\).

**End of proof.**
