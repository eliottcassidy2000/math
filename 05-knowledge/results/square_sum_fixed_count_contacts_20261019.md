# Fixed-count contact atlases and polynomial-height finite gluing

**Status: PROVED**, with targeted **FINITE-EXACT** controls and an
independent proof audit passed. For every fixed positive number of requested neighbors,
the two bounded-section atlases extend to a finite atlas producing that many
distinct neighbors. The constants may depend on this fixed number. No bound
polynomial in the number of requested neighbors, novelty claim, or optimized
exponent is asserted.

## 1. The fixed-count theorem

A primitive positive Pythagorean triple `T=(a,b,C)` is **admissible** when
`C²>4ab`. Its spread is `f_T=2ab/C² in (0,1/2)`, and exchanged legs represent
the same primitive petal.

**Theorem 1.** For every fixed integer `r>=1`, effective constants
`A_r,B_r>0` and a deterministic exact algorithm exist with the following
property. Given any admissible primitive `T` and either prescribed port
`epsilon in {-1,+1}`, the algorithm returns `r` admissible primitive petals
`U_1,...,U_r`, pairwise distinct and all different from `T`, and positive
rational certificates `rho_i` such that

\[
 \rho_i^2=1+\varepsilon f_T-\varepsilon f_{U_i},
 \qquad C_i\le B_r C^{A_r}\quad(1\le i\le r).             \tag{1}
\]

For fixed `r`, it uses a bounded number of rational field operations and
comparisons on `O_r(log C)`-bit operands. Output and certificate lengths
are `O_r(log C)`, with polynomial bit time in `log C`. Here and below the
constants may grow arbitrarily with `r`.

We inherit, rather than repeat, the proofs of:

- [Universal signed Pythagorean contacts](square_sum_ppt_bridges_20261019.md):
  the quartic/cubic maps, non-torsion of the translating point, unique
  primitive petal for a spread, and common-center graph transfer.
- [Exact rational recurrence budget](square_sum_contact_recurrence_bound_20261019.md):
  the invariant differential, monotone admissible arcs, circle pigeonhole
  lemma, and rational period upper bound.
- [Plus-port bounded atlas](square_sum_bounded_section_contacts_20261019.md)
  and [minus-port bounded atlas](square_sum_minus_bounded_section_contacts_20261019.md):
  low-cone normalization, the cusp sections for `r=1`, ordinary exceptional
  cases, and the bounded-field-operation height argument.

## 2. All fixed odd cusp sections have distinct first-order displacement

Use the inherited normalization

\[
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad
 1<t<\cot(5\pi/24)<4/3.                                \tag{2}
\]

For reduced `t=n/m`, one has `n²+m²=delta C`, with `delta in {1,2}`, and
`H(t)<=sqrt(2C)`. The rational leg swap `t -> (t+1)/(t-1)` chooses this
low representative. The function `f` is strictly increasing on `(1,5/4)`,
and `f(5/4)=720/1681<1/2`.

For the two port families put

\[
 k_+(t)=\frac{t^2+2t-1}{t^2+1},\qquad
 k_-(t)=\frac{-t^2+2t+1}{t^2+1},\qquad K=k_\varepsilon^2,
\]

\[
 E_t:Y^2=(X-K)(X^2-1),\qquad P=(0,k_\varepsilon).
                                                               \tag{3}
\]

For plus, let `Q` be the image of `(u,z)=(t,t²+1)`; for minus, let `Q_-`
be the image of the **negative-sheet** point `(-t,-(t²+1))`. For every odd
`d=2j+1>=3`, define

\[
 R_d^+=dQ+(d-1)P=Q+jS_+,\qquad S_+=2(Q+P),
\]
\[
 R_d^-=dQ_-=Q_-+jS_-,\qquad S_-=2Q_-.                 \tag{4}
\]

Let `(u_d,z_d)` be their inverse quartic coordinates, and set
`v_d^+=u_d`, `v_d^-=-u_d`. These are rational functions of `t`.

**Lemma 2.** For each fixed odd `d>=3`, these rational functions are regular
at `t=1` and, writing `h=t-1`,

\[
 v_d^\varepsilon(1+h)=1+dh+O(h^2),\qquad
 z_d(1+h)=2\varepsilon+O(h).                           \tag{5}
\]

**Proof.** The cubic fiber at `h=0` is nodal, but the identity and all points
used in this local calculation are in its smooth locus. Set
`B_+=Q+P` and `B_-=Q_-`. Their exact coordinates from the two predecessor
notes have the common first-order behavior

\[
 X(B_\varepsilon)=-1+O(h^2),\qquad
 Y(B_\varepsilon)=-2h+O(h^2),\qquad K=1+O(h).
\]

Indeed both X-coordinates equal `(t²-4t+1)/(t²+1)`; their Y-coordinates
are respectively `-8t²(t-1)/(t²+1)²` and `-8t(t-1)/(t²+1)²`.
In the doubling formula the slope is

\[
 \lambda=\frac{3X^2-2KX-1}{2Y}=-h^{-1}+O(1).
\]

It follows that

\[
 X(S_\varepsilon)=h^{-2}+O(h^{-1}),\qquad
 Y(S_\varepsilon)=h^{-3}+O(h^{-2}),\qquad
 w(S_\varepsilon):=-X/Y=-h+O(h^2).                    \tag{6}
\]

The formal parameter `w=-X/Y` at the smooth identity satisfies
`[j]w=jw+O(w²)`. This formal group law is defined over `Q[[h]]` even
though the special fiber is nodal: it uses only the smooth identity chart.
Thus `w([j]S_epsilon)=-jh+O(h²)`.

The invariant differential satisfies both

\[
 \omega=\frac{dX}{2Y}=-\frac{du}{z},\qquad
 \omega=(1+O(w))\,dw \quad\hbox{near the identity}.     \tag{7}
\]

Here the differentials are relative to the fixed fiber. Translation
preserves `omega`, so, at either smooth anchor `Q_epsilon`,

\[
 u(Q_\varepsilon+A)
   =u(Q_\varepsilon)-z(Q_\varepsilon)w(A)+O(w(A)^2).
                                                               \tag{8}
\]

The anchor values are `u=epsilon*t` and `z=epsilon*(t²+1)` for these
chosen sheets. Substituting (6) into (8) gives
`u(R_d)=epsilon*(1+(2j+1)h)+O(h²)`, proving (5). Translation by a formal
identity point preserves the smooth anchor chart. At the special fiber
the plus anchor is `(X,Y)=(3,4)` and the minus anchor is `(-1,0)`; both
have `X(X-K)!=0`. Thus the ordinary quartic inverse is also regular,
including its z-coordinate. QED.

This argument proves (5) for every fixed odd `d`, rather than inferring it
from the first several computations. The formal-law and smooth-reduction
facts are the same ones already used and audited in Section 6 of the
plus-port atlas; no new rank or specialization claim is needed.

### An effective common cusp cutoff

For the finite list `d=3,5,...,2r+1`, Lemma 2 gives an effective rational
`eta_r in (0,1/16]`, common to both ports, such that throughout
`0<h<=eta_r`,

\[
 1<t<v_3^\varepsilon(t)<v_5^\varepsilon(t)<\cdots
     <v_{2r+1}^\varepsilon(t)<5/4.                     \tag{9}
\]

All inverse-coordinate denominators are nonzero there. Here is an exact
way to obtain the cutoff, so no unspecified asymptotic remainder is used
as an algorithm. Compute and reduce the finite rational functions

- `(v_3-t)/h`;
- `(v_{d+2}-v_d)/h` for consecutive d;
- `5/4-v_{2r+1}`;
- the denominators of every `v_d` and `z_d`.

Each of the first three types has positive value at `h=0`, and each
denominator has nonzero value there. For every resulting numerator or
denominator polynomial `p(h)`, change its sign if necessary so `p(0)>0`.
With `M_p=sum_{i>=1}|p_i|`, the choice

\[
 0<\eta\le\min\left(1,\frac{p(0)}{2M_p}\right)         \tag{10}
\]

makes `p(h)>=p(0)/2>0`; omit the second bound when `M_p=0`.
Take the minimum of these finitely many rational bounds and `1/16`.
For a positive rational function, its numerator and denominator have the
same sign at zero, so this procedure proves its strict positivity too.

By monotonicity of `f`, (9) gives `r` different admissible spreads, all
strictly greater than the anchor spread. The quartic equation supplies

\[
 \left(\frac{z_d}{v_d^2+1}\right)^2
       =1+\varepsilon f(t)-\varepsilon f(v_d).          \tag{11}
\]

The radicands lie in `(1/2,3/2)`, so taking the absolute rational roots
gives the required positive certificates.

## 3. Fixed-count recurrence on the complement

Set `delta_r=f(1+eta_r)>0`. On the rest of the low cone, `f_0>delta_r`.
For **both** ports use the positive-z anchor `(epsilon*t,t²+1)` from the
recurrence theorem; this sheet differs from the minus cusp sheet.

The inherited invariant lengths to the two spread boundaries are

\[
 L_- =\frac14\int_0^{f_0}
 \frac{ds}{\sqrt{1-s^2}\sqrt{1+\varepsilon(f_0-s)}},\quad
 L_+ =\frac14\int_{f_0}^{1/2}
 \frac{ds}{\sqrt{1-s^2}\sqrt{1+\varepsilon(f_0-s)}}.
                                                               \tag{12}
\]

Choose the direction with greater **spread width**. Its width is at least
`1/4`, and the displayed integrand is strictly greater than `1/5`.
Therefore its open, strictly spread-monotone arc has length `ell>1/20`.
As inherited, `ell<Omega/2`, where `Omega` is the real identity-component
period. The translating point `A=2P` is non-torsion and lies in that
component.

**Lemma 3 (r-point signed pigeonhole).** Such an arc contains
`Q+sigma*n*A,...,Q+sigma*r*n*A`, for a sign `sigma` and an integer

\[
 1\le n\le\left\lfloor\frac{r\Omega}{\ell}\right\rfloor.
                                                               \tag{13}
\]

These points have pairwise distinct spreads, all different from the anchor.

**Proof.** Apply the inherited circular-pigeonhole lemma to an arc of
length `ell/r` beginning at Q, with its chosen orientation. It supplies a
nonzero signed displacement `alpha` with `0<alpha<ell/r` and the bound
(13). The first r multiples have displacements `alpha,...,r*alpha`, all
strictly between zero and ell. There is no wraparound since `ell<Omega/2`.
Strict spread monotonicity proves the final assertion. QED.

For a completely rational uniform bound, choose `J_r>=0` with
`2^{J_r}delta_r>=5/2` and put

\[
 W_r=\frac57(4+J_r),\qquad
 N_r=\lfloor20rW_r\rfloor,\qquad M_r=2rN_r.             \tag{14}
\]

The inherited period bound gives `Omega<W_r` for every anchor in this
complement, for either port (`D<=5/2`). Equations (12)--(13) then give
`n<=N_r`. Consequently the finite candidate list

\[
 Q+qP,\qquad q\in2\mathbb Z,\quad0<|q|\le M_r           \tag{15}
\]

contains at least r successful candidates with distinct spreads. Its
indices include `q=2*sigma*j*n` for `1<=j<=r`.

Scan (15) with exact rational group operations and the inherited ordinary
identity/inverse/doubling branches. Skip unusable inverse charts, and retain
only positive-z inverse points with `v=epsilon*u>1`, `0<f(v)<1/2`, spread
different from the anchor and from all retained spreads. The guaranteed
monotone-arc points are in ordinary charts, so the scan must obtain r.
There is no need to compute the real period, identify n, or test an
irrational interval endpoint. Nor must the actual returned list be monotone.

## 4. Primitive normalization, height, and forbidden sets

For any accepted reduced parameter `v=n'/m'>1`, divide
`(n'²-m'²,2n'm',n'²+m'²)` by 2 when both components are odd, otherwise by 1,
and reorder its legs. This preserves spread and produces a primitive
positive PPT. The inherited uniqueness-by-spread lemma converts all the
strict spread distinctions into distinct unordered primitive petals.

For fixed r, Sections 2--3 form a finite rational-section atlas with a
bounded number of field operations and exceptional cases. Its initial
parameter has `O(log C)` bits. Applying the predecessor's elementary
height inequalities

\[
 H(xy),H(x/y)\le H(x)H(y),\qquad H(x+y)\le2H(x)H(y)
\]

through this bounded computation gives effective (r-dependent) polynomial
numerical heights for all coordinates, normalized hypotenuses, and square
certificates. Equivalently, homogenize the finitely many rational functions
and use `H(t)<=sqrt(2C)`. Taking maxima over both ports proves (1).
The conservative schoolbook bound `O_r((1+log C)^3)` time and
`O_r(1+log C)` working bits carries over. This is a fixed-r statement;
it supplies no uniform polynomial bound in r.

**Corollary 4 (height-independent forbidden-set avoidance).** Given fewer
than r forbidden primitive petals of any heights, at least one of the r
outputs in Theorem 1 lies outside that set. Its hypotenuse still satisfies
`C'<=B_r C^{A_r}`, independently of the forbidden petals' heights.

Generate r outputs, compare primitive petals or exact spreads, and filter.
The extra reading/comparison cost depends on the encoded forbidden list's
bit length; it is not claimed to depend only on `log C`.

## 5. Fixed-length AP-paired fans and odd square-sum cycles

**Corollary 5.** Fix `k>=1`. From any admissible primitive seed T of
hypotenuse C, one can construct a chain of k distinct petals beginning at T
with contacts `(T_i,+)~(T_{i+1},-)`, whose common-center integer realization
has polynomially bounded numerical labels in C. More generally, from an AP
seed already at integer scale g, the realization contains an integer-square
scaled copy of that seed and has polynomially bounded labels in its size.

**Proof.** For `k>=2` use Theorem 1 with `r=k` at each last petal's plus
port, forbidding the earlier petals. There are always fewer than k forbidden
petals, so Corollary 4 supplies the next distinct petal. The number of steps
is fixed. Absorb `B_k` into a larger effective integer exponent `a_k>=2`,
so every such extension has `C_{i+1}<=C_i^{a_k}`. Then

\[
 C_i\le C^{a_k^{i-1}},\qquad
 S_k=\sum_{j=0}^{k-1}a_k^j,\qquad
 H=g\operatorname{lcm}(C_1,\ldots,C_k)\le gC^{S_k}.
                                                               \tag{16}
\]

For `k=1` take `S_1=1`. The inherited graph transfer gives center `M=2H²`
and outer labels `v_i^epsilon=2H²(1+2epsilon*f_i)`, each strictly below
`4H²<=4g² C^{2S_k}`. Their explicit square-sum cycle is

\[
 M,v_1^-,v_1^+,v_2^-,v_2^+,\ldots,v_k^-,v_k^+,M.        \tag{17}
\]

All labels are positive and distinct. Each outer pair with M is an AP
square-sum triangle. The center is adjacent to every outer vertex, so this
also realizes the AP-paired fan `K_1 join P_{2k}` as a subgraph.

The original seed labels are
`2g²(C²-4ab,C²,C²+4ab)`. Their new copy is multiplied by the integer square
`(H/(gC))²`. If N is their largest original label, then
`g²C²<=N/2`, so all resulting labels are at most `4(N/2)^{S_k}`.
The original center and labels need not remain fixed. QED.

The same fixed-size argument also realizes all prescribed edges of any
fixed rooted tree of petals (additional contacts may occur) when each
parent-child edge specifies opposite ports: choose a new petal
from the specified parent port and forbid all earlier petals. Only one
parent constraint is imposed at each insertion. This observation does not
supply simultaneous contact constraints to multiple existing petals.

No prescribed global spread direction, increasing-gap property, graceful
rank labeling, fixed-center extension, interval-label coverage, or
Hamiltonicity of the whole interval square-sum graph follows here. The
polynomial bounds concern fixed-length sparse constructions only.

## 6. Targeted exact controls

The companion checker derives `R_d` by rational-function group operations
for `d=3,5,7` on both ports, checks the ordinary inverse quartic and square
identities, and verifies the formal identity parameter's leading coefficient
and `v_d'(1)=d`. It also tests exact rational cusp probes and compact
complement probes returning three pairwise different neighbors. These are
implementation controls, not the proof for every fixed r; Lemmas 2--3
supply the uniform argument. The predecessors' complete censuses are not
repeated. Normal and optimized runs have byte-identical output, passing
242 exact rational and 38 exact symbolic checks. Eighteen cusp outputs
cover three rational inputs on both ports with d=3,5,7; eight compact
port/input probes each return three distinct neighbors.

The exact d=5,7 Taylor controls, beyond Lemma 2's required first derivative,
are

\[
 v_5^+(1+h)=1+5h+O(h^4),\qquad
 v_7^+(1+h)=1+7h-14h^3+O(h^4),
\]
\[
 v_5^-(1+h)=1+5h+20h^2+80h^3+O(h^4),\qquad
 v_7^-(1+h)=1+7h+42h^2+238h^3+O(h^4).
\]

These finite controls are not used to infer any growth law in d.

Run, with SymPy installed:

```bash
python3 04-computation/experiments/square_sum_fixed_count_contacts_20261019.py
python3 -O 04-computation/experiments/square_sum_fixed_count_contacts_20261019.py
```

The connection contract is: an admitted primitive seed and prescribed
mixed-sign port; r pairwise distinct new admissible opposite-port neighbors;
a finite list of cusp sections plus bounded signed even translates; exact
rational square certificates and polynomial numerical height. Necessary
sidecars are fixed r, the effective common cusp cutoff, the minus-port sheet
switch, duplicate-spread filtering, primitive normalization, and the
r-dependent bounds. Forbidden-set height independence does not erase the
cost of reading that set.
