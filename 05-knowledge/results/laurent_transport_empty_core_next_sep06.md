# Signed transport: two exact boundaries of the mixed-product reduction

Status: **EXACT REFUTATIONS / PROVED ANALYTICALLY + FINITE-EXACT +
INDEPENDENTLY AUDITED** for the
two stronger abstract implications below. The actual all-height A2B3
beta-skip sign and actual-versus-virtual signed transport remain **OPEN**.
These are auxiliary-carrier counterexamples, not actual Laurent-support
counterexamples. No theorem ID or external priority is claimed.

The new obstruction is precise: the two incoming differential identities,
even together with all three full auxiliary polynomials having simple
negative roots, do not imply the remaining skip sign. A separate minimal
degree-four example shows that generic strict interlacing does not supply
the missing mixed-product inequality either. The first counterexample
nevertheless retains the stronger full-response inequality `Q_model<V<0`:
a positive skip term can be absorbed by the already-proved negative term.

## 1. Inheritance and the next boundary

This continuation read the incoming notes at shared HEAD `6dd59c9c41`:

- [Midpoint transport](overnight7_20260906_laurent_midpoint_transport.md)
  retains the full negative auxiliary support and decomposes the beta
  midpoint into one hit class and two skip classes.
- [Carried quartic family](overnight7_20260906_laurent_quartic_carry.md)
  closes the specific `(-27,2g-27,3g-27)` family. No fixed-degree family is
  rederived here; the preceding cubic family was also already closed.
- [Alpha completion](overnight8_20260906_alpha_completion.md) proves the
  all-A2 beta-hit term `W<0` at every first-row root. In its A2B3
  specialization only a coupled mixed-product coefficient remains as a
  sufficient sign target.

The closest proved mechanism is
[THM-4440, signed duplication SOS](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md)
applied to that same real-rooted carrier. The canonical hostile remains
the actual `h=1,g=5,rho=-2` example: replacing a factor by its contiguous
neighbor loses the original zero coefficient. The corrected near miss is
the loss of real-rootedness after truncating a full auxiliary PF sequence.
The least-used relevant sidecar is
[local factorial interlacing](nc2_channel_contiguous_overnight_hexagon_sep05.md),
which already distinguishes local Euler moves from an unjustified global
root-order transfer.

The six live objects used here were: full auxiliary support; the same
first-row zero; completed alpha parity; the beta hit/skip decomposition;
the two contiguous differential identities; and the retained negative
square margin. Targeted searches for `K_B`, coupled mixed products,
coupling-only signs, beta/contiguous interlacing, and the exact named root
tuple found the incoming reduction and the local interlacing sidecar,
but no existing claim covering either witness below. This is a repository
overlap audit, not a literature novelty assertion.

In the incoming A2B3 notation, put

```text
q=x+h, g=x+3h+1=q+2h+1, k=2h+1,
H_B=(1+u)^g K_B, H_C=(1+u)^g K_C, H_D=(1+u)^g K_D.
```

The actual composition kernels satisfy

```text
K_B'=(2/3)u[(q+2)K_C+u K_C'],
(q+1)K_D+u K_D'=2K_C+(3/2)u K_C'.                 (1)
```

At negative phase `t`, the remaining sufficient target is

```text
[u^k]H_B=0  ==>  [u^(2k-2)]H_C H_D >=0.          (2)
```

The first test asks whether (1), the exact common binomial factor, and
the three real-rooted kernels suffice for (2). They do not.

## 2. A full-support, Euler-coupled counterexample

Fix `q=6,h=3,x=3,g=13,k=7`. Define the following polynomials in `v`:

```text
B(v)=(v-1)(v-2)(v-3)(v-5)(v-9)(v-12)
    =v^6-32v^5+380v^4-2110v^3+5739v^2-7218v+3240,

C(v)=v^5-30v^4+(2280/7)v^3-(3165/2)v^2
                   +(17217/5)v-10827/4,

D(v)=v^5-28v^4+(25080/91)v^3-(12660/11)v^2
                   +1913v-10827/14.
```

For every `s>0`, set

```text
K_B(s,u)=s^6 B(u^2/s),
K_C(s,u)=s^5 C(u^2/s),
K_D(s,u)=s^5 D(u^2/s),
H_i(s,u)=(1+u)^13 K_i(s,u).                       (3)
```

Writing `B(v)=sum b_i v^i`, the coefficients were defined by

```text
c_i=3(i+1)b_(i+1)/(q+2+2i),
d_i=(2+3i)c_i/(q+1+2i),       0<=i<=5.            (4)
```

Thus (1) holds identically in both `s` and `u`; it is not merely checked
at the eventual first root. The formal phase is `t=-s`.

All roots of B,C,D are positive and simple. For C and D, the following
five disjoint rational sign-change intervals exhaust their degree five:

| Kernel | Five intervals in increasing order |
|---|---|
| C | `(19/10,2)`, `(14/5,29/10)`, `(47/10,24/5)`, `(17/2,43/5)`, `(59/5,119/10)` |
| D | `(1/2,3/5)`, `(31/10,16/5)`, `(22/5,9/2)`, `(81/10,41/5)`, `(58/5,117/10)` |

The verifier evaluates both endpoints using rational arithmetic. Degree
exhaustion proves completeness and simplicity without approximate roots
or a root-isolation oracle. Consequently all three even kernels in (3)
have simple real nonzero roots, and all three H_i are real-rooted. The
common factor `(1+u)^13` is retained, with its allowed repeated root.

Equivalently the full positive-coefficient polynomials are

```text
G_B(t)=product_(a in {1,2,3,5,9,12})(1+a t),
G_C(t)=1+30t+(2280/7)t^2+(3165/2)t^3
                       +(17217/5)t^4+(10827/4)t^5,
G_D(t)=1+28t+(25080/91)t^2+(12660/11)t^3
                       +1913t^4+(10827/14)t^5.
```

Each has all roots negative and simple. The full Laurent sequences are
`t^-3 G_B`, `t^-3 G_C`, and `t^-3 G_D`; none has been truncated before
convolution. What has been removed is the specific binomial composition
coefficient law, not its full-support root geometry.

The original coefficient constraint is exact:

```text
[u^7]H_B(s,u)=26s^3 p(s),
p(s)=213840s^3-357291s^2+63129s-1055.              (5)
```

Signs at `0,1/10,1,3/2` are successively `-,+,-,+`. Hence p has exactly
three distinct positive roots, one in each intervening interval. Let
`s_*` be its unique root in `(1,3/2)`.

For the mixed response define

```text
[u^12]H_C(s,u)H_D(s,u)=s^4 T(s),
T(s)=sum_(i=0)^6 [v^i](C D) binom(26,12-2i)s^(6-i).
```

Exact polynomial division gives

```text
T(s) mod p(s) = N(s)/920602755072000,

N(s)=-3799751043113696471769929961 s^2
       +762374596098079635101930499 s
        -12900924131028983563588405.             (6)
```

Every coefficient of `N(1+w)` is strictly negative. Therefore
`T(s_*)<0`, and (3)--(6) prove

```text
[u^7]H_B(s_*,u)=0,
[u^12]H_C(s_*,u)H_D(s_*,u)<0.                    (7)
```

This is an exact counterexample to the proposed coupling-only extension
of (2), even with the original zero, full auxiliary support, common
binomial factor, and all three real-rooted kernels preserved. No
minimality in q, h, or coefficient height is claimed.

### What information the witness loses

C even strictly interlaces B in the squared variable. D does not: it
has two roots in `(3,5)` and none in `(2,3)`. In particular

```text
C(1)D(1)=-3320830339/25480<0,
C(3)D(3)=-11201247/40040<0.
```

This is a concrete failure of aligned contiguous root geometry. It also
prevents an automatic real-rooted perturbation of `B(v)^2` through
`B(v)^2-lambda v C(v)D(v)` for all sufficiently small positive lambda:
near either of those simple B roots the double root splits away from
the real axis. These extra conclusions are local observations about the
witness, not an assertion that interlacing alone repairs (2).

The actual `F_18,F_17,F_16` composition triple satisfies the same (1),
but it is a different triple. Its defining recurrence and exact doubled
composition identity carry additional information absent from (4).

### The negative margin survives this hostile

Keep the full binomial alpha sequence at g=13 and construct the formal
completed response by the same hit/skip formula. Denote it `Q_model`,
because these auxiliary beta coefficients have no asserted actual
multinomial-fibre interpretation. At `t=-s`, define

```text
W(s)=s^-6 [u^14]H_B(s,u)^2,
R_skip(s)=-2s^-5 [u^12]H_C(s,u)H_D(s,u)
         =-2s^-1 T(s),
Q_model(s)=W(s)+R_skip(s).
```

The virtual response V is the literal `O^2 star B_raw^2` at `t=-s`,
where `O_j=binom(13,2j+1)` and `B_raw=t^-3 G_B`. The exact raw Laurent
path includes the lower doubled index `j=-1`; both W and R_skip have a
nonzero contribution there. The binomial parity identity
`binom(26,2+2j)=(O^2)_j+(E^2)_(j+1)` is checked at every retained index.

At `s_*`, (6) says `R_skip>0`. Nevertheless exact remainders for
`sW`, `sQ_model`, `sV`, and `s(Q_model-V)` all have strictly negative
coefficients after shifting `s=1+w`. Thus this same witness proves

```text
R_skip(s_*)>0,        Q_model(s_*)<V(s_*)<0.       (8)
```

The full integer numerator/denominator certificates are in the frozen
output. This retains the distinction already made by the incoming note:
skip negativity is sufficient, whereas the necessary-and-sufficient
payment condition is `R_skip<-W`. The first failed implication is (2);
neither actual two-rung noncancellation nor even formal full-response
negativity is refuted by this witness.

## 3. Strict generic interlacing is also insufficient

A second cheap test challenges an attempted repair that would discard
the common binomial and even-kernel structure in favor of generic
interlacing. Put

```text
F(u)=(1-3u)(1-2u)(1+u)(1+u/4)
    =1-(15/4)u+(25/4)u^3+(3/2)u^4,
J(u)=u(1-5u/2)(1+3u/4)
    =u-(7/4)u^2-(15/8)u^3.
```

All roots are simple, and their exact order is

```text
-4 < -4/3 < -1 < 0 < 1/3 < 2/5 < 1/2,
 F     J     F   J    F      J      F.
```

Thus J strictly interlaces F, vanishes at the origin, and even may serve
as both factors of the proposed mixed product. But

```text
[u^2]F=0,             [u^4]J^2=-11/16<0.         (9)
```

This refutes the generic implication

```text
F real-rooted, J strictly interlaces F, J(0)=0,
[u^k]F=0  ==> [u^(2k)]J^2>=0.
```

Degree four is minimal for this square obstruction with
`deg J=deg F-1`, `J(0)=0`, and `0<k<deg F`: at degree two the only
candidate coefficient is the square of the linear coefficient of J;
at degree three the two candidates are respectively the squares of its
linear and leading coefficients. Degree one has no interior k.

Even a fully real-rooted completion pencil does not repair the generic
claim. Strict interlacing gives real-rooted `F+aJ` for every real a:
evaluate at the three J roots and use the unchanged quartic leading
term at both infinities to force four real roots. Consequently

```text
F^2-lambda J^2=(F-sqrt(lambda)J)(F+sqrt(lambda)J)
```

is real-rooted for every `lambda>=0`, but its middle coefficient is

```text
[u^4](F^2-lambda J^2)=-351/8+11lambda/16,
```

which becomes positive for `lambda>702/11`. Real-rootedness of a
completion pencil does not preserve the original coefficient constraint.
This generic quartic is not of the actual coupled composition form.

## 4. A proved same-zero lowering operation and the next test

The incoming all-A2 theorem `W<0` remains intact. The same-root negative
margin, together with the full contiguous path data, is the strongest
useful surviving mechanism. Both new hostiles are consistent with the
actual signed-transport conjecture.

The first identity in (1) gives one positive operation without requiring
the composition coefficient law. For `g=q+k`, define

```text
L(u)=(1+u)^(g-1)(K_B(u)+(2/3)uK_C(u)).
```

Then the following coefficient identity is **PROVED**:

```text
k[u^k](1+u)^gK_B = g[u^(k-1)]L.                 (10)
```

Indeed, differentiate `(1+u)^gK_B` and extract degree `k-1`.
The contribution involving `K_B'` is `(2/3)` times

```text
[u^(k-2)](1+u)^g((q+2)K_C+uK_C')
 = (q+k)[u^(k-2)]H_C
       -g[u^(k-3)](1+u)^(g-1)K_C
 = g[u^(k-2)](1+u)^(g-1)K_C,
```

where the last step uses exactly `g=q+k`. The other derivative term is
`g[u^(k-1)](1+u)^(g-1)K_B`, proving (10). In particular this operation
retains the same original zero while lowering its coefficient index.

There is a **CONDITIONAL** same-constraint SOS consumer. If `uK_C`
interlaces `K_B`, their real linear pencil is real-rooted, and so is L.
For the present even-kernel types,

```text
L(0)=K_B(0)!=0,
deg L=g+2q-1,
0<k-1<g+2q-1,             (g+2q-1)-(k-1)=3q>0.
```

Here `k=2h+1>=3`, and the degree is exact since `K_B` has degree `2q`
while `uK_C` has degree `2q-1`. Thus THM-4440 and (10) imply

```text
[u^k]H_B=0  ==> [u^(2k-2)]L^2<0.               (11)
```

The interlacing hypothesis is explicit; (1) and the three individual
real-rootedness statements do not automatically supply it. In fact our
Section 2 hostile already has this interlacing. It therefore has both
negative square margins, while its mixed skip coefficient is negative.
Equation (11) is a retained operation and a candidate ingredient for a
magnitude estimate; it does not repair skip negativity by itself.

The incoming [all-degree second-coefficient margin](signed_degree5_empty_core_next_sep06_uniform.md)
has a direct conditional consumer here. At `h=1`, the lowered index is
`k-1=2` and the exact degree of L is `n=3x+5`. Whenever the interlacing
hypothesis above holds, write `L=L(0)product_i(1+r_i u)`. Its nonzero real
roots give positive `e_2(r_i^2)`, so the incoming sharp finite-degree
constant `c_(3x+5)` strengthens (11) to

```text
-[u^4]L^2 >= c_(3x+5) L(0)^2 e_2(r_i^2).
```

The sharp uniform constant `(13-8sqrt(2))/3>1/3` also applies strictly.
This is a quantitative connection between the two concurrent lanes;
the required real-rootedness and same-zero predicate are kept explicit.
It is not another closure of the already-proved h=1 return family.

The next structural test should retain all of

```text
K_B=u^(2q)F_(3q)(t/u^2),
K_C=u^(2q-2)F_(3q-1)(t/u^2),
K_D=u^(2q-2)F_(3q-2)(t/u^2),
F_n=F_(n-1)+tF_(n-3),
H_i=(1+u)^(q+2h+1)K_i,
[u^(2h+1)]H_B=0.
```

A precise target remains the coefficient inequality (2) on that exact
composition family, or a quantitative bound paying a positive skip
term against the original and lowered square margins. A generic
interlacer inequality is already excluded by (9). A useful first
algebraic attempt is a quadratic coefficient identity involving L and
the step-one/step-three recurrence, with its residual term kept
explicitly; neither a fixed-h characteristic escalation nor another
broad sign bank addresses the lost coordinate identified here. No such
mixed-product identity is claimed in this note.

| Connection | Preserved data | Lost data / needed sidecar |
|---|---|---|
| Actual path kernels to the Euler-coupled PF model | Full support, the original zero, common binomial factor, all three real-rooted kernels, both differential identities | The exact composition coefficient law and aligned contiguous roots; hostile (7) |
| Coupled carriers to generic interlacers | Simple real roots, strict interlacing, a zero at the origin, original coefficient zero | Even-kernel/binomial geometry; minimal hostile (9) |
| Hit/skip signs to the full response | Exact lower carry and both parity classes | A positive skip requires its magnitude, not only its sign; survivor (8) |
| Source-paper replacement method to this probe | A smallest useful witness retaining the proposed invariants | The analogy is methodological only; no geometric theorem is transported |

### Incoming weighted-midpoint invariant

The subsequently integrated [weighted midpoint theorem](nc2_weighted_midpoint_overnight_hexagon_sep05.md)
proves a stronger predicate for the actual individual path factors: the
full midpoint-deletion pencil is PF for every nonnegative weight. In the
squared-root variable of Section2 this gives the real-rooted pencil
`w B(v)^2-2v C(v)D(v)`. Root and `orthogonal_returns` read the incoming
proof and compared its exact predicate with this note.

Our Euler-coupled hostile fails that added predicate already at the B-roots
v=1 and3: `C(v)D(v)<0`, so the defect `-2v C(v)D(v)` is positive. The
incoming square-pencil degeneration lemma then excludes real-rootedness
for arbitrarily large positive w. This pinpoints an actual-path invariant
that the Euler/PF abstraction lost. The hostile remains a valid refutation
of its stated weaker implication.

Conversely the incoming cubic Hadamard hostile retains the whole weighted
pencil but does not retain our two Euler identities and common binomial
carrier. Neither hostile settles the model preserving both collections
of predicates. A cheap next test can therefore retain the weighted pencil,
the Euler relations, the common binomial factor and the original joint
zero together. The exact composition recurrence remains available as the
fully specified target. Individual-factor signs concern their own roots;
the zero-preserving lowering above acts at the original joint constraint,
so transferring one conclusion to the other still requires a proof.

## 5. Verification, controls, and reproduction

The standalone [source](../../04-computation/laurent_transport_empty_core_next_sep06.py)
uses standard-library Fraction arithmetic only. Its universe is exactly
the q6/h3 hostile, the strict quartic hostile, and two actual controls:
the same-constraint `h=1,x=1,g=5,t=-2` row and the actual
`F_18,F_17,F_16` midpoint identity. The first actual control gives
`mixed=4591`, `sV=-1420`, `sW=-3428`, and `sR_skip=-9182` at s=2.

Independent constructions inside the verifier compare literal ordinary
carrier multiplication against full Laurent sequence convolution and
Hadamard multiplication. Root brackets, quotient reductions, and shifted
coefficient signs are exact. No numerical root or finite extrapolation
is a proof input. The simple degree counts provide the root-completeness
certificate. Every gate remains active under optimization.

```text
python3 -B 04-computation/laurent_transport_empty_core_next_sep06.py
python3 -B -O 04-computation/laurent_transport_empty_core_next_sep06.py
```

Both modes pass **126 explicit gates** and have byte-identical output.
The frozen [output](laurent_transport_empty_core_next_sep06.out) contains
all five exact remainder certificates and the declared root brackets.
Semantic digest:

```text
19f4aea7d102744370e21a4c50c7111fdff8cdf4019fc827c650c95e7910a08a
```

Raw LF SHA-256:

```text
dd2b50e96f6baf2abeeb15dc5c0383759eca00f8ac3815c341598ce75b6646ea  source
d450c7976facaf6895fceace0807320705e0dc2605842136bc84752dd8ee79ae  output
```

The probe changes no canon ID, shared navigation, or Git state. It
records a proof-method boundary and its stronger surviving response;
it does not promote the open actual Laurent theorem.

**Independent root audit: PASS.** Root checked the two hostile proofs,
all scalings and degree-exhaustion arguments, the distinction between
formal and actual responses, strict-interlacer minimality, and the full
completion pencil. Root independently derived (10) using coefficient
Euler and checked the conditional scope of (11).

Additional independent final text audit by `certificate_audit`: **PASS**.
The referee reconstructed both Euler identities, the first-row polynomial,
the mixed coefficient and all formal response remainders without importing
the producer; audited the strict quartic interlacer and minimality; and
checked the zero-preserving lowering and its conditional h=1 sharp-margin
consumer. The final 126-gate replay and both raw hashes agree.
