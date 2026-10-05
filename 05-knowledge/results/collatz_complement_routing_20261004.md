# A disjoint complement family by composing before testing source size

2026-10-04. **PROVED:** guarded smaller-dependency families, a positive-density
extension of the specified bank, a general mod19 refinement, and a whole
infinite family excluded from this schema. **FINITE-EXACT:** the declared
controls. **OPEN:** universal coverage and independence of all the constructed
families. No historical-priority claim.

Artifacts: [program](../../04-computation/experiments/collatz_complement_routing_20261004.py)
and [saved output](collatz_complement_routing_20261004.out).

## 1. Inheritance and the changed operation

The closest mechanisms are the ternary sibling rule in
[binary/ternary fusion](collatz_binary_ternary_guard_fusion_20261004.md) and
the unbounded [reset-at-least-three switch](checked_switch_phase19_20261004.md).
The [sibling grammar](creative_sibling_20260925.md) requires retaining
intermediate choices; the [frontier compiler](frontier_family_compiler_20261004.md)
retains the original source in a composed size comparison. These operations
are inherited. Here a sibling rule is pulled through a reset-two prefix and
then composed with a second common-future switch.

Canonical hostiles are7,27,703. The corrected near miss is discarding a
composition because its first auxiliary child exceeds the original request.
The least-used sidecar is the coefficient of the whole transport. The live
board is **original source / composed child / exact word / ternary guard /
disjoint selector class / mod19 phase**. Anchor: an infinite complement
family. Niche: precision versus child size. Wildcard: phase-preserving
refinement. This does not alter the Collatz map, unlike the pairing
modifications of [THM-4475, price of provable descent](../../01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md).

Targeted constant/formula searches found no earlier explicit copy. No novelty
is claimed for the inverse identities. The proved extension is relative to
the sixteen selected binary rows, the entire single-stage ternary bank, and
direct/reset-at-least-three tests, not every historical proof grammar.

## 2. Compose first, then compare with the original source

Write $U(n)=\operatorname{oddpart}(3n+1)$ and $S(x)=4x+1$. Initially
restrict $n=123\pmod{128}$, fixing valuations $(1,2,1,2)$. Then

\[
 x=U^2(n)=\frac{9n+5}{8},\qquad J=U(x)=\frac{27n+23}{16}.
\]

The first four affine multipliers are $3/2,9/8,27/16,81/64$ and all carries
are positive, so all four source values exceed the original $n$.

If $3^\ell\mid S^k(x)+1$, form

\[
 m=2^\ell\frac{S^k(x)+1}{3^\ell}-1,\qquad h=(m-1)/2.
\]

The intermediate $m$ follows $1^\ell$ to $S^k(x)$, whose next valuation is
$2k+1\ge3$. Composing the inherited reset switch gives the actual diagram

\[
 n\xrightarrow{(1,2,1)}J
 \xleftarrow{(1^{\ell-1},2,2k-1)}h.                 \tag{1}
\]

We do not require $m<n$ and do not introduce an unproved certificate for$m$.
The final child in(1) is the only home-certificate premise.

The exact ternary guard is

\[
 4^k(27n+23)+16=0\pmod{3^{\ell+1}}.                \tag{2}
\]

The size condition below requires $\ell\ge3$. Modulo27, (2) becomes
$4^k=4\pmod{27}$, equivalent to $k=1\pmod9$ (the nine powers are distinct).
For any such$k$, define

\[
 g_k=\frac{23\,4^k+16}{27},\quad M=3^{\ell-2},
 \quad n=-g_k4^{-k}\pmod M.
\]

Together with $n=123\pmod{128}$ this gives one arithmetic progression.
On it,

\[
 h=\frac{2^{\ell+2k-4}n+2^{\ell-4}g_k-M}{M}
   =\rho(n+c_k)-1,\quad
 \rho=\frac{2^{\ell+2k-4}}{3^{\ell-2}},\quad
 c_k=\frac{23+16/4^k}{27}\le1.                     \tag{3}
\]

At the exceptional smallest case $\ell=3$, $g_1=4$ makes the fractional
power-of-two notation integral. Choose the least positive$\ell$ with

\[
 (3/2)^\ell>(9/16)4^k.                            \tag{4}
\]

It exists for every$k\ge1$ and is at least3. At$k=1,\ell=2$ equality
would give$h=n$, explaining the strict boundary. Formula(3) gives
$h\le\rho(n+1)-1<n$.

Positivity and oddness need a separate guard argument: $(S^k(x)+1)/3^\ell$
is positive even. Hence$m\ge2^{\ell+1}-1$ and$h\ge2^\ell-1>0$; the
exact inverse-one and reset identities make$h$ odd. This proves(1) and
the strict smaller dependency at every parameter height. The uncomposed
coefficient is$2\rho$ and can exceed1; index19 at the former single-stage
depth65 explicitly has$m>n>h$.

Every$k=1+9q$ therefore supplies a valid family, with no irrational-rotation
argument. They are not asserted independent or all new. Index1 is redundant
with the old ternary index0.

## 3. An explicit positive-density extension

For$k=10$, the optimized depth is$\ell=33$. The whole family is

\[
\begin{aligned}
 n_t&=27727075633746555+79062194724345216t,\\
 h_t&=25270565367447551+72057594037927936t,\qquad t\ge0.
\end{aligned}                                    \tag{5}
\]

Its child word is$1^{32},2,19$ and
$h=(2^{49}n-138123117816363)/3^{31}$, with slope
$562949953421312/617673396283947<1$.

Every member misses the sixteen binary rows: initial run1, reset2 and
next valuation1 select only the$(s,e)=(1,1)$ row, which asks for fourth
valuation6, whereas ours is2. No other row has this initial run/debt type.
The first four values grow and the first reset is2, also excluding direct
and reset-at-least-three source rules.

Write the old ternary guard as$n=R(j)\pmod{3^{E(j)}}$, where$E(j)$ is
least with$(3/2)^{E(j)}>4^j$. The new ternary class is
$549446197252887\pmod{3^{31}}$. Its exact residues rule out every old
ancestor of depth at most31:

|j|E(j)|new residue at that depth|old R(j)|
|---:|---:|---:|---:|
|0|1|0|2|
|1|4|9|40|
|2|7|252|273|
|3|11|87732|94109|
|4|14|2744937|2503585|
|5|18|261025263|62804493|
|6|21|5684912109|5439587969|
|7|24|120748797342|221027314255|
|8|28|403178333823|7257210008829|
|9|31|549446197252887|1814302502207|

All other old guards have$j\ge10$, $E(10)=35$, and
$E(j+1)\ge E(j)+3$. The last inequality follows because, by minimality,
$(3/2)^{E(j)+2}<(27/8)4^j<4^{j+1}$.
A compatible guard occupies fraction$3^{31-E(j)}$ of(5); otherwise its
intersection is empty. The entire infinite old union therefore occupies
at most

\[
 \sum_{j\ge10}3^{31-E(j)}
 \le\frac{3^{-4}}{1-3^{-3}}=\frac1{78}.              \tag{6}
\]

This is a relative **natural-density** statement. The inherited positive
child bound permits only$O(\log X)$ eligible old indices for integer sources
at most$X$. Counting each intersection progression contributes at most one
endpoint error. These total$O(\log X)$ errors vanish after division by$X$;
the geometric tail then proves density existence by finite truncation.
The binary modulus128 is coprime to every ternary modulus, so it does not
change the relative fractions.

Define the new selector class as(5) minus all old selected domains.
Membership is effective using that finite source-height bound. The class
is disjoint by definition and has odd-relative density at least

\[
 \frac{77}{78}\frac1{64\,3^{31}}
 =\frac{77}{3083425594249463424}>0.                 \tag{7}
\]

Thus infinitely many previously unselected sources acquire a strictly smaller
dependency. This is neither a nonconvergence density nor an unconditional
home certificate for every free parameter.

## 4. General initial runs and an all-length height obstruction

The composition extends to any actual prefix$(1^r,2)$, $r\ge1$. Now

\[
 x=\frac{3^{r+1}(n+1)}{2^{r+2}}-\frac12.
\]

Choose$k=1\pmod{3^{r+1}}$ and put

\[
 C=\frac{2^{r+1}(4^k-4)}{3^{r+2}},\quad
 M=3^{\ell-r-1},\quad
 \lambda=\frac{2^{\ell+2k-r-3}}{3^{\ell-r-1}}.
\]

Take the least$\ell\ge r+2$ with$\lambda<1$. The exact ternary guard is
$4^k(n+1)=C\pmod M$. Its necessity follows from the valuation
$v_3(4^k-4)=1+v_3(k-1)$, with$k=1$ handled by zero carry.
Together with the exact binary prefix guard, CRT gives a nonempty source
progression. The final child satisfies

\[
 h=\lambda(n+1-C/4^k)-1<n,                         \tag{8}
\]

because$C\ge0$; positivity follows from the same inverse-integrality
argument. If$a=v_2(3x+1)$, the actual common words are

\[
 (1^r,2,a),\qquad (1^{\ell-1},2,a+2k-2).           \tag{9}
\]

No condition that the intermediate$m=2h+1$ be smaller is imposed.
Minimal depths for$(r,k)=(1,10),(2,28),(3,82)$ are checked explicitly.

This also gives a complete finite applicability search for a supplied
source within this schema. Positivity gives$h+1\ge2^\ell$; since$h<n$,
$2^\ell<n+1$. The allowed depths are finite. From$\lambda<1$ and$3^b<4^b$,
also$2k<\ell-r+1$. Enumerate the finitely eligible$k=1\pmod{3^{r+1}}$,
testing their minimum-depth guards. A missed minimum guard cannot become
true at greater precision because those guards are nested.

**PROVED obstruction to universality of this schema.** The whole family

\[
 n=2^{2j+1}-1,\qquad j\ge1
\]

has an initial run$r=2j$ and first reset2. Every construction above requires
$\ell\ge r+2$, but$h+1\ge2^\ell>n+1$, contradicting$h<n$. Thus no index
or greater depth in this schema handles those supplied sources. This
includes7 and is an all-length obstruction, not a finite search failure.
Other constructions may handle them. A universal family dispatcher still
needs an operation outside this one.

## 5. A general phase-preserving refinement

For every general construction,$k=1\pmod9$, so$4^k=4\pmod{19}$.
The ratio$2/3=7\pmod{19}$ has order3. Direct substitution yields

\[
 h+1=(2/3)^{\ell-r-1}(n+1)\pmod{19}.               \tag{10}
\]

For example,$S^k(x)+1=4x+2=
3^{r+1}(n+1)/2^r\pmod{19}$. Choose$d\in\{0,1,2\}$ with
$\ell+d-r-1=0\pmod3$. Increasing depth by$d$ refines the source guard
by$3^d$, and on that refined class

\[
 h_{\rm new}+1=(2/3)^d(h+1),\qquad h_{\rm new}=n\pmod{19}. \tag{11}
\]

The congruence is preserved by construction while size decreases over the
integers. It is not a residue-only proof. Phase$n=-1$ is fixed at any depth;
the stated depth congruence preserves every phase.

For(5), the refinement is$t=6\pmod9$, $\ell=35$. Writing$t=6+9u$ gives

\[
\begin{aligned}
 n_u&=502100243979817851+711559752519106944u,\\
 h^{(35)}_u&=203384946486673407+288230376151711744u.
\end{aligned}                                    \tag{12}
\]

Both constant terms are0 mod19 and both periods16 mod19. Every phase
occurs and is retained. An independent coefficient check is
$h^{(35)}=(2^{51}n-D)/3^{33}$, where
$D=19\cdot191624181720273$ and
$2^{51}-3^{33}=-19\cdot174066355414225$.
This child has first reset2 followed by19, whereas its binary-bank row
asks6; phase preservation is not preservation of a selector class.

There is a useful exact recursive partition, inherited from index0.
For any resulting child let$e=v_3(h+1)$ and
$g_e=2^e(h+1)/3^e-1$. For$e>0$, these are$e$ strict, positive odd
predecessor decreases, with word$1^e$ from$g_e$ to$h$.
At maximal depth,$g_e\bmod3\in\{0,1\}$, so that specific rule stops.
For$e=0$ leave$h$ unchanged. These depth classes partition the parameter
space; they do not claim all resulting states enter another known class.
On(5),$g_e+1=7^{e+1}(n+1)\pmod{19}$, so the phase returns for$e=2\pmod3$.
On(12), it returns for$e=0\pmod3$.

## 6. Reproduction and exact scope

    python 04-computation/experiments/collatz_complement_routing_20261004.py
    python -O 04-computation/experiments/collatz_complement_routing_20261004.py

Controls comprise126 actual-word replays (indices1..181 step9, six parameter
values through1000000); general phase refinements for21 indices; three
general-run pairs at three heights;512 members of(5), all dispatched to
the new disjoint class;32 comparisons against the inherited selector;
512 maximal predecessor depths; all19 phases of(12); and24 Mersenne
applicability controls. These finite counts are not the proofs of the
infinite statements.

Four AST demonstrations explicitly obtain the smaller-child premise by a
bounded route search of at most5000 odd steps, then splice the supplied
certificate to the original source. Independent source replay is performed
only afterward. There is no free-parameter home inference. Root cannot
occur early on either displayed route in(1), since its join$J>n>1$.
Six malformed/guard/source cases are rejected, and7,27,703 remain pending
under the specialized new selector.

The actual interface returns a smaller obligation with both exact words.
An unavailable child remains unavailable. The [fair scheduler](fair_frontier_extension_20261004.md)
can retain this information without abandoning original-source exploration.
The next problem is to dispatch the remaining complement with other
well-founded constructions, including the proved hard Mersenne family,
rather than extrapolate from the successful members.
