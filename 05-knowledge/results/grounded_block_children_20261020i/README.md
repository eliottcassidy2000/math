# Stable reset blocks, the pending 31 source, and grounded child factories

20 October 2026. **PROVED:** an exact invariant interval for discarded
debt, its block-composition and variable-last-letter laws, source-preserving
paid rays, and an explicitly grounded family with unbounded reset count.
**FINITE-EXACT:** all175 previous controller-pending sources are
accepted and grounded by supplied words:146 use the unit-cap monoid and29
use a separate cap3 with their actual source gaps. Universal controller coverage,
arbitrary-ray child grounding, and signed ROOT-energy estimates remain
**OPEN**. No literature-priority claim is made.

Files: [checker](check.py), [supplied words](inputs.json),
[receipt](receipt.json). Production discovers no orbits or ROOT suffixes.

## 1. Inheritance and the operation being added

The [previous terminal controller](../debt_annulus_20261020h/README.md)
left 175 of the 5,000 odd sources from 3 through 10,001 pending. Its
terminal grammar allowed a single debt-producing reset followed by twos.
It grounded 76 additional named sources, with 63 as its useful ray.
The incoming [source63 separator theorem](../source63_separator_debt_20261021/SOURCE63_SEPARATOR_DEBT.md)
adds a separator-normalized cone and an all-depth generator architecture.
Its separator is an actual word only where the source guard permits it.
The [native tail recursion](../collatz_native_tail_recursion_20261013c.md)
already supplies the fixed-child exponent lift; the present result
transports a new block-payment certificate through that inherited mechanism.

The canonical hostile is still a later unpaid descent, such as the
inherited source55 example: an arbitrary descending word does not inherit
a discarded-debt estimate by forgetting its earlier resets. The corrected
near miss is a finite bank of fixed words being mistaken for universal
coverage. The least-used sidecar here is the **incoming debt**, considered
as a variable rather than evaluated only at zero.

The live board is **actual source / ordered block / incoming debt /
ternary phase / retained child / first ROOT and three clocks**. The
anchor is the pending source31. The niche is an invariant interval under
block composition. The wildcard is an unbounded reset architecture whose
source and final phase change with the architecture.

Connection contract: a complete native receipt is mapped to its affine
debt action, retaining the source and word separately. The preserved
predicate is that every incoming debt below one stays below one. Positive
emissions, native guards, endpoints and ROOT are not recovered from that
scalar action; they must be authenticated from the retained word and child
certificate. The inverse child factory changes the source explicitly; it
does not replace the unresolved child of an old source by a neighbour.

## 2. The exact unit-debt block monoid

For a word $w=(a_1,\ldots,a_m)$ of positive exact odd valuations, retain

$$
F_w(n)-1=\lambda_w(n-1)+P_w-D_w,
\qquad F_w(n)=\frac{3^m n+B_w}{2^{A_w}},\quad A_w=\sum a_i.
$$

The inherited nonnegative emissions satisfy

$$
D_{uv}=\lambda_vD_u+D_v,\qquad
\lambda_{uv}=\lambda_v\lambda_u.                       \tag{1}
$$

For a single letter, $\lambda_a=3/2^a$, while $D_1=D_2=0$
and $D_a=1-4/2^a$ for $a\ge3$.

**Interval theorem (iff).** A word $v$ maps every incoming debt
$d\in[0,1)$ to $\lambda_vd+D_v\in[0,1)$ exactly when

$$
\boxed{\lambda_v+D_v\le1.}                            \tag{2}
$$

Sufficiency follows from $\lambda_v>0$ and $d<1$. Necessity follows by
taking $d\uparrow1$; if the sum exceeds one, a rational $d<1$
already fails. The empty word is the identity. The satisfying words form
a monoid under actual ordered concatenation: composing two interval maps
again preserves the interval. This is an exact maximal class for the
stated invariant interval, not a claim that every paid word must satisfy it.

Every single letter $a\ge2$ satisfies (2), so arbitrarily many consecutive
resets are admitted. Every word over $\{1,2\}$ with
$3^{|v|}\le2^{A_v}$ also satisfies it, since its discarded debt is zero.
Such blocks can be freely interspersed with the reset letters **provided
the complete resulting word is native at the actual source**. For example
$(1,2,2)$ has multiplier $27/32$ and zero discarded debt.

Order matters: $(1,2,3)$ satisfies (2), but its rearrangement $(2,3,1)$
does not. The criterion is not a count of resets or a residue class of the
valuation cost.

**Last-letter theorem (iff).** For any head $b$ and any $a\ge2$,

$$
ba\text{ satisfies (2)}\quad\Longleftrightarrow\quad
\lambda_b+D_b\le\frac43.                             \tag{3}
$$

At $a=2$ the left-hand sum is $3(\lambda_b+D_b)/4$.
At $a\ge3$ it is
$1+[3(\lambda_b+D_b)-4]/2^a$. Thus increasing the actual final
valuation never destroys a stable block. The equality boundary is
included, and the domain $a\ge2$ is essential.

## 3. Source-specific payment and whole native rays

Let $u$ be a cone prefix, meaning
$\lambda_u\ge1$ and $D_u\le\lambda_u-1$. Let the nonempty terminal
$v$ satisfy (2). If the complete **actual** word $uv$ takes a positive
odd source $n>1$ to an odd $h<n$, then its positive envelope is paid:

$$
D_{uv}<n-h.                                         \tag{4}
$$

Indeed the complete word has positive raw carry, so
$\lambda_v\lambda_u<h/n<1$. Hence, using (1) and (2),

$$
D_{uv}\le\lambda_v(\lambda_u-1)+1-\lambda_v
       =\lambda_v\lambda_u+1-2\lambda_v<2\le n-h.
$$

This includes an empty cone prefix. It requires an actual descending
endpoint, but not that this is the first descent. The interval theorem
alone supplies neither a source guard nor an endpoint or ROOT certificate.

Now fix such a supplied word $w$, of length $m$, cost $A$, last valuation
$a_0\ge2$, carrier $(P,Q,B)=(3^m,2^A,B)$, and actual endpoint $h$.
For every integer $t\ge0$, set

$$
n_t=n+Qt,\qquad
a(t)=a_0+v_2(h+Pt),\qquad
h_t=\operatorname{oddpart}(h+Pt).                    \tag{5}
$$

The first $m-1$ valuations stay fixed: their complete native guard needs
only the corresponding cost plus one odd-terminal bit, at most $A$.
The final numerator is exactly $2^{a_0}(h+Pt)$, proving (5).
Since the original descending word has $P/Q<h/n<1$, one has
$h_t\le h+Pt<n+Qt$. The same cone prefix remains fixed, and (3)
retains terminal stability. Thus the entire coarse ray is paid.

This is a paid-child theorem. It is not ROOT for all $h_t$, and for a
general supplied word it is not a first-descent theorem for every lift.
Earlier proper coefficients can be contracting even if the original
seed's first descent happened later.

## 4. The exact finite pending closure

The input contains the complete old list of 175 pending sources, each
with its actual first-exit word. The consumer checks every proper
checkpoint against the original source, computes the latest eligible
cone/stable-block cut, and replays a separately supplied strict ROOT suffix
at each accepted child. The 146 accepted sources use 95 distinct children;
all 95 suffixes are retained. Their discovery cap was 10,000 odd steps
per child; the longest supplied suffix has 68 letters. Production performs
only replay and exact inverse source recovery.

| Disjoint group | Count |
|---|---:|
| Previously pending sources |175|
| New consecutive-$\ge2$ terminal cases |117|
| Additional mixed stable blocks |29|
| Pending under the unit-cap controller |29|

Thus the inherited controller union increases from 4,825 to 4,971 of its
declared 5,000 sources. The least unit-cap residual source is927. The complete
29-element unit-cap remainder is retained in the receipt; Section7 supplies
a stronger source-gap controller for exactly these29.

Examples beyond consecutive resets include source487 with terminal
$(3,1,2,6)$ and source2047 with terminal $(5,1,2,2,3)$.
The new Mersenne controller anchors in the declared universe are exactly
31, 127 and 2047. Their paid rays (5) contain only their respective
Mersenne anchor: if $E\ge A$, $2^E-1$ has residue $-1$ modulo $2^A$;
for $E<A$ it lies below $2^A$, where the retained seed is the sole
possible representative. This prevents a finite Mersenne observation
from being promoted to an all-exponent result.

## 5. Grounded factories retain the child, denominator, and source

Take a fixed head $b$ of length $l$, cost $S$, and carrier
$(3^l,2^S,B_b)$. Supply a cone cut inside it, such that the remaining
terminal head satisfies (3). Also supply a strict ROOT certificate at a
positive odd child $h>1$ with $3\nmid h$. Put

$$
P=3^{l+1},\quad C=3B_b+2^S,\qquad
N(a)=\frac{h2^{S+a}-C}{P}.                           \tag{6}
$$

Because $C\equiv2^S\pmod3$, the congruence

$$
2^a\equiv C(h2^S)^{-1}\pmod P                       \tag{7}
$$

has a unique exponent class modulo $H=2\cdot3^l$. Indeed 2 generates
the units modulo every power of three: parity selects the residue modulo
three and each successive ternary digit has a unique lift. Choose the
least exponent in this class above both the stated minimum $a\ge2$ and
the strict positivity/descent cutoff

$$
h2^{S+a}>C+Ph.                                      \tag{8}
$$

Call it $a_*$. For every $j\ge0$, the source $N(a_*+Hj)$ is an odd
integer greater than $h$, and its actual word is $b$ followed by that
last valuation. Integrality plus the odd terminal identity supplies the
full native guard, not just a rational affine equality. Here is a direct
induction. For a complete word of length m with first valuation a_1,
its carry is 3^(m−1) modulo 2^a_1. The endpoint identity therefore forces
3n+1 divisible by 2^a_1. Substitute n_1=(3n+1)/2^a_1 into the remaining
carrier equation. If a remaining letter exists, its carry is odd and its
right side even, forcing n_1 odd; at the final letter the supplied h is
odd. Thus each division has exactly the stated valuation, and the same
induction continues. The original source is odd by the whole identity
because its carry is odd (also for an empty head, C=1). Forward evolution
starts positive, so every intermediate state is positive.

The prefix cannot hit ROOT early: actual odd evolution fixes 1 forever,
whereas the authenticated endpoint is $h>1$. Appending the supplied
first-hit suffix at $h$ therefore gives a strict ROOT word at the same
source $N(a_*+Hj)$. Statement (4) pays the entire prefix to $h$.

For each of the 146 finite new records, take its last letter as the minimum
in this construction. Then $a_*$ is precisely that seed's actual final
valuation and $N(a_*)$ is precisely the original seed. Consequently there
are 146 explicitly anchored infinite grounded factories; intersections
between different factories are not counted as new disjoint coverage.
They lie inside the corresponding paid rays (5).

For any query modulus $M\ge1$, read the same symbolic source by

$$
N(a)\bmod M=
\frac{h\,[2^{S+a}\bmod PM]-C}{P}\pmod M.             \tag{9}
$$

The numerator is divisible by $P$. Reduction modulo $M$ before division
is invalid when the discarded multiple is not divisible by $PM$.
The query keeps $a$, the source index, head, child, and denominator.
The checker compares the compatible mod15/mod35 registers with their
joint mod105 value; it does not choose those residues independently.

## 6. Source31 with arbitrarily many actual reset blocks

The first-exit word at31 has 35 letters and cost56, taking31 to23.
Write it as $u(4,3)$, where $u$ has 33 letters and cost49. Its exact cone
data are

$$
\lambda_u=\frac{5559060566555523}{562949953421312},\qquad
D_u=\frac{209651139}{134217728}<\lambda_u-1.
$$

The terminal $R=(4,3)$ satisfies

$$
\lambda_R=\frac9{128},\quad D_R=\frac{25}{32},\qquad
\lambda_{R^s}=(9/128)^s,\quad
D_{R^s}=\frac{100}{119}\bigl(1-(9/128)^s\bigr).       \tag{10}
$$

The old single-reset controller had no eligible cut for31. The new
stable terminal does. The actual child23 has supplied suffix
$(1,1,5,4)$, giving31 a strict ROOT word of39 odd steps, valuation cost67,
and106 ordinary unaccelerated steps.

For **every** $r\ge1$, use the head and terminal architecture

$$
b_r=u(4,3)^{r-1}4,\qquad w_{r,a}=u(4,3)^{r-1}(4,a),
\quad a\ge3.                                       \tag{11}
$$

Every terminal letter in (11) is at least two, so the monoid theorem
proves stable terminal debt for every $r,a$. Apply (6)-(8) with the same
child23 and its same supplied suffix. Here

$$
l_r=32+2r,\quad S_r=46+7r,\quad P_r=3^{33+2r},
\quad H_r=2\cdot3^{32+2r}.                           \tag{12}
$$

This yields grounded sources $N_{r,j}$ for every $r\ge1,j\ge0$.
For $r=1$, the least phase is $a_*=3$, so $N_{1,0}=31$ literally:

$$
N_{1,j}=
\frac{23\,2^{56+(2\cdot3^{34})j}-106346764803351611}{3^{35}}.
                                                               \tag{13}
$$

For $r=2$, the least final valuation is 265513516350064983;
it is stored symbolically, never expanded into its enormous source.
The all-$r$ theorem comes from (2), (3), and the exact exponent phase,
not from finite agreement at these two values.

The ROOT clocks for (11), including the four-letter child suffix, are

$$
\text{odd steps}=37+2r,\quad
\text{valuation cost}=57+7r+a,\quad
\text{ordinary steps}=94+9r+a.                       \tag{14}
$$

Only at $r=1$ do we assert the displayed coarse ray's first odd descent:
all proper prefixes of the seed word have multiplier greater than one.
For $r>1$, intervening contracting blocks can descend earlier. The theorem
retains original-source payment at23 and first ROOT, without asserting
a common first-descent location. One source is never claimed to realize
all architectures $r$.

For fixed $r$, the source registers have minimal index period six. Since
$23$ is coprime to $105P_r$, equality of two registers is equivalent to

$$
2^{H_r(j-j')}=1\pmod{105P_r}.
$$

The order is $4\cdot3^{33+2r}$, whose quotient by its gcd with $H_r$
is six. For $r=1$, the exact mod105 table at indices0 through5 is

$$
(31,3,101,73,66,38).
$$

Its mod15 and mod35 projections are compatible views of the same source.
The six phases are not six unrelated CRT choices or six ROOT residues.

## 7. Paying a larger terminal cap with the actual source gap

For an exact integer c≥1, replace (2) by

$$
\lambda_v+D_v/c\le1.                                \tag{15}
$$

This is precisely preservation of the incoming interval [0,c). With
an ordinary cone prefix u and a complete actual descent n→h, put
μ=λ_v. Positive raw carry again gives μλ_u<h/n<1. Therefore

$$
D_{uv}\le\mu(\lambda_u-1)+c(1-\mu)
       <(c+1)(1-\mu)<c+1.                           \tag{16}
$$

Consequently the **original-source** gap n−h≥c+1 pays the prefix.
The gap is an essential additional premise. CapPaid retains source,
word, cut and c; audit_cap checks all four, every native valuation, and
the actual endpoint before accepting. It is separate from the unit-cap
factory API rather than silently broadening (3).

At c=3, every one of the29 unit-cap residual first-exit words satisfies
(15) at a cone cut. Their minimum actual source drop is50, so (16)
pays them all. They use18 additional supplied strict child ROOT suffixes,
bringing the complete used child bank to113. These18 were found with the
same10000-step cap; their longest supplied suffix has55 letters. Production
replays supplied words only. Thus the union of the inherited controller,
146 unit-cap records and29 cap3 records covers **exactly the declared5000
odd sources3 through10001**. It does not prove an all-source assertion.

The incoming [THM-4735 fixed-threshold obstruction](../../../01-canon/theorems/THM-4735-valuation-dependent-debt-cones-optimal-terminal-splits-and-the-debt-two-obstruction-to-universal-cone-coverage.md)
is compatible with this result. Its source4575 has first drop110 but

$$
D=2377272981729/1099511627776>2.
$$

A fixed debt-two sufficient criterion misses it. The cap3 criterion instead
retains the actual gap110; it does not deny the debt-two obstruction or
claim that an optimized split can lower the already-fixed total debt.
The incoming larger census is a cited dependency boundary, not rerun here.

For any terminal head b and a≥3 the exact completion expression is

$$
\lambda_{ba}+D_{ba}/c
 =1/c+\frac{3\lambda_b+(3D_b-4)/c}{2^a};           \tag{17}
$$

at a=2 it is 3(λ_b+D_b/c)/4. For c≥1, once this expression is≤1,
increasing a preserves that inequality: from a2 to a3 the value becomes
half the old value plus1/(2c); thereafter it approaches1/c along a line.
Therefore the cap3 paid seeds also lift through (5). Their actual drops
are at least n−h+(Q−P)t, so the original cap budget remains paid.
These coarse rays are **paid**, not automatically grounded.

For c>1, any fixed terminal head eventually satisfies (17) by taking a
large enough last valuation. Together with the exact phase (7) and a
larger source cutoff n≥h+c+1, this yields fixed-child grounded inverse
families for arbitrary supplied heads and retained ternary-unit h>1.
This is a source-changing completion operation. It does not ground an
original small source whose last valuation is already fixed. The concrete
symbolic factory implementation in this packet remains the sharper146
unit-cap factories; no new generalized-factory API is implied.

## 8. Tests, remaining obligation, and scope

The checker consumes all175 supplied first-exit records and113 strict ROOT
words, validates the immutable old pending list, and recovers each ROOT
source independently by reverse integer divisions. It checks the 146
anchored factories, six literal coarse lifts each, symbolic large-index
readers, exact old/new controller counts, and the complete unit-cap residual list plus its29 cap3 certificates.

Independent suffix formulas cover every word of length0 through5 over
valuations1 through4. Stable-word composition is exhaustively tested for
stable words through length3. All ternary-unit exponent phases through
precision6 are independently enumerated. The repeated architecture is
checked at $r=1,2,3,4,8,16$, always as a finite control for its all-$r$
proof. Type, wrong-child, changed-head, changed-phase, ROOT-padding and
materialization-cap hostiles stay active under optimization.

```text
python -B 05-knowledge/results/grounded_block_children_20261020i/check.py
python -B -O 05-knowledge/results/grounded_block_children_20261020i/check.py --verify
```

The declared5000-source controller list is now complete. The arbitrary
source complement is not bounded by this finite test set. The
new algebra permits composable reset blocks and certifies actual grounded
families, but supplies no universal guard selector or proof that every
orbit reaches one of them. This positive discarded-debt bound is also
distinct from the signed, clocked terminal-energy and annulus obligations.

The sibling-interaction and [temporal pair-block ledger](../owned_pair_blocks_20261020i/README.md)
retain different signed Gram terms and arrival cuts. The present discarded
positive-envelope debt is not their potential; none of the grounding or
cap inequalities supplies their still-open uniform signed-energy bound.

Normal and independent optimized `--verify` runs agree at **54,400 checks**.
LF receipt SHA256:
`f07e35ba635943e89c66f58ba3beeaf6c03cbc839c2aff0e90122bd3238ba6da`.
The parent independently reviewed the full note and production interfaces,
including the cap3 original-source gap proof, with no mathematical correction.

Independent peer `tiling_fold` completed the full proof/API audit and a
separate optimized-versus-saved replay: **PASS**, the same54,400 checks
and receipt hash. It also independently confirmed that all113 supplied
child ROOT suffixes are used. The earlier saved-output mismatch was the
superseded22,440-check receipt before the declared final refresh, not a
mathematical failure. All four files are frozen after this audit footer.
