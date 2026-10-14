# Squaring residues: exact fibres, signed observers, and false anchor collisions

**PROVED:** the elementary fibre formulas, marked signed dynamics, and complete
false squared-anchor classification below. **FINITE-EXACT:** the declared
residue, word and source controls. **OPEN:** universal Collatz completion.
This packet adds no source coverage and discovers no ROOT word.

Artifacts: [checker](../../04-computation/experiments/residue_power_defects_20261014d.py)
and [saved JSON](residue_power_defects_20261014d.json).

## 1. Inheritance and the question that survives the quotient

The closest previous mechanism is the
[marked square-norm guard compiler](collatz_square_norm_transfer_20261007f.md).
It already transports binary and ternary source guards at unchanged precision
when their square-root branches are retained. We do not claim that mechanism
anew. The [signed basin observer](collatz_basin_observer_20261009.md) retains
the distinguished terminal, rather than identifying different signed cycles.
The new incoming anchor algebra is
[THM-4640, repair algebra at minus one](../../01-canon/theorems/THM-4640-the-repair-algebra-at-minus-one-collisions-are-increasing-two-three-representations-ladders-and-gersonides-carries.md).
Its header explicitly says that independent audit is pending. Only the short
numerator calculation reproduced below is used here; none of its census,
coverage trend or external finiteness claims is a dependency.

The live concepts are source residue, pushforward multiplicity, sign,
ordered word, anchor value and authenticated terminal. The anchor is unchanged
Collatz source guards; the niche is the exact failure of a power quotient;
the wildcard is the prime-power counting operator in the
[prime-race companion](prime_power_bias_20261014d.md). That companion owns the
counts involving 1473 and 26861. Those numbers are not Collatz states or word
lengths in this packet.

The initial correction is elementary but necessary: squaring modulo 3 sends
$0$ to $0$, and sends each of $1,2$ to $1$. Thus it is false on all integers
that squaring produces residue 1. The exceptional prime 3 must be retained.
The cheap hostile to dropping sign is $3$ versus $-3$; the more consequential
hostile below is an actual source 67 and its proposed child 33.

## 2. Power maps redistribute counts; they do not correct the source

For a finite modulus $m$ and source weights $w_x$, the exact power pushforward is

\[
 (P_d w)_y=\sum_{x^d\equiv y\pmod m}w_x.                    \tag{1}
\]

It preserves total weight. Empty fibres contribute zero; a target with two
preimages receives both weights. Assigning one unit of weight to each distinct
image instead would silently discard multiplicity. The implementation retains
all source residues, including nonunits.

For an odd prime $p$, write $\chi_p$ for the quadratic character, with
$\chi_p(0)=0$. Under uniform counting weights, the square fibre has size

\[
 \#\{x:x^2=y\}=
 \begin{cases}1,&y=0,\\1+\chi_p(y),&y\ne0.\end{cases}       \tag{2}
\]

Indeed a nonzero square has exactly the two distinct roots $x,-x$; a nonsquare
has none. More generally the standard cyclic description of the finite-field
unit group gives fibre size $\gcd(2^j,p-1)$ on the image of $x\mapsto x^{2^j}$.
There are $(p-1)/\gcd(2^j,p-1)$ nonzero image residues. If $p-1=2^s m$ with
$m$ odd, repeated squaring eventually maps onto the subgroup of order $m$;
it is not generally a map to the single residue 1.

For example the square pushforwards modulo 3 and 4 are respectively
$(1,2,0)$ and $(2,2,0,0)$. On integers there is also the precise predicate

\[
                 n^2\equiv1\pmod {24}\iff\gcd(n,6)=1.       \tag{3}
\]

An odd square is 1 modulo 8, and a square is 1 modulo 3 exactly for a 3-unit.
Conversely the displayed square congruence forces both conditions. Thus every
prime $p\ge5$ passes (3), but so does the composite 35; the prime 3 instead has
$3^2\equiv9\pmod {24}$. This is a coprimality test, not a primality test or a
convergence certificate.

There is a distinct idempotent residue selector. Choose a marked nonresidue
$\nu$ modulo $p$ and set

\[
 C_\nu(x)=\begin{cases}\nu x,&\chi_p(x)=-1,\\x,&\chi_p(x)\ge0.\end{cases} \tag{4}
\]

It fixes zero and every square; each nonzero square has two preimages under
this selector. If $p\equiv3\pmod4$, one may take $\nu=-1$. Then $C_{-1}$
chooses the square from each pair $\{x,-x\}$, and

\[
                  C_{-1}(x)=(x^2)^{(p+1)/4}\pmod p.          \tag{5}
\]

The right side is the unique square square-root of $x^2$: the square subgroup
has odd order and $2(p+1)/4\equiv1\pmod{(p-1)/2}$. If $p\equiv1\pmod4$,
negation preserves the quadratic character and cannot repair a nonresidue.
For instance both $2$ and $-2$ are nonsquares modulo 5. Formula (4) remains
valid with another marked nonresidue, but changes the source residue. It does
not preserve an arbitrary affine Collatz guard.

## 3. Sign is necessary for deterministic signed odd dynamics

Let $U_+(n)=(3n+1)/2^{v_2(3n+1)}$ for nonzero signed odd integers, using the
valuation of the absolute numerator. The sign is invariant. For $x>0$ odd,

\[
 U_+(-x)=-U_-(x),\qquad U_-(x)=\frac{3x-1}{2^{v_2(3x-1)}}. \tag{6}
\]

This conjugates the negative component to a different positive map. It does
not identify that component with the positive $3x+1$ dynamics.

**Sharp square-quotient obstruction.** For every positive odd $x>1$,

\[
                    U_+(x)^2\ne U_+(-x)^2.                 \tag{7}
\]

If the magnitudes were the same positive odd $y$, then
$3x+1=2^a y$ and $3x-1=2^b y$. Their difference gives
$2=(2^a-2^b)y$. Therefore $y=1$, $a=2$, $b=1$ and $x=1$.
Conversely $x=1$ does give equal squares. In particular the identical square
9 has next squares 25 and 1 when its original sources are 3 and $-3$.
There is no deterministic squaring quotient of the full signed map away
from this single magnitude.

The lossless repair is simply the marked record $(n^2,\operatorname{sgn}n)$.
The checker uses an exact integer-square test, reconstructs the signed odd
source, and requires equality with the supplied source before taking a step.
It also consumes positive-source ROOT words separately: $3\xrightarrow{(1,4)}1$
is accepted, ROOT itself has the empty word, while a certificate at another
source or a padded ROOT word is rejected. A signed record at $-1$ is not a
positive ROOT certificate.

There is a second small loss: every actual signed step with valuation $a$
has endpoint

\[
                         U_+(n)\equiv(-1)^a\pmod3.         \tag{8}
\]

So the endpoint residue retains the parity of the valuation; its square is
always 1 modulo 3 and forgets that parity. This is true *after an actual
Collatz step*, because $3n+1$ is a 3-unit. It does not justify sending a
source multiple of 3 to residue 1 merely by squaring it.

The same branch loss appears at the two anchors. For $E\ge3$, both
$(2^E-1)^2$ and $(2^E+1)^2$ are 1 modulo $2^{E+1}$, while their first
$3x+1$ valuations are respectively 1 and 2. The square observer merges the
$-1$ and $+1$ neighbourhoods that the inherited Mersenne rules distinguish.

## 4. Complete classification of false squared-anchor collisions

For a positive valuation word $w$, write $F_w(x)=(3^\ell x+B)/2^A$.
Call it reduced when its first letter is at least 2. If its proper-prefix
partial sums minus one are $s=(e_1<\cdots<e_m)$, direct composition gives

\[
 F_w(-1)=\frac{K(s)}{2^{A-1}},\qquad
 K(s)=-3^m+\sum_{j=1}^m3^{m-j}2^{e_j}.                    \tag{9}
\]

The numerator is odd. Removing a possible leading exponent 1 leaves its
value unchanged. For a nonempty canonical sequence, $e_j\ge j+1$, and hence

\[
 K(s)\ge3^{m+1}-2^{m+2}.
\]

This is positive for $m\ge1$, and at least 11 for $m\ge2$. Consequently
$K=-1$ occurs exactly at $s=()$ and $(1)$; $K=1$ occurs exactly at $(2)$
and $(1,2)$. These elementary cases suffice for the following complete result.

**Theorem.** If two reduced positive words satisfy
$F_u(-1)^2=F_v(-1)^2$ but $F_u(-1)\ne F_v(-1)$, they have a common total
$A\ge4$, and one word is in each of the following two lists:

\[
 \begin{array}{c|c}
 F(-1)=-2^{1-A}&(A),\ (2,A-2)\\
 F(-1)=+2^{1-A}&(3,A-3),\ (2,1,A-3).
 \end{array}                                             \tag{10}
\]

Conversely every cross-pair in (10) is such a false collision. Equality of
squared fractions with odd numerators forces equal denominators and equal
absolute numerators. Opposite signs then force the only negative possibility,
$K=-1$, and the classification of $K=1$ above proves completeness.
Thus retaining the sign alongside the squared anchor value exactly repairs
the equality test. It does not recover the word: genuine anchor collisions
still exist and the ordered native word must still be retained.

These are rational anchor evaluations, not actual signed orbits starting at
$-1$. The inherited uniform-switch compiler uses equality at that anchor to
prove an affine identity; it must not receive a squared equality instead.

Here is the precise failure in its source interface. Choose $A\ge4$ and $r\ge1$.
The falsely paired words $(A)$ and $(3,A-3)$ have length difference one.
For the prescribed child $h=(n+1)/2-1$, the complete paths obey

\[
 F_{1^{r-1},3,A-3}(h)-F_{1^r,A}(n)=2^{2-A}.              \tag{11}
\]

Their slopes agree, but the retained anchor signs leave the nonzero constant
on the right. It lies strictly between 0 and 1. Therefore the two endpoints
cannot both be integers: their native guards are disjoint under this child
map. This is stronger than a failed numerical sample.

For every source word $w$ of total cost $C$, its exact valuation cylinder is
$n\equiv(2^C-B)(3^{|w|})^{-1}\pmod{2^{C+1}}$.
The hostile constructor selects a positive member above $2^C$ of the source
cylinder for $w=(1^r,A)$, replays that source, and displays the nonintegral
partner endpoint. The source bound excludes an earlier ROOT, since an inverse
prefix from ROOT would have source less than its power-of-two denominator.
No source is replaced merely to make the partner guard hold.

At $A=4,r=1$, source 67 has the actual word $(1,4)$ ending at 19.
Its proposed child 33 has formal partner value $77/4$ under $(3,1)$;
already $v_2(3\cdot33+1)=2$, not 3. A Mersenne example is
$63\xrightarrow{1^5,4}91$: proposed child 31 has formal partner endpoint
$365/4$ under $1^4,3,1$. Neither squared equality supplies a common future.

## 5. Preserved predicate, loss, and the next obligation

| Object and map | Preserved | Lost or required |
|---|---|---|
| Residue weights under $x\mapsto x^d$ | Exact total mass with all fibre multiplicities | Original residue within a fibre |
| Nonresidue selector (4) | Membership in the chosen square image | The unchanged original source |
| Signed odd source to $n^2$ | Magnitude | Sign, next state and terminal component |
| Signed source to marked square | Full source and actual step | No computational speedup or new ROOT proof |
| Anchor value to its square | Absolute rational value | Sign; precisely the false pairs (10) |
| Anchor value with sign and retained word | Correct affine collision test | Native source and child ROOT guards remain separate |

The root-square branch repair is classical; the specific contribution is its
complete failure family at the current anchor compiler, with an exact
source-owned native-guard contradiction. This applies the repository's
**retain local incidence through counting** and **verify the consequence,
not merely the model's assumptions** research moves. No new general
meta-pattern is proposed.

The next constructive requirement remains an actual source-native common
future and, when completion is claimed, authenticated terminal evidence.
Neither a constant square residue nor an idempotent correction operation
supplies those missing obligations. The prime-power Möbius inversion in the
companion is a valid count-recovery operation; it is not an orbit repair.

## 6. Reproduction

Run `python -B -X utf8 04-computation/experiments/residue_power_defects_20261014d.py`
and repeat with `-O`. The explicit universe comprises odd primes through 101,
six square-iteration levels, all integers from $-1001$ through 1001 for (3),
all odd magnitudes below 2000 for the signed observer, all 32,767 reduced
positive words of total at most 16, and 135 native false-switch family controls.
The word census has 104 ordered false pairs, exactly eight per total $A=4,\ldots,16$.
Strict source, sign, perfect-square, valuation, prime, weight and terminal
types have positive and malformed controls. The large check count includes
repeated elementary validations, not that many independent mathematical cases.

Normal and optimized execution agree with the saved JSON on **614,213 exact
checks**. LF-normalized SHA256:
`fc191ae6dcdc8b67089e70014a4befde66367b277d9ada37b630557925283d20`.
