# Collatz passage flags: compress expansions, retain residues, and separate the return clocks

2026-11-26. **PROVED**, with the cited inputs and exact rational tests below: an all-height bound on positive cycles admitting bounded expanding-packet descriptions, several explicit exclusions, and a large representative-free extension of the mod-16 descent criterion. **OPEN:** arbitrary packet counts and dictionaries, divergent orbits, and Erdős 1135 in full. No submission to the bounty or statement about independently verifying the attached paper.

## Reading and inheritance

The supplied [Uniform bounds for planar polynomial limit cycles (September 24)](https://github.com/openai/math/tree/main/preprints/uniform-bounds-for-planar-polynomial-limit-cycles-September-24-2026) *claims* uniform boundedness for planar polynomial vector fields. The transferable feature is the proof architecture, not a theorem about discrete maps. Sections 2–5 retain competing scales, exact identities, and the original domain before isolating zeros; Sections 9–11 separately need a uniformly finite passage description and use exact endpoint matching to count cycles. Real asymptotic flatness by itself is explicitly **not** the separation hypothesis in the paper.

The repo already made the [continuous-tube and adelic-clock comparison](../collatz_cycles_tubes_debt_walk_openai_20261006.md), so this is not another tube construction. We use the exact native guard inherited from the [successor/holonomy compiler](../collatz_successor_holonomy_20261014b.md), the [maximal expanding run decoder](../collatz_maximal_run_reset_20261019.md), the signed anchors of the [mixed return compiler](../collatz_mixed_return_compiler_20261003.md), and the difference between accumulated gap and endpoint term from the [marked log-gap note](../collatz_log_gap_20261014a.md). The repo's [finite-lasso audit](../collatz_run_classification_scope_20261007e.md) is the hostile control: bounded denominators of rational states do not force repetition.

[Conjectures.io's Erdős 1135 page](https://conjectures.io/problems/erdos1135-erdos-1135) links the [published Lean piece](https://github.com/conjectures-io/conjectures-contribution/blob/main/contributions/erdos-1135/48c2d77b3b8320ea78c4be696dedd9ab4a2557a3b6a6234ecdc073a69658520b/script.lean). It uses the Terras affine formula on shortcut parity cylinders, then induction from actual descent, with three surviving mod-16 cells. The extension in §4 does not require enumerating representative cells of larger moduli. It is a mathematical proof with exact Python arithmetic, not yet a Lean contribution.

## 1. An exact packet can grow in real size only as fast as it consumes its binary guard

Write $U(x)=(3x+1)/2^{v_2(3x+1)}$ on positive odd integers. A nonempty word $w$ of $a$ valuations with total $b$ has real affine carrier
\[
 F_w(x)=(3^a x+B_w)/2^b,\qquad f_w=B_w/(2^b-3^a).
\]
Call a repetition $w^m$ of an expanding word ($\lambda_w=3^a/2^b>1$) an *expanding packet*. Individual valuations $e\ge2$ are the other allowed passages; their number and exponents are **not** bounded. This is a cyclic parsing of the actual valuation word into expanding packets and individual contracting steps. The letters **inside** an expanding packet may include arbitrary rises and falls.

In lowest terms $f_w=p/q$, both $p,q$ odd with $q>0$. Then
\[
 q(F_w^m(x)-f_w)=\lambda_w^m(qx-p).
\]
The exact native guard is $2^{bm+1}\mid qx-p$. Necessity follows by oddness of the endpoint and the coprime numerator $3^{am}$. For sufficiency, the guard first gives an odd integer endpoint. Apply the inverse maps from it: all intermediate rationals have odd numerators and odd denominators, hence are 2-adic odd. Starting from the given integer at the other end, a forward step therefore divides an integer by exactly the required power of two and is again an odd integer. This verifies every letter. An expanding anchor is negative, so a given positive $x$ cannot sustain infinitely many repeats of it.

Pick constants $C\ge\max_w|f_w|$ and $L\ge\max_w\log_2 q_w$, and let $Y(x)=\log_2(x+C)+L$. On a legal expansion packet, $bm\le Y(x)$ and
\[
 x'+C\le\lambda_w^m(x+C).
\]
Put $\beta=\max_w a/b$, and choose a rational $\Gamma>\beta\log_2 3$ (we use $317\beta/200$). The number of odd steps in the packet is at most $\beta Y(x)$, and $Y(x')\le\Gamma Y(x)$. If $x\ge C$, the ordinary log-height increase is at most $(\Gamma-1)Y(x)+1$. Thus a near-critical expanding return can be a *cheaper* passage than a pure valuation-one run: the coefficient is $a\log_2(3)/b$, not its raw number of steps.

## 2. A stitched length bound, with no bound on the contracting words

Use $N=2^{71}$, from [Bařina's verified convergence range](https://link.springer.com/article/10.1007/s11227-025-07337-0), and assume $C\le N$. In a nontrivial cycle all vertices are at least $N$. An individual step $e\ge2$ loses at least $\mu=83/200$ in base-two log height, because $(3x+1)/(4x)\le(3+1/N)/4<2^{-\mu}$. The same argument applies to any valid prefix staying above $N$ and ending at or above its original source. The actual trajectory, not a constructed off-source return, supplies the telescoping height inequality.

At the beginning of a packet or a contraction, cut the cycle as close to its minimum $n_*$ as follows. If the minimum occurs inside a packet, it is in the first repetition: at every fixed phase the states increase across repeats by the rotated expanding affine map. A prefix of that one copy has positive additive carry, so the start of its packet is at most $2^{b_*}n_*$ where $b_*=\max b_w$. If the minimum is already a boundary no loss is needed. Hence an admissible initial budget is
\[
 H=b_*+L+1+\log_2 n_*.
\]
For a first-descent-free prefix starting at its retained source $x$ itself, one may use $Y(x)$ instead. There is no bound on the lengths or the exponents in the contraction-only gaps.

If there are at most $r$ expanding packets and $k$ odd steps in total, summing expansion heights and paying for all contractive steps by the log height which expands gives the **packet bound**
\[
 k\le
 \left(\beta+\frac{\Gamma-1}{\mu}\right)H
 \sum_{j=0}^{r-1}\Gamma^j+\frac r\mu .
 \tag{*}
\]
Indeed the height at the start of packet $i$ is at most $\Gamma^{i-1}H$, as all intervening contractions lower the height function, and the total contraction loss is at most the total expansion gain. In a cycle the total change is zero. This is a marked-source bound, not an assumption of a random or independent valuation process.

Bounded $r$ is not ordinary bounded local minima. For $w=(1,2)$, the anchor is $-5$ and $2\cdot8^m-5$ follows $w^m$ exactly, with a rise and fall in each copy. Its boundary states are $-5+2\cdot9^j8^{m-j}$. The packet's fixed dictionary can compress arbitrarily many local extrema. This hostile case is why the result is not just a rephrasing of the classical bounded-$m$-cycle theorem.

## 3. Ruling out cyclic returns in the prepared grammar

For a cycle with $k$ odd steps and total valuation $E$,
\[
 0<\Lambda=E\log2-k\log3
 =\sum_x\log(1+1/(3x))<k/(3n_*).
\]
[Simons–de Weger, Lemma 12](https://deweger.net/papers/%5B35%5DSidW-3n%2B1-ActaArith%5B2005%5D.pdf) gives the cycle-context lower bound $\Lambda>\exp[-13.3(0.46057+\log k)]$. Thus $n_*<2^8k^{14.3}$. This is an imported arithmetic separation estimate, not an invocation of the planar paper's complex-domain theorem. Inserting it in (*) yields a computable upper bound on $k$ for each fixed library and packet count, since $k/\log k$ tends to infinity. Then $E<2k$ on a nontrivial cycle, so only finitely many words remain to test. Each word has at most one fixed point, which must pass the exact native guard and replay. This yields finite isolation inside the specified grammar.

For a specified small library one can do better. For each phase of an expanding packet, successive states grow by at least $\lambda_w$; summation of reciprocals across repeats costs at most $a_w\lambda_w/(\lambda_w-1)$ times $1/n_*$. For each contraction-only stretch, sum backwards from its smallest vertex to get less than $5/n_*$, as the individual ratio is below $4/5$. In a cycle there are at most $r$ such stretches. With $M=\max_w a_w\lambda_w/(\lambda_w-1)$ and $r\ge1$ this gives the optional bound $n_*<2^8r(M+5)k^{13.3}$. The code uses the better of the two bounds.

### Exact lower barrier and exclusions

Writing $\alpha=\log_2 3$, a cycle has
\[
 0<E/k-\alpha<1/(3N\log2).
\]
The replay bounds these logarithms by rational series with explicit remainders, then walks down adjacent Farey endpoints until the mediant first enters the slightly enlarged closed interval. Farey adjacency guarantees that no smaller denominator lies strictly between the two guards. The first relaxed admissible fraction is **114,208,327,604 / 72,057,431,991**. In particular the odd period of any hypothetical nontrivial cycle is at least the denominator. This is a modest sourced barrier, not a world-record claim; here it is used in a different grammar.

For the 32-unit library, allow **every** expanding word with total valuation $b\le32$ in the indicated slope range: no need to enumerate the library. In general $B_w\le 2^{b-a}(3^a-2^a)<3^b$ and $q_w\le3^b$, so take $C=3^{32}$ and $L=\lceil 32\cdot317/200\rceil$. At the lower period barrier the right side of (*) is already too small for each entry below. Dividing the right side by $k$ makes it decreasing from the barrier onward; each upper height used is a positive constant plus a positive multiple of $\log k$.

| Allowed expanding packets around the cycle | Excluded if there are at most |
|---|---:|
| Total valuation at most 32 per motif, any expanding slope | 37 packets |
| Same, also $a/b\le2/3$ | 289 packets |
| Same, also $a/b\le7/11$ | 1,658 packets |
| Rotations of the negative $-5$ anchor word only | 292 packets |
| Rotations of the negative $-17$ anchor word only | 1,676 packets |
| Rotations of the $-17$ word or the inherited length-22 expanding reset | 1,670 packets |

An unrestricted number of positive contracting steps remains allowed between packets. The negative cycles supply legal positive **prefixes** via their native guards; they are not themselves positive cycles. The theorem says a nontrivial positive cyclic return must escape these compressed descriptions, not that it cannot exist.

## 4. A direct upgrade to the published mod-16 descent route

Use the **shortcut** map $T(x)=x/2$ on even numbers, $(3x+1)/2$ on odds. Take a native parity path and let $p$ be the *first* time its multiplier falls below 1; let $j$ be its odd count. Every earlier prefix has multiplier at least 1, so every odd additive contribution, when normalized by the full multiplier, is at most $1/3$: at an odd step taken after $i$ shortcut steps it is $2^i/(3^{j_i+1})\le1/3$. Therefore
\[
 T^p(n)\le(3^j/2^p)(n+j/3).
\]
If the first slope drop does **not** yield strict source descent (necessarily $j>0$) from $n\ge N$, necessarily $0<p/j-\alpha\le1/(3N\log2)$, exactly the one-sided Farey interval used for cycles. Thus $j\ge72,057,431,991$. Since $p>\alpha j$, the certificate verifies $p\ge114,208,327,604$. The least counterexample, if any, cannot have a first multiplier drop within **114,208,327,603** shortcut steps: any actual descent would contradict minimality. This tests every native residue, not just one representative. It is a statement about coefficient drift before descent, not an upper bound on stopping time.

This is also why a uniform finite number of ordinary steps cannot substitute for the paper's prepared passage words. The integers $2^{L+1}-1$ have arbitrarily long initial expanding blocks. The useful open target is **source-preserving**: prove a bounded number of admissible packets pays an actual descent (or prove entry into such a grammar by a separate argument). Constructing large sources that follow a chosen word, or checking a finite family of bounded packet counts, is not universal entry.

## Replay

- **python3 bounds.py**: exact rational logarithm and Farey certificates, packet budgets, and 32 native signed-anchor controls. No floating comparisons enter the certificates.
- **python3 first_crossing.py**: checks the mod-16 examples and a finite exact cross-check of representative-free clock thresholds through shortcut depth 100,000. The full Farey result does not depend on this scan.
