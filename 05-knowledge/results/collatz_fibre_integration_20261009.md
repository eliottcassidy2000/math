# What a Mersenne deletion pays: two clock defects and terminal obligations

**PROVED:** the source-authenticated defect ledger, its composition law, and conditional terminal transport below. **FINITE-EXACT:** the declared small rooted universe and label accounting. **INHERITED FINITE-EXACT:** the 1,949-member giant fan of THM-4605. **OPEN:** ROOT for that fan and universal Collatz.

Artifacts: [checker](../../04-computation/experiments/collatz_fibre_integration_20261009.py) and [output](collatz_fibre_integration_20261009.out).

## 1. Inheritance and the changed target

The governing source is [THM-4605, Mersenne exits by pair-chain absorption](../../01-canon/theorems/THM-4605-mersenne-exits-by-pair-chain-absorption-the-t23-fan-escapes-at-depth-1.9e7.md), especially statements 2 and 8, with the proof in [Mersenne-line barriers, §4.1](mersenne_line_barriers_20261007.md). [MISTAKE-588](../../01-canon/MISTAKES.md) repairs the noninterval absorbed set, the stopping-time scope, and the trivial-cycle exception. These corrections take precedence over our earlier [parameter-23 packet](collatz_residual_t23_20261007f.md).

That old packet emitted 24 ungrounded Mersenne children. The later computation joins them into a class of **1,949 exponents**, down to $99,708,991,741$. It supplies additional genuine receipts, not a ROOT suffix. Its escape set below the old bottom is

\[
\{1,\ldots,1910\}\cup\{1929,1930,1931,1932,1935,1936\}.
\]

There is no interval between its minimum and maximum that may safely replace this set. The separate level-2 cluster is merely unabsorbed through its declared horizon; it is not proved permanently separate. The new checker verifies the 1,949 distinct labels and the missing offset 1911. It does **not** rerun the 19,000,765-step absorption or the level-2 experiment.

The anchor here is terminal grounding, the niche is the exact relation between the two orbit clocks, and the wildcard is composing receipts with nonzero clock defects. The retained coordinates are the actual source, actual child, two prefix words, common endpoint, both clock defects, and ROOT status. The inherited [clock-holonomy compiler](collatz_clock_holonomy_20261007.md) supplies the seam-composition mechanism; the present ledger identifies what that seam changes in a Mersenne fibre.

## 2. The invariant preserved by the old deletion grammar

Use the Terras map

\[
T(n)=\begin{cases}(3n+1)/2,&n\text{ odd},\\n/2,&n\text{ even},\end{cases}
\qquad M_E=2^E-1.
\]

After $E-1$ initial odd steps, $M_E$ is at $x_E=2\cdot3^{E-1}-1$. For an exponent whose orbit reaches ROOT, define

\[
s(E)=\sigma_T(M_E)-(E-1),\qquad
o(E)=\text{number of odd Terras steps before first ROOT}.
\]

The exact equality classes of $(s,o)$ are the same as those of $(\sigma_T(M_E)-E,o)$; the first coordinate differs by 1.

THM-4605 concerns absorption of the normalized pair chain comparing $x_E$ and $x_K$ at equal Terras times. At a merge value at least 3, such an absorption preserves both $s$ and $o$. Conversely, **when the two ROOT times are finite**, equality of both invariants gives absorption at the shared value 8 no later than residual time $s-3$. Equality of one coordinate alone is not the theorem. The converse does not establish finite ROOT time.

The initial ascent accounts for the entire clock difference $E-K$. Thus a long sequence of these deletions may substantially lower the exponent while leaving the remaining stopping-time problem exactly unchanged. It consolidates equivalent obligations and their reusable suffixes. It does not, by itself, cross to an independently grounded fibre.

## 3. A computable defect for every authenticated pre-ROOT join

Let $E>K\ge2$. A receipt consists of actual finite parity words $u,v$ with

\[
T^{|u|}(M_E)=z=T^{|v|}(M_K),\qquad z\ge3,
\]

and no earlier ROOT in either prefix. Let $o(u)$ count the 1-bits of $u$. Define

\[
\boxed{\delta(u,v)=\bigl(|u|-|v|-(E-K),\ o(u)-o(v)\bigr).} \tag{1}
\]

This pair is computable from the checked prefixes even when both ROOT obligations remain open. If either orbit reaches ROOT, both do, and cancellation of their common suffix gives

\[
\boxed{\delta(u,v)=\bigl(s(E)-s(K),\ o(E)-o(K)\bigr).} \tag{2}
\]

When both times are infinite, (1) remains meaningful; (2) is not an operation on $\infty-\infty$.

For odd-step valuation words $a,b$, the same formula is

\[
\delta=\bigl(\operatorname{cost}(a)-\operatorname{cost}(b)-(E-K),\ |a|-|b|\bigr),
\]

because a valuation $a_i$ expands into the Terras parity block $1\,0^{a_i-1}$. Our old compensated deletion heads have equal full odd-word lengths and valuation-cost difference $E-K$, so their defect is $(0,0)$.

**A decisive distinction.** The integers $7=M_3$ and $3=M_2$ have the actual join

\[
7\xrightarrow{(1,1,2,3)}5\xleftarrow{(1)}3.
\]

Its defect is $(7-1-1,4-1)=(5,3)$. Thus an arbitrary smaller-Mersenne common-future receipt need not preserve the normalized deletion fibre. Attaching the supplied child suffix $(4)$ gives the actual ROOT proof $(1,1,2,3,4)$ of 7. This example is a scope boundary, not newly discovered convergence.

Accordingly, the inherited phrase “deletion routes cannot ground the giant children” must be read with its stated **normalized absorption / zero-defect** hypothesis. It does not rule out a different authenticated Mersenne-to-Mersenne join with nonzero defect, a non-Mersenne child, or a direct symbolic ROOT proof. Any such proposed exit still owes its exact source guards and terminal evidence.

## 4. Defects add under actual clock alignment

Suppose one receipt joins $M_E$ to $M_K$, using prefix lengths $a,b$, and another joins $M_K$ to $M_J$, using lengths $c,d$. Compare the two prefixes of the **same actual middle orbit**. Extend the shorter prefix to time $\max(b,c)$. The composed lengths are

\[
\left(a+\max(c-b,0),\ d+\max(b-c,0)\right). \tag{3}
\]

The extra parity block is already supplied by the longer middle prefix. Hence the new receipt is authenticated without discovering a new orbit segment. Its endpoint is still pre-ROOT: both supplied middle prefixes end before ROOT. Equation (3) gives addition of the time difference. Counting the same shared middle parity block gives addition of the odd difference. Subtracting the telescoping exponent differences yields

\[
\boxed{\delta(r_1\circ r_2)=\delta(r_1)+\delta(r_2).} \tag{4}
\]

The checker implements this as `compose`. It validates the middle source and every parity seam. `discharge` consumes a supplied strict first-hit child word, checks its prefix against the receipt, and replaces that prefix by the source prefix. It performs no ROOT search. A cycle of assumed obligations has no terminal word to supply to this API.

The production checker has a declared literal exponent cap of 4096. Its finite control discovery is confined to `main()`. No giant integer or giant ROOT word is materialized here.

## 5. An explicit necessary clock bill for the giant fan

One Terras step decreases a positive state by at most a factor of 2. Thus, whenever $s(E)$ is finite,

\[
s(E)\ge\log_2 x_E>(E-1)\log_2 3>\frac{19(E-1)}{12},
\]

where the final strict inequality follows from the exact integer check $3^{12}>2^{19}$.

Consider a route of arbitrary authenticated pre-ROOT receipts from $M_E$ to an independently grounded bank whose maximum residual time is $S$. From (2) and (4), the **total first defect coordinate** must satisfy

\[
\boxed{\sum_i\delta_1(r_i)\ge
\left\lfloor\frac{19(E-1)}{12}\right\rfloor+1-S.} \tag{5}
\]

This is a necessary clock-offset condition, not a sufficient ROOT criterion and not a lower bound on running time or proof-description length. A compressed symbolic rule can carry a large clock offset.

The independently recorded bank in THM-4605 has $K\le12800$ and $S=97982$. For the deepest new fan member $E=99,708,991,741$, (5) requires a net time defect of at least

\[
\boxed{157,872,472,274.}
\]

Arbitrarily many zero-defect deletions contribute zero toward this bill. A single authenticated ROOT proof for any member of the giant class would still ground the entire class by suffix transport. Neither the class size nor this bound supplies that proof.

## 6. What the no-meeting lemma does and does not exclude

THM-4605 statement 2 has the explicit height condition

\[
E\ge(1+\log_3 4)N+D+2.
\]

Under it, absence of normalized pair-chain absorption through $N$ rules out **all** meetings of the two post-run orbits at times $m,m'\le N$, including asynchronous ones. At these heights the coefficient comparison prevents a distinct-time accidental equality. Thus merely changing the clock alignment inside the same bounded window cannot rescue a negative result satisfying this condition.

This is not an all-time no-meeting theorem. It leaves a longer window, a different child, a different source representation, and a later nonzero-defect receipt to investigate. The tiny $M_3/M_2$ example lies outside the large-height hypothesis and is not a counterexample to it. The finite level-2 negative must likewise retain its exact horizon and deletion range.

## 7. Exact controls and the next target

The bounded control universe is $2\le E\le96$, with a declared discovery cap of 100,000 Terras steps per source. All 95 sources finish; their words are then treated as supplied evidence by the production APIs. There are 22 distinct invariant pairs, the largest containing 14 exponents. All 4,465 smaller-exponent pairs produce authenticated pre-ROOT joins: 327 have zero defect and 4,138 cross fibres. In every pair, the two defect coordinates agree with the separately replayed ROOT invariants, and suffix transport reproduces the exact source word. Another 435 compositions verify the actual middle-prefix alignment.

These counts are a control of the ledger, not new giant-source coverage. The zero-defect quotient compresses repeated terminal obligations. A terminal outside that quotient still needs a supplied ROOT word or a proved rule that actually changes the defect and reaches grounded evidence.

The next focused target is therefore a **source-preserving exit across these fibres**, with a complete clock ledger and a genuine terminal. The [grounded ordinary-source factory](collatz_grounded_child_ports_20261007f.md) supplies such terminals only when the original source passes its exact membership equation. Its [fixed-head S-unit boundary](collatz_terminal_head_finiteness_20261007f.md) explains why a finite fixed-middle-head library cannot silently close every Mersenne source.

```text
python -B 04-computation/experiments/collatz_fibre_integration_20261009.py
python -B -O 04-computation/experiments/collatz_fibre_integration_20261009.py
```

The checks use explicit failures rather than assertions. Malformed types, corrupted parity words, wrong child suffixes, ROOT-padding joins, mismatched seams, and precision-cap overruns are rejected.

Normal and optimized output agree: **40,094,204 primitive exact checks** (including each authenticated parity step). LF SHA256: `f43503ddb400b79ef37f74a70776ee4fdfe30edb37554728cfb750c2917b0cb1`.

Independent peer proof/API audit and optimized replay passed with the same saved output. The audit includes normalized versus raw clocks, both composition defects, strict suffix transport, and the conditional interpretation of the large clock bill.
