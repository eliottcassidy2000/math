# Exact source floors under guarded refinement and finite deadlines

**PROVED:** exact positive-threshold decisions, uniform pointwise intervals, the sharp largest source at each threshold, and conditional floor transport through verified prefixes. **FINITE-EXACT:** the bounded controls below. **OPEN:** an independently justified positive floor for every supplied source. All ROOT receipts stop on their first visit to 1.

Artifacts: [script](../../04-computation/experiments/collatz_floor_transport_deadlines_20261005.py) and [output](collatz_floor_transport_deadlines_20261005.out).

## 1. Inheritance and the source coordinate

The weight is inherited from [adaptive mixture flow, P2](collatz_adaptive_mixture_flow_20261005.md), not introduced here. The closest underused operation is a pathwise deadline in a fair-bit extractor. [THM-2160, dyadic checksum under the critical-run deadline](../../01-canon/theorems/THM-2160-dyadic-checksum-extracts-a-fair-bit-under-the-critical-run-deadline.md), [THM-2225, dyadic critical-run extractors](../../01-canon/theorems/THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection.md), and [THM-2253, online dyadic contrast](../../01-canon/theorems/THM-2253-online-dyadic-contrast-tournament-extractor.md) separate two obligations: preserve the relevant measure under an involution, and prove that the entire observation episode stops within its promised deadline.

The transferable part here is a **whole-episode resource bound**, together with a retained source and guarded word. Bernoulli exchangeability does not itself preserve Collatz guards. The corrected near miss is “each stage terminates, therefore there is one positive floor through arbitrarily many stages.” The explicit hostile below changes the source through inverse siblings and drives its weight to zero.

The concept board is **fixed source / endpoint floor / actual valuation word / aggregate counter bill / deadline / zero-support boundary**. The least-used sidecar is the original floor and cumulative bill, instead of repeatedly compressing them into one new scalar. This package supplies no new source basin. The companion [finite threshold compiler](collatz_weight_threshold_receipts_20261005.md) enumerates the whole finite superlevel, whereas the API here decides the threshold for one supplied source by bounded forward steps.

## 2. Hidden counters give a sharp source-size bound

Put \(U(n)=\operatorname{oddpart}(3n+1)\) on positive odd integers, with no outgoing first-hit edge at ROOT 1. If the actual strict ROOT word of \(n>1\) is \(a=(a_1,\ldots,a_\tau)\), define

\[
 L=\tau-1,\qquad K=\sum_i\left\lfloor\frac{a_i-1}{2}\right\rfloor,
 \qquad N=L+K,\qquad
 W(n)=w(L,K)=\frac{2K!(L+1)!}{(L+K+2)!}.
 \tag{1}
\]

The inherited definition sets \(W(n)=0\) on unrooted sources and \(W(1)=1\). No finite counters are assigned to an unrooted source. The last nonroot predecessor of 1 has an even valuation at least 4, so every nonroot receipt has \(K\ge1\). In particular \(\tau\le N\).

**F1 — PROVED.** Every rooted nonroot source satisfies

\[
 W(n)\le\frac{2}{(N+1)(N+2)},\qquad
 n\le\frac{4^{N+1}-1}{3}=S^N(1),\quad S(x)=4x+1.
 \tag{2}
\]

For the first inequality,

\[
 w(L,K)=\frac{2}{(N+2)\binom{N+1}{K}},
\]

and \(1\le K\le N\) implies \(\binom{N+1}{K}\ge N+1\). For the second, let \(A=\sum_i a_i\). Each \(a_i\le2\lfloor(a_i-1)/2\rfloor+2\), hence \(A\le2(N+1)\). The exact affine ROOT identity is

\[
 2^A=3^\tau n+C,\qquad C\ge1.
\]

Consequently \(3n+1\le4^{N+1}\). Both bounds in (2) are simultaneously attained by \(S^N(1)\): for \(N\ge1\) its strict word has the single valuation \(2N+2\), and its counters are \((0,N)\). ROOT supplies the \(N=0\) case separately. The source maximum is sharp; this does not assert that the bound \(\tau\le N\) is attained at every \(N\).

For rational \(\epsilon=p/q\in(0,1]\), define

\[
 B=\max\{b\in\mathbb N:\epsilon(b+1)(b+2)\le2\}.
 \tag{3}
\]

It is computed with integer arithmetic, not a floating square root:

\[
 M=\left\lfloor\frac{2q}{p}\right\rfloor,\qquad
 B=\left\lfloor\frac{\operatorname{isqrt}(1+4M)-3}{2}\right\rfloor.
 \tag{4}
\]

**F2 — PROVED.** The set \(\{n>0\text{ odd}:W(n)\ge\epsilon\}\) is finite, and its largest integer is exactly

\[
 \boxed{\ \frac{4^{B+1}-1}{3}\ }.
 \tag{5}
\]

Indeed (2) excludes larger sources, including unrooted sources because their weight is zero. The root-ray witness in (2) belongs to the superlevel by (3). Every nonroot member reaches ROOT in at most \(B\) odd steps. The witness at the largest source has just one step when \(B\ge1\); these are different extremal questions.

## 3. A total threshold decision and a pointwise enclosure

**F3 — PROVED.** After following an actual source for at most \(T\) odd steps, stop immediately if ROOT is seen and evaluate (1) exactly. If ROOT has not appeared after those steps, then

\[
 0\le W(n)\le\frac{2}{(T+2)(T+3)}.
 \tag{6}
\]

If the source eventually reaches ROOT, its total rank is at least \(T+1\), hence \(N\ge T+1\); apply (2). If it never reaches ROOT, the definition gives weight zero. This is an unconditional, uniform pointwise interval and an always-terminating approximation algorithm. To reach absolute error at most \(\epsilon\), choose an integer \(T\) satisfying \(\epsilon(T+2)(T+3)\ge2\). The number of odd steps is \(O(\epsilon^{-1/2})\); this statement does not bound intermediate integer bit lengths or bit-operation cost.

Positive-threshold **comparison**, including equality, is also total. Follow at most \(B\) actual steps, where \(B\) is (3). If ROOT appears, compare its exact weight with \(\epsilon\). Otherwise (6) gives

\[
 W(n)\le\frac{2}{(B+2)(B+3)}<\epsilon.
 \tag{7}
\]

Thus a negative threshold decision does not classify the source as unrooted; it may reach ROOT later with smaller positive weight. Positive decisions return the actual strict ROOT word. The endpoint case \(B=0\) needs no edge: ROOT qualifies, and every other source is rejected.

The exported function threshold_receipt(source, epsilon) implements this decision. The function pointwise_interval(source, steps) implements (6). They accept exact positive odd integer sources and exact rational thresholds; floats and booleans are rejected. The returned reached_root flag reports only whether the bounded replay reached1; false is not a nonconvergence verdict. The helper word_weight is only a formula for a given valuation tuple, not an independent validity certificate; exported guarded operations replay the tuple at the retained source.

The missing global step remains an initial positive lower bound for every source. Computable arbitrarily narrow intervals with lower endpoint zero do not supply it. Conversely, an independently proved \(\epsilon\le W(n)\) turns the finite deadline into a ROOT receipt for that exact integer.

## 4. Conditional transport through a verified prefix

Suppose an actual guarded prefix sends \(x\) to a nonroot endpoint \(y\), using \(r\) odd steps and

\[
 k=\sum_{\text{prefix}}\left\lfloor\frac{a_i-1}{2}\right\rfloor.
\]

If \(W(y)\ge\epsilon>0\), its suffix is rooted, so its counters \((L,K)\) exist. They satisfy \(K\ge1\), \(N=L+K\le B\), and necessarily \(\epsilon\le1/3\). The full source has counters \((L+r,K+k)\). Write \((z)_j=z(z+1)\cdots(z+j-1)\), with \((z)_0=1\).

**F4 — PROVED, conditional on the endpoint floor and actual prefix.**

\[
 \frac{W(x)}{W(y)}
  =\frac{(K+1)_k(L+2)_r}{(N+3)_{r+k}}
  \ge\frac{(k+1)!(r+1)!}{(B+3)_{r+k}}.
\]

Consequently

\[
 \boxed{\ W(x)\ge
 \epsilon\,\frac{(k+1)!(r+1)!}{(B+3)_{r+k}}\ }.
 \tag{8}
\]

The equality is factorial cancellation; the inequality uses \(K\ge1\), \(L\ge0\), \(N\le B\). It retains the actual endpoint, not just an endpoint residue or a formally matching counter. If the prefix ends at ROOT, replay already gives an unconditional completed receipt, and its exact weight replaces (8). Prefixes that continue beyond ROOT are rejected.

The formal finite improvement is

\[
 \min_{\substack{K\ge1,\ L+K\le B\\w(L,K)\ge\epsilon}}
      w(L+r,K+k).
 \tag{9}
\]

It is a safe bound over a finite superset of realizable counters. It is not asserted to be the optimal bound over actual integer sources.

The structured record also supplies a better deadline directly: retain the endpoint floor and the prefix word, then ROOT is reached from \(x\) within \(r+B\) odd steps. There is no need to discard that information and recompute a potentially much larger deadline from the smaller scalar (8).

## 5. The entire refinement episode must have a bounded bill

If an episode has a verified aggregate increment \(r+k\le D\), then every output source in that episode obeys the common conditional floor

\[
 \boxed{
 \epsilon\,
 \frac{(\lfloor D/2\rfloor+1)!(\lceil D/2\rceil+1)!}
      {(B+3)_D}\ }.
 \tag{10}
\]

For fixed \(r+k\), the factorial product is minimized when the two arguments differ by at most one. These minimized bounds decrease as the total increment increases: adding to a smaller argument multiplies the expression by at most \((d+2)/(B+d+3)<1\). This proves (10) for every bill at most \(D\).

For a completed unary valuation code \(0^{a-1}1\), the bit length is \(a\), whereas that valuation contributes

\[
 1+\left\lfloor\frac{a-1}{2}\right\rfloor
 =\left\lceil\frac a2\right\rceil\le a
\]

to \(r+k\). Thus a finite deadline for the **whole completed code** yields a sufficient budget \(D\). A rule that merely makes each next stage finite provides no uniform total \(D\).

Retaining the original floor and aggregate bill can be materially stronger than repeated scalar recompression. At \(\epsilon=1/3\), two steps with total increment \((r,k)=(2,0)\) give \(1/10\) in (8); applying (8) separately twice and rebuilding \(B\) from the intermediate scalar gives \(1/15\). This example is a comparison of valid envelopes; it is not a claim that every formal intermediate counter state is realized.

There are two different meanings of refinement:

1. **Refining information about the same integer.** A proved lower bound remains true. Keep it, or take the maximum of independently valid lower bounds. No mass is lost merely by learning more digits.
2. **Changing the integer by adding an actual inverse prefix.** The new source may have a smaller weight. Bound the cumulative bill and use (8)-(10), or retain the resulting vanishing floor explicitly.

The second operation has an exact hostile. For \(j\ge0\),

\[
 n_j=S^j(3),\qquad
 n_j\xrightarrow{\,1+2j\,}5\xrightarrow{\,4\,}1,
 \qquad
 W(n_j)=w(1,j+1)=\frac{4}{(j+2)(j+3)(j+4)}\longrightarrow0.
 \tag{11}
\]

Every individual stage and receipt is finite. There is no positive lower bound uniform in this changing-source sequence.

## 6. The fair-coin symmetry is not an arithmetic guard symmetry

The words \((1,2)\) and \((2,1)\) have unary codes 101 and 011, with equal counts of zeros and ones. Their affine maps are nevertheless

\[
 F_{12}(x)=\frac{9x+5}{8},\qquad
 F_{21}(x)=\frac{9x+7}{8}.
\]

An integer inverse at a fixed endpoint \(y\) requires respectively \(y\equiv4\pmod9\) and \(y\equiv2\pmod9\). No integer endpoint supports both. Positive actual controls are \(11\to17\to13\) and \(9\to7\to11\). A composition-preserving Bernoulli permutation therefore cannot simply be used to exchange these inverse histories while retaining their arithmetic endpoint.

The exact transfer here has source **a verified endpoint floor plus guarded prefix**, target **a lower bound and finite ROOT deadline at the original source**, map **counter addition with factorial transport**, preserved predicate **actual first-hit suffix and source identity**, lost information in the scalar bound **exact suffix counters and order**, and required sidecar **the original source, valuation word, and endpoint-floor proof**. No fair-coin construction is claimed to furnish the initial Collatz floor.

## 7. Reproduction and finite scope

Run from the repository root:

~~~text
python 04-computation/experiments/collatz_floor_transport_deadlines_20261005.py
python -O 04-computation/experiments/collatz_floor_transport_deadlines_20261005.py
~~~

The verifier uses exact integers and Fraction, explicit checks that survive optimization, and independently completed finite control routes. The universe is all 512 positive odd sources below 1024 (control cap 512 odd steps), eight thresholds, seven pointwise horizons, prefix cuts up to four actual steps at nonroot endpoints, and 100 formal transport envelopes. There are 4,096 threshold controls, 3,584 pointwise intervals, 1,966 guarded prefixes, seven unbounded-family samples, and 12 malformed-type or invalid-guard controls. The theoretical infinite claims are proved above; finite success is not used to infer universal convergence.

Normal and optimized runs agree: **32,229 exact checks**. The saved output records all eight sharp threshold extrema, the guard-exchange hostile, aggregate versus scalar floors, and the changing-source examples. The polynomial-dual companion can provide a truthful atom floor from certified moment intervals; its truth premise must be retained. The finite compiler companion can turn such a floor into a source-labelled receipt. Neither interface manufactures positivity from unverified moment data or from approximation alone.
