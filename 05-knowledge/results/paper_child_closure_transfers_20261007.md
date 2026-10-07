# Exact child retraction, retained guards, and receiver payment

**Status:** PROVED elementary statements below; CONDITIONAL transport of supplied ROOT receipts; FINITE-EXACT controls. The four supplied paper headlines are accepted premises, as requested. Their mechanisms motivate the interfaces; none is used to assert Collatz convergence or to supply a missing ROOT receipt.

## 1. Inheritance and the four paper interfaces

The closest arithmetic mechanism is the guarded negative-five lift in [the collision compiler](collatz_collision_dp_20261007.md), summarized in [the reset-two synthesis](collatz_reset2_synthesis_20261007.md). It already gives

\[
G(m)=\frac{8m-5}{9},\qquad U^{(1,2)}(G(m))=m,
\]

on its positive odd native domain. Here \(U(n)=(3n+1)/2^{v_2(3n+1)}\), and a word records **actual valuations**, in chronological order. The old root-preserving routing interface and [finite graph grounding](paper_local_global_grounding_20261007.md) keep a supplied terminal receipt separate from an ungrounded common future. The present addition is a total maximal-prefix decoder, a finite-chain partition, and a parameter-space receiver compiler. It does not rebrand the inherited negative-five recurrence as new.

The canonical hostile is the signed cycle \(-5\to-7\to-5\) with word \((1,2)\): the shift \(n+5\) vanishes there, so no finite valuation budget exists. The corrected near miss is to call a decrease from a changed core a decrease from the original source. The least-used sidecars retained here are the original source, the exact block depth, and the nonunit part of a receiver congruence.

| Supplied source and inspected interface | Retained mechanism | Concrete transfer and limit |
|---|---|---|
| [Hilbert–Smith, paper20](C:/Users/Eliott/Downloads/paper-20.pdf), 46 pages; pp.3–5,19–20,39–40; Figure1 p.5 inspected visually | A fixed integral lattice is chosen before the refinement parameter; rational equality alone does not bound its denominators. The compact source remains marked. | For one fixed positive core \(b\), the nonzero integer \(b+5\) bounds inverse depth by its ternary valuation. Changing \(b\) changes this bound; there is no all-source depth bound. |
| [Grothendieck homotopy hypothesis, paper19](C:/Users/Eliott/Downloads/paper-19.pdf), 26 pages; pp.1–10,16–18; exact-gluing diagram p.18 inspected visually | Retractions alone are insufficient; exact boundary data and their joining witnesses are retained. Empty sets of witnesses are allowed. | Prefix expansion and retraction give a bijection of strict ROOT-receipt sets, including when both are empty. A formally reversible prefix therefore creates no terminal witness. |
| [Integer Lipschitz heights with two-arc boundary, paper18](C:/Users/Eliott/Downloads/paper-18.pdf), 62 pages; pp.4–8,12–19 | Exact cuts retain internal pins/equality constraints; rare-event comparisons retain their conditioning. The absolute height needs an anchor in addition to local differences. | Receiver guards are intersected with the existing prefix, including their ternary gcd obstruction. Source comparison retains the original integer anchor. We import no stochastic independence claim. |
| [Balanced six-vertex GFF](C:/Users/Eliott/Downloads/The-Gaussian-free-field-limit-of-the-balanced-six-vertex-model-with-variance-multiplier-1-over-arcsin-c-over-2-September-23-2026.pdf), 93 pages; pp.13–16,22–26,30–31 | Fixed-word formulas retain ordered marks; bounds for a fixed word are explicitly not uniform over unbounded word length. Sparse boundary modes require separate control. | A fixed finite bank of forward receivers has an exact depth limit for source payment. Deeper prefixes require an adaptive receiver or a different authenticated rule. This is proved arithmetically below, not inferred from the field-limit theorem. |

These are construction-level transfers with different mathematical domains. In particular, no homotopy equivalence, spin law, or GFF normalization is being identified with the integer Collatz map.

## 2. A total maximal-\(12\) decoder

Write \(12=(1,2)\). For every positive odd integer \(n\), define

\[
q(n)=\left\lfloor\frac{v_2(n+5)-1}{3}\right\rfloor,
\qquad b(n)=\frac{9^{q(n)}(n+5)}{8^{q(n)}}-5.
\tag{1}
\]

**Theorem 1.** The maximal initial repeated \(12\)-word of \(n\) has exactly \(q(n)\) blocks and endpoint \(b(n)\). Every checkpoint after the source in this prefix is strictly larger than \(n\). The endpoint is a positive odd integer and \(q(b(n))=0\).

**Proof.** One actual \(12\)-block occurs exactly when \(n\equiv11\pmod {16}\), equivalently \(v_2(n+5)\ge4\). Its two maps give

\[
F(n)=\frac{9n+5}{8},\qquad F(n)+5=\frac98(n+5).
\]

The shifted binary valuation decreases by exactly three. A block remains legal while that valuation is at least four, proving (1) and maximality. Its first step increases the current state by \((n+1)/2\); the full block increases it by \((n+5)/8\). Thus every checkpoint is above the original source. In particular, no ROOT can occur inside a nonempty stripped prefix. ∎

The canonical core set is

\[
\mathcal B=\{b>0\text{ odd}:v_2(b+5)\in\{1,2,3\}\}
=\{\text{positive odd integers not congruent to }11\pmod {16}\}.
\]

For a fixed core \(b\), inverse prefixes are

\[
E_r(b)=\frac{8^r(b+5)}{9^r}-5,
\qquad0\le r\le R(b):=\left\lfloor\frac{v_3(b+5)}2\right\rfloor.
\tag{2}
\]

**Corollary 2.** All positive odd integers partition into the finite chains
\(\{E_r(b):0\le r\le R(b)\}\), indexed by \(b\in\mathcal B\). The decoder recovers the unique pair \((b,r)\).

Indeed, the divisibility in (2) is necessary and sufficient for integrality. The quotient \((b+5)/9^r\) is positive even, so the child is positive odd (at least eleven if \(r>0\)). Its shifted binary valuation is \(v_2(b+5)+3r\), which recovers exactly \(r\). The two formulas are inverse. Partial inverse composition satisfies \(E_r(E_s(b))=E_{r+s}(b)\) whenever the joint divisibility guard holds.

This is a partition by finite actual paths, **not** a proof that any core reaches ROOT. There is no uniform bound on chain length: \(b_r=2\cdot9^r-5\) is a core and \(E_r(b_r)=2\cdot8^r-5\). The negative fixed shift \(b=-5\) lies outside the positive domain and defeats precisely the nonzero-lattice argument.

There is also no fixed binary precision that recovers every depth. For any modulus \(2^s\), sufficiently large consecutive \(r\) give
\(2\cdot8^r-5\equiv2\cdot8^{r+1}-5\pmod {2^s}\), while their decoded depths differ.

## 3. Exact cylinder and receiver-guard transport

The entire native \(12^r\) cylinder has a particularly small parameterization:

\[
n=2\cdot8^r t-5
\quad\xrightarrow{12^r}\quad
m=2\cdot9^r t-5.
\tag{3}
\]

For \(r\ge1\), the positive domain is exactly the integers \(t\ge1\). At \(r=0\), positivity requires \(t\ge3\). The prefix is maximal exactly when \(8\nmid t\), since
\(q(n)=r+\lfloor v_2(t)/3\rfloor\). Equivalently the unmaximized prefix guard is the single cylinder

\[
n\equiv-5\pmod {2\cdot8^r}.
\]

Suppose a receiver requires its input \(m\equiv a\pmod M\). On the **already retained** prefix (3), this is precisely

\[
2\cdot9^r t\equiv a+5\pmod M.
\tag{4}
\]

Let \(g=\gcd(2\cdot9^r,M)\). There is no compatible parameter if \(g\nmid a+5\). Otherwise (4) gives one residue class modulo \(M/g\), obtained by inverting \((2\cdot9^r)/g\). The positive parameter cut and any maximality requirement remain part of the domain. This works for mixed binary/ternary moduli; it does not silently invert a nonunit.

For example, after two blocks the parent is always \(4\pmod9\). A receiver requiring \(0\pmod9\) is impossible, whereas the guard \(4\pmod9\) adds no restriction. Counting the unconditional density of the first guard would miss this exact incompatibility.

With normalized odd Haar measure (or the corresponding natural densities),

\[
\Pr(q\ge r)=8^{-r},\qquad
\Pr(q=r)=\frac7{8^{r+1}}.
\tag{5}
\]

This follows directly from the nested dyadic cylinders, whose tails have density \(8^{-r}\to0\). It is a coding statistic, not a claim that different children are independent or that the cores are grounded. On the 4096 normalized odd residues at precision \(2\cdot8^4\), the exact counts for depths \(0,1,2,3,\ge4\) are \(3584,448,56,7,1\).

## 4. The original source must survive a change of receiver

Let an independently authenticated receiver send \(m\) to \((Pm+B)/Q\), with \(P,Q>0\) integers and integer \(B\). It may be an actual forward word or a certified common-future dependency; that distinction and its guard remain external premises. Composition with (3) gives

\[
z=\frac{P9^r n+5P(9^r-8^r)+B8^r}{Q8^r}.
\tag{6}
\]

Thus strict payment against the original source is **equivalent** to

\[
(Q8^r-P9^r)n>5P(9^r-8^r)+B8^r.
\tag{7}
\]

The implementation transports this full carrier and verifies actual supplied receiver words separately. It does not authenticate a receiver merely from its affine coefficients.

If the receiver is an actual nonempty positive Collatz valuation word, then \(P=3^\ell,Q=2^A,B>0\). In this special case, payment requires \(P9^r<Q8^r\). Consequently a fixed finite bank of such forward receivers cannot pay at arbitrarily large forced-prefix depths: beyond a finite depth, **none of its words can pay any source on which it is actual**. This is only a barrier for that fixed placement and that fixed finite bank. It does not prohibit adaptive words, earlier checkpoints, or guarded common-future receivers with negative intercepts. Formula (7), rather than this special sign conclusion, applies to those receivers.

An exact budget, without logarithmic rounding, is
\[
R_w=\max\{r\ge0:P9^r<Q8^r\},
\]
with value \(-1\) for an empty set. The routine `forward_depth_budget` computes it by successive integer multiplications; termination follows from \(9/8>1\). For a finite bank the maximum of its individual budgets is a necessary source-payment bound. Within this bound the native guard and the strict height inequality (7) still need checking; the budget is not a sufficient payment criterion.

The one-letter receiver \((2)\) illustrates the sharp boundary. At depth two, its transported map is \((243n+319)/256\), so payment is exactly \(13n>319\). At depth three it is \((2187n+3767)/2048\), which never pays a positive source. Concretely,

\[
3067\xrightarrow{12^3}4369\xrightarrow{2}3277.
\]

Although \(3277<4369\), it is greater than the original \(3067\). This is an actual-word hostile, not merely a formal slope example. Returning a child to its larger core alone likewise supplies no numerical descent.

## 5. The requested Mersenne children and strict ROOT receipts

For odd \(K\ge5\), put \(m=2^{K-2}-1\). Then \(v_2(m+5)=2\), so \(m\) is a canonical core. Every legal requested child

\[
h_r=\frac{8^r(2^{K-2}+4)}{9^r}-5
\]

has maximal prefix exactly \(12^r\) and returns to \(m\). Its pair \((K,r)\) is recovered as follows: apply (1), check that \(b+1\) is a power of two of odd exponent at least three, and recover \(K=2+\log_2(b+1)\). Thus the child itself, as an exact positive integer, retains both coordinates. A residue-only observer need not.

The inherited legality bound is

\[
1\le r\le\left\lfloor\frac{1+v_3(K-4)}2\right\rfloor.
\]

It follows from \(m+5=4(2^{K-4}+1)\) and the odd-exponent lifting formula \(v_3(2^e+1)=1+v_3(e)\). The present package checks this formula but does not claim it as a new result.

For any positive core \(b\), let \(\mathcal R(b)\) be the set of supplied, exactly replayed, first-hit ROOT words at \(b\). Prefixing \(12^r\) is a bijection

\[
\mathcal R(b)\longleftrightarrow\mathcal R(E_r(b)),
\]

whose inverse removes the maximal forced prefix. This is valid even if both sets are empty. The first-hit odd rank increases by precisely \(2r\) on expansion, when a receipt exists. No discovery algorithm or convergence premise is hidden in this statement. All legal children of one Mersenne core have exactly the same ROOT obligation, so their number must not be counted as independent terminal evidence.

The APIs `expand_root` and `retract_root` consume an explicitly supplied receipt, validate all actual valuations, and reject any step after ROOT. `attach_receiver` records both the core and original source, and exposes separate `paid` and `reached_root` fields. The parent’s deeper receiver construction can use (4), (6), and (7); its own guard and terminal proof remain necessary. This note does not claim the newly decoded child family is fully grounded.

## 6. Exact controls and reproduction

[Program](../../04-computation/experiments/paper_child_closure_transfers_20261007.py) and [saved output](paper_child_closure_transfers_20261007.out).

```text
python -B -X utf8 04-computation/experiments/paper_child_closure_transfers_20261007.py
python -B -O -X utf8 04-computation/experiments/paper_child_closure_transfers_20261007.py
```

The exact finite universe comprises all 10,000 positive odd integers below20,000 for an independent literal maximal-prefix decoder; 574 affine prefix pairs; 76,320 receiver-congruence comparisons over twelve moduli including nonunits; every legal requested child for odd \(K=5,7,\ldots,605\) (112 pairs); and 252 supplied ROOT transports from canonical cores below512. Only the test fixture computes those finite ROOT words, with an explicit 2000-step cap. Production splicing never searches for them.

Additional controls cover the finite dyadic depth partition, signed affine receiver composition, the original-source hostile, finite-precision ambiguity, illegal ternary inverses, ROOT padding, and boolean/float/container type aliases. All assertions are explicit checks and remain active under `-O`. The script has no source-bit cap in the exact decoder; materializing a receipt with \(r\) blocks naturally uses a word of length \(2r\).

Both modes reproduce **140,447 exact checks** and the saved normalized-LF SHA256 `b44b8d31d153d7fe5d81e7e3c508eefd12f80d3e5159458fd14b45ea52f922d1`.

The constructive survivor is a total lossless change of boundary and a compact exact guard pullback. The remaining obligation is to provide a paid or grounded receiver at the decoded core **with the original source still marked**; paper limit theorems and finite inverse depth do not supply that missing proof.
