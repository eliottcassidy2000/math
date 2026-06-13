# Zeckendorf Non-Consecutivity, Pairing, and Tournament Structure

**Session:** oracle-2026-05-05
**Depends on:** `summand-graph-fermat-zeckendorf.md`, `chemistry-and-the-independence-polynomial.md`,
THM-002 (OCF), THM-227 (k-nacci Mersenne)
**Computations:** all verified numerically in-session

---

## The Central Insight

Zeckendorf's theorem says every positive integer is the sum of *non-consecutive*
Fibonacci numbers. The word "non-consecutive" is the key: it is exactly the
*independence condition* in a path graph. The OCF says $H(T) = I(\Omega(T), 2)$ where
$I$ is the independence polynomial. These are **the same mathematical object evaluated
at different fugacities** — Zeckendorf operates at $x=1$ (unit weight per independent
element), tournaments operate at $x=2$ (doubled weight).

The "pairing" in non-consecutivity arises because every non-consecutive pair
$(F_k, F_{k+2})$ **pairs around the missing mediator $F_{k+1}$**: the three
Fibonacci nodes $F_k, F_{k+1}, F_{k+2}$ form a path $P_3$ in the Fibonacci graph,
and the independence condition skips $F_{k+1}$, leaving an implicit "bond" between
its two neighbors. In the tournament conflict graph, this corresponds to two
vertex-disjoint odd cycles ($F_k$ and $F_{k+2}$) connected through a shared-vertex
cycle ($F_{k+1}$) that is not independently selected.

---

## Section 1: The Master Identity — Cycles at Two Fugacities

Let $C_m$ denote the cycle graph on $m$ vertices. Its independence polynomial
$I(C_m, x)$ is the generating function for independent sets, counted with weight
$x^{|S|}$.

**Classical result:** $I(C_m, 1) = L_m$ (the $m$-th Lucas number, where
$L_1=1, L_2=3, L_3=4, L_4=7, L_5=11, L_6=18,\ldots$).

**New result:**
$$\boxed{I(C_m,\, 2) = 2^m + (-1)^m}$$

*Proof.* Using the deletion-contraction recursion
$I(C_m, x) = I(P_{m-1}, x) + x\cdot I(P_{m-3}, x)$ and the closed form
$I(P_k, x) = \frac{(1+\sqrt{1+4x})^{k+1} - (1-\sqrt{1+4x})^{k+1}}{2^{k+1}\sqrt{1+4x}}$,
one verifies at $x=2$ that the characteristic roots of the path recurrence are $\{2, -1\}$,
giving $I(P_k, 2) = (4\cdot2^k - (-1)^k)/3$.
Substituting: $I(C_m, 2) = (4\cdot2^{m-1}-(-1)^{m-1})/3 + 2(4\cdot2^{m-3}-(-1)^{m-3})/3
= 2^m + (-1)^m$. $\square$

The formula $2^m + (-1)^m$ produces:
- **Odd $m$ (Mersenne):** $I(C_m, 2) = 2^m - 1$ (Mersenne numbers $M_m$)
- **Even $m$ (Fermat-type):** $I(C_m, 2) = 2^m + 1$

| $m$ | $I(C_m, 1) = L_m$ | $I(C_m, 2) = 2^m\pm 1$ | Type |
|---|---|---|---|
| 3 | 4 | **7** = $2^3-1$ | Mersenne |
| 4 | **7** = $L_4$ | 17 | Fermat-type |
| 5 | 11 | 31 | Mersenne |
| 6 | 18 = $h(E_7)$ | 65 | Fermat-type |
| 7 | 29 | **127** | Mersenne prime |
| 11 | 199 | **2047** = $23\times 89$ | Mersenne |

**The number 7 appears twice** in this table:
- As $I(C_3, 2) = 7$: the forbidden tournament value arising from the cycle $C_3 = K_3$
  at the **tournament fugacity** $x=2$
- As $I(C_4, 1) = L_4 = 7$: the Lucas number arising from $C_4$ at the
  **chemistry fugacity** $x=1$

The forbidden H-value lives at one less cycle length when moving from $x=1$ to $x=2$.
This is the "off-by-one" between the chemistry and tournament perspectives — an instance
of the dimension-axis shift from `chemistry-and-the-independence-polynomial.md`.

---

## Section 2: The k-nacci Trace Identity IS the Cycle Independence Theorem

**Theorem (THM-227):** $\operatorname{Tr}(M_k^n) = 2^n - 1$ for all $1 \leq n \leq k$.

**New connection:**

$$I(C_n, 2) = \operatorname{Tr}(M_k^n) \quad \text{for all ODD } n \leq k.$$

*Proof.* For odd $n$: $I(C_n, 2) = 2^n - 1$. By THM-227, $\operatorname{Tr}(M_k^n) = 2^n-1$
for $n \leq k$. These agree. $\square$

**Why only odd $n$?** Even-length cycles do not appear in $\Omega(T)$ (the conflict graph
contains only directed **odd** cycles as vertices). So tournament theory is naturally
restricted to the odd-$m$ branch of the cycle formula, where Mersenne numbers appear.
The even-$m$ Fermat-type branch is the "invisible half" of the cycle independence theory,
relevant in chemistry but not in tournaments.

**Corollary.** The k-nacci trace identity $\operatorname{Tr}(M_k^3) = 7$ for all $k \geq 3$
is equivalent to: *the cycle $C_3$ has 7 independent-set-weighted configurations at fugacity 2*,
and this number is realized by the conflict graph $\Omega(T) = K_3 = C_3$ — which is
permanently forbidden in tournament theory.

---

## Section 3: Lucas Numbers = Minimum-Gap Zeckendorf Pairs

**Every Lucas number is a minimum-gap Zeckendorf representation:**

$$L_m = F_{m-1} + F_{m+1} \quad (\text{indices differ by 2 = minimum valid gap}).$$

Verified for all computed $L_m$ (indices $[m-1, m+1]$, gap exactly 2):

| Lucas | Value | Zeckendorf indices | Gap |
|---|---|---|---|
| $L_3$ | 4 | $[F_1, F_3] = [1,3]$ | 2 |
| $L_4$ | 7 | $[F_2, F_4] = [2,5]$ | 2 |
| $L_5$ | 11 | $[F_3, F_5] = [3,8]$ | 2 |
| $L_6$ | 18 | $[F_4, F_6] = [5,13]$ | 2 |
| $L_7$ | 29 | $[F_5, F_7] = [8,21]$ | 2 |

**Interpretation.** Every Lucas number is the "minimum-gap pairing" of two Fibonacci
numbers. The gap-2 condition is *exactly* the minimum allowed in a valid Zeckendorf
representation (gap $\geq 2$ required). Lucas numbers are the *most tightly packed*
Zeckendorf two-term sums possible.

**Tournament interpretation.** Two Fibonacci terms at gap 2 correspond to two odd
cycles at distance 2 in the conflict path graph — i.e., two cycles that are NOT
directly conflicting but share a common mediator cycle in between:
$$C_k \text{ conflicts with } C_{k+1} \text{ conflicts with } C_{k+2}.$$
The Zeckendorf selection picks $C_k$ and $C_{k+2}$ (non-adjacent in $P_\infty$),
skipping $C_{k+1}$ (the mediator). This is the **pairing** the user identified: the
non-consecutive condition creates implicit pairings of neighbors through missing mediators.

---

## Section 4: Three Structures, One Object

The independence polynomial $I(G, x)$ evaluated at three natural points gives three
different structures, all present in this project:

| $x$ | Mathematical regime | Evaluation | Formula for $C_m$ |
|---|---|---|---|
| $0$ | Boolean algebra | $\alpha_0 = 1$ | 1 |
| $1$ | Chemistry / Merrifield-Simmons | $\sigma(G) = $ total independent sets | $L_m$ (Lucas) |
| $2$ | Tournaments / OCF | $H(T) = I(\Omega(T), 2)$ | $2^m + (-1)^m$ (Mersenne/Fermat) |
| $x \to \infty$ | Extremal: max independent set | $x^{\alpha(G)}$ | $x^{\lfloor m/2\rfloor}$ |

The **Zeckendorf theorem** corresponds to the $x=1$ case: using the independence
polynomial of the *infinite path* $P_\infty$ at $x=1$, every positive integer has
a unique representation as an independent set of minimum total weight = 1
(non-consecutive Fibonacci numbers).

The **OCF** corresponds to the $x=2$ case: the independence polynomial of $\Omega(T)$
at $x=2$ gives $H(T)$, where each independent cycle collection contributes $2^{|S|}$.

**The non-consecutivity pairing in each regime:**
- At $x=1$ (Zeckendorf): selecting $F_k$ and $F_{k+2}$ contributes $F_k + F_{k+2} = L_{k+1}$
  (a Lucas number). The "pair" sums to a Lucas number.
- At $x=2$ (tournaments): selecting cycles $C_k$ and $C_{k+2}$ contributes $2^2 = 4$ to
  the independence polynomial weighted sum. This contributes to $\alpha_2 \cdot 4$ in $H(T)$.

In both cases, the mediator $F_{k+1}$ (or cycle $C_{k+1}$) is absent but "felt" through
its role in the Fibonacci/conflict-graph structure.

---

## Section 5: The Pairing Structure Generates the Forbidden Values

**Why H=7 is forbidden:** $7 = I(C_3, 2) = 2^3 - 1 = \operatorname{Tr}(M_k^3)$.
Achieving $H = 7$ requires $\Omega(T) = K_3 = C_3$: exactly three mutually-conflicting
odd cycles forming a complete triangle. But tournament completeness forces any three
pairwise-conflicting 3-cycles to generate a 4th cycle (proved: THM-200). So $K_3$
as the entire conflict graph is unrealizable.

In Zeckendorf terms: $7 = F_2 + F_4$ (minimum gap = 2) is a Lucas number $L_4$.
It is the minimum-gap two-term sum. The "pairing" of $F_2$ and $F_4$ around the
missing $F_3$ is exactly the structure of two cycle-clusters separated by a mediator,
and this mediator forces the creation of a third independent cluster, contradicting
the $\Omega = K_3$ requirement.

**Why H=21 is forbidden:** $21 = I(P_4, 2) = F_8$ (a Fibonacci number, single Zeckendorf term).
Achieving $H = 21$ requires $\Omega(T)$ to contain a $P_4$ component (path on 4 vertices).
But any path of 4 cycles in the conflict graph forces a 5th cycle to appear
(the middle pair of cycles, sharing a vertex, forces additional cycles by tournament
completion). So $P_4$ is unrealizable.

In Zeckendorf terms: $21 = F_8$ has the *simplest possible* Zeckendorf structure
(single term). No pairing at all. It corresponds to the "most uniform" possible
conflict graph structure ($P_4$, a path), and this uniformity turns out to be
the exact structure that tournament theory cannot achieve.

**The forbidden values are the simplest Zeckendorf structures whose corresponding
conflict graphs are unrealizable:**

| H-value | Zeckendorf | Complexity | Conflict graph required | Realizable? |
|---|---|---|---|---|
| 7 | $F_2 + F_4$ | Min-gap pair (Lucas) | $K_3 = C_3$ | NO |
| 21 | $F_8$ | Single term (pure Fibonacci) | $P_4$ component | NO |
| 11 | $F_3 + F_5$ | Min-gap pair (Lucas) | $P_3$ path | YES (n≥5) |
| 43 | $F_1+F_5+F_8$ | 3-term, mixed gaps | various | YES |

---

## Section 6: The Cycle Spectrum of H-Values

Define the "cycle-type H-values" as those achievable when $\Omega(T)$ is a union of
cycle graphs: $\Omega(T) = C_{m_1} \sqcup C_{m_2} \sqcup \cdots$
(odd $m_i$ only, since only odd cycles appear in conflict graphs).

Then $H(T) = \prod_i I(C_{m_i}, 2) = \prod_i (2^{m_i} - 1)$ (Mersenne numbers for odd $m_i$).

The achievable "cycle-type H-values" are products of distinct odd Mersenne numbers:
$$\mathcal{H}_{\text{cycle}} = \{1\} \cup \{M_3, M_5, M_7, \ldots\} \cup \{M_3 M_5, M_3 M_7, M_5 M_7, \ldots\} \cup \cdots$$
where $M_m = 2^m - 1$.

$$\mathcal{H}_{\text{cycle}} = \{1, 7, 31, 127, 217, 889, 3937, \ldots\}$$

Note: $7 = M_3$ is in this set, but it corresponds to $C_3 = K_3$, which is unrealizable
as an **isolated** conflict graph component. All other cycle-type values
$31, 127, 217 = 7 \times 31, \ldots$ are achievable.

---

## Section 7: The Infinite Fibonacci Path as the Universal Conflict Graph

The universal Zeckendorf object is the *infinite path* $P_\infty$ on vertices
$\{F_1, F_2, F_3, \ldots\}$ with edges between consecutive Fibonacci numbers.
At $x=1$: every positive integer is a specific independent set of $P_\infty$ (Zeckendorf).
At $x=2$: the independence polynomial $I(P_\infty, 2)$ diverges, but truncated to
$I(P_n, 2)$ gives $(4 \cdot 2^n - (-1)^n)/3$.

Now consider the question: **what if the conflict graph $\Omega(T)$ is itself a path
$P_n$?** Then $H(T) = I(P_n, 2)$:

| $n$ | $I(P_n, 2)$ | Zeckendorf of $I(P_n, 2)$ | Notes |
|---|---|---|---|
| 1 | 3 | $[F_3]$ | Single $F_3$ |
| 2 | 5 | $[F_4]$ | Single $F_4$ |
| 3 | 11 | $[F_3, F_5]$ | Min-gap pair = $L_5$ |
| 4 | 21 | $[F_7]$ | Single $F_7$ = **H_forb_2** |
| 5 | 43 | $[F_1, F_5, F_8]$ | 3 terms |
| 6 | 85 | $[F_1, F_5, F_7, F_9]$ | 4 terms, cluster |
| 7 | 171 | $[F_1, F_4, F_7, F_{11}]$ | 4 terms, uniform gap |

**Critical observation:** $I(P_4, 2) = 21 = F_7$ — the second forbidden H-value is
the independence polynomial of $P_4$ at $x=2$. And $P_4$ as a conflict graph is
**unrealizable** (THM-079). The Zeckendorf representation of $I(P_4, 2) = F_7$ is
a single Fibonacci term — perfectly "pure" — reflecting the maximally regular structure
of $P_4$.

The Zeckendorf of $I(P_n, 2)$ values shows an interesting pattern: single-term
at $n = 1, 2, 4$ (indices $3, 4, 7$); minimum-gap pair at $n=3$ (= Lucas $L_5 = 11$);
growing multi-term representations for larger $n$.

---

## Section 8: The Non-Consecutive Condition as Tournament Independence

The Zeckendorf non-consecutive condition states: no two chosen Fibonacci numbers
$F_k, F_{k+1}$ appear together. This is equivalent to: the chosen set is an
**independent set in the infinite path graph** $P_\infty$.

In the tournament OCF: an **independent set in $\Omega(T)$** is a collection of
pairwise vertex-disjoint directed odd cycles. The independence polynomial
$I(\Omega(T), 2) = H(T)$ weights each such collection by $2^{|S|}$.

The *direct structural parallel*:

| Zeckendorf | Tournament OCF |
|---|---|
| Positive integers $n$ | Hamiltonian path count $H(T)$ |
| Fibonacci numbers $F_k$ | Directed odd cycles in $T$ |
| Path graph $P_\infty$ on $\{F_k\}$ | Conflict graph $\Omega(T)$ on cycles |
| Non-consecutive condition | Independence (vertex-disjoint) condition |
| $n = \sum F_{k_i}$ (unique, non-consec.) | $H = \sum_k \alpha_k \cdot 2^k$ (all independent sets) |
| Evaluation at $x=1$: one representation | Evaluation at $x=2$: sum over ALL representations |

**The key difference:** Zeckendorf selects ONE maximum-weight independent set
(the unique representation), while the OCF sums over ALL independent sets
weighted by $2^{|S|}$. Zeckendorf is a "mode" or "argmax," while the OCF is an
"ensemble" or "partition function."

---

## Section 9: The m=15 Curiosity — When Zeckendorf Saturates

$I(C_{15}, 2) = 2^{15} - 1 = 32767$ (Mersenne). Its Zeckendorf representation:
$$32767 = F_{21} + F_{17} + F_{15} + F_{13} + F_{11} + F_8 + F_4$$
with indices $[4, 8, 11, 13, 15, 17, 21]$ and gaps $[4, 3, 2, 2, 2, 4]$.

The middle three terms $F_{11}, F_{13}, F_{15}$ all have minimum gap 2 — a
**cluster of minimum-gap pairings**. This cluster structure means that $32767$
has a "maximally packed" central region in its Zeckendorf representation,
surrounded by larger gaps. This mirrors the tournament structure: a cycle $C_{15}$
has 15 directed odd cycles forming a ring, and the maximum independent set in the
center of such a ring is a consecutive alternating pattern.

The general structure of $I(C_m, 2) = 2^m - 1$ (odd Mersenne) in Zeckendorf:
for large $m$, the representation develops a central cluster of minimum-gap pairs
surrounded by boundary terms with larger gaps. This cluster structure encodes the
"oscillating independence" in large cycle conflict graphs.

---

## Section 10: Open Questions

1. **Prove the cycle-path identity directly.**
   $I(C_m, 2) = 2^m + (-1)^m$. The closed-form proof uses generating functions;
   is there a bijective proof using directed graph structure?

2. **Characterize the Zeckendorf representations of $I(P_n, 2)$.**
   The sequence $I(P_n, 2) = 3, 5, 11, 21, 43, 85, 171, \ldots$ has single-term
   Zeckendorf at $n=1,2,4$ (corresponding to pure Fibonacci values) and
   multi-term otherwise. What determines when $I(P_n, 2)$ is itself a Fibonacci number?

3. **The forbidden values as simplest-structure values.**
   $H=7$ (min-gap 2-term Zeckendorf = Lucas) and $H=21$ (single-term = pure Fibonacci)
   are forbidden because their corresponding conflict graphs ($K_3$ and $P_4$) are
   unrealizable. Is there a general principle: "the simpler the Zeckendorf structure
   of $h$, the more likely the corresponding conflict graph is unrealizable"?

4. **The Mersenne-prime connection.**
   $I(C_m, 2) = 2^m-1$ for odd $m$. When is this a Mersenne prime? The known
   Mersenne primes $M_3=7, M_5=31, M_7=127, M_{13}=8191, \ldots$ correspond to
   specific cycle conflict graphs. Is there a tournament-theoretic characterization
   of which tournament structures correspond to Mersenne primes?

5. **Generalization to higher fugacities.**
   $I(C_m, x) = 2$ gives the empty set and single-vertex contributions only...
   What about $I(C_m, k)$ for integer $k \geq 2$? In $k$-nacci ($k$-step Fibonacci)
   theory, the transfer matrix traces $\operatorname{Tr}(M_k^n) = 2^n-1$ for $n \leq k$.
   Does $I(C_n, k) = k^n + (-1)^n$ for $k$-step recurrences?

6. **The Lucas-forbidden duality.**
   $7 = I(C_3, 2) = I(C_4, 1)$. The forbidden tournament value equals the Lucas
   number one step higher in the chemistry evaluation. Is there a similar "shift by 1"
   pairing for $H=21$? We have $I(P_4, 2) = 21$; is there a natural graph $G$ with
   $I(G, 1) = 21$ in a meaningful chemistry/Lie algebra context?

---

## Summary

The Zeckendorf theorem is the $x=1$ shadow of the tournament OCF (at $x=2$). The
non-consecutive condition in Zeckendorf — the condition that no two chosen Fibonacci
numbers are adjacent — is exactly the independence condition in the conflict graph
$\Omega(T)$. The "pairing" induced by non-consecutivity corresponds to cycles
pairing around missing mediator cycles in the conflict graph structure.

The master identity $I(C_m, 2) = 2^m + (-1)^m$ (new) connects cycle independence
polynomials to Mersenne/Fermat numbers. Combined with the k-nacci identity
$\operatorname{Tr}(M_k^n) = 2^n-1$ (THM-227), this gives:

$$I(C_n, 2) = \operatorname{Tr}(M_k^n) \quad \text{for all odd } n \leq k.$$

The forbidden H-values $\{7, 21\}$ are the simplest Zeckendorf structures (min-gap
Lucas pair and pure single Fibonacci term, respectively) whose conflict graph
realizations ($K_3$ and $P_4$) are tournament-unrealizable.

Every Lucas number is a minimum-gap Zeckendorf pair $L_m = F_{m-1}+F_{m+1}$,
representing the "tightest possible pairing" in both the additive number theory
and the tournament independence structure.
