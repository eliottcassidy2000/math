# Order fourteen: complete two-shore exclusion

October 21, 2026. **FINITE-EXACT:** no strong self-converse tournament of
order fourteen has black Hamilton-path-orbit deficit six. All 2,048 necessary
normalized sources have eight supplied Hamilton words with distinct full-score
strings, each failing the blue reverse-complement condition. The certificate
consumer performs no path search or automorphism enumeration.

Together with the [order-twelve exclusion](../deficit_twelve_20261021c/README.md)
and the [combined all-order cutoff](../deficit_six_centers_20261021b/COMBINED_CUTOFF.md),
these two new packets alone leave

\[
\boxed{15,18,20,21,24,25,28,35}.
\]

Those coefficients remain **OPEN in this packet**; a separate simultaneous
result may remove twenty-five. Here
$R_6=z^5+z^6+z^7+O(z^{15})$, while its degree bound35 and the order18
recurrence's guaranteed onset48 are unchanged.

## 1. Exact normalized universe

Use the same literal source definitions as the order-twelve packet:
$\Delta=h-b$, where $h$ counts directed Hamilton-word orbits under the full
source automorphism group and $b$ counts the blue orbits whose positional
reversal is an anti-automorphism.

For a marked Hamilton cycle with successor $R$, write
$m=|\operatorname{Aut}(T)\cap\langle R\rangle|$ and $q=14/m$.
The tournament automorphism group has odd order. The marked-cycle cuts occupy
$q$ path orbits, at most $\gcd(2,q)\le2$ of them blue. Under $\Delta=6$,
therefore $q\le8$. The odd divisors of fourteen are1 and7, so necessarily

\[
\boxed{q=2,\qquad m=7,\qquad R^2\in\operatorname{Aut}(T)}. \tag{1}
\]

No vertex-transitivity assumption is used. This short reduction uses the
inherited [exact cycle-cut source](../paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md),
not an unverified higher-order census.

Label the cycle $0\to1\to\cdots\to13\to0$. The shift by two has thirteen
orbits on the91 unordered pairs, each of size seven. The prescribed cycle
fixes two of these orbits, leaving eleven independent orientation bits:

\[
2^{13-2}=2^{11}=2048.                                    \tag{2}
\]

Every normalized frame contains the actual directed cycle and is strong.
Every candidate in (1) occurs after relabelling its cycle. Different masks
give distinct labelled sources, although isomorphic sources can occur more
than once. This is a complete necessary-frame census, not an ordinary
isomorphism-class census or an exact deficit histogram.

For an independent source reader write $e_i=2i,o_i=2i+1$ with indices in
$\mathbb Z/7$. Each shore is a cyclic tournament; choosing one sign from
each pair $\{j,-j\}$ for $j=1,2,3$ gives three bits on each shore. Cross
arcs satisfy

\[
e_i\to o_j\quad\Longleftrightarrow\quad f(j-i)=1.
\]

The cycle requires $f(0)=1,f(6)=0$, leaving five cross bits. These are the
same eleven bits as (2). The exact edge-orbit mask order has internal even
bits0,2,4; internal odd bits8,9,10; and the free cross values $f(1),\ldots,f(5)$
at bits1,3,5,6,7. The producer sweeps edge orbits; the consumer builds these
two shore connection sets and the cross difference function directly.

If $k=\sum f$, the two full scores are $3+k$ and $10-k$, each repeated seven
times. Their sum is13 and they are distinct because $2k=7$ is impossible.
Thus all frames pass the necessary complementary-score histogram. Passing
does not certify self-converseness; all2,048 frames are deliberately included.

## 2. Eight actual black orbits, with source and target retained

The [score-port lemma](../deficit_twelve_20261021c/README.md) applies unchanged:
full-score strings are invariant under the full automorphism group, and a blue
word must satisfy

\[
d^+(v_i)+d^+(v_{13-i})=13\quad\text{for every }i.
\]

Eight native Hamilton words with distinct score strings, each violating this
condition, give eight distinct black path orbits. The complete bank supplies
eight for **every one** of the2,048 frames. Hence $\Delta\ge8$ whenever the
source is self-converse, and no order-fourteen deficit-six center remains.

The producer visits147,456 score-prefix states and finishes in about0.26s
in the recorded local run. This timing is not a proof premise. It records
16,384 literal words and their full targets: retain the thirteen path arcs
and reverse all other source arcs. The compact saved `inputs.json` contains
the complete311,304-byte producer stream under lossless zlib/base85 encoding.

The consumer reconstructs each source independently, verifies the marked
cycle and shift-two automorphism, validates every word and target, and checks
the sufficient score predicates. It additionally tests the positional
reversal directly on every saved word, independently of the score-based
blackness test. Source masks, vertex labels and targets remain literal; no
isomorphism quotient replaces them. The known six-vertex deficit-six center
still has exactly six distinct detectable black score strings and rejects
the eight-port request. Wrong source, wrong target, malformed path, duplicate
word, negative index and type-alias controls reject.

The generic `score_ports` helper is imported from the frozen order-twelve
checker, whose bytes and receipt are pinned. Only its source/word/target
validation APIs are used; the order-twelve census is not rerun by this checker.
The new `audit_record` binds each received port packet to the actual
fourteen-vertex source and eleven-bit mask.

## 3. Reproduction and scope

```text
python -B 05-knowledge/results/deficit_fourteen_20261021c/check.py --verify
python -B -O 05-knowledge/results/deficit_fourteen_20261021c/check.py --verify
```

Both runs pass24,601 always-active gates and match `receipt.json`, LF SHA-256

```text
2b150fd9d93ea708cd201aef0eaa41fa5d5205daed5b1dba41757cf0da3c329d
```

Independent root proof/API review and optimized replay both passed, with all
24,601 gates and the exact saved receipt hash above.

The independent producer can be rebuilt with a C++17 compiler. Place its
executable and temporary raw stream outside the repository:

```text
g++ -std=c++17 -O3 05-knowledge/results/deficit_fourteen_20261021c/census.cpp -o <scratch>/census
<scratch>/census 2048 <scratch>/bank.bin
python -B 05-knowledge/results/deficit_fourteen_20261021c/pack_bank.py <scratch>/bank.bin
python -B 05-knowledge/results/deficit_fourteen_20261021c/check.py --verify
```

This is an additional completed finite base using the same invariant, not an
all-order theorem that score strings always distinguish enough paths. The
regular Paley7 hostile from the inherited packet prevents that promotion.
No claim about any remaining center coefficient or universal metagraph
connectivity follows.
