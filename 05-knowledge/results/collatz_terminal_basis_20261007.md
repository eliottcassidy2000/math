# A minimum terminal basis for the frozen finite Collatz atlas

**Status:** PROVED finite-graph statements and cylinder lemma; FINITE-EXACT census and certificate discharge. The statement about every positive input remains OPEN. This package certifies exactly the declared finite universe and retains reusable receipts.

**Artifacts:** [script](../../04-computation/experiments/collatz_terminal_basis_20261007.py), [exact output](collatz_terminal_basis_20261007.out).

## 1. Inheritance and the distinction being tested

Use the accelerated odd map

\[
U(n)=\frac{3n+1}{2^{v_2(3n+1)}}
\]

on positive odd integers. ROOT is the first occurrence of 1; its certificate is the empty word. A nonempty word executed after ROOT is rejected.

The source universe is all **16,384 odd inputs from 1 through 32,767**. The frozen routing rules are the eight-macro adaptive selector, the two inverse families, and every paid sibling at its final inspected frontier from [uncovered join routes](collatz_uncovered_join_routes_20261007.md). No learned rule or oracle certificate is added during routing. Each requested source gets one call to the frozen selector, including sources another call might already certify.

The closest proved mechanisms are [retained observation union](adaptive_observation_union_20261004.md) and the [finite seed kernel](finite_seed_kernel_20261005.md): checked actual edges can be shared, and ROOT reachability is a least fixed point. The canonical hostile is an ungrounded dependency cycle: mutually conditional receipts do not supply a ROOT witness. The earlier near miss is counting a smaller dependency as a completed source while its child is still unproved. The useful sidecar here is the exact partial word and its endpoint, including observations attached to an unsuccessful routing attempt.

There are three different objects:

1. A **strict-size dependency DAG** stores each source's advertised smaller children.
2. A **literal observation graph** retains every checked odd edge in all actual source words and advertised child words.
3. A **terminal receipt** supplies an independently verified first-hit word from a previously missing endpoint to ROOT.

The distinction matters: no new oracle is needed to pass from the first object to the second.

## 2. Frozen finite results

| Object | Exact count |
|---|---:|
| B8 ROOT or REDUCED requests | 16,141 |
| B8 PENDING requests | 243 |
| No-route leaves after the two inverse rules and final-frontier inspection | 212 |
| ROOT-only closure using just the strict-smaller dependency DAG | 10,474 |
| Submitted literal edge occurrences retained from the routing records | 60,322 |
| Distinct initially observed edges | 25,663 |
| Initial graph vertices, including ROOT | 25,701 |
| Requested sources already connected to ROOT in this graph | 16,258 |
| Requested sources still unresolved | 126 |
| Distinct missing-edge terminals demanded by those 126 sources | **37** |
| Ungrounded observed cycles | 0 |

Every smaller child is itself in the declared source universe. Conditional on the 212 no-route leaves being ROOT, the strict-size DAG proves all requested sources by ascending induction. In that particular implication system every leaf needs a seed: it has no outgoing deduction. This is a conditional statement, not 212 silently accepted certificates.

Retaining literal observations reduces the actual missing obligations to 37 terminals. For example, the source 4,591 has the observed 26-letter prefix

```
1,1,1,2,1,1,2,1,1,1,1,2,1,1,1,2,1,1,1,2,1,2,1,1,1,1
```

ending at 2,717,873. That terminal is shared by 32 requested inputs. Its supplied 35-letter ROOT word introduces only eight previously unobserved edges; the other edges already occur in retained records.

This comparison is against a graph that collects observations from **every** requested source. It is not a claim of identical selector calls or runtime to the previous ascending closure, which skips already-certified inputs.

## 3. Explicit assumption discharge and its exact minimum

For a finite authenticated functional graph, each unresolved path either reaches a missing outgoing edge or remains in an ungrounded cycle. The latter never proves itself. In this census every unresolved requested path ends at one of 37 missing terminals; no cycle occurs.

Let their set be \(S\). The program first freezes \(S\) from the finite observed graph. Only then does the explicitly labelled experimental routine `discover_selected` explore an orbit, with cap 10,000 odd steps, for each selected terminal. It succeeds in every case. The 37 results are now literal data in `TERMINAL_WORDS` and are repeated in the output's `TERMINAL_RECEIPTS_JSON` record. Every receipt is independently replayed and has no earlier ROOT.

The production function `complete(graph, supplied_terminal_words)` performs no orbit discovery. It authenticates the retained edges, verifies the supplied terminal words, inserts their actual edges with provenance, and compiles every ROOT-connected first-hit word by increasing proved graph rank. Omitting a needed receipt leaves sources unresolved. Supplying an empty nonroot receipt, an incorrect exponent, a forged edge, or ROOT padding is rejected.

**PROVED finite conditional theorem.** If each member of \(S\) has the supplied ROOT receipt, all requested inputs have a compiled ROOT receipt. Here all 37 premises have been discharged by explicit words.

The result is stronger than greedy irredundancy. Independent replay verifies

\[
\operatorname{Orbit}_{\text{first ROOT}}(s)\cap S=\{s\}
\qquad(s\in S).
\]

Thus no actual first-hit word, even one starting outside the original graph, can contain two different members of \(S\): after its first terminal, its deterministic future is already fixed and avoids the others. Each missing terminal's outgoing edge is indispensable to the actual ROOT path of at least one requested source. If completion is restricted to inserting literal first-hit ROOT words followed by rooted edge closure, at least **37 such words** must be supplied. The displayed basis attains the bound.

This is a minimum in that proof representation. It is not a lower bound against symbolic identities, a packet that contains many distinct point receipts, new common-future rules, other proof systems, or computational runtime. The general fixed-component theorem is inherited from the finite seed kernel; the new census also checks the terminal-future condition that keeps the minimum valid when certificate interior edges are imported.

Two independent finite controls agree:

- A demand-greedy order chooses the terminal serving the most remaining requested inputs, breaking ties by integer size. It imports all 37 and produces the same final graph as simultaneous insertion.
- Deleting any one terminal receipt leaves exactly its original demand set unresolved. All 37 deletion tests pass.

The smallest requested representative of each unresolved component is

```
4591, 9663, 10087, 12447, 13551, 15519, 16551, 17439, 17647,
17691, 18943, 20095, 20895, 21231, 21351, 21787, 25243, 25927,
26047, 26623, 26863, 27135, 27327, 27367, 28207, 28671, 28783,
28911, 28923, 30439, 30879, 31167, 31215, 31743, 31903, 31911,
32767.
```

Its ROOT word is obtained by concatenating its retained prefix with the supplied terminal suffix, without another orbit search. These representatives can seed a separate symbolic guard compiler.

## 4. What the completed proof stores

The 37 explicitly selected discoveries total **1,623 odd steps**; their longest word has 89 letters. They add **435 distinct edges** to the existing graph.

The final proof has **26,098 edges and 26,099 vertices**, all ROOT-connected. It returns all 16,384 requested first-hit words; their longest length is 113. Storing these routes separately would repeat **559,560 odd edges**. Shared suffixes are preserved in the finite proof graph, not approximated by their lengths or residues.

An independent replay of every requested certificate reconstructs exactly the final edge dictionary. Hence the final graph is precisely the union of the canonical first-hit trajectories of this universe, with no unrelated edges. In a literal-edge representation, the **435 missing edge observations** are therefore also necessary and sufficient. This edge count differs from the 37 terminal-certificate queries and from the 1,623 steps spent independently discovering the selected tails.

The output supplies the full terminal words, component demand sets, greedy provenance, and SHA256 digests of the initial graph, final graph, and every requested source/word pair. Rebuilding the graph from the frozen routing records and these words preserves every requested route. `complete` returns the complete certificate dictionary for downstream reuse.

## 5. A useful redundancy control

The inherited uncapped reset switch applies when

\[
K=v_2(n+1)\ge2,\quad t=(n+1)/2^K,\quad
r=v_2(3^Kt-1)\ge2.
\]

Its smaller child is \((n-1)/2\). The source word \(1^{K-1},r+1\) and child word \(1^{K-2},2,r-1\) have a common future; at child 1 the latter is replaced by its empty first-hit word. This is a known reset identity, not a new one.

There are 4,096 guarded requested sources, including 22 previously no-route inputs. Authenticating and inserting **all** their source and child words adds **zero distinct edges** to the retained graph. The entire edge dictionary is unchanged. This is a finite stopping result: a locally stronger dispatcher can reorganize already-retained evidence without closing another missing terminal. Its separate all-height usefulness is not denied by this census.

## 6. All-height reuse and the remaining obligation

Let \(w\) be any nonempty actual first-hit ROOT word for \(s>1\), with affine data

\[
F_w(n)=\frac{Pn+B}{Q},\qquad P=3^{|w|},\quad Q=2^{\sum w}.
\]

Since \(Ps+B=Q\), one has \(Q>P\). For every integer \(t\ge0\),

\[
n=s+2Qt,\qquad F_w(n)=1+2Pt<n.
\]

This is an exact native valuation cylinder. At each prefix of cost \(A_i\), the change from its value at \(s\) is \(2^{A-A_i+1}3^it\), an even integer. Every intermediate value stays odd and positive, so all listed valuations remain exact. Proper prefix values exceed 1 because they do at \(t=0\) and increase with \(t\). Thus the same stored word supplies an all-height smaller-endpoint receipt. It need not be the earliest descent.

The program checks this symbolic lemma on four parameters for each of the 37 stored words. For \(t>0\), the endpoint is generally not ROOT: its future remains a separate obligation. Replacing a finite terminal basis by finitely many universally quantified rules does not prove that those rules cover every source or recursively end at the known seeds. [The companion terminal lifts](collatz_terminal_lifts_20261007.md) study stronger guards from the same authenticated data.

The precise global gap is still a coverage or well-founded completion theorem outside the declared universe. The finite basis has no missing premises; extending the universe can introduce new terminal components. Within this universe, adding more observation scheduling alone cannot improve a proof that already contains every requested ROOT path.

## 7. Reproduction and validity controls

From repository root:

```
python -B 04-computation/experiments/collatz_terminal_basis_20261007.py
python -B -O 04-computation/experiments/collatz_terminal_basis_20261007.py
```

Both runs and the saved output agree: **4,329,251 explicit checks**, LF-normalized SHA256 `9b840e84f730a9cb41805edb3e5cb6824d7bcb5a6a17e9994a2eebd136d9e6b6`.

There are no import-time computations, default file writes, or network calls. The only orbit discovery is the explicitly counted 37-terminal reproduction control, bounded at 10,000 steps. Production completion consumes supplied receipts. All other route checks authenticate already supplied words or frozen selector records.

Controls include independent reverse closure versus forward graph walks, every requested word's literal replay, exact final-route-union equality, the 37 deletion tests, greedy/simultaneous agreement, atomic rejection of corrupt words, strict integer and ROOT typing, and an abstract ungrounded-cycle closure hostile. That abstract cycle is not asserted to be a legal Collatz cycle. A discovery cap expiry returns pending, not nonconvergence.

The independent API audit exposed a Python numeric-alias boundary in raw retained graph targets: `True` compared equal to 1, and `5.0` to 5. The replayed result was still mathematically correct, but the stated exact-integer contract was too weak. Retained targets now receive an explicit positive-odd exact-integer check before edge comparison; both hostile examples are permanent regressions.
