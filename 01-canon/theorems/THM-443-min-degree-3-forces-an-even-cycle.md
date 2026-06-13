# THM-443 — (EVEN-cycle lemma) Every finite graph with minimum degree ≥3 contains an even cycle (length ≥4); even-cycle-free ⟺ cactus of odd cycles ⟺ min degree ≤2

> **⚠ CORRECTION (MISTAKE-064, S709).** This theorem is the **EVEN-cycle** statement (length `2k`
> = an even number `≥4`), which is **classical and settled-TRUE**. It is **NOT Erdős Problem 64**.
> Erdős Problem 64 is the **Erdős–Gyárfás conjecture**: min degree ≥3 ⟹ a cycle of length a **POWER
> OF TWO** (`2^k`: 4,8,16,…), which is **OPEN** — treated honestly in **THM-444** (S709). S708
> wrongly conflated the two; this file is rescoped to the (correct) even-cycle lemma only.

**Status:** PROVED (classical even-cycle result; clean longest-path proof recorded) + VERIFIED
(exhaustive over all graphs on ≤7 vertices). The theorem is **classical/folklore**; the new content is
the proof record, the characterization, and the directed-analog (Pfaffian) note. **Does NOT address
Erdős 64** (see THM-444).
**Source:** opus-2026-06-08-S708 (rescoped S709). Connects to CLAUDE.md cut⊕cycle / even graphs `E_n`,
THM-174/S707 (Pfaffian — applies to EVEN dicycles, not power-of-2).

## The theorem

> **THM-443.** Let `G` be a finite graph with `δ(G) ≥ 3`. Then `G` contains a cycle of length `2k`
> for some `k ≥ 2` (an even cycle of length `≥ 4`).

**Proof (longest path + pigeonhole).** `δ≥3` forces `|V|≥4`, so `G` has an edge; take a **longest
path** `P = v_0 v_1 ⋯ v_l`. By maximality every neighbour of the endpoint `v_0` lies on `P` (else
extend), so `N(v_0) ⊆ {v_1,…,v_l}`, and `|N(v_0)| = deg(v_0) ≥ 3`. Let `v_0` have neighbours at
positions `p_1 < p_2 < p_3` (all `≥1`). Three positions, two parities — by pigeonhole some `p_i < p_j`
have **the same parity**. Then
```
   v_0 , v_{p_i} , v_{p_i+1} , … , v_{p_j} , v_0
```
is a cycle of length `(p_j − p_i) + 2`, which is **even** (`p_j−p_i` even) and `≥ 4`
(`p_j−p_i ≥ 2`). ∎

## The characterization (contrapositive)

> **A finite graph has no even cycle of length `≥4` ⟺ every block is an edge or an odd cycle (a
> "cactus of odd cycles") ⟺ it has a vertex of degree `≤2`** (in fact `δ ≤ 2`).

**Proof.** A 2-connected graph that is **not** a single cycle contains a **theta** subgraph (two
vertices joined by 3 internally-disjoint paths of lengths `a,b,c`), whose three cycles have lengths
`a+b, b+c, a+c` summing to the even `2(a+b+c)` — so an even number of them are odd, hence at least one
is **even** (verified all `a,b,c≤5`). Thus an even-cycle-free 2-connected block is a single cycle, and
(being even-cycle-free) an **odd** cycle. So every block is an edge or an odd cycle; the graph is a
cactus of odd cycles, whose leaf-blocks contribute a vertex of degree `≤2`. ∎

## Verification

Exhaustive over the **graph atlas** (all 1252 graphs on `≤7` vertices, `…s708.py`): of the **173**
graphs with `δ≥3`, **every one** has an even cycle (`0` exceptions); the even-cycle-free graphs that
contain a cycle have maximum min-degree `2`. Witnesses: `W₅` (even 4-cycle via same-parity endpoint
neighbours), Petersen (3-regular girth 5 — even cycle present via the theta/`E>V` block).

## Repo bridges

- **(cut ⊕ cycle / even graphs `E_n`).** An even cycle is a specific element of the GF(2) cycle space;
  the theorem is a *realizability* statement — `δ≥3` (i.e. `|E| ≥ 3|V|/2`, large cycle space) forces an
  even cycle to be realised, not merely an even-degree subgraph (the `E_n` objects).
- **(directed analog = Pfaffian, S707/THM-174).** The **directed** even-cycle problem — does a digraph
  have a directed cycle of even length — is the heart of **Pólya's permanent problem**: a digraph has
  **no even directed cycle ⟺ it admits a Pfaffian orientation ⟺ `per(A) = ±det(A)`** (poly-time,
  Robertson–Seymour–Thomas / McCuaig). In a **tournament**, a strong component of size `≥4` has an
  even directed cycle (Moon: strong tournaments are vertex-pancyclic, cycle lengths `3..n`), while a
  transitive tournament has none. (Verified n=3..6: strong → even dicycle at `n≥4`.) So the undirected
  `δ≥3` theorem and last session's Pfaffian thread are two faces of the even-cycle / parity question.
- **(covering / parity theme).** Both this and the covering-systems exploration are "parity-covering"
  problems; see the reflection for the LRC-as-danger-covering lens.

## Scope / honesty

- The theorem is **classical** (longest-path argument is standard); this entry records a clean proof,
  the block characterization, an exhaustive `n≤7` check, and — the genuinely repo-relevant part — the
  **directed-even-cycle ↔ Pfaffian-orientation** bridge to S707 and the cut⊕cycle framing. No new
  extremal bound is claimed.

**Artifacts:** `04-computation/even_cycle_mindeg3_and_covering_s708.py` (+`.out`). Reflection
`07-reflections/even-cycles-covering-systems-and-the-parity-covering-lens-s708.md`. Relates: **HYP-2313**.
Builds on CLAUDE.md (cut⊕cycle / `E_n`), THM-174/S707 (Pfaffian), Moon (pancyclic), Pólya/RST/McCuaig.
