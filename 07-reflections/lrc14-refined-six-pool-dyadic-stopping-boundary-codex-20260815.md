# Refined six-pool and dyadic-sample stopping boundary

**Status:** FINITE-EXACT stopping package; no theorem ID, no new LRC(14)
row terminal, and no physical-cover claim.

## Inheritance pass

- **Closest proved mechanism:** THM-3366, `THM-3366-all-sector-complement-clock-completion.md`, compiles an `r`-clock cover of the unsupported open cells into a global cover by at most `7+r` clocks in the positive-aligned sectors.  Its exact current compositions leave `4,056 / 200,069,517,203` k=2 rows/occurrences and `1,897 / 2,548,893,834` k=3 rows/occurrences.
- **Canonical hostile:** a six-clock completion reaches thirteen global clocks, exactly the open LRC(14) boundary; it does not inherit the cited theorem through twelve nonzero speeds.
- **Corrected near miss:** a block on a large modulus need not have size at most `ceil(Q/7)` when it is a lower-order pullback.  At `Q=3822`, an order-three block has size `1274`, not `546`.
- **Least-used sidecar:** the aligned multiplier parity word, together with the quotient-tail numerator residues that would record forced collisions with a completion set.

## Exact refined six-pool composition

The companion reconstructs the current k=2 post-located-`d=6` ledger and k=3 post-one-spike ledger from their exact source programs.  It intersects the actual `(body,D)` keys with the already proved pool-14 completions of size at most five, then examines the surviving rows for six-clock completions.  Raw support-row counts are never subtracted.

The split is unexpectedly disjoint:

| sector | body is the only found six-set | a nonbody six-set only | no pool-14 six-set | body and nonbody |
|---|---:|---:|---:|---:|
| k=2 rows | 2,968 | 382 | 706 | 0 |
| k=2 occurrences | 199,075,988,191 | 133,049,669 | 860,479,343 | 0 |
| k=3 rows | 1,823 | 20 | 54 | 0 |
| k=3 occurrences | 2,547,058,578 | 313,120 | 1,522,136 | 0 |

Representative exact controls are:

- k=2 tautological: `F=(1,2,3,4,8,12)`, `L=D=336`, `|S_D|=160`, completion `F`;
- k=2 nontrivial: `F=(1,2,3,4,7,12)`, `L=1176`, `D=588`, `|S_D|=358`, completion `(1,2,6,7,8,10)`;
- k=2 no-six boundary: `F=(1,2,4,6,9,10)`, `L=2520`, `D=1260`, `|S_D|=646`; its first pool-14 cover has seven clocks `(1,2,3,5,8,9,10)`;
- k=3 tautological: `F=(1,3,4,5,6,10)`, `L=D=840`, `|S_D|=366`, completion `F`;
- k=3 nontrivial: `F=(1,2,5,7,8,10)`, `L=3920`, `D=1960`, `|S_D|=1110`, completion `(1,4,5,7,8,12)`;
- k=3 no-six boundary: `F=(2,3,5,7,10,12)`, `L=5880`, `D=2940`, `|S_D|=1574`; its first pool-14 cover has seven clocks `(1,5,6,7,8,10,12)`.

There is no new compiler terminal.  A hypothetical quotient cover plus six complement clocks gives a global cover by at most thirteen clocks.  The current `(body,D)` key retains neither a forced equality between a tail clock and a completion clock nor the tail numerator residues from which such an equality could be checked.  Calling the 402 nonbody completions terminals would therefore assume LRC(14).

The nonbody rows nevertheless suggest a future operation: regard a nontrivial completion as a body-substitution edge `F -> C`.  Reapplying the critical chart would require the actual quotient tail clocks and the new ruler `14 lcm(C)`.  The present row quotient has discarded exactly that numerator/alignment sidecar, so no descent claim is made here.

## Typed dyadic sample compiler

For a hypothetical row

```text
Y_D(A) subset union_i D_(b_i),
```

the half sample is exact when every aligned multiplier is odd.  Then `1/2` lies in `R_A`, and for every `r in S_D`,

```text
x=(r+1/2)/D=(2r+1)/(2D).
```

Thus the quotient tails cover the supported sheets of the literal half-twist modulus `D`.  An auxiliary half-twist block cover of the unsupported sheets would complete a global literal half cover.

The paired quarter sample is exact when no aligned multiplier is divisible by four.  Both `1/4` and `3/4` lie in `R_A`, and

```text
(r+1/4)/D=(4r+1)/(4D),
(r+3/4)/D=(4r+3)/(4D).
```

These are sheets `2r` and `2r+1` of the literal half-twist modulus `2D`.  The compiler retains both child bits.  Projecting them to one base-support bit would forget the selected dyadic coset, precisely the loss proved by THM-3435.

The proved rank-support gates leave no hit:

- For k=2, only two post-THM3366 rows have `D` outside the THM-3416 rank-six divisor support.  Both have `D=3822`, occurrence weight `641`, and unsupported sheet counts `1530` and `1560`.  Exhaustion of every transverse residue modulo `2D` gives maximum half-block size `1274`, attained by residue `2548` of quotient order three.  Hence one auxiliary block cannot work.  Every `2D` quarter modulus is already inside the rank-six support.
- For k=3, every surviving `D` and every `2D` is already in the relevant THM-3415 rank-five and THM-3416 rank-six divisor supports.  There is also no odd `D` outside THM-3434's rank-seven support.  The support theorems therefore yield no contradiction before an auxiliary search is needed.

These are zero parity-conditional hits as well as zero unconditional hits.  Independently, parity cannot be inferred from the row key: `A=(1,4)` in k=2 and `A=(1,2,4)` in k=3 meet neither sufficient chart.  The first chart fails because not all entries are odd; the second fails because `4` is divisible by four.

## Connection and loss ledger

| field | exact content |
|---|---|
| source | current refined THM-3366 `(body,D)` survivors |
| target | a pool-14 six-completion, or literal half-twist sheets on `D` / `2D` |
| map | strict open-cell cover; evaluation at `u=1/2`; paired evaluation at `u=1/4,3/4` |
| preserved | exact row key, occurrence weight, strict endpoints, every sampled sheet, and both quarter-child orientations |
| destroyed upstream | aligned multiplier parity, tail numerator/unit residue, forced equality with an auxiliary clock, and cross-scale phase compatibility |
| needed sidecar | parity word for `A` plus quotient-tail residues modulo the completion and new-body rulers |
| cheapest decisive hostile | six clocks reach the open thirteen-clock boundary; `D=3822` defeats the only k=2 rank-support opportunity; `(1,4)` defeats chart exhaustiveness |

## Reproduction

```bash
python3 -B 04-computation/lrc14_refined_six_pool_dyadic_stopping_audit_20260815.py
python3 -B -O 04-computation/lrc14_refined_six_pool_dyadic_stopping_audit_20260815.py
```

Normal and optimized transcripts are required to match the stored output byte for byte.  The script pins the exact THM-3366, THM-3415, THM-3416, THM-3434 and refined-composition dependency hashes.  The LF-normalized hashes are

```text
script  e4f57a245ca8c077968c37f6576721162313c9f127b764567b04678f46102d83
output  eddf158869ca046ab1c67ebcf6f47316d3fe02d38aaaf8f8811e91d739e948c8
```
