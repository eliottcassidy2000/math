# Independent audit: clock-owned runner event certificates

**Verdict: PASS, with the two identified issues corrected in the frozen source.**

The exact signed-event certificate, unit-grid dichotomy, forward no-reciprocal-edge lemma, prescribed-cycle realization equivalence, explicit two-sided hostile, and quadratic-residue hostile family are correct in their stated positive, distinct, common-start integer-speed scope. The packet does **not** prove general LRC, close LRC(14), establish a smallest hostile, or claim novelty.

The audit froze and independently checked these revised producer files:

| File | SHA-256 |
|---|---|
| `source/README.md` | `cf49309a88e3113843a0b7cfc03ea3b15ec8fa789d98c85e314e8df3e9e7fdc3` |
| `source/verify.py` | `94e0c6b1a3bf5721c2cfffb24e74d5bad5895335c5eebdc26368c9337f5c35e0` |
| `source/certificate.json` | `d02d1800e54f2dd2b150ca1ee4d658d1df5fe625d1f89a6b0f7ea4bb3bf520c7` |

The initial source is preserved separately for lineage. It contained a literal strict-extrema wording error in §6 and relied on `assert` for all verification checks. The root corrected the ratio interval to the closed interval `[1/R,R]`, replaced all producer checks with explicit exceptions, and added the optimization replay requirement. The mathematical conclusions and certificate bytes did not change. The audit inspected the exact diff and checks the revised hashes before execution.

## Reproduction

Python 3.9+; standard library only. No repository, network, or external package is required.

```sh
python3 /path/to/session19-runner-events-audit/audit.py > /tmp/audit-normal.json
python3 -O /path/to/session19-runner-events-audit/audit.py > /tmp/audit-optimized.json
cmp /tmp/audit-normal.json /tmp/audit-optimized.json
```

Both runs must report `PASS` and produce identical bytes. The checker itself contains no `assert` statements. It does not import producer computations. It separately runs the producer in both modes and compares its raw stdout with the frozen certificate, then independently reconstructs every certificate field from physical tooth intervals and direct integer modular occupancy. Replaying the producer alone is not treated as independent validation.

`--source-dir PATH` accepts another copy of the frozen three producer files; their hashes must match. `--certificate PATH` audits a substitute JSON certificate after checking the canonical source pins, which permits explicit negative controls without modifying frozen files.

## Independent proof review

### 1. Tied signed events really decide existence

For `n=k+1>=3`, an individual runner cannot be both an entry and an exit owner at the same time. Its forbidden intervals are open. Thus every exit owner was bad immediately before the event but is safe at the event, and every entry owner is safe at the event but bad immediately after it. All nonowners retain their state. This gives, simultaneously rather than in an arbitrary tie order,

`N(e)=N^- - |X(e)|=N^+ - |I(e)|`.

There is no event at zero, and the starting open cell has occupancy `k`. The finite periodic event partition has constant occupancy in each open cell. Any good cell has safe boundary events; isolated good events are retained. Consequently an event with occupancy zero is equivalent to existence of a lonely point. The AP rows expose why entry/exit ties cannot be broken sequentially: both neighboring cells have occupancy one while the tied point has occupancy zero.

The independent checker constructs event times as `(nj+1)/(nv)` and `(nj+n-1)/(nv)`, determines entry/exit sets by exact phases, and computes each neighboring cell's occupancy directly. This is separate from the producer's running-debt sweep.

### 2. The grid dichotomy includes composite clocks correctly

With `q=2k+1=2n-1`, `1/q<1/n<2/q`. A unit speed cannot have residue zero at a nonzero clock numerator, even when the numerator is not a unit. The only bad residues are therefore `1` and `q-1`.

Every unit speed excludes exactly one inverse antipodal time class. There are `k` such nonzero classes. Repeated speed sign classes imply fewer than `k` excluded classes, so a grid witness exists. A complete set of `k` distinct unit sign classes forces `phi(q)=q-1`, hence primality. For odd composite `q`, at least one nonzero residue is a nonunit, so `phi(q)/2<k`; the all-unit sector is discharged already. Nothing here discharges rows containing nonunit speeds.

For the complete signed system, at `t_i=a_i/q` only owner `i` is unsafe. Its phase is `1/q`, its tooth center is `m_i=(a_i v_i-1)/q`, and its exact repair endpoints are `(m_i±1/n)/v_i`. This verifies both stated duration formulas. Reflection sends the two endpoints at `t_i` to the two endpoints at `1-t_i` with their signs interchanged; it adds no new tests. A negative backward endpoint, possible for speed one, is interpreted on the periodic real time line and is checked without silently changing the signed travel direction.

### 3. Three graphs must remain distinct

At an owner clock, write the target phase as `J+r/q`, with `2<=r<=q-2`. Its first forward tooth is centered at `J+1`; its first backward tooth is centered at `J`.

The physical first forward entry is `(J+1-1/n)/v_j`; subtracting the owner clock yields `(n(q-r)-q)/(nqv_j)`. For backward travel the physical entry is `(J+1/n)/v_j`, yielding backward duration `(nr-q)/(nqv_j)`. Comparing signed entry times with the actual owner endpoint produces the stated corridor inequalities.

The target lies inside that first tooth precisely when **both** strict bounds hold. The lower equality is entry and the upper equality is exit; either equality is safe. Actual endpoint occupancy can instead occur in a later tooth, so first-tooth blockers are a subset of actual endpoint blockers, which are a subset of corridor blockers. Both containments can be proper:

- In `{11,19,50}`, owner 11's backward endpoint is `1/4`. Targets enter danger in the corridor but none is bad at the endpoint.
- In `{1,18}`, owner 1's forward clock is `1/5`, first target entry is `11/54`, and repair endpoint is `1/3`. Target 18 is bad there in a later tooth, not the first following tooth.

The explicit strict-boundary fixtures use `q=5`, owner 2, and endpoint `2/3`: target 1 sits exactly at first entry, so there is no corridor edge; target 11 sits exactly at first exit, so the corridor edge exists but endpoint occupancy does not.

### 4. No reciprocal forward edges

Let `a=b_ij`, `b=b_ji`. Then `a,b>=2`, `ab=1 mod q`, and hence `ab>=q+1=2n`. Also `(a-2)(b-2)>=0`, so `a+b<=ab/2+2`. Therefore

`(na-q)(nb-q) >= nab/2-q >= n^2-q = k^2`.

Equality requires both `ab=2n` and one of `a,b` to be 2, which is exactly `{a,b}={2,n}`. The proof uses only odd `q`, unit inverse pairs and the specified ranges, not primality. Two strict forward corridor edges would multiply to the contradictory inequality `k^2>c_ij c_ji`. In particular the actual forward endpoint graph, being a subgraph, also has no reciprocal edges. No corresponding claim is made about the backward graph or a union of the two graphs.

A finite acyclic forward graph has a sink and hence a lonely repair endpoint. It need not be acyclic, and no graph completion to a tournament is permitted.

### 5. The cycle equivalence has both necessity and sufficiency

The residue labels are nonzero and distinct up to sign, as required for their inverses and positive coefficients. They are fixed labels of a genuine cycle, not arbitrary edge weights satisfying only a modular product relation. Writing `r_s=(h_(s+1)h_s^-1 mod q)` and `c_s=n(q-r_s)-q`, every `c_s>=1`.

Necessity follows by multiplying the strict lower inequalities: all speeds cancel and give `product c_s<k^ell`. Product equality is excluded, not a limiting realization.

For sufficiency, set `L_s=c_s/k`, `U_s=(c_s+2q)/k`. The hypotheses give `product L_s<1<product U_s`, with every `L_s>0`. The producer's continuous interpolation argument correctly gives a closed positive weight cycle with all strict bounds. Openness then permits residue-constrained integer approximation at arbitrarily large scales. Every weight is positive, so every designated speed becomes arbitrarily large. Distinct residue labels guarantee distinct integer speeds. Filling unused signed classes does not change the designated clock phases or thresholds.

The audit supplies a separate **exact rational construction**, avoiding numerical roots or an unquantified rounding step. While choosing ratios successively, maintain remaining product `P`. For the next ratio choose a rational point in

`(max(L_s,P/product remaining U), min(U_s,P/product remaining L))`.

This interval is nonempty by the current product feasibility inequalities. Choose its midpoint, divide `P` by it, and at the last edge use `P` itself. All ratios are strictly interior and multiply exactly to one. Clearing denominators gives positive integers `W_s` with

`c_s W_s < k W_(s+1) < (c_s+2q) W_s`.

Now set `v_s(T)=q T W_s+h_s`, with canonical residues `1<=h_s<q`. Each lower and upper slack is an affine function `A T+B` with positive integer slope. The explicit threshold `T0=max(1, floor(-B/A)+1)` over all slacks guarantees every inequality for **all** integers `T>=T0`. This produces positive, distinct, residue-correct lifts and actual first-following-tooth endpoint blockers, not merely corridor entries. The checker verifies these affine slopes and representative values `T0,T0+1,T0+37`; the affine proof, not the finite samples, proves the all-height conclusion.

Example: `q=7`, residue cycle `(1,2,3)`, coefficients `(13,1,1)` has product `13<27`. The audit obtains weights `(18,120,47)`, and every `T>=1` gives the first-tooth cycle via speeds `(126T+1,840T+2,329T+3)`.

### 6. Graph cycles are not chronological paths

An edge reaches the target's first entry time `E_ij`, not the target's own inverse clock `t_j`. The next graph arrow is anchored at that other clock. For the hostile row, `4->14` travels from `3/11` to `23/84`, whereas the next owner clock is `4/11`. A graph cycle therefore does not produce a cyclic chronological event chain. The source explicitly preserves this missing reanchoring theorem rather than assuming it.

### 7. The QR family defeats the method while actual maximum loneliness tends to 1/2

For prime `q=3 mod 8`, the supplementary quadratic-residue laws give `(-1/q)=(2/q)=-1`, hence `(-2/q)=1`. Thus the nonzero squares `H` form a complete signed system and are closed under `-2` and `4`. For `q>=11` both multipliers give targets different from the owner.

The forward target `-2h` has `r=q-2`, coefficient 1, and first-tooth bounds `1<k(v_j/v_i)<4k+3`. The backward target `4h` has coefficient `2k+3`, yielding `2k+3<(3k+2)(v_j/v_i)<6k+5`.

With `R=U/L`, speed ratios lie in `[1/R,R]`. Let `A=(3k+2)/(2k+3)`. The exact comparisons needed are `A<k`, `A<(4k+3)/k`, and `A<(6k+5)/(3k+2)`, all true for `k>=2`. For example, the last cross-product difference is `3k^2+16k+11>0`. The **strict** spread bound `R<A` supplies strict tooth bounds even for the extremal ratio pairs. The corrected closed-interval wording is therefore the right argument.

For integer `M` sufficiently large, all `qM+h` are positive and distinct and their ratio spread tends to one, so the hypothesis holds. At the actual time `t=1/(2U)`, all phases are in `[L/(2U),1/2]`; the attained minimum distance is exactly `L/(2U)`, because the minimum-speed runner attains it. Since `A<(k+1)/2`, this is a lonely witness. If `Lambda(V)=max_t min_v ||vt||` denotes actual maximum loneliness, then

`L/(2U) <= Lambda(V_M) <= 1/2`.

Hence `Lambda(V_M)->1/2` by squeezing. This is not merely a claim about the repair endpoints or a finite sample. For `q=11`, the spread condition is equivalent to `M>25/11`, so every integer `M>=3` works. A fixed prime already supplies infinitely many rows; no theorem on infinitude of the prime progression is required.

## Computational controls and limits

The exact output records the authoritative counts. Controls include:

- All 2,006 producer inverse-pair cases, extended independently to 50,532 pairs over all odd `q=5,...,501`, including composite moduli and equality classification.
- All AP controls `k=2,...,13`, plus the hostile's full event tape: 934 direct event checks.
- All signed systems and height choices in `{0,1,3}` at `q=5,7`: 252 rows.
- Eleven composite all-unit rows for composite odd `q<=51`.
- Every directed cycle modulo rotation with no antipodal repetition at `q=5,7,11`: 1,480 cycles. Exactly 88 have strict product deficiency and receive rational first-tooth lift constructions; 20 have equality and 1,372 surplus.
- The complete bundled three-runner and five-runner certificates, 18 bundled QR rows, and 28 additional QR rows at primes through 251.
- Explicit entry-equality, exit-equality, corridor-only, later-tooth, and noncomposing-clock controls.
- Four syntactically valid corrupted certificates are rejected: tied occupancy, endpoint blocker, substituted arrival clock, and a false margin of `1/2`. The revised producer also rejects duplicate speeds in both normal and optimized modes.

Finite checks are controls for the written proofs, not universal LRC evidence. The audit does not certify priority, literature completeness, unrelated inherited repository claims, or generalized arbitrary-phase/arbitrary-real-speed statements.

## Artifacts

- `audit.py`: independent explicit-failure checker.
- `audit-normal.json`, `audit-optimized.json`: identical exact audit outputs.
- `source/`: final producer files, pinned above.
- `initial-source/`: pre-correction lineage only, not the accepted revision.
- `producer-normal.json`, `producer-optimized.json`: revised producer byte replays.
- `external-negative-controls.json`: subprocess exit results for corrupted-certificate checks.
- `manifest.json`: portable, relative-path SHA-256 manifest and verdict.

All audit work was staged outside the repository. The auditor made no producer-file or Git writes.
