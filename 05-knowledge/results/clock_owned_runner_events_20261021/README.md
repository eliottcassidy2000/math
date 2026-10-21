# Clock-owned runner events: oriented repair graphs and their exact obstruction

Status: PROVED elementary certificate/obstruction lemmas below; FINITE-EXACT replay controls. No general LRC proof, no LRC(14) closure, no novelty/priority claim.

Session 19, 2026-10-21. Based on repository `de0463bc8ddf`. This packet is staged outside the repository and makes no Git changes.

## Useful result first

At the critical regular-polygon clock `q=2k+1`, unit speeds either have an immediate grid witness or form a complete signed residue system. In the second case there is an intrinsic, physically timed **oriented graph**, not generally a tournament: an edge records another runner entering danger before the clock's unique unsafe owner can leave danger. Reciprocal edges are impossible by an exact modular-inverse inequality. A sink is a sufficient lonely-time certificate.

A cycle has an exact multiplicative invoice. For a prescribed distinct residue cycle, product deficiency is necessary and sufficient for realizing that cycle at arbitrarily large integer heights. The realization can ensure that its blockers are unsafe AT the repair endpoints, with no extra tooth wrap, not merely somewhere earlier along the corridor.

Neither no-2-cycles nor two-sided repair closes LRC. The explicit common-start row `{4,5,9,12,14}` defeats all ten inverse-clock repair endpoints, while `t=1/24` is lonely. More strongly, quadratic-residue rows for every prime `q=3 mod 8`, `q>=11`, defeat all these two-sided endpoints at large height even though their actual loneliness tends to `1/2`.

## Inheritance and concept board

Closest proved mechanism:
- `THM-357-lonely-runner-endpoint-protection-trichotomy.md`: finite open forbidden intervals, protected endpoints, and the strict/boundary/full-cover distinction. Section 1 below is a signed-event implementation of that mechanism, not a new claim that endpoints decide the problem.
- `THM-524-lrc-binding-pair-reduction-and-the-regions-reframe.md`: actual lower-envelope switches and the requirement that every other runner remain clear. Its on-grid and order-tournament losses motivate keeping time and owner here.

Canonical hostile:
- Tight AP `{1,...,k}` has an isolated good point at `1/(k+1)` where an exit and an entry tie. Testing only open cells misses it.
- `THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate.md`: fixed-bank blindness is method failure, not LRC failure. This packet adds a different adaptive-repair blind family, with high loneliness rather than hidden difficulty.

Corrected near miss:
- The repository's `q<=25` claim is false; `26*{1,...,12} union {339}` is its hostile. No bounded clock bank is promoted here.
- A forbidden-set complement need not contain an open interval: boundary-only witnesses must remain.

Least-used relevant sidecar:
- `THM-817-return-satellite-max-speed-cell-decomposition.md`: retain actual cell and left/right endpoint owners instead of reconstructing them from anonymous pair data. We retain the source clock AND target first-entry time on every arrow.

Historical prompt only, not proof dependency:
- `07-reflections/lrc-event-media-tournaments-holes-gates-certificates-s539.md`: unanchored event/gate media mix good and bad states. Its future-priority versus present-occupancy distinction is implemented exactly below.

Live board: tied signed event tape; critical-grid occupancy; timed first-entry graph; multiplicative residue-cycle invoice; endpoint-to-clock arrival debt; high-loneliness hostile family. None may discard another's load-bearing coordinate.

## 1. Exact signed event tape

Let `V={v_1,...,v_k}` be distinct positive integers, `k>=2`, `n=k+1`, and `delta=1/n`. All phases start at zero. This is an integer-speed statement; arbitrary real speeds and arbitrary starting phases are not being identified with it.

Runner `v` is bad in the open teeth

`(m/v-delta/v, m/v+delta/v)`.

Its right endpoint `(nm+1)/(nv)` is an exit event; its left endpoint `(nm-1)/(nv)` is an entry event. Work modulo one, but order the events by the actual time cut at zero. At a tied event time `e`, retain the full exit-owner set `X(e)` and entry-owner set `I(e)`.

If `N^-` and `N^+` count bad runners immediately before and after `e`, then

`N(e) = N^- - |X(e)| = N^+ - |I(e)|`,

`N^+ = N^- - |X(e)| + |I(e)|`.

These identities hold because each forbidden tooth is open; no event owner is bad at its own endpoint. Initially, just after zero, every runner is bad, so `N=k`. Thus a finite signed tape reconstructs the complete occupancy and tests loneliness exactly: a witness exists if and only if some event has `N(e)=0`. A good open component has good boundary events, so no separate midpoint witness search is required for existence. Midpoints remain useful as an independent replay.

For `{1,...,k}`, at `e=1/n`, `X={1}`, `I={k}`, and `N^-=N^+=1`, while `N(e)=0`. This is the canonical reason not to split a tie into an arbitrary ordered pair of events.

Ordering all distinct events creates only a transitive comparison. Discarding event signs/owners in favor of that order is not a new LRC invariant.

## 2. The critical polygon clock

Put `q=2k+1=2n-1` and assume each `v_i` is a unit modulo `q`. At any nonzero grid time `a/q`, a runner is bad exactly when

`a v_i = +1 or -1 (mod q)`.

Indeed `1/q < 1/n < 2/q`, and unit times are not required: multiplication by a unit speed makes residue zero impossible whenever `a!=0 mod q`.

There are `k` antipodal nonzero time classes `{a,-a}`. Each speed excludes exactly one such class, the inverse pair `{v_i^{-1},-v_i^{-1}}`.

**Grid certificate.** If two speeds are equal up to sign modulo `q`, some time class is excluded by no speed, giving a lonely grid point. This is sufficient, not necessary outside the unit hypothesis.

The only remaining unit case is that the speeds form a complete system up to sign. This forces `q` prime: otherwise there are only `phi(q)/2<k` unit sign classes. Composite `q` is therefore already discharged within the all-unit sector, not globally.

For each owner `v_i`, define

`a_i = v_i^{-1} mod q`, with `1<=a_i<=q-1`,

`t_i = a_i/q`,

`m_i = (a_i v_i-1)/q`.

At `t_i` exactly owner `i` is bad, at phase `1/q`. The two endpoints of that owner's current forbidden tooth are

`t_i^+ = t_i + k/(n q v_i) = (n m_i+1)/(n v_i)`,

`t_i^- = t_i - (3k+2)/(n q v_i) = (n m_i-1)/(n v_i)`.

At the antipodal clock `1-t_i`, these endpoints are reflected. Thus testing these `2k` endpoints also covers both endpoints at all `2k` nonzero critical grid clocks. Reflection creates no additional candidates.

## 3. Intrinsic timed forward graph

For `j!=i`, let

`r_ij = a_i v_j mod q`, so `2<=r_ij<=q-2`,

`b_ij = q-r_ij`,

`c_ij = n b_ij-q >=1`.

Starting at `t_i`, runner `j` is initially good. Its first entry into danger occurs at

`E_ij = t_i + c_ij/(n q v_j)`.

Define the forward corridor graph by

`i -> j` iff `E_ij < t_i^+`, equivalently `k v_j > c_ij v_i`.

The vertices are clock-owned repair obligations, not unmarked runners. The observable is actual first-entry order versus actual repair time. Equality is NOT an edge: at equality both relevant runners are at allowed boundary points. Missing edges are honest missing comparisons. A sink has `t_i^+` lonely.

This graph is generally not a tournament. Completing missing edges changes the certificate and has no license here.

### No reciprocal edges, for every k in scope

The reverse values satisfy `b_ij b_ji=1 mod q` and both are at least two. Write `a=b_ij`, `b=b_ji`. Then `ab>=q+1=2n` and `(a-2)(b-2)>=0`, hence `a+b<=ab/2+2`. Therefore

`(na-q)(nb-q)`
`= n^2 ab - n q(a+b)+q^2`
`>= n ab/2-q`
`>= n^2-q = k^2`.

Equality holds exactly for `{a,b}={2,n}`. Reciprocal graph edges would instead multiply to `k^2>c_ij c_ji`, a contradiction.

Consequently this is an oriented graph. If it is acyclic it has a sink and proves loneliness for the given row. No claim that it is always acyclic is made.

### Corridor versus endpoint graph

For a fixed row, entering danger before owner repair does not imply remaining bad at owner repair. The exact condition for runner `j` to be bad in the FIRST following tooth at `t_i^+` is

`c_ij v_i < k v_j < (c_ij+2q)v_i`.

The actual endpoint graph allows later bad teeth as well and requires the exact integer wrap sidecar. It is a subgraph of the corridor graph. Example: at the backward endpoint of owner 11 in `{11,19,50}`, corridor entries occur, but every runner has recovered by the endpoint `1/4`.

## 4. Cycle invoice and an all-height realization theorem

Fix a prime `q=2k+1` and distinct residue labels `h_0,...,h_(ell-1)`, with `2<=ell<=k` and no two equal up to sign. For the cyclic next-edge, set

`r_s = h_(s+1) h_s^{-1} mod q`,

`b_s=q-r_s`, `c_s=n b_s-q`.

Indices are cyclic. Residue consistency gives

`product b_s = (-1)^ell (mod q)`.

**Theorem.** The following are equivalent:

1. There are positive integer speed lifts `v_s=h_s mod q` realizing every cyclic edge as a forward corridor blocker.
2. `product c_s < k^ell`.
3. There are arbitrarily large positive integer speed lifts in those residues realizing every cyclic edge as an ACTUAL forward-endpoint blocker, all in the first following tooth.

Other residue classes, if needed, may be filled arbitrarily to complete the critical signed system; the designated cycle persists. This does not imply that the added vertices also lack sinks.

Proof of `1=>2`: multiply `k v_(s+1)>c_s v_s`; all speed factors cancel.

Proof of `2=>3`: for each edge set `L_s=c_s/k`, `U_s=(c_s+2q)/k`. We have `product L_s<1`, while every `U_s>1`, so `product U_s>1`. Interpolate each ratio linearly from `L_s` to `U_s` with a common parameter. By continuity there is a parameter strictly between zero and one for which the product of ratios equals one. Choose positive weights `w_0=1`, `w_(s+1)=ratio_s*w_s`. They close because their ratio product is one, and they satisfy all strict lower AND upper bounds.

For arbitrarily large `T`, choose positive integer `v_s=h_s mod q` within `q` of `T w_s`. Their ratios converge to the chosen ratios, so all strict inequalities eventually hold. Distinct residue classes ensure distinct speeds. The upper bounds keep each blocker inside its first following tooth, ruling out the wrap ambiguity. Finally, `3=>1` is immediate.

This is a constructive existence proof by open inequalities and residue-lattice approximation, not a finite census and not an LRC counterexample theorem.

### Arrival-clock debt remains

A graph cycle is NOT a chronological cycle of physical events. Edge `i->j` reaches runner `j` at `E_ij`, generally not at its own inverse clock `t_j`. The next edge is defined at `t_j`, so graph edges do not compose as actual time evolution without an arrival/reanchoring theorem.

For `{4,5,9,12,14}`, the edge `4->14` starts at `3/11` and reaches runner 14's first danger entry at `23/84`, not runner 14's own clock `4/11`. Thus the concrete mismatch is already visible in the smallest hostile recorded here. No minimality claim is made.

## 5. Backward repair and an explicit two-sided hostile

At `t_i`, backward travel first brings runner `j` into danger after duration

`(n r_ij-q)/(n q v_j)`.

Thus the backward corridor edge is

`(3k+2)v_j > (n r_ij-q)v_i`.

The corresponding first-tooth endpoint condition is

`(n r_ij-q)v_i < (3k+2)v_j < (n r_ij+q)v_i`.

For `{11,19,50}`, every forward repair endpoint is blocked. Its forward corridor graph is a 3-cycle. Backward repair rescues the row, e.g. owner 50 gives `27/200`, and owner 11 gives `1/4` (the latter after other runners enter and leave danger during backward travel).

But `V={4,5,9,12,14}`, `k=5`, `q=11`, defeats both endpoints for every owner:

| Owner | Clock | Forward endpoint | An actual blocker | Backward endpoint | An actual blocker |
|---:|---:|---:|---:|---:|---:|
| 4 | 3/11 | 7/24 | 14 | 5/24 | 5 |
| 5 | 9/11 | 5/6 | 12 | 23/30 | 4 |
| 9 | 5/11 | 25/54 | 4 | 23/54 | 5 |
| 12 | 1/11 | 7/72 | 9 | 5/72 | 14 |
| 14 | 4/11 | 31/84 | 5 | 29/84 | 9 |

The forward graph is the actual directed 5-cycle `4->14->5->12->9->4`. Nevertheless `t=1/24` is lonely: its phases are `1/6,5/24,3/8,1/2,7/12`. This is failure of the repair strategy only.

## 6. Two-sided repair can miss rows whose loneliness tends to 1/2

**Theorem.** Let `q>=11` be a prime with `q=3 mod 8`, and `k=(q-1)/2`. Let `H` be the nonzero quadratic residues modulo `q`, represented in `{1,...,q-1}`. Put

`V_M={q M+h : h in H}`.

Write `L=min V_M`, `U=max V_M`. If

`U/L < (3k+2)/(2k+3)`,

then every one of the `2k` inverse-clock repair endpoints `t_i^+,t_i^-` is unsafe. This occurs for all sufficiently large `M`. Yet the valid witness `t=1/(2U)` has loneliness `L/(2U)`, which tends to `1/2`.

Proof. Since `q=3 mod 8`, `-1` and `2` are nonsquares, while `-2` is a square. Thus `H` is a complete signed residue system and is closed under multiplication by `-2` and `4`.

For a forward endpoint, choose the unique target whose residue is `-2` times the owner's. Then `r=q-2`, `c=1`, and the exact first-tooth endpoint condition is

`1 < k (v_j/v_i) < 4k+3`.

For a backward endpoint, choose the target with residue four times the owner's. Then `r=4`, and the exact condition is

`2k+3 < (3k+2)(v_j/v_i) < 6k+5`.

All speed ratios lie in the closed interval `[1/R,R]` with `R=U/L`; these non-strict outer bounds suffice because the stated spread bound is strict. For `k>=2`, the bound `R<(3k+2)/(2k+3)` implies both displayed pairs of inequalities: it is smaller than `k`, `(4k+3)/k`, and `(6k+5)/(3k+2)`. The targets differ from the owners because neither `-2` nor `4` is one modulo `q>=11`.

Finally every phase at `t=1/(2U)` lies in `[L/(2U),1/2]`, so its distance is at least `L/(2U)`. The spread bound is below `(k+1)/2` for `k>=2`, making this a lonely witness. Its margin tends to `1/2` as `M` grows. QED.

For `q=11`, `H={1,3,4,5,9}` and every `M>=3` satisfies the displayed spread bound. Thus the fixed five-runner family already gives infinitely many explicit method hostiles. No assertion about infinitude of the prime progression is needed.

Mechanism: the proposed repair samples one tooth per owner at clocks separated from zero, while the near-equal large-speed row has an easy common safe interval at time of order `1/M`. The lost coordinate is the accessible time scale, not hardness or low loneliness.

## 7. Bridge/loss ledger

| Source | Target/map | Preserved predicate | Destroyed information | Required sidecar / decisive test |
|---|---|---|---|---|
| Full bad teeth | Signed tied event tape | Exact lonely-point existence | Nothing if times, signs, owner sets retained | Tight AP simultaneous exit/entry |
| Critical grid | Antipodal excluded-class map | Grid witness existence in unit sector | Off-grid time and velocity magnitude | Explicit `2k` repair endpoint formulas |
| Timed repair obligations | Forward corridor graph | Sink implies lonely endpoint | How long target stays bad; actual arrival if arrows compose | Entry/exit times, tooth index, test `{11,19,50}` |
| Corridor cycles | Residue weights `c_s` | Product deficiency iff cycle realizable at some integer lifts | Specific height ratios and actual row | Strict interval realization, source clocks retained |
| Actual endpoint graph | Unweighted oriented graph | A sink implies endpoint certificate | Margins, multiplicities, clock scale | QR high-loneliness family refutes universal sink |
| Clock-owned cycle | Chronological event path | No automatic transfer | `E_ij` is not `t_j` | `3/11 -> 23/84` versus next clock `4/11` |

## 8. Reproduction and audit scope

Run from any directory:

`python3 /path/to/session19-runner-events/verify.py > /tmp/session19-runner-events.json`

The script uses only Python's standard library and no repository paths. Its output should equal the bundled `certificate.json`. Repeating with `python3 -O` must produce the same bytes; all checks use explicit exceptions and remain active under optimization.

Controls:
- 2,006 admissible modular inverse pairs across every odd `q=5,...,101`, checked exactly (the output records the actual count, authoritative if changed).
- Tight AP rows for `k=2,...,13`; direct Fraction occupancy and integer residue tests independently agree with the tied-event sweep.
- Composite-unit grid positive control `{1,2,4,5}` at `q=9`.
- Three-runner forward-only obstruction, with successful backward branches.
- Five-runner two-sided obstruction, plus independently recovered actual witness `1/24`.
- QR families at `q=11,19,43,59,67,83` and `M=10,100,1000`, every endpoint checked exactly, with explicit `-2`/`4` first-tooth blocker identities and an actual high-margin witness.

These are theorem controls, not an exhaustive LRC search. No minimal-height or smallest-k claim is made for the hostile. `probe.py` is a disposable discovery probe, not a proof dependency; do not publish it as a minimality certificate.

Targeted inbox peek was unavailable because the checkout had no `.machine-id`; none was invented. All worktree files remained untouched.
