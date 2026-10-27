# Session74 inheritance: occurrences of11 and usable sidecars

Read-only recovery; frozen source pin `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`. `SOURCES.json` supplies exact repository paths, Git blob IDs, SHA256s, byte sizes and local snapshots for15 load-bearing sources. Connector search lagged at4c3ba61; every load-bearing source was reread at the actual pin. Search was bounded, not an exhaustive census of textual11. No repository edit or push was made. One transient exec transport disconnect resolved on retry.

## 1. The user's exact{2,3,11}: prime doubling escape

**PROVED, inherited.** Put `B_m=((2m−1)^2,(2m+1)^2]`. Integers sharing B_m with a nontrivial integer multiple are exactly`{2,3,4,10,11,12}`; prime restriction gives`{2,3,11}`. Doubling is decisive. For m>=3,
`2((2m−1)^2+1)−(2m+1)^2=4m²−12m+3>0`.

Sources:
- `05-knowledge/results/procgen_numerology_20260923_odd_square_brackets.md`, §2; first commit`4e8f39478a1b5b3cddd3c3374d3512f05324569f`.
- `05-knowledge/results/procgen_brackets_20260924_pairings_transitions.md`, §1; generalized in commit`a754aac33ebcc827b6dc3907491aec5800a210fd`.
- `05-knowledge/results/decoder_prime_square_20260925.md`, §2; precise integer/endpoint convention, commit`5332bfc90def62d20529728b73d8eafa1ce55591`.

**Already proved—do not reclaim as new:** for real multiplier r>1, integer threshold`g(m)=(2m+1)^2/((2m−1)^2+1)` differs from open-real-interval threshold`h(m)=((2m+1)/(2m−1))²`. Define`rho(n)=(2m(n)+1)^2/n`; same-bracket set is`rho>=r`. Prime set is{2,3,11} precisely for`25/13<r<=25/11`. Complete integer-multiplier pairs(k,p) are`(2,2),(3,2),(4,2),(2,3),(3,3),(2,11)`. Strict upper endpoint deletes(3,3), not any prime. The13 skip-free primes and complete internal Collatz-step list are already in §1.

**Connection contract:** n ->(bracket,rho,prime tag) preserves same-bracket multiplication; loses orbit order and carry. Restoration for dynamics requires actual native transitions/ordered carries. Cheapest hostiles:13 enters at25/13;11 exits above25/11;26 in bracket3 fails doubling. Every-fifth prime-gap pulse is local residue density, not a reset: old strong law fails at(841,961], prime853, distance12. Pairing maps with identical brackets/primes include explicit density-zero divergent mutants. Bracket-only facts cannot distinguish Collatz within that class. Arithmetic lane was alerted directly to inherited generalization.

## 2. Working11 sidecar: square-sum attachment port

`05-knowledge/results/decoder_prime_square_20260925.md`, §§1,4–5: **PROVED+AUDITED** Q24 has no Hamiltonian path; every Q25 HP uses`11−25−24`, sums36 and49. Suppress25 to get a HP of`Q24+marked edge{11,24}` using that edge; conversely subdivide. Because11+24=35 is nonsquare, forgetting artificial-edge status leaves the arithmetic category.

**Current correction** in`01-canon/MISTAKES.md`, heading *2026-09-25 Square-sum degree-two forcing forgot the endpoint coordinate*: graph-degree2 does not imply path-degree2. Minimal hostile: triangle. Internal hostile: Q23 endpoint22 omits22−14. Repaired Q24 proof certifies second endpoint in{2,9,11,22}, saturates5, then forces a proper11-vertex cycle. It does not revive the retracted unguarded rule.

**PROVED all-N halving:** with Q_N^(c) edge sums c times squares,`Q_(2N)^(c)[even]/2 ≅ Q_N^(3−c)` for c1,2. Halving toggles square class; quadrupling restores it. Preserves induced-even adjacency; loses odd vertices and attachment/coverage obligations. Hostile1+3=4 but2+6=8 nonsquare; positive4+12=16. Unequal dyadic layers can meet only when their smaller valuation is even; not sufficient (1+2=3).

**Replayed here:** exact pinned original checker normally and with`python3 -O`; outputs byte-identical PASS (`decoder-replay-normal.out`, `decoder-replay-optimized.out`, hashes in SOURCES). Covers both sheets N<=128; independent adjacency through256; doubling through10000 plus proved tail; all23 second-endpoint cases atQ24; independent DFS14..25 counts`0,1,1,1,0,0,0,0,0,3,0,10`. Q25 endpoints are18 paired with`{2,3,4,8,9,10,11,13,22,23}`.

## 3. Later inheritance prevents redoing solved work

`05-knowledge/results/session40_marked_ports_square_sum_20261023/square-sum/README.md`, commit`14df2c7f019d3294692e8e9cd540109710ad6f1f`, accepted by`audit-second/SECOND_AUDIT.md`, already proves all-size edge-marked birth polynomial

`H_(G+v)=sum_a x_va E_a(G)+sum_(a<b) x_va x_vb F_ab(G)`.

E marks an old HP endpoint. F counts spanning two-path forests with a,b in DIFFERENT components; singleton components count once. AtQ25, all oldQ24 E vanish but`F_(11,24)=10`. HP counts/endpoints alone lose all ten births. It also proves endpoint+cotree parity decoding, separately gated by degree/connectivity/orientation. Finite minimum global edge-membership banks at25..32 are`0,1,3,4,1,1,5,6`; not arbitrary-information minima. AtQ25, individually active edges1−24 and4−5 NEVER coexist: an exact hostile to independent marginal gluing.

**Critical supersession:** README §7's coarse single-edge frontier is superseded by`LATE_INHERITANCE.md` and`audit-fourth/FOURTH_AUDIT.md`, integrating commit`d45865c3323af6d1b3ce44525d6d4789b22778d7`. Every edge in a Hamiltonian-cycle graph is includable and omittable in some HP (cycle cutting/chord insertion). With CITED Gerbicz all-N square-sum cycles, every fixed square-sum edgea,b is optional for EVERY`N>=max(32,a,b)`. Full external cycle construction was not re-audited here. Do not reopen eventual single-edge activity.

**Actionable still-open question:** transport specified endpoints AND simultaneous marked edges / realized path-difference hypergraphs / multi-component forest ownership through an explicit interval construction. Cheapest test: incompatibleQ25 pair above; thenQ26 endpoint ambiguities resolved by3−6 membership. Prior protected-edge affine-lift obstruction482−3362 and conditional Hall absorber must be inherited before deriving an induction. Tournament transfer needs actual vertex alignment, detached-cycle exclusion and separate mod8 guards. This is a proposed next research question, not an established transfer.

## 4. Other11s: preserve distinct mechanisms

### Golden11
`01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`: standard Collatz C-period5 no-adjacent-ones words map by`kappa=sum w_j phi^(n−1−j)` to`Z[phi]/(phi^n−1)`. Atn5 this isF11 withphi->4; cyclic shift becomes×4. Size11 isL5, not Terras period11. Nonzero cycles are quadratic-residue cosets; affine group<+1,×4>=<+1,×3> because<3>=<4>. **PROVED:** equivariance,n5 labels, affine group equality. **FINITE-EXACT:** all-odd-n bijection only throughn22 in theorem. Quotient forgets positive integrality and arithmetic fixed-point denominator. Atn10, period5 lifts to`(2kappa5,0)`, reversing residue colour. Repair MISTAKE-584, commit`904f01fcba3`, original`d9f3744328b`. All-period extension requires proof, not extrapolation.

### Terras period11
`01-canon/theorems/THM-4484-free-and-sporadic-cycles.md`: wordw shape(L,k) has fixed point`c_w/(2^L−3^k)`; free iff clock divides shift. Shape(11,7), clock−139, has1 integral necklace among30. `THM-4591-integer-cycles-of-the-terras-map-to-period-301993-are-the-five-pythagorean-cycles.md` gives FINITE-EXACT two-sign periods{1,2,3,11} through301993. Projection to(L,k) loses ordered numeratorc_w; other29 necklaces are cheapest hostile. Not prime bank{2,3,11} or golden period5. MISTAKE-584 repairs negative-factor monotonicity and3x+23 sign. Unmapped Golay/PSL11/etc links remain numerology.

### Finite prime bank{2,3,11}
`01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`, promoted`32a4a306cc6`: **PROVED+AUDITED** no`c log n+f(finite polynomial valuation vector)+bounded g`,c>0, is globally eventually nonincreasing for shortcut/accelerated/max-rise Collatz. {2,3,11}273-step witness is one control; theorem covers every finite prime bank. Expanding rational negative-cycle shadows have positive CRT lifts with identical endpoint features and unbounded height ratio. Missing coordinate: unbounded within-fibre information. Witness changes with requested length; no divergent positive orbit is constructed.

### Mills11
`05-knowledge/results/seam_mills_20260925.md`, §5, commit`8b5088d66e0548820f84d6ec6c14290abd71ee45`: only prime pairq=p³+3 is(2,11), by parity. Distinct from bracket escape. Five cubic prime children11,13,17,19,23 of2 show parent loses branch. Exact gap/fractional carry is needed. Hostilep2,theta9/10: droppingtheta³ reports prime23 while true floor24. Quintic counterpart37=2^5+5. No universal Collatz decoder follows.

### Camion11
`01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md`:341=11*31 fails because31∤D(11);1023=3*341 passes because multiplier3 sends local residues into universal triangle class3. **PROVED** quasi-period/CRT laws require prime powers and valid thresholds. MISTAKE-651, correction`3596c5d0b19`: connected graphs were a falsely tested hostile; modp² law fails only at early indices. Genuine rooted-count hostile changes period. Tower continuation`153d17c527a`: frozen valuation gives passing prime-power initial segment. Local-residue map preserves pass/fail with sidecars, loses original tournament structure. No LRC consequence.

### Eleven spectra: count only
`01-canon/theorems/THM-1455-why-eleven-spectra.md` explicitly scopes11 as count`1+1+1+3+3+1+1`, not prime11. Seven-vertexc3=35+8k uses four-subset parity; full spectrum needsc1. Treat as anti-numerology hostile, not bridge to prime11 or LRC eleven-cores.

## 5. Recommended board and disposition

Board: bracket capacity / marked forest joins / golden cyclic-word quotient / ordered carry divisibility / residue CRT thresholds. Closest proved mechanism: quotient plus missing-coordinate sidecar. Canonical hostile:Q24->Q25 (zero HPs yet ten forest joins). Corrected near miss: endpoint-blind degree2 forcing. Least-used sidecar: actual component ownership with joint incidence. This pass is complete; it asserts no new theorem or novelty.
