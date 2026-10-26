# Goldbach: independent scoped audit

**PASS: FINITE-EXACT table and native controls; CITED/source-inspected formal reduction.**
No claim of Goldbach, infinite nonhalting, minimum state count, or independent Lean kernel verification follows.

## 1. The objects and directions are correct

The literal table in the saved pinned `GoldbachTm/Tm25/TuringMachine25.lean` has 25 nonhalting state labels, binary alphabet, 50 state/read slots, 49 ordinary instructions, and the unique undefined branch `(24,0)`. The specified semantics write first, then move; blank input is state0 with all tape cells zero. “Halts” means encountering an undefined transition, not visiting a control state that also has a continuing branch.

Independently parsing the table line-by-line and comparing every slot to the wiki encoding verifies the bijection

```
Lean 0..24 → A D E G H I V J W K L M N U O P Q R X Y F S B C T
```

It preserves the read label, write symbol, movement, successor state, initial state and halt branch. Extending this map by identity on tape and head gives a genuine full-configuration transition-system isomorphism. Four deliberate changes to the initial instruction—successor, write, direction and halt—are all rejected. This is the precise isomorphism available here; the separate statement about strong Goldbach is a global halting equivalence.

The saved source's `halt_lemma'` and `halt_lemma_rev` give opposite implications between an even counterexample and blank-tape halting. The definition of primality is ordinary `Nat.Prime`; even targets are `2*n+4`. `never_halt` explicitly assumes nonhalting. Root imports reach `Tm25.Content`, and the saved local Tm25 dependency cone is complete. A no-placeholder text scan is only a static check, not an axiom audit. [Pinned Leng source](https://github.com/lengyijun/goldbach_tm/tree/6cd292062516a0a14ea1b34f2ab75154cae7ab1e)

All 19 saved Leng files match a freshly fetched pinned Git tree by Git blob ID, and the transition file also matches a fresh file fetch. Upstream CI run26210149415 reports the exact SHA `6cd292062516a0a14ea1b34f2ab75154cae7ab1e`, completed/success on May21,2026; job77118861892 and its `leanprover/lean-action@v1` step also report success. This is independently re-observed upstream metadata, not a local build or inspection of kernel output. [Exact CI run](https://github.com/lengyijun/goldbach_tm/actions/runs/26210149415)

## 2. Lineage and the one correction

The original gist was created May8,2016 and explicitly lists states000..026. Rows022 and024 are identical. Leng's dated Oct17,2024 comment describes the 26-state simplification and 25-state compression. Thus the wiki's 2016/25-state row combines historical stages; it does not date the 25-state formalization to2016. [Original gist and dated comment](https://gist.github.com/anonymous/a64213f391339236c2fe31f8749a0df6)

The current wiki also lists a 21-state,41-instruction candidate with an independent-verification caveat. This defeats any unqualified statement that25 is the current smallest claimed size; neither entry supplies a minimum-state theorem. [Cryptids index, observed October26,2026](https://wiki.bbchallenge.org/wiki/Cryptids)

At DeepMind pin `6c28d8524fa008c6fb5aadb759df74a55ce96423`, the primary source has `candidate_fewer_states`, `candidate_equivalent`, and the arithmetic reduction endpoint. Its semantic wrapper proves equivalence to the project-defined `ConjectureMachines.GoldbachConjecture`; this is **not a Mathlib-owned declaration**. `Cryptids/Propositions.lean` defines that proposition using Mathlib `Nat.Prime`. The producer's attribution was corrected, with its seals refreshed. The mathematical statement was already correct. [Pinned defining source](https://github.com/google-deepmind/cryptids/blob/6c28d8524fa008c6fb5aadb759df74a55ce96423/Cryptids/Propositions.lean)

The DeepMind README's toolchain and standard-axiom/comparator statements remain author-reported claims inspected in source. This audit does not build its dependencies, run the comparator, or independently verify the entire21-state development. [Pinned README](https://github.com/google-deepmind/cryptids/blob/6c28d8524fa008c6fb5aadb759df74a55ce96423/README.md)

## 3. Independent native controls

The new checker imports only the Python standard library and reads producer files as data. It verifies the source inventory and seals before the following controls:

* All50 table slots and the exact state renaming, including the sole halt branch.
* 30,000 native microsteps of independently represented full configurations in both labellings. Their configurations agree at every step under the renaming. The unary macro-prefix agrees with the producer's recorded prefix; first four stage clocks are4→45,6→132,8→301,10→596.
* A fresh Eratosthenes sieve checks249 even targets4..500, comparing prime-indicator coefficient counts with direct split counts. This is a finite dictionary check, not a proof of the conjecture.
* All2,025 ordered odd-prime pairs through200 obey the inverse `(p,q)=(2E-O,O-E)`, swap reflection `O↦3E-O`, and diagonal test. Disjoint pair fibres with identical marginals demonstrate why separately projected witnesses cannot be joined.
* Two configurations with control24 and different read symbols have opposite immediate halting behavior. Control state alone is not a sufficient quotient.

The producer's 2,000,000-step execution is integrity-checked, but **not replayed in full** by this audit. Its later stage clocks beyond this independent30,000-step prefix are not newly established here. The separate LRC cutoff control was inspected in the producer checker but not independently rerun: it is not needed for either arithmetic-machine interface.

## 4. Padding: full proof and guarded tests

Given an ordinary old instruction `(q,b)↦(t,w,d)`, with `d∈{-1,+1}`, introduce fresh numeric states25 and26. Replace the instruction by `(25,w,d)`. At state25 on either bit `c`, copy `c`, move right, enter26. At26 on either bit `c`, copy `c`, move left, enter `t`. No other transition changes.

After the replacement's first step, the tape and head equal the old successor except for control25. The two following moves cancel, both writes preserve the cells they see, and the final control is`t`. This is true even for a left-going original instruction, where the return path can revisit the originally written cell. Hence each selected old step becomes exactly3 steps with the identical old macroconfiguration; all unselected steps remain1 step.

The guards are explicit: the chosen branch is ordinary; the two states are fresh and are not any old target; both bit branches exist; writes copy; bounce moves cancel; return target is the old target; no halt branch is modified. Every fresh-state segment exits in2 more steps, so there is no new internal halting or infinite detour. Induction on old steps yields equivalence of finite halting; decomposition of padded runs into these bounded segments yields the converse. Thus halting/nonhalting is preserved, but microstep clocks are not.

Independent tests cover **all49 ordinary instruction choices**, both directions,2,352 guarded five-cell local patterns and2,352 opposite-read controls, translated to heads-11,0,17 with distant tape sentinels. The chosen `(0,0)` padding reaches both fresh states from the blank tape. A coupled10,000-old-step run matches10,020 padded steps exactly at old boundaries, with10 replaced visits. A separately authored toy that actually halts does so on attempt2 versus4, confirming positive halting behavior as well as continuing-prefix controls. These finite controls corroborate the elementary general argument; they do not replace it.

Twenty-five therefore cannot be an invariant of the halting predicate even after requiring reachable extra states. This constructed27-state machine is not being identified with the historical2016 machine.

## 5. Historical atlas scope

The S62c broadcast literally advertises25+ cross-field connections, but its selected companions do not provide25 individually proved maps: the Weil script labels a15-connection dictionary, the tropical output advertises20+, and its mirror/Gromov–Witten/cluster-character claims lack the needed typed constructions in the inspected material. This is provenance for an atlas request, not a proof count.

HYP-2936 explicitly identifies21 carriers and marks itself a proof-route/synthesis rather than an LRC proof. HYP-2218/2219 provide an honest exact linear pair carrier with an OPEN synthesis scope. The June29 atlas proposes DFA minimization as future work; it does not supply continuation equivalence. The historical inputs are independently compared with repository pin `6a8ad5d3bf01a232da4de41bcbff68390a9337ee` in `../REPO_PROVENANCE.json`: 12 Goldbach-lane inputs are byte-identical and four differ only in trailing newlines; the separate inherited multiperfect snapshot is byte-identical. No historical theorem or topology claim is promoted merely because its file is preserved.

## 6. Boundary

Neither the finite trace nor the number of atlas cards decides Goldbach. A certified Busy Beaver **runtime** upper bound in the identical model would make the finite simulation complete; none is supplied. Known champions, raw state count, finite missing branches, or a bounded successful search are not substitutes. The only remaining formal-verification gap identified here is the explicitly unperformed local Lean build and endpoint axiom audit; it does not invalidate the finite table certificate or elementary padding proof.
