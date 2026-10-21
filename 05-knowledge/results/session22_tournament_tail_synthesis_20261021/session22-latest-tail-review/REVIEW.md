# Final bounded review of incoming linked tail gates and prefix reserves

2026-10-21. PASS, with the scope distinctions below. Incoming mathematical work is attributed to `tail_geometry_20261021c` and `uniform_prefix_20261021c`, introduced at incoming commit `359ddd20`. This review contributes proof checking, a separate support-normal arithmetic consumer, and fresh checker replays, not the incoming theorems or native discovery. Observed checkout HEAD was `336b76f137a6bc4843a9b1cb702132ad41574d58`; file hashes, rather than a moving shared checkout name alone, identify the reviewed inputs.

## PROVED: linked majorant, optional absence, and exact cluster bookkeeping

For the positive definite Q and its polarization beta, g-f=lambda^-1 h and g+f=lambda(2w+h). Therefore Q(g)-Q(f)=beta(g-f,g+f)/2=Q(h)+beta(w,h)=delta(w,h). Polarizing this exact quadratic identity gives the stated kappa identity. Both channels share a sign and remain attached to their actual same-depth current node; the coordinate change does not modify native dynamics.

The destination lemma is sound for any disjoint assignment, with optional absence: independently randomizing signs of entire terminal sums makes the expected Q of their total equal the sum of terminal Q values. Completing absent generators with independent signs adds nonnegative Q diagonals. Each complete signed sum is bounded by the global signed maximum. Thus dropping the negative f norm provides an upper bound, not an equality or a new physical energy. Current-port interactions are separately bounded by sum_i max_t |kappa(c_t,b_i)|. Empty marked ports remain permissible destinations. Convexity applies to positive Q on the projected planar zonotope, not to indefinite delta on four-space. The sorted-edge implementation covers arbitrary directions, including repeated/collinear directions and zero generators.

Cluster maxima include absence, signs and destinations. Within-cluster terms are optimized jointly; cross-cluster unordered pairs receive one absolute bill each. Merging clusters cannot increase the resulting bound at a fixed state. The exact two-block solver's cases are exhaustive: neither, either, common terminal with all four sign choices, or distinct terminals. For distinct terminals the top two support destinations in each list suffice, since at most one best destination is forbidden by collision. Greedy pairing is only a choice of disjoint beneficial pairs, not a maximum-weight matching claim. The consumer rejects repeated actual-node owners and recomputes each claimed local gain.

The scalar hostile b=((-5,0),(8,0)) has delta=-16, so absence beats the two full-sign endpoints. Two pending copies give projected lambda2 bound 0; after one arrives the same freshly evaluated bound becomes16. This is not a soundness defect: an older bound remains valid for every later cut. A running minimum of retained, authenticated all-future certificates is nonincreasing; including the inherited monotone U preserves this. No claim of individual angular-bound monotonicity follows. Earlier cuts require their own certificates. The minimum selects a complete bound and never adds overlapping credits. For a fixed ownership ledger, Delta=A+R'+J' and A+M=153N give R'<=B-A-J'; B<=153N gives precisely R'<=M-J'.

## FINITE-EXACT arithmetic; CITED native provenance

The separate `audit.py` imports no incoming checker. Its support-normal/infinitesimal-tangent vertex construction does not sort angles. On the frozen pending blocks it independently recovers these lambda1 values:

- N4096,D105: 421276 + 74092 + 56476 = 551844; inherited U557164.
- N8192,D110: 793715 + 152203 + 172963 = 1118881; inherited U1138407.
- N32768,D130: 2838302 + 649176 + 654736 = 4142214; inherited U4241482.
- N7,D12: lambda1 value1153; incoming disjoint-pair certificate improves this to1152, still above budget1071.

All positive rows have zero aggregate pending old vectors; this is not a statement that every original old source has arrived. The independently recovered angular gains are5320,19526,99268. No new N or new tail cut was sampled. The earlier-cut predicates and tail source/edge authentication are CITED from the hash-pinned `future_tail_blocks_20261021b/receipt.json`, not independently rerouted in this review. Current ports, aggregate blocks, saved source witnesses and receipt consistency are checked by the incoming consumer; a hash is not itself a proof of provenance. Combining the new arithmetic with the existing earlier-cut certificates preserves only the same finite all-horizon scale conclusions.

## PROVED: moving potential debit and excursion summary

Psi=Delta+sum_pending max(0,delta(b)) is at least Delta. At depth0 each literal source has its own node: old-only delta is0, new-only delta is its original Q norm, and initially terminal sources contribute the same quantity to Delta. Hence Psi0=153N conditional on the inherited exact annulus norm identity; this review does not newly prove that all-N charge identity.

Group every pending node by its actual successor with its actual clock sign. For a terminal successor the delta increment is delta(c+v)-delta(c); for a pending successor the new positive part is delta(v)_+. Subtract each predecessor's positive part exactly once. Summing yields the exact stated debit Psi_(D+1)-Psi_D. Both-channel-zero erasure is harmless; total-charge-zero alone is not. This is a stronger moving potential, not an all-future certificate just because it is below budget at one cut.

For a finite debit word, T is total and P the maximum of0 and all prefix sums. Concatenation is (T1+T2,max(P1,T1+P2)); associativity follows by expanding to max(P1,T1+P2,T1+T2+P3). The identity(0,0) is an identity on valid summaries, whose domain satisfies P>=max(0,T), not arbitrary unrestricted pairs. A block with initial reserve B pays every cut iff P<=B, and ends with reserve B-T. Nonpositive net T alone misses positive internal excursions. Infinite-horizon coverage requires finite consecutive blocks with unbounded endpoints and the payment condition at every block, or an authenticated all-future stopping certificate. Existence of these blocks uniformly in N remains OPEN.

## Fresh native replay, strictly bounded

Both normal and optimized runs of the incoming prefix checker independently route all393216 original primitive addresses at N4096 to exactly cuts17 and18. They recover Delta166652, positive moving part268827, Psi435479 at17; Delta184903, moving254791, Psi439694 at18. Global increase4215 is positive, while both are below626688. Terminal811 local debit2004 uses actual children1081(clock2) and2163(clock3). This refutes pointwise local/global monotonicity of Psi, not153N coverage.

The saved discovery file contains a0..150 run. The checker reads its saved values and checks their excursion arithmetic, but only routes original sources independently to17 and18. No discovery sweep was rerun here. The maximumPsi=626688 over0..150 remains a CITED producer observation plus saved-data arithmetic, not a fresh full-range native replay.

## Replay artifacts and stopping boundary

Run from the workspace root:

    python3 -B session22-latest-tail-review/audit.py math
    python3 -B -O session22-latest-tail-review/audit.py math
    python3 -B math/05-knowledge/results/tail_geometry_20261021c/check.py
    python3 -B -O math/05-knowledge/results/tail_geometry_20261021c/check.py
    python3 -B math/05-knowledge/results/uniform_prefix_20261021c/check.py
    python3 -B -O math/05-knowledge/results/uniform_prefix_20261021c/check.py

All checks use explicit exceptions and remain enabled under optimization. Replays omit --write and do not alter incoming artifacts. Fresh stdout receipts and stderr captures are alongside this review; the saved incoming full receipt hashes are recorded in source_pins.json. This bounded review performs no Git or remote write, no new benchmark, no scale search, and no discovery full sweep. Uniform153N coverage, uniform stopping/payment, and arbitrary Collatz entry remain OPEN. No mathematical blocker was found.
