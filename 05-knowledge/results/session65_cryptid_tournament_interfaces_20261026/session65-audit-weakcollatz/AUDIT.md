# Session 65 independent audit: MB orbit, 3-adic cylinders, and weak Collatz

**Verdict: PASS in the stated scope.** The elementary bridge is correct. The
producer's final explicit-gate verifier replays byte-identically in normal and
optimized Python. This audit also independently implements the native five-state
machine; it does not rerun downloaded tinycryptids or Rust software.

Audit target: `../session65-weakcollatz/{REPORT.md,verify_bridge.py,
bridge_results.json,bridge_results_optimized.json,verification_receipt.json}`.
Only this audit directory was written. `session65-repo` was read-only throughout;
no commits, pushes, ref changes, or source edits were made by this auditor.

## 1. Arithmetic and native clock — FINITE-EXACT / VERIFIED

For the partial map

    g(n)=(5n+18+4r)/3, r=n mod 3 in {0,1}; halt for r=2,

the orbit of zero is exactly

    0, 6, 16, 34, 64, 114, 196, 334, 564, 946, 1584,
    2646, 4416, 7366, 12284, HALT.

The native costs, including the final halt branch, are

    15, 73, 277, 907, 2757, 7957, 22777, 64407, 180307,
    504027, 1403967, 3906393, 10861903, 30196527, 24576.

They sum to **47,176,870**. The supplied positive list has only one
multiply-perfect member, **6**. The only displayed live positive successor
identity involving the divisor sum is `sigma(946)=1584=g(946)`. Prime-power
factorization in the independent checker reproduces every producer abundancy
ratio; the producer instead enumerates divisors.

`native_mb.c` is a newly authored, fixed-table C interpreter, separate from the
producer's tinycryptids interpreter. It makes every native transition, verifies
the entire tape at all 15 macro checkpoints, and verifies the complete final
pattern `...01(H0)1(001)^4095 10...`. It returns **4098 ones**, a nonzero tape
span of **12289**, and halts at exactly **47,176,870 transitions**. An additional
121 fresh native configurations `C(n)`, `0<=n<=120`, check all three source macro
formulas and full endpoint tapes. Both C `-O2` and `-O0 -DNDEBUG` builds produce
identical JSON. Explicit checks do not depend on `assert` or `NDEBUG`.

The native checkpoint durations are predictions from the cited formulas, not
independently discovered symbolic formulas. Direct transition execution verifies
the predictions and that no earlier native halt occurred. This proves the
listed run and finite cases, not the general macro formulas or BB maximality.

## 2. Binary cylinders and exact coding — PROVED

Work in `Z_3`. Each inverse branch

    h_r(y)=(3y-18-4r)/5, r in {0,1},

is a bijection onto `r+3 Z_3` and multiplies distances by `1/3`.
For a length-m word `w`, composition of its inverse branches maps `Z_3` onto
one residue class modulo `3^m`. Distinct words give disjoint classes: at their
first differing position, forward images have different residues modulo 3.
Consequently exactly `2^m` classes survive m successful steps.

Nested inverse cylinders for an infinite binary word have a unique limit

    pi(r)=-sum_{j>=0}(18+4r_j)3^j/5^(j+1).

The intersection is the infinite-survivor set K, and `g pi = pi shift`.
The itinerary gives the inverse of pi. If codes first differ at j, then
`v_3(pi(s)-pi(r))=j`; thus pi is an isometry from binary sequence space with
distance `3^(-j)` to K. In particular K is a compact Cantor set carrying the
one-sided full shift. It is not literally the standard ternary-digit Cantor set.

The independent checker uses a forward-composition path rather than the
producer's reversed Fraction inverse recursion: set `B_0=0` and
`B_(j+1)=5 B_j+(18+4r_j)3^j`. Then

    g^m(x)=(5^m x+B_m)/3^m,
    x=-B_m/5^m modulo 3^m.

All residues modulo `3^m`, for `1<=m<=10`, agree exactly with this construction.
The periodic extension of every checked word is `-B_m/(5^m-3^m)`; 2046 such
negative rational points were checked exactly, including their itineraries and
returns. Generally the same formula proves a periodic extension for every
finite word. Constant words give the integer fixed points -9 and -11.

## 3. Measure, dimension, and integer density — PROVED

The m-step survivor set has Haar measure `2^m/3^m`, so K has Haar measure zero.

**Hausdorff dimension, supplied explicitly here:** let `s=log(2)/log(3)`.
The `2^m` depth-m cylinders cover K by sets of diameter `3^(-m)`, which proves
the upper bound `dim_H K <= s`. Give each depth-m cylinder Bernoulli mass
`2^(-m)=(3^(-m))^s`. Any nonempty subset of K of nonzero diameter `3^(-m)` lies
in one depth-m cylinder, so its mass is at most its diameter to the power s.
Singletons have mass zero. Hence every countable cover of K has sum of s-powers
of diameters at least one. This proves the lower bound and therefore

    dim_H(K)=log(2)/log(3).

This dimension statement was absent from the audited REPORT.md, rather than
incorrect there. It may be added with the preceding proof.

Let E be the nonhalting nonnegative integer seeds. For every fixed m, E lies
in exactly `2^m` permitted residue classes modulo `3^m`. The upper natural
density of E is at most `(2/3)^m` for every m, hence is zero. Thus halting has
natural density one. This argument uses integer residue counting, not merely
Haar measure zero.

For integer `N>=1`, choose the smallest `m>=0` with `3^m>=N`. Each permitted
class contains at most one integer of `[0,N)`, so

    # (E intersect [0,N)) <= 2^m <= 2 N^(log(2)/log(3)).

No conclusion that E is empty follows. The native clock is not bounded by the
number of macro steps. Nor do the density limits justify an exchange with a
mean runtime or an expectation limit.

## 4. Tournament orientation and loss ledger — PROVED

Use only the intrinsic unit-part definition from
`01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md`:
`x -> y` iff the 3-adic unit part of `y-x` is 1 modulo 3.

For x,y with the same admissible m-step itinerary,

    g^m(y)-g^m(x)=(5/3)^m (y-x).

Thus each successful step reverses the orientation. For two infinite codes
first differing at j, the common inverse prefix contributes `(3/5)^j`.
The remaining branch difference is `s_j-r_j` modulo 3. Therefore

    unit_3(pi(s)-pi(r))=(-1)^j(s_j-r_j) modulo 3.

Complementing every odd-position digit makes lexicographic order orient every
edge forwards. The tournament induced by K is transitive. The same holds for
the tournament of depth-m safe leaves (or one representative per leaf).
The independent checker verifies 174251 pair orientations at all depths 1..9,
including the first-difference valuations. This finite count supports, but is
not the proof of, the general statement.

The map from integers to the 3-adic carrier preserves safe itineraries and
guarded affine differences. Retaining only the tournament loses integer
membership, designated seed, branch cylinders, native runtime, and the halt
exit. Transitivity is about this pairwise relation, not the temporal successor
graph: K contains fixed points and many periodic points. The bridge cannot
establish E is empty and gives no implication for LRC or multiply-perfect search.

## 5. External sources and precise formalization boundary — CITED

- [Michel, Behavior of Busy Beavers](https://bbchallenge.org/~pascal.michel/beh),
  read 2026-10-26, independently confirms the table, configurations, formulas,
  orbit, final ones count, and total. The page is marked September 2026 and
  credits BB maximality to bbchallenge's 2024 proof. Its Collatz-like section
  still lists the general-input problem as open. These are known-source facts,
  not a new Busy Beaver claim.
- [Aaronson, The Busy Beaver Frontier, section 5.6](https://scottaaronson.com/bb.pdf)
  states Conjecture 20 as eventual halting for every natural seed, and credits
  Michel with its equivalence to the machine halting on every finite input.
  For a nonnegative live input, g(n)>n; thus no nonnegative cycle is possible,
  while nonhalting implies unbounded growth. This is not the ordinary weak
  Collatz cycle assertion.
- [tinycryptids README at the exact pin](https://github.com/calr3/tinycryptids/blob/3c43f97768d4ea06d70e0eb154ba88060bc36dec/weakcollatz/README.md)
  and [proof guide](https://github.com/calr3/tinycryptids/blob/3c43f97768d4ea06d70e0eb154ba88060bc36dec/weakcollatz/PROOF.md)
  were independently read through the GitHub connector. The 14-state binary
  and 10-state ternary designs detect a nontrivial positive Collatz cycle, not
  an arbitrary divergent orbit. Completeness uses the smallest cycle maximum;
  a round trapped below its threshold would imply an earlier cycle. The
  `halts_iff` source declarations have exactly this existential cycle scope.
- All **33** source-manifest paths and hashes were checked against GitHub's
  tree at the pinned commit; all **30** locally downloaded blobs have matching
  Git object hashes. The three documentation paths not locally present are
  correctly labeled as connector-only by the producer. `remote_git_pin.json`
  preserves the independent remote association, not just local hash agreement.
- The union of the two targeted Lean import closures has **15** local files and
  no raw-word matches for `sorry`, `admit`, `axiom`, `unsafe`, or `native_decide`.
  Both main files are imported by the root source. This is only a source check;
  absence of those tokens does not replace elaboration, kernel checking, or
  an evaluated axiom audit. Neither `lean` nor `lake` is available in PATH.
- Additional reproduction limitation: seven other modules imported by the
  downloaded `TinyCryptids.lean` are absent from this partial local snapshot.
  A full `lake build` / root `Check.lean` run would require those sources as
  well as the toolchain. See `normal.json` for their names. No Lean build or
  `#print axioms` run is claimed by either this audit or the producer.
- [The current Cryptids page](https://wiki.bbchallenge.org/wiki/Cryptids) was
  independently read. It distinguishes the weak-Collatz machines, the
  23-state ternary odd-multiply-perfect machine, and the 37-state binary
  odd-perfect machine, with their distinct verification labels. The request
  for historical revision 9619 returned an inaccessible-page error in the web
  tool; the historical revision number is therefore not independently
  authenticated here. The current page supports the substantive scope.
- The pinned local Ligocki excerpt contains the complete named
  `test_weak_collatz_14` rule certificate with the machine string and no
  `Admit` in that test. It is not the full Rust validator implementation and
  was neither compiled nor run here. Source inspection is not an independent
  global proof of that validator's soundness.

The producer's external Python/C outputs are expressly identified as
transcriptions of observed runs, not raw logs. This audit checked source,
receipt consistency, and independent arithmetic/native evidence; it did not
authenticate those transcriptions by rerunning external software. The trusted
task handoff reports the runs occurred. No external-run blocker was triggered
because no such rerun was attempted.

## 6. Replay and sealed evidence

Run from this audit directory:

```sh
gcc -std=c11 -O2 -Wall -Wextra -Werror native_mb.c -o native_mb
./native_mb > native_mb.json
python3 -I independent_audit.py > normal.json
python3 -I -O independent_audit.py > optimized.json
cmp normal.json optimized.json
python3 -I ../session65-weakcollatz/verify_bridge.py > producer_replay_normal.json
python3 -I -O ../session65-weakcollatz/verify_bridge.py > producer_replay_optimized.json
cmp producer_replay_normal.json producer_replay_optimized.json
cmp producer_replay_normal.json ../session65-weakcollatz/bridge_results.json
gcc -std=c11 -O0 -DNDEBUG -Wall -Wextra -Werror native_mb.c -o native_mb_unoptimized
./native_mb_unoptimized > native_mb_unoptimized.json
cmp native_mb.json native_mb_unoptimized.json
```

Observed: every command completed with status 0. Python is 3.12.14; GCC is
Debian 14.2.0-19. Target receipt hashes match current bytes. `manifest.json`
seals target inputs and audit deliverables, excluding itself and generated
executables. No all-input halting, ordinary Collatz convergence, new BB
maximality, independent Lean build, or external novelty claim is made.
