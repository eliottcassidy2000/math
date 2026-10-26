# Signed modular tables, exact carries, and the 14/15 split

**Session 66, 2026-10-26.** Repository source pin:
`63f52936f8dd8592d4f2fee869748baad2e1997e`.
**Status:** elementary propositions proved below; FINITE-EXACT controls in
`check.py`; independent hostile audit pending. No novelty claim. No publication
or change to the pinned repository was made by this lane.

## Findings first

1. The recovered precise alternation is **between successive inverse-fibre
   exponents**, not automatically between arbitrary odd moduli. If a signed
   odd Collatz edge has target `q`, raising its valuation by one changes
   `(n,b)` to `(2n+b,-b)`, with the same integer target. The sheet cannot be
   discarded when continuing forward.
2. On odd source residues modulo `210=14*15`, with the sheet retained, this
   lift is a permutation with **seven 6-cycles and fourteen 12-cycles**. The
   source transformation `z=3n+b mod315` makes it ordinary doubling on
   residues prime to three. The factor five controls the 6/12 split; seven
   distinguishes orbits but adds no period because its order three is already
   present in the mod-nine coordinate.
3. The signed pair-orbit defect machinery from Session 65 gives minimum
   orientation defect **21** for this permutation. If antipodal vertex pairs
   are required to be actual clone modules, the exact minimum instead is
   **77**. The lost term is `4*14=56`: the weighted quotient's even-cycle
   obstruction. This is an explicit, quantified clone-merger failure, not
   an identification of this doubling map with arithmetic Frobenius.
4. Neither a residue table nor its clone quotient determines a native
   valuation word, positive-integer orbit, terminal certificate, or descent.
   Even `1` and `211`, equal modulo both 14 and 15, have plus-sheet valuations
   two and one and odd successors `1` and `317`.

## 1. Inheritance, corrections, and live concepts

The seven exact input blobs and SHA-256 hashes are in `SOURCE_MANIFEST.json`.
All source references below mean those bytes at the pin, not an unpinned
current branch. The source titles identify separate files even where theorem
numbers overlap elsewhere in the repository.

- **Closest proved mechanism:** `signed_collatz_pair_lifts_20261014d.md`,
  §§1–3, gives the complete two-sheet inverse fibre, exact valuation, the
  alternating lift, and the false-seam example `U_+(3)=U_-(7)=5`.
- **Prior modular rows:** `arithmetic_braids_20260917_collatz.md`, §§1–5,
  and `collatz_mod6_20260917_row_braid_typing.md`, §§1–4, already give
  the three source/image rows, the same-sheet `4n+b` lift, the distinction
  between mod-seven row information and the mod-three gate, and ordered
  carry cylinders. These are not new here.
- **Signed tournament action:** THM-4754,
  `three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md`,
  supplies the least nonzero ternary digit orientation and its doubling
  anti-action. Its same-branch hypotheses must survive any forward use.
- **Defect machinery:** Session65 `session65-fifteens/GENERAL_PERMUTATION.md`
  and `REPORT.md` give the signed pair-orbit enumerator and its Frobenius
  specialization. We use the elementary permutation theorem, not the
  arithmetic genus identification. The argument needed here is reproved.
- **Clone boundary:** THM-4463,
  `homogeneous-tournament-contacts-and-pinned-scale-floor.md`, defines
  clones by exterior rows, not endpoint equality. In a transitive triple,
  the endpoint-deleted zero-contact relation is not transitive. We therefore
  specify an actual disjoint partition into homogeneous pairs below rather
  than taking an unjustified equivalence closure.

The canonical hostile is a shared endpoint with different sheet continuation.
The corrected near miss is integrality versus *exact* valuation: one extra
binary bit is required. The least-used sidecar is the signed twist of the
internal edge in a merged pair.

**Concept board:** inverse-fibre clock / modular carry / exact binary guard /
common CRT sheet / signed orientation defect / clone-module quotient.
The anchor is the signed Collatz interface; the niche is the 14/15 split;
the wildcard is the extra defect forced by clone modules. Each bridge below
states what those coordinates preserve and what they destroy.

The targeted inheritance pass did **not** locate a literal theorem saying
that all odd-modulus multiplication tables alternately equal the two Collatz
graphs. The exact interpretation below is therefore stated explicitly;
it does not silently attribute a broader table theorem to the user or canon.

## 2. The alternating sign is a genuine inverse-fibre law

For `b in {−1,+1}` and a nonzero odd integer `n`, put

`a_b(n)=v_2(3n+b)`, `U_b(n)=(3n+b)/2^{a_b(n)}`.

The valuation uses the absolute numerator, so signed sources are permitted.
Since `3n+b` is nonzero and even, `a_b(n)>=1`.

Fix a nonzero odd target `q` with `3` not dividing `q`. For each `a>=1`
there is one sheet `b_a in {−1,+1}` with `b_a=2^a q mod3`. Its source is

`n_a=(2^a q-b_a)/3`.

It is odd, and `3n_a+b_a=2^a q` proves the exact valuation. Conversely,
every signed integer edge ending at `q` has this form. No such edge enters a
target divisible by three. Since `2=-1 mod3`,

`b_{a+1}=-b_a`, `n_{a+1}=2n_a+b_a`.                 (1)

Thus the inverse-port lift and its square are

`L(n,b)=(2n+b,-b)`, `L²(n,b)=(4n+b,b)`.           (2)

More generally, the ordered lift has its exact carry

`L^j(n,b)=(2^j n+b(2^j-(-1)^j)/3, (-1)^j b)`.    (3)

It preserves the integer endpoint and raises the edge's valuation by `j`.
It is **not** a forward orbit conjugacy. The inherited hostile is

`U_+(3)=U_-(7)=5`, but `U_+(5)=1`, `U_-(5)=7`.

Negation is a different exact operation:
`U_{-b}(-n)=-U_b(n)`. It preserves forward time and the valuation; it is
not reversal of graph arrows. These two meanings of a sign change must not
be merged.

## 3. An arbitrary odd-modulus table needs both carry and guard

Let `m>=3` be odd, choose `r,s in {0,...,m-1}`, and prescribe `a>=1` and
`b=±1`. The modular branch is

`s = 2^{-a}(3r+b) modm`.

Define the **signed integer carry**, which may be negative, by

`c=(3r+b-2^a s)/m`.

Write a prospective integer source as `n=r+mt`. Then the complete lift is

`y=s+m(3t+c)/2^a`.                              (4)

It is an exact odd-to-odd `U_b` edge of valuation `a` if and only if

`3t+c = (1-s)2^a mod2^{a+1}`.                    (5)

**Proof.** Equation (4) is simply `3n+b=2^a y`. Integrality requires
`2^a | 3t+c`; after division, oddness of `y` requires
`(3t+c)/2^a = 1-s mod2`. Combining them gives (5), and conversely (5)
makes `y` odd and forces exact valuation. This also forces `n` odd.
For positive `n`, both sheets have positive `y`; signed heights retain their
sign instead of being replaced by a positive representative. QED.

Equivalently the exact source guard is

`n = 3^{-1}(2^a-b) mod2^{a+1}`.                   (6)

The inverse exists because three is odd. CRT joins (6) to *every* residue
`n=r modm`, with infinitely many positive representatives. Thus an odd
modulus never determines the valuation. In fact `U_b` does not descend to
a function on `Z/m`: take `r=0` and sources satisfying respectively `a=1`
and `a=2`. Their target residues are `b/2` and `b/4`, distinct modulo every
odd `m>1`. A residue graph retaining *all* labelled branches is legitimate;
a deterministic quotient pretending the guard is known is not.

For an exact successor modulo `210`, source residue modulo `105` together
with (6) suffices: the target modulo105 follows by invertible halving and
its parity is one. The combined source period is `105*2^{a+1}=210*2^a`.
The bare 14/15 table carries only the first parity bit, not this whole guard.

**Hostiles.** `1` and `211` have identical residues modulo14 and15, but
`U_+(1)=1, a=2` and `U_+(211)=317, a=1`. Dropping the top bit of (6)
also fails systematically: if `n` has exact valuation `a`, then
`n+m2^a` has valuation strictly larger than `a`, while it satisfies the
same coarse source congruences modulo `m` and `2^a`.

## 4. The {3,5}/{7,2} split as an exact signed-port permutation

Use the finite state space

`P={(n mod210,b): n odd, b=±1}`; `|P|=210`.

The map `L` of (2) preserves this space and is bijective: modulo105 its
linear coefficient two is a unit, and an odd residue modulo210 is the
unique odd lift of its residue modulo105. Define

`Phi(n,b)=3n+b mod315`.

This is a bijection from `P` to the 210 residues modulo315 prime to three.
Indeed the residue modulo3 determines `b`, after which division by three
recovers `n mod105`, hence its unique odd lift modulo210. Directly,

`Phi(L(n,b))=3(2n+b)-b=2(3n+b)=2 Phi(n,b) mod315`. (7)

**Why 315 instead of 210?** Recording the signed carry `3n+b` uses one
additional factor of three, while known oddness removes the factor two.
This is a consequence of the map, not a substitution of unrelated moduli.

Since `315=9*5*7`, CRT and the elementary orders

`ord_9(2)=6`, `ord_5(2)=4`, `ord_7(2)=3`

give the full cycle classification. Every `Phi` coordinate is a unit
modulo9, so its period contains six. If `5|Phi`, its mod-five component
is fixed and the total period is six. Otherwise it contributes four and
the period is twelve. The mod-seven component has period one or three,
both already dividing six.

There are `6*1*7=42` states in the first case and `6*4*7=168` in the second:

`cycle type(L on P)=6^7 12^14`.                    (8)

For an actual inverse port, `Phi=2^a q` modulo315. Thus `5|Phi` is exactly
`5|q`. This retains target divisibility and a finite valuation-clock phase;
it does not determine the target's subsequent orbit.

### The two smaller tables and the common-sheet seam

On `(n mod14 odd,b)` the cycle type is `2^1 6^2`; on `(n mod15,b)` it is
`6^1 12^2`. Their exact joining is the **fibre product over the same b**:

`P = P_14 ×_{ {−1,+1} } P_15`.

The two integer coordinates join by ordinary CRT. Independently choosing a
sheet in each factor gives `14*30=420` states, of which half have mismatched
sheets and do not encode one signed source. Under the joint lift the sheet
equality or inequality is preserved; it is not repaired automatically.

This 210-vertex signed-port carrier is distinct from a divisor-complement
carrier built from the same prime set; no defect count transfers by the shared
number 210 alone.

This proves a useful role for the split `{3,5}` versus `{7,2}`: 15 retains a
ternary carry and the order-four five-clock, while 14 retains oddness and
the seven-clock. It does not privilege this split as a universal arithmetic
classification. The detailed preservation statement is (7)–(8).

## 5. Ternary version and the distinction between native and optimized defect

For any `k>=0`, signed ports `(n mod3^k,b)` map bijectively by `3n+b`
onto `(Z/3^{k+1})^×`, again intertwining `L` with doubling. The elementary
lifting law `ord_{3^{k+1}}(2)=2*3^k` shows this is one full cycle.

Restrict THM-4754's tournament to those units: orient `x->y` when the first
nonzero ternary digit of `y-x` is one. Multiplication by two negates that
digit, so `L` reverses every arc; its square preserves every arc. This is a
native **anti-automorphism**, not an ordinary automorphism.

For a vertex permutation `sigma`, define `D_sigma(T)` as the number of
unordered pairs whose orientations are not preserved by sigma. The inherited
signed pair-orbit proof is short: in each orbit of unordered pairs, the
parity of edge failures is even if the endpoints return individually and
odd if they return swapped. Every allowed failure pattern has exactly two
orientation lifts. The swapped-endpoint orbits occur exactly once per even
vertex cycle, namely its antipodal pairs. Consequently

`min_T D_sigma(T)=e(sigma)`,                       (9)

where `e` counts even vertex cycles. The full enumerator is

`prod_O [(1+z)^{|O|}+epsilon_O(1-z)^{|O|}]`,

with `epsilon=-1` exactly on those antipodal pair orbits. In particular,
for the one ternary signed-port cycle the minimum is **one**, whereas its
native tournament has defect **all** `binom(2*3^k,2)` pairs.

This is the precise limitation of importing the Frobenius-defect language:
the optimization in (9) is over *all* tournaments. It does not supply a
preferred arithmetic orientation. The Session65 arithmetic specialization
also needs its actual Frobenius involution and geometric theorem. `L` is
not asserted to be that Frobenius.

## 6. Clone-merging preserves defect only with its weights and internal twist

This supplies a quantitative connection to actual clone machinery rather
than calling equal Collatz endpoints clones.

Let sigma have only even vertex cycles. In each cycle of length `2r`, form
the `r` disjoint antipodal pairs `{v_i,v_{i+r}}`. Require every such pair
to be a **homogeneous pair module** of a tournament `T`: each exterior
vertex sees the two members identically. This is a genuine allowed family
of tournament substitutions. It need not be the maximal canonical clone
partition of THM-4463; no maximality or uniqueness claim is needed.

The quotient tournament `Q` has one vertex per pair. Sigma induces a
permutation `pi` whose cycle lengths are the original lengths divided by
two. Every quotient arc represents four original cross-pair arcs. Their
orientation failures agree, while the internal edges are disjoint. Thus

`D_sigma(T)=4 D_pi(Q)+I_sigma(T)`.                (10)

The internal pair edges in each original vertex cycle form one orbit with
negative twist: after `r` applications sigma swaps the endpoints. Therefore
that orbit has an odd failure count, whose minimum is one. Its edge
orientations are freely selectable independently of the quotient. If `c`
is the number of vertex cycles of sigma, (9) and (10) prove the sharp formula

`min_{antipodal-pair modules} D_sigma(T)`
`   = c + 4 e(pi)`
`   = c + 4 #{sigma cycles whose length is divisible by4}`. (11)

**Attainment:** choose an orientation of `Q` attaining its permutation
minimum by the signed pair-orbit construction, and independently select
one internal failure in each negative pair orbit. All homogeneous-pair
conditions still hold, since internal orientations do not affect exterior
rows. This proves both the lower bound and its attainability.

### The new 210-state invoice

By (8), `c=21`. The antipodal-pair quotient has cycle type `3^7 6^14`,
so `e(pi)=14`. Hence

`unrestricted minimum =21`,
`antipodal-clone-module minimum =21+4*14=77`.      (12)

The missing 56 is a **weighted quotient defect**, and the remaining 21 is
the internal signed twist. Neither may be thrown away at the merger.
The smallest strict hostile is a single four-cycle: unrestricted minimum
one, but its two antipodal clone pairs force minimum five (four cross arcs
plus one internal failure). A six-cycle instead permits minimum one because
its three quotient blocks admit an invariant cyclic tournament. This isolates
the actual boundary: divisibility of the cycle length by four, not its
being even alone or the sign of the entire permutation.

No claim is made that the native Sierpinski orientation or an actual Collatz
graph has these antipodal clone modules. They are a precise constrained
orientation model whose cost is computed, with a hostile showing why
arbitrary merging cannot preserve the unrestricted answer.

## 7. Connection ledger and proof boundary

| Source → target | Map / preserved predicate | Information lost / needed sidecar |
|---|---|---|
| Integer inverse edge → next inverse port | (1); same exact integer target | Sheet and valuation change; continuation needs its own sheet certificate |
| Exact odd edge → residue table | (4)–(6); exact edge iff guard holds | Carry and binary guard cannot be discarded |
| 14/15 tables → 210 signed states | CRT with common sheet | Independent sheets are spurious; height and valuation still absent |
| Signed ports → doubling modulo315 | (7); complete finite permutation | Only a quotient of the infinite inverse fibre, not integer orbit recurrence |
| Permutation → signed pair-orbit defect | (9); minimum over all orientations | Arithmetic orientation, primality, genus, and descent are not retained |
| Antipodal clone modules → quotient | (10); exact weighted defect | Internal negative twists and block weights must remain |

The complete inverse fibre can have full modular support while its common
target enters a different basin. A cycle of the finite lift changes the
integer source by a nonzero carry and raises its valuation; it is not an
integer Collatz cycle. Finite CRT closure does not prove a terminating
choice of guards. Collatz convergence, integer-cycle classification,
LRC(14), and unrestricted tournament connectivity remain OPEN.

## 8. Exact controls and reproduction

Run from this packet directory:

```
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
```

The checker uses only Python integer arithmetic and explicit exceptions,
not optimization-disabled assertions. Its declared universes are:

- both sheets; every odd endpoint `q` in `[-101,101]` prime to three;
  valuations 1 through20;
- every odd modulus3 through101, every residue, both sheets, and exact
  valuation1 through8; independent CRT construction, the signed-carry
  identity, lifted height, exact guard, and one missing-top-bit hostile;
- complete signed-port permutations at 14, 15, 210, their cycle census,
  CRT reconstruction, and the `3n+b mod315` intertwiner;
- all pair comparisons in the native ternary unit tournament for
  `k=0,...,4`, and the complete signed-port clock at each level;
- complete module-constrained tournament enumeration for vertex-cycle
  types `(2),(4),(6),(8),(10),(4,2),(4,4),(6,2)`, comparing direct edge
  failures with (10) and checking the sharp minimum (11);
- inherited false-sheet seam, the `1/211` guard hostile, and SHA-256
  verification of every archived source blob.

Both executions passed and are byte-identical: **261,902 exact checks**,
including 35,978 complete module-constrained orientations. The infinite
statements follow from the proofs, not extrapolation of these finite runs. The 210-state **unrestricted** and constrained minima follow
from the proved formulas and its complete cycle census; all `2^21945`
tournaments are emphatically not enumerated. No Lean/Coq build is claimed.
