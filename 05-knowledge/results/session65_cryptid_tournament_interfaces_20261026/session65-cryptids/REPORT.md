# Three total machines: exact guarded macros and a tournament tape observer

2026-10-26. Research at repository base `6a8ad5d3bf01a232da4de41bcbff68390a9337ee`.
No repository files were edited by this subtask. Elementary proofs below are
self-contained derivations from the literal transition tables; their finite
controls are separate. This is not a claim of priority for the macros.

## 1. Identification and status boundary

The requested names are interpreted as the canonical examples in the
[wiki's pinned nonhalting zoology, revision 8776](https://wiki.bbchallenge.org/w/index.php?title=Non-halting_Turing_machine&oldid=8776):

| Name used here | Exact 4-state / 2-symbol machine |
|---|---|
| Spaghetti sandwich | `1RB1LC_1LA1RD_1RA0LC_1LB0RD` |
| Bouncer + bell | `1RB0LA_0RC1LA_1RD0RA_0LB1RB` |
| Typical irregular bell | `1RB0RC_1LC1RA_1RA1LD_0LC0LA` |

The zoology treats spaghetti and irregular bells as qualitative classifications;
it does not supply a complete asymptotic proof for these particular examples.
The third is one example, not a uniquely defined meaning of “irregular bell.”

**PROVED directly:** all eight entries of each machine are defined and target
one of A,B,C,D. Each therefore runs forever from every tape, not merely blank.
Their difficult questions concern morphology, recurrence and growth, not
halting. In particular they should not be advertised as undecided Cryptids in
the [Cryptids page's blank-tape-halting sense](https://wiki.bbchallenge.org/wiki/Cryptids).
The generic machine-page “Undecided” status is not mathematical evidence
against this elementary totality observation.

[Primary discussion by Yongyi Chen, 21 Oct 2024](https://groups.google.com/g/busy-beaver-discuss/c/lcr_6buFz_8/m/e2BiAob8AgAJ)
reports a bounded search on a list including the sandwich: no period up to
10^11 after 8 trillion steps; the suggested behavioral equivalence with
`1RB0LA_1RC1LA_1LB1RD_1LC0RD` is explicitly tentative. These are historical
computational reports, not our replays or proofs of nonperiodicity.

The relevant positive-control literature is more precise:
[Bouncer revision 8027](https://wiki.bbchallenge.org/w/index.php?title=Bouncer&oldid=8027)
gives an inductive shift/bounce analysis of `1RB0LB_1LA0RA`;
[Bell revision 6789](https://wiki.bbchallenge.org/w/index.php?title=Bell&oldid=6789)
gives example transcripts, distinguishing ordinary and cubic bells. Neither
example is silently substituted for the requested 4-state bouncer+bell.

The primary Coq source
[`Bouncers.v` at 1c23203761602c199a7ac5ea6f7947fb80851d93](https://github.com/meithecatte/busycoq/blob/1c23203761602c199a7ac5ea6f7947fb80851d93/verify/Bouncers.v)
contains `shiftrule_left_exec`, `steps_some`, and
`verify_bouncer_correct`: successful certificates yield native positive
progress and then nonhalting. Its `Repeat` constructor denotes positive
repetitions, while alignment and subsumption preserve membership by one-way
inclusion, not language equality.
The pinned file was fetched and compared equal to the current fetched file;
we did not build Coq, audit all imports/axioms, or produce certificates for
our three machines. Thus formal-verification claims here remain **CITED**.

## 2. Native convention and all-parameter macro theorem

A configuration is `(S,h,q)`, with finite `S subset Z` the 1-cells, head
position `h`, and control `q`. Symbol 0 is everywhere outside S. Reading
`1RB` means write 1, move right one cell, enter B. Time counts native
instructions; no halting-instruction convention is needed. Unmentioned cells
are unchanged. The statements also hold on arbitrary tapes whenever the
specified finite run exists; finite support guarantees every selected run
ends. Counts refer to *maximal* matching runs.

### 2.1 Sandwich: complete section states A,B

Let `d(A)=-1`, `d(B)=+1`. At state `q in {A,B}`:

* If the current cell is 0, write 1, move to `h-d(q)`, and switch A/B.
  Time 1; changed bits 1; change in number of ones +1.
* If the current cell is 1, let k>=0 be the number of consecutive 1-cells
  at `h+d(q),...,h+k*d(q)`, so `h+(k+1)*d(q)` is 0. Erase those k cells,
  put 1 at that terminal zero, keep the original head cell 1, and return
  to state q at `h+k*d(q)`. Time k+2; changed bits k+1; ones change 1-k.

**Proof.** A1 enters C one cell left; C1 erases and continues left until
C0 writes 1 and returns right to A. B1 similarly enters D right; D1 erases
right until D0 writes 1 and returns left to B. The other two section
branches are the literal A0/B0 entries. This proves every k by induction.
The k=0 branch moves no net distance but changes the adjacent zero to 1;
it must not be discarded as an idle step. All macro clocks are positive.
Starting from blank A, the section is closed and exhaustive. This is a
complete variable-word acceleration, not a bounded-dimensional recurrence
or proof of the “sandwich” morphology.

### 2.2 Bouncer+bell: complete section state A

If A reads 1, let k>=1 count the consecutive 1-cells starting at h and
continuing left. Erase all k; return A at h-k after k steps. Changed bits k;
ones change -k.

If A reads 0, let k>=0 be maximal such that the right tail starting at
h+1 begins with `(001)^k`. Write 1 at h and replace those k copies by
`(011)^k`. Put `j=h+1+3k`. Exactly one following guard holds:

| Tail beginning at j | Remaining edits | Final head/state | Total native time | Changed bits | Ones change |
|---|---|---|---:|---:|---:|
| `1` | none | j-1, A | 3k+2 | k+1 | k+1 |
| `01` | erase cell j+1 | j+2, A | 3k+3 | k+2 | k |
| `000` | write 1 at j+1 | j, A | 3k+5 | k+2 | k+2 |

**Proof.** After A0, the head is B at h+1. In B, the four prefix cases
are `1`, `01`, `000`, `001`. The first takes B1 back to A; the second
uses B0,C1; the third uses B0,C0,D0,B1. The last uses B0,C0,D1, transforms
001 to 011, and returns B three cells right. Induct over its repetitions.
The four prefixes are disjoint and exhaustive; an infinite repeated 001
suffix is excluded by finite support. No coordinate was deduced from a
plotted envelope. In particular the guard `000` cannot be weakened to `00`.

### 2.3 Irregular bell: complete section states A,C

At A the cases are:

| Guard | Edits | Final head/state | Time |
|---|---|---|---:|
| current 1 | erase h | h+1,C | 1 |
| current 0, right 1 | write 1 at h | h+2,A | 2 |
| current 0, right 0 | write 1 at both h,h+1 | h,C | 2 |

At C, repeatedly pass left across pairs 01 with the head initially on
the pair's 1. Let k>=0 be maximal, and j=h-2k. The traversed `(01)^k`
block is unchanged. The two remaining cases are:

| Guard | Edit | Final head/state | Time |
|---|---|---|---:|
| cell j is 0 | write 1 at j | j+1,A | 2k+1 |
| cells j-1,j are 11 | erase j-1 | j-2,A | 2k+2 |

**Proof.** A0 uses B0 or B1 according to its right symbol. A1 is literal.
The two-instruction C1,D0 loop crosses a 01 pair left without changing it.
If that loop stops, either C0 writes 1 and enters A right, or C1,D1 erases
its left neighbor and enters A left. Finite support forces a stop. This
proves closure and clocks without claiming asymptotic regularity.

These macro systems are exact on their stated sections. They preserve the
full tape, head and control; replacing the retained word by a few counts
would be a different map needing a new proof of closure.

## 3. A precise, noncosmetic tournament connection

### 3.1 Isometric tape encoding

Let `G0` be the transitive tournament 0->1->2, 0->2. Let `G1` be the
triangle obtained by reversing only 0->2. For a binary word of length m set

    E(w)=G_(w1) join ... join G_(wm),

where all arcs between blocks point toward the later block. For fixed m,

    d_iso(E(u),E(v)) = sum_i |u_i-v_i|.                 (T1)

**PROVED using inherited THM-2183.** Each factor has equal order 3;
d_iso(G0,G1)=1 since they are nonisomorphic and differ by one literal arc.
Iterate the order-join metric product theorem. Its distance is the unmerged
isomorphism quotient, with no complement identification. It permits arbitrary
factor tournaments and requires equal corresponding orders; both hypotheses
are met. There is no claim of the same formula for arbitrary cyclic quotient
substitutions; THM-2195 explicitly excludes that universal extension.

Dependency pins, both at the base commit above:

* `01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md`,
  SHA256 `00c34a935dc73191edbd7088e1bd9ce2029af23c9d89c9889e1af5b5425499bc`.
* Scope hostile: `01-canon/theorems/THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md`,
  SHA256 `60ddaa62c9d6e594b994e2e9c7ab2600edf9c763146c950d5ebf009ab025a006`.

**Decoding even after zero-block mergers.** The strong-component order of
E(w) is canonical. Every 1 contributes one C3 strong component; every 0
contributes three consecutive singleton components. Each intervening run of
3r singletons therefore decodes as r zeros. This proves injectivity on words
without pretending that each original zero block remains a maximal module.
Adjacent zero blocks form a larger transitive interval; a 0->1 write can
split it, and a 1->0 write can merge neighboring transitive intervals.
The literal labels and cuts remain available to execute a write.

### 3.2 Dynamics preserved only with its sidecar

Choose a fixed finite tape window containing every read/write/head position
of the finite trajectory under consideration. Include every cell in that
window, including endpoint zeros. The record

    (E(tape word), window origin, head index, machine control, program)

recovers its native configuration within that window. Each native write
changes zero or one distinguished internal arc; the head move updates a
marker. Thus the native map has an exact realization as a **controlled,
marked, local arc rewrite**, with read guards and exit-from-window checks.
No single fixed finite window contains an arbitrary infinite execution.
Padding two compared words identically preserves (T1), but deleting a
zero prefix loses the physical origin unless retained separately.

This is not an intrinsic autonomous tournament dynamics and **not** the
repository's K_P operation, which reverses all non-path arcs. The sharp
n>=8 closed endpoint-four three-black theorem cannot be invoked without
producing actual directed Hamiltonian words, successive-source blackness,
and the correct literal K target. Those obligations are absent here.
The useful transfer is a faithful static carrier, exact edit metric, and
explicit clone split/merge examples—not a black-core connectivity proof.

### 3.3 Scalar collapse and clock hostile

Let H count directed Hamiltonian paths and c3 count directed triangles.
Order joins force every Hamiltonian path to traverse the blocks in order,
and no directed triangle crosses blocks. Hence, directly,

    c3(E(w))=#ones(w),     H(E(w))=3^(#ones(w)).         (T2)

The macro tables therefore give exact c3 increments and H ratios 3^Delta.
Nevertheless 01 and 10 have equal `(c3,H)=(1,3)` while their tournaments
are at distance 2. At head 0 in bouncer+bell state A, 01 takes two steps
to 11 with head 0, whereas 10 takes one erasing step to 00 with head -1.
Thus these scalar invariants do not even determine the next macro.

A sharper metric/clock separation comes from the irregular machine. Take
C on the rightmost cell of `0(01)^k`. Its sweep changes only the initial
zero to 1. Thus the input and output tournament tapes have distance 1,
independently of k, but the native clock is 2k+1. Include a padding cell
for the final head when k=0. This all-k example disproves any bound on
native macro time depending solely on tournament edit distance.

## 4. Reproducible exact controls and limitations

Run from this directory:

    python3 -B check.py > normal.json
    python3 -B -O check.py > optimized.json
    cmp normal.json optimized.json

`check.py` independently interprets raw transition strings and compares its
native outputs to the algebraic macro implementation. It tests every binary
word of length n=0..10, every head position -1..n, and all section states
(A/B, A, A/C respectively), with zeros elsewhere. It checks literal tape,
head, control, native time, changed-bit count and ones increment. Additional
controls use every k=0..200 for all bouncer+bell exits and both irregular
sweep exits: **113,645 total macro cases**. No blank-orbit reachability
filter is imposed on this universe. A separate blank trajectory passes
10,000 native steps for each machine, replaying every macro.

The tournament check exhausts all pairs of words of length<=2 and all
vertex permutations up to order6, comparing unlabeled distance with Hamming;
it also counts all Hamiltonian paths and cyclic triples. The all-length
claims rest on the proofs, not this small computation. Explicit negative
controls reject dropping the 000 guard and forgetting the machine state.
The strong-component decoder additionally checks every word through length8
after reversing all vertex labels (511 words); a literal pair across a zero
block boundary is verified to lose its twin relation under a one-bit write.
Normal and optimized outputs agree byte-for-byte.

No infinite-orbit classification, nonperiodicity theorem, Coq build, new
cryptid halting resolution, or tournament K-simulation was established.
The next nontrivial dynamical obligation would be a smaller **closed**
word grammar for a marked blank orbit, with exact clock and maximal-run
boundary data. Reducing the symbolic tape to its population or envelope
without a continuation proof has already failed the hostiles above.
