# Independent adversarial audit of opus-2026-11-05-S81 (THM-4884, HYP-9477, definability-lenses note)

Auditor: independent (own code, nothing in `C:/tmp/wt-s81` or the main repository modified). Audited state: worktree
`C:/tmp/wt-s81` at `97a5da3458` (includes the author's corrections `886693f53a` and the Session 200 merge
`97a5da3458` that the coordinator asked to be included). All audit files are under
`C:/Users/Eliott/AppData/Local/Temp/claude/800365894db7d24c/f973fbc2-2adf-4ee3-a3b5-b87fe066eb0a/scratchpad/s81_audit/`
(below `$SA`).

## VERDICT: PASS WITH CORRECTIONS

Every numerical and structural claim of THM-4884 (1), (2), (2b), (3) in the *current* text reproduces with my own code
(labelled-tournament enumeration, isomorphism-class census to order 8, direct re-reading of the `P_d(1)` lines, subset
DP on the witnesses, direct-tournament Hankel blocks, an independent O/Z/J profile replay). The proof of (1) is correct.
The corrections concern (i) an error that was present in the committed canon text for 7 minutes and is not logged,
(ii) prior-art credit, (iii) several loose or over-stated sentences in the reasoning sections, and (iv) an arbitrary
cutoff in HYP-9477 (A).

---

## Majors

1. **A false value was committed to canon and repaired without a ledger entry.** The pre-audit checkpoint `841dc8e7f6`
   (THM-4884 (2), note section 4, HYP-9477 (A)/(B), both INDEX rows) stated `m(41) = 3`, "m(h) >= 4 for every
   attained h in (63, 625)" and a row "m = 3 for h in {25, 41, ...}". The `P_d(1)` data give `m(41) = 6`, `m = 3` only
   at 25, and `m >= 6` above 63 (minimum 6 at 67). The author corrected all of this in `886693f53a` (15:57:37, seven
   minutes after the checkpoint) and the corrected numbers are exactly what I recompute (section A2 below). Required:
   record the slip (a mis-read of the `P_d(1)` index) in `01-canon/MISTAKES.md` or at least in THM-4884's status line,
   which currently gives no sign that statement (2) was ever different; the pre-correction INDEX rows ("at least four
   classes above 63") were also rewritten silently. Nothing else in the packet depended on the wrong values.

2. **Prior-art credit for THM-4884 (1) is incomplete.** The mechanism of (1) is already in canon twice: THM-4746 (2)
   proves `s(p) = (2/p) (mod p)` "by the rotation action" (orbits of size 1 or p on labelled tournaments on `Z_p`,
   fixed points = circulants), and THM-4747 (4) sums over *all labelled tournaments on `Z_p`* under the same rotation
   action (fixed pairs = circulant + arithmetic cycle) to recover Wilson's theorem. The inline argument "the rotation
   acts freely on their Hamiltonian paths (the first vertex moves), so each circulant count is a multiple of p" is
   LEM-003 (Aut(T) acts freely on Hamiltonian paths; `|Aut| | H`), which THM-4747 lists as its dependency. THM-4884 (1)
   is the H-graded version of the same two-line argument; it is correct and worth recording (THM-4747's own status
   says the same of its statements 1-4), but the status line must cite THM-4747 (4) and LEM-003 explicitly beside
   THM-4746, and the word "new" should not attach to the mechanism. (The note's own text, section 3, already calls it
   "the spectrum's version of the Korselt mechanism of THM-4746".)

## Minors

1. **HYP-9477 (A): the cutoff 38 is arbitrary.** Below 625 the multiplicity function omits, besides 9, 13, 21, 22, 34,
   37, also 41, 44, 50, 51, 55, 57, 60, 62, 64, 67 and 79 (my recomputation; every value 0..38 other than the six is
   taken, the maximum below 625 is 327 at h = 615). The argument given for (A) ("m(h) >= 73 on [375, 623] and the
   counts grow") applies equally to every omitted value below 73. Either conjecture the full list of omitted values
   below the floor 73, or explain why only the values up to 38 are singled out.

2. **Note section 6 misquotes THM-4748's crossing formula.** THM-4748 line 62 reads
   `H = 1 + 2^(p+r+1) + 2^(q+r+1) - [q = 0]*2^p - [p = 0]*2^q + 2^(p+q+1)*3^r`; the note drops the two boundary
   corrections `-[q=0]2^p - [p=0]2^q`.

3. **Note section 2 overstates Session 200.** "is the reachable set ... existentially definable inside its lattice by
   equations of bounded size - which Session 200 answered positively for C1/C2 words (O(log n) variables)": Session
   200's systems are a *target-dependent family* of quadratic systems of size O(log n) (its own header: "This is a
   target-dependent family of finite systems, not one fixed polynomial defining an unbounded integer predicate"); that
   is not a definition "by equations of bounded size". Say "by a target-dependent quadratic system of size O(log n)".
   Also, the height-bounded lifting is typed "KNOWN" in the section 1 table while its packet is "PROVED conditional on
   the frozen physical C1 source/adjoint classification"; add the conditionality.

4. **Section 0 / section 1 table: "every integer is aU - bD on a one-apex source" and "the signed functional is onto
   Z".** THM-4882 (3): for fixed coprime (a, b) the range `S(a, b)` is all odd integers (a even, b odd) or all even
   integers (a, b odd), never all of Z. The statement is true only with the slope chosen by parity; say so.

5. **THM-4884 (2b) / note 4b: "matching the counts of Session 200 and of its auditor" at U <= 50,000.** Session 200's
   auditor enumerated only U <= 5,000 (126,910) and *inspected* the author's 50,001 receipts ("The larger stated
   universe was reproduced by the author, not independently re-enumerated by this audit"). My own replay now supplies
   the independent 3,816,297 (section A7); cite that rather than the Session 200 auditor for the 50,000 figure.

6. **(2b) "the range (5001, 50001] is new".** Session 200's census already reports 24,971 represented targets through
   50,001, i.e. exactly 29 missed odd values; what is new in S81 is the identification of those 29 with Session 197's
   list (Session 200's packet does not list them: no "2647" anywhere in its PROOF/README/VALIDATION files). Word the
   novelty as the identification, not the range.

7. **Prior art for the small multiplicities of (2).** Session24's packet
   (`actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`, an inventory "complete for H <= 73",
   with the remark "H=19,35,39,63 have no self-converse classes at any order") and the S50 note's section 3 (the
   uniquely-attained-implies-self-converse parity argument, and the per-order lists `h_patterns.out` lines 65-69: at
   n = 5 the uniquely attained values are 1, 11, 13) are the earlier sources for everything in (2) below 73 and for
   observation (a); cite them beside THM-4744 (3).

8. **THM-4880 is not cross-referenced in (3).** The ribbon lane (THM-4880 (A), (E)) has a rank-14 representation for
   span-2 words and proves "47 is not a span-2 count; every realiser of 47 has width 3". S81's crossing chains are a
   different sub-language, but rank 14 and the loss of 47 recur; the theorem should relate the two (and say they are
   different languages: span-2 words contain shared-endpoint pairs of span <= 2, crossing chains do not).

9. **Two different lower bounds for the minimal complete rank.** Note sections 0/5: "between 5 and 17"; HYP-9477 (D):
   "at least 10 and at most 17". Both are reasoning, but state once that 5 is "no complete candidate of rank <= 4 is
   known" and 10 is the conjectured bound from the failure of the rank-10 nested chains.

10. **Q3's "program output" is not archived.** The note says "The full coefficient tables for n <= 7 are in the
    program's output" but no output of `spectrum_polynomial_congruence_20261105.cpp` exists in the worktree or in
    `Temp/s81`. My tables are in `$SA/spectrum_poly_audit.out`; commit the author's.

11. **Pre-registration Q1 is a prediction about data already in the repository** (the `P_d(1)` lines are from
    2026-10-21). Git can certify only that the prereg commit `ab4dd938fb` (15:43:57) precedes the census outputs
    (15:44:05-06, runs of 4.7 s / 0.3 s / 0.7 s - section A5), not that the author had not read the `P_d(1)` lines.
    The two refuted Q1 clauses ("contains 41", "every multiplicity 1..10 occurs") are the only evidence of blindness;
    the note should say this rather than present Q1 as a pre-registration of the same strength as Q2.

12. **Minor wording in (2b)**: "profiles are neither words nor isomorphism classes" is Session 200's caveat (its
    README: "Profile multiplicity is not word multiplicity or graph-isomorphism multiplicity"); attribute the quote.
    The quoted sentence in section 3 ("a residue fiber can have many physical labelled objects while an exact scalar
    fiber is empty") is verbatim from `session200/sharp-dyadic-fibers/PROOF.md` line 105 - fine.

## Upgrades (optional, supported by this audit's computations)

1. (3)'s ranks can be typed as verified from both sides: my direct-tournament Hankel blocks (no automaton, prefixes and
   suffixes of length <= 5, H by subset DP) give rank exactly 14, 14, 10 modulo two primes (section A6), and the
   automaton pairing gives the same numbers as upper bounds.
2. (2) is now independently reproduced for every h < 75 by a class census of orders <= 8 (section A3); the status line
   may say that the inherited single-enumeration caveat of THM-4744 (the 625 file was never re-enumerated by the S50
   audit) applies only to 375 <= h < 625.
3. (2b)'s profile counts and the pf table are independently reproduced (section A7); the 29 missed values through
   50,001 coincide with Session 197's list exactly (plus n = 1, which is not odd >= 3).
4. Observation: at n = 7 the non-circulant regular tournament (H = 171, |Aut| = 3) has labelled count 1680 = 7!/3,
   divisible by 7 as the theorem requires; the circulant classes have 720 = 7!/7 (H = 175, six circulants on Z_7) and
   240 = 7!/21 (H = 189: the two circulants S = {1,2,4}, {3,5,6} are the Paley tournament and its converse, which is
   isomorphic to it since -1 is a non-residue mod 7). This explains the residues 6 and 2 and could be stated.

---

## A. Reproductions (commands and outputs)

Compilers: `g++ -O2`; `python3` 3.x with numpy 2.3.2 (the repository's rank script needs numpy; `python3 -I` hides the
user site-packages, so it was run without `-I` on dumps I generated myself).

### A1. THM-4884 (1): spectrum polynomial, invariant tournaments, circulants (own code)

`$SA/spectrum_poly_audit.cpp` enumerates all `2^C(n,2)` labelled tournaments, counts Hamiltonian paths by its own
subset DP (bitmask adjacency), tests invariance under `i -> i + n/p` for every prime `p | n`, and lists the circulants.

```
$ for n in 2 3 4 5 6 7; do ./spectrum_poly_audit.exe $n; done
n = 3: 8 labelled tournaments. P_3 coefficients h:count 1:6 3:2
p = 3 (rotation by 1): 2 invariant tournaments; F coefficients h:count 3:2
  residues mod 3 where P != 0 mod 3: h=3:P=2(2 mod 3),F=2
  congruence P_3 == F_{3,3} mod 3: HOLDS
n = 4: 64 labelled tournaments. P_4 coefficients h:count 1:24 3:16 5:24
p = 2 (rotation by 2): 0 invariant tournaments ... congruence P_4 == F_{4,2} mod 2: HOLDS
n = 5: 1024 labelled tournaments. P_5 coefficients h:count 1:120 3:120 5:240 9:240 11:120 13:120 15:64
p = 5 (rotation by 1): 4 invariant tournaments; F coefficients h:count 15:4
  residues mod 5 where P != 0 mod 5: h=15:P=64(4 mod 5),F=4
  congruence P_5 == F_{5,5} mod 5: HOLDS
circulants on Z_5: S={34}:H=15 S={24}:H=15 S={13}:H=15 S={12}:H=15
n = 6: 32768 labelled tournaments. P_6 coefficients h:count 1:720 3:960 5:2160 9:2960 11:1440 13:1440 15:2208 17:1440
  19:1440 23:2880 25:1440 27:480 29:2880 31:1440 33:2640 37:3600 41:720 43:1440 45:480
p = 2 (rotation by 3): 0 invariant tournaments ... congruence P_6 == F_{6,2} mod 2: HOLDS
p = 3 (rotation by 2): 32 invariant tournaments; F coefficients h:count 9:8 45:24
  residues mod 3 where P != 0 mod 3: h=9:P=2960(2 mod 3),F=8
  congruence P_6 == F_{6,3} mod 3: HOLDS
n = 7: 2097152 labelled tournaments. P_7 coefficients ... 171:1680 175:720 189:240
p = 7 (rotation by 1): 8 invariant tournaments; F coefficients h:count 175:6 189:2
  residues mod 7 where P != 0 mod 7: h=175:P=720(6 mod 7),F=6 h=189:P=240(2 mod 7),F=2
  congruence P_7 == F_{7,7} mod 7: HOLDS
circulants on Z_7: S={456}:H=175 S={356}:H=189 S={246}:H=175 S={236}:H=175 S={145}:H=175 S={135}:H=175 S={124}:H=189 S={123}:H=175
```
(full tables in `$SA/spectrum_poly_audit.out`). Every statement of (1) checks: the congruence for all (n, p) with
p | n <= 7; n = 5: only h = 15 escapes (64 = 4 mod 5, four circulants); n = 7: only 175 (720 = 6 mod 7) and 189
(240 = 2 mod 7); the two circulants with 189 are S = {1,2,4} (Paley) and S = {3,5,6} (its converse); the other six
have 175; n = 6, p = 3: thirty-two invariant tournaments, residue 2 at h = 9 only (F has 9:8, 45:24). The even cases
(p = 2) are trivial: the rotation by n/2 reverses the arcs {i, i+n/2}, so no tournament is fixed and P_n is even
coefficientwise (which also follows from the converse involution on labelled tournaments).

**Proof check.** Orbits of the order-p rotation subgroup have size 1 or p and H is a relabelling invariant, so
P_n = F_{n,p} (mod p) coefficientwise; for n = p the fixed tournaments are exactly the circulants (arc i -> j iff j - i
in a half-system S); the rotation is an automorphism of a circulant and a nontrivial automorphism fixes no Hamiltonian
path (the first vertex moves) - this is LEM-003 - so p | H(circulant); hence for p not dividing h, N(p, h) = (number of
circulants with count h) = 0 (mod p). Correct. The wording "the rotation acts freely on the Hamiltonian paths of a
circulant" is precise (free action of the cyclic group Z_p, orbits of size p).

### A2. THM-4884 (2): the `P_d(1)` lines re-read (own parser)

```
$ python3 -I multiplicity_audit.py
375 file: 186 values, d = 3..373; 625 file: 125 values, d = 375..623
125 explicit polynomials P_d in the 375 file cross-checked against P_d(1): mismatches 0
m = 0: [7, 21]
m = 1: [3, 5, 11, 13]
m = 2: [9, 17, 19, 35, 39, 47, 63]
m = 3: [25]
m = 4: [15, 23, 31, 43, 55]
m = 5: [27, 61]
min m(h) for attained 63 < h < 625: 6 at [67]
all attained h in (63,625) have m >= 6: True ; >= 4: True
largest h with m(h) <= 3: 63 ; largest with m(h) <= 5: 63
min m(h) on [375, 623]: 73 at [395] ; >= 73: True
multiplicity values NOT taken below 625, up to 40: [9, 13, 21, 22, 34, 37]
every value 0..38 other than 9,13,21,22,34,37 taken: True
max multiplicity below 625: 327 at [615]
values 39..80 not taken: [41, 44, 50, 51, 55, 57, 60, 62, 64, 67, 79]
table h<=101: 3:1 5:1 7:0 9:2 11:1 13:1 15:4 17:2 19:2 21:0 23:4 25:3 27:5 29:6 31:4 33:8 35:2 37:7 39:2 41:6 43:4 45:12
  47:2 49:8 51:12 53:10 55:4 57:7 59:12 61:5 63:2 65:14 67:6 69:18 71:11 73:19 75:12 77:16 79:11 81:26 83:11 85:10 87:11
  89:23 91:29 93:15 95:17 97:23 99:24 101:18
m(41) = 6 ; m(15) = 4 ; m(25) = 3 ; m(67) = 6 ; m(61) = 5 ; m(63) = 2
Q1: U within [3,100]: True ; |U| <= 12: True ; contains 3,5,11,13: True ; contains 41: False ; every element > 15 prime or 3*prime: vacuous
Q1: every multiplicity 1..10 occurs: False ; missing: [9]
```
Every number in the corrected statement (2) and in the note's section 4 table (including the h <= 101 table) agrees.
Definition check: `strong_ears_625.py` builds `inv = {m: Counter() for m in range(3, 13)}`, loads orders 3..9 from the
class census and orders 10-12 from the ears, and prints `P_d(1) = sum over m of inv[m][d]`: P_d(1) is the number of
strong classes **of all orders** with H = d, as THM-4884 uses it; the 375 file's explicit polynomials (e.g.
`P_125 = 1z^7 + 10z^8 + 2z^9 + 1z^10`) confirm that orders below 9 are included. Completeness below 625 is THM-4744 (3)
(f(13) = 625, audited by S50 for d < 375; the 625 file was reproduced by the author only - inherited caveat).

### A3. Independent class census to order 8 (own code; exact for h < f(9) = 75)

`$SA/class_census_audit.cpp`: vertex extension + canonical form (colour refinement, permutations within colour
classes, lexicographic minimum), subset-DP H, strong test, self-converse test, |Aut|.

```
$ ./class_census_audit.exe 8
order 1..8: classes 1, 1, 2, 4, 12, 56, 456, 6880 (all as expected); strong 1, 0, 1, 1, 6, 35, 353, 6008 (as expected)
order 3: H=3 self-converse=1 |Aut|=3          order 4: H=5 self-converse=1 |Aut|=1
order 5 strong classes: H=15 (sc, |Aut|=5), H=11 (sc), H=13 (sc), H=9 (sc), H=9 (sc), H=15 (sc, |Aut|=3)
m(h) over strong classes of orders <= 8: 3:1 5:1 9:2 11:1 13:1 15:4 17:2 19:2 23:4 25:3 27:5 29:6 31:4 33:8 35:2 37:7
  39:2 41:6 43:4 45:12 47:2 49:8 51:12 53:10 55:4 57:7 59:12 61:5 63:2 65:14 67:6 69:18 71:11 73:19 (then 75:9 - orders
  <= 8 only; the P_d(1) value 12 adds E(9) = 3)
odd h < 75 with m(h) = 0: 7 21
self-converse counts: h=3:1 5:1 11:1 13:1 (the four unique classes are self-converse); 9: 2 of 2; 17: 2 of 2;
  19, 35, 39, 47, 63: 0 of 2 (converse pairs)
```
Every m(h) for h < 75 agrees with the `P_d(1)` lines (34 values), the four unique classes are self-converse, and the
order-5 realisations of 11 and 13 are each a single self-converse class (in fact all six strong 5-classes are
self-converse). Also checked by brute force over the 1024 labelled 5-tournaments: the H = 13 class has exactly 2
Hamiltonian cycles (so it is not a C1 readout, as (2b) says), the H = 11 and both H = 9 classes have 1.

### A4. THM-4884 (3): witnesses (own subset DP)

```
$ python3 -I witness_dp_audit.py
claimed 43: m=7 arcs=[(2,0),(4,2),(6,2),(6,4)] -> H=43 OK; shared heads [2], shared tails [6]; cutwidth 2
claimed 47: m=7 arcs=[(6,0),(6,1)] -> H=47 OK; shared heads [], shared tails [6]; cutwidth 2
claimed 63: m=8 arcs=[(3,0),(6,3),(7,3)] -> H=63 OK; shared heads [3], shared tails []; cutwidth 2
claimed 79: m=8 arcs=[(7,0),(7,2)] -> H=79 OK; shared heads [], shared tails [7]; cutwidth 2
claimed 103: m=9 arcs=[(3,0),(6,3),(8,3)] -> H=103 OK; shared heads [3], shared tails []; cutwidth 2
claimed 113: m=9 arcs=[(2,0),(6,2),(6,3),(8,6)] -> H=113 OK; shared heads [], shared tails [6]; cutwidth 2
claimed 235: m=10 arcs=[(6,0),(6,4),(8,6),(9,6)] -> H=235 OK; shared heads [6], shared tails [6]; cutwidth 2
S79 file (Temp/s79/minfa_decoded.txt): all seven entries match the note's transcription (m and arcs)
```
The S79 minimal-arc file comes from a `dedupe = 0` run to X = 511 (`census_k0_w2_X511_minfa5.txt` header), so the arc
minima are exact, as the note says.

### A5. THM-4884 (3): censuses, dumps, ranks (author's engine recompiled; rank script rerun)

```
$ g++ -O2 -o spine.exe /c/tmp/wt-s81/04-computation/experiments/spine_grammar_sublanguages_20261105.cpp
$ SPINE_MODE=1 ./spine.exe 0 2 2047 64 0 301   # real 4.68 s
G(k=0, w=2), X = 2047 ...: 64 automaton states, 2726467 DFS nodes ...; odd values attained directly: 1013 of 1024
with ordered-join products: 1014 of 1024; 10 misses; first misses: 7 21 43 47 63 79 103 113 235 709
$ SPINE_MODE=2 ./spine.exe 0 2 2047 64 0 301   # real 0.26 s
... 64 automaton states ...; odd values attained directly: 994 of 1024
with ordered-join products: 998 of 1024; 26 misses; first misses: 7 21 37 43 47 63 67 79 103 113 127 139 163 173 235
  431 491 519 559 563 709 751 797 1373 1427 1499
$ SPINE_MODE=3 ./spine.exe 0 2 2047 64 0 301   # real 0.71 s
... 58 automaton states ...; 440 of 1024; 584 misses; first misses: 7 11 13 19 21 23 35 39 41 ...
diff against Temp/s81/census_mode{1,2,3}_X2047.txt: identical (all three)
$ for m in 0 1 2 3; do SPINE_MODE=$m ./spine.exe 0 2 10 -1 > dump_mode${m}_audit.txt; done
STATES 83 / 64 / 64 / 58; diff against Temp/s81/dump_mode{0,1,2,3}.txt: identical (all four)
$ python3 spine_grammar_ranks_mod_20261104.py dump_mode{m}_audit.txt --empty1
mode 0: Hankel rank = 17 (two primes agree); seams c=0: 1, c=1: 3, c=2: 13      (calibration = THM-4878)
mode 1: Hankel rank = 14; seam c=2: 10 (50 cut states, reachable dim 25)
mode 2: Hankel rank = 14; seam c=2: 10
mode 3: Hankel rank = 10; seam c=2: 6 (44 cut states, reachable dim 20)
```
All numbers of the note's section 5 table and of THM-4884 (3) reproduce (states, ranks, seams, attained counts, miss
lists: 24 misses besides 7, 21 in mode 2; 582 in mode 3; first mode-3 misses 11, 13, 19, 23, 35).

**Letter filter (read, `letterAllowed`).** Mode 1: `b > 1 || popcount(T) > 1` rejected - "no double events" exactly as
described (a vertex opening two slots is a vertex beaten by two later vertices = shared head; a vertex terminating
two slots beats two earlier vertices = shared tail; T = 1, b = 1 at one vertex is allowed and is not a double event).
Mode 2 additionally rejects `c == 2 && T == 2`; mode 3 rejects `c == 2 && T == 1`. Slot 0 is the older slot: the
engine renumbers "surviving old slots, then the new ones" in creation order after every letter, so with two pending
slots bit 0 is the earlier-opened one. Hence mode 2 = only the older slot may close = first-in-first-out = heads
h1 < h2 and tails t1 < t2 (crossing); mode 3 = last-in-first-out = nested. The description is correct. The pruning
devices stay sound under the filter: the closer bound takes the minimum over *allowed* closing letters with b = 0, and
in every mode the first future tail of any completion in the sub-language is such a letter (mode 2: the older slot,
mode 3: the newer), so the induced prefix-plus-one-vertex tournament is among the candidates; sink runs, dedupe and the
ordered-join closure are mode-independent (a disconnected sub-language word splits at an empty cut into sub-language
words).

### A6. Independent Hankel-rank lower bounds without the automaton (own code)

`$SA/hankel_block_audit.cpp` decodes letter words directly into chain tournaments (letters indexed per c exactly as
`lettersFor(c)` does, mode filters applied, double arcs rejected), computes H(T_uv) by subset DP, and takes the rank of
the block [H(T_uv)] over prefixes |u| <= Lp and suffix letter strings |v| <= Ls modulo 2^31-1 and 2^31-19 (empty word
= 1; the series value of a word with pending slots is H of its induced tournament, which is what the automaton's
"complete" functional computes).

```
$ ./hankel_block_audit.exe m 4 4 :  mode 1: rank 13; mode 2: 13; mode 3: 9     (110/89/89 prefixes)
$ ./hankel_block_audit.exe m 5 5 :  mode 1: 428 prefixes x 3906 suffixes, 638160 nonzero entries, rank 14 = 14
                                    mode 2: 286 x 1365, rank 14 = 14
                                    mode 3: 286 x 1365, rank 10 = 10      (38 s in all)
```
So rank >= 14, 14, 10 from the tournaments themselves; the automaton pairing (A5) gives <= 14, 14, 10. The ranks are
exact.

### A7. (2b): Session 200's profile census replayed (author's script, then own code)

```
$ python3 profile_vs_class_multiplicity_20261105.py 624 strong_ears_375.out strong_ears_625.out   # 0.20 s
profiles with U <= 624: 5757; pf = 0: [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421,
  449, 503, 533, 569, 583]; pf = 1: [3, 11]; m = 1 & pf >= 2: [5]; pf = 1 & m >= 2: []; pf > m: [5, 9, 17, 19, 27, 35,
  39, 53, 61, 63]   (output identical to the committed .out apart from CRLF)
$ python3 profile_vs_class_multiplicity_20261105.py 5000     # 0.40 s: 126910 profiles; 29 zeros [7 ... 2647]; pf=1: [3, 11]
$ python3 profile_vs_class_multiplicity_20261105.py 50000    # 8.4 s: 3816297 profiles; the same 29 zeros; pf=1: [3, 11]
$ ./ozj_profiles_audit.exe 624 / 5000 / 50000     (own BFS over O/Z/J from (1,0,2,1), readout F = D - E + U)
B = 624: 5757 distinct profiles; B = 5000: 126910; B = 50000: 3816297
odd n <= 50001 with pf(n) = 0: 1 7 13 21 41 43 55 67 77 83 163 173 231 233 247 251 277 299 421 449 503 533 569 583 637
  821 1067 1339 1861 2647   (n = 1 aside: exactly the 29 values of session197/lowcycle/C1_COMPLETE.md lines 84-86)
pf(n) = 1: 3 11 ; pf table to 101 identical to the script's
```
Checks: tO = (1,0,4,1) and tJ = (3,1,3,1) both have F = 5 (the collision at the strong 4-tournament, Session 200's
audit REPORT.md: "the two profiles at target5 close to isomorphic order-four tournaments"); 13 has pf = 0 and its
unique strong class has two Hamiltonian cycles (A3). F >= U + 1 (D - E >= d >= 1), so U <= 50,000 covers every F <=
50,001. The pf < m list (23, 25, 29, 33, 37, 41, 43, 45, 49, 51, 55, 57, 59) and the table h : pf : m to 63 in note
4b agree with the script's table. Session 200's PROOF.md section 6 states 3,816,297 profiles and singleton fibers 3 and
11; its auditor enumerated U <= 5,000 (126,910) independently and inspected the 50,001 receipts (24,971 represented
targets; 25,000 - 24,971 = 29).

Section 3's new paragraph on `sharp-dyadic-fibers` is accurate: Theorem 1 there says that at block length k = 2^r - 2
every odd residue mod 2^r occurs an odd number of times among the 2^k O/Z words and at every larger length a positive
even number of times, the critical coefficient being the parity of the number of perfect matchings of the square of a
path on 2t vertices (= F_(t+1), Fibonacci); the quoted sentence is verbatim (PROOF.md line 105).

### A8. Pre-registration and commit order

```
56f673ca07 15:39:36 reserve THM-4884 / HYP-9477
ab4dd938fb 15:43:57 pre-register Q1-Q3 + the two programs (prereg md, congruence cpp, sublanguages cpp)
Temp/s81 outputs: dump_mode0 15:44:05.12, dump_mode1 .26, census_mode2 .38, dump_mode2 .42, dump_mode3 .53,
                  census_mode3 15:44:06.06, census_mode1 15:44:09.70
841dc8e7f6 15:50:39 checkpoint (note, THM, HYP, index rows, ledger; Q1/Q2 outcomes)
dea60f8efb 15:53:19 witnesses listed;  6e3f07dc41 15:54:36 index repair
886693f53a 15:57:37 multiplicity table corrected (m = 3 only at 25; m(41) = 6; minimum 6 at 67)
97a5da3458 16:10:25 Session 200 carmichael-fibers merged (4b, (2b), script + .out)
```
The runs take 4.7 / 0.3 / 0.7 s, so the outputs written at 15:44:05-06 were produced after the 15:43:57 commit:
Q2's predictions are certified as prior to the runs. Recorded outcomes match the data: Q2 "mode 1 complete" REFUTED,
"mode 2 complete" REFUTED, "mode 3 misses >= 10" CONFIRMED (582), "rank(2) < rank(1)" REFUTED (14 = 14),
"rank(1) <= 17" holds, "rank(3) < rank(1)" CONFIRMED (10 < 14). Q1: containment, size and the 3, 5, 11, 13 clause
CONFIRMED; "41" REFUTED (m(41) = 6, the corrected value); "prime or 3 x prime above 15" vacuous; "every multiplicity
1..10 occurs" REFUTED (9). Q3 is a check, not a prediction (correctly labelled). See Minor 11 for Q1's weakness.

### A9. Status typing, reasoning sections, Session 200 description

- (1) PROVED: correct (A1), with the credit issue of Major 2. (2), (2b), (3) FINITE-EXACT: appropriate; (2) rests on
  THM-4744 (3)'s completeness (inherited status), (3) on the S79 engine (audited) plus the modular ranks (now exact by
  A6). KNOWN/CITED items: THM-4746 (2), THM-4744 (3) converse pairing, THM-4878 ranks/censuses (17, 83 states, seam 13,
  1022/1024, G(1,0) 15 states, rank 4, 763 of 1024 - all present in THM-4878 / S79 dumps), THM-4882 (4) barrier,
  THM-4716/4719 ladders, THM-4862/4865 (C2 digit steps, tribonacci), Session 193 modular_sync, Session 198 modular
  seeds, Session 194 module barrier: all exist as cited.
- Sections 1, 2, 6, 7 make no claim about the paper's correctness ("taken as stated, not verified"; the Mazur and HTP
  consequences are reported as the paper's claims). Session 200 is described accurately except for Minors 3 and 5-6.
  The Session 200 packets and their times (compiler 15:38:06-15:39:12, lifting 15:41:49, sharp-dyadic 15:45:02,
  carmichael-fibers 15:52:12, signed parity-class coverage + fiber transports 15:56:25, root-weighted height 15:58:34)
  match the note's scope paragraph ("seven packets pushed between 15:38 and 15:58").
- Index rows and the PROBLEM-LEDGER line type the packet "PROVED / FINITE-EXACT / reasoning; audit PENDING" -
  consistent with the files; after this audit set AUDIT to PASS WITH CORRECTIONS and apply the items above.
- No MISTAKE entry exists for the 2026-11-05 S81 session (the only `S81` hits in MISTAKES.md are the July session of the
  same label).

### A10. Prior-art search (untruncated `git grep` over the worktree, excluding .out/.json/HISTORICAL)

Patterns: orbit counting / Burnside / "fixed by the rotation" / "spectrum polynomial" / `x^H` / `N(p, h)` / "divisible
by p"; "unique class" / "exactly one strong" / "uniquely attained" / "m(h) = 1" / Carmichael+spectrum / Ford / "omitted
multiplicit"; "sub-language" / "crossing chain" / "nested chain" / FIFO / LIFO / "rank 14" / "rank 10" / "double event".
Findings: THM-4746 (2) and THM-4747 (4) (Major 2); Session24 packet and S50 note / `h_patterns.out` (Minor 7); THM-4880
(Minor 8); no earlier statement of the coefficientwise congruence P_n = F_{n,p} or of the all-order unique-class set, and
no earlier sub-language census of G(0,2) (the ribbon lane's span-2 language is different). The S79 audit's brute-force
Hankel blocks (G(0,2): 477 x 320, rank 17) are the model for A6.

## Files produced

`$SA/spectrum_poly_audit.cpp/.exe/.out`, `class_census_audit.cpp/.exe/.out`, `multiplicity_audit.py/.out`,
`witness_dp_audit.py/.out`, `hankel_block_audit.cpp/.exe`, `hankel_block_audit_44.out`, `hankel_block_audit_55.out`,
`spine.exe` (the author's engine), `census_mode{1,2,3}_X2047_audit.txt` + `_time.txt`, `dump_mode{0,1,2,3}_audit.txt`,
`ranks_audit.out`, `ozj_profiles_audit.cpp/.exe/.out`, `profile_script_{624,5000,50000}.out`.
