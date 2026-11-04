# The ternary parent and C2: small digits are not inherited (digit 1 dies at the grandchild), the end-count lemma behind it, and where a ternary parent must come from

**Session:** opus-2026-11-04-S77. Worktree `C:/tmp/wt-s77`. Canon: THM-4873. Pre-registration: HYP-9466.

**Owner prompt (verbatim).** "start another session on the ternary parent for full spectrum coverage"

## 0. The problem and what is already known

**The target.** Chan-Pak Conjecture 6.2 (TCS 1015 (2024) 114776, arXiv:2308.10214): every odd number except 7 and 21 is
the Hamiltonian path count of some tournament.

**The C2 law.** For a tournament B and a mark a, add x with x -> a only (every other vertex beats x), then y beaten by
all of B with y -> x. The child has

    H = h(B) + 2 B_a(B) = 3 h(B) + 2 sigma_a,    sigma_a = B_a - h,  B_a = sum_{S subset V - a} h(B[S]) s_a(B[V - S]).

This is `multiport_defect_grammar_20261103/construction/RESULT.md` law (4). klein's THM-4865 states it at ANY vertex:
the C2 child is the height-two member of his sole-exit tower, and the sole-exit condition is only the joint guard.

**The ternary parent.** Every odd N has a unique odd predecessor under 3h, 3h + 2, 3h + 4, chosen by N mod 3. A
*ternary parent* of N is a realiser of that predecessor carrying a mark of the needed digit 0, 1 or 2.
- If every C2 child again carried marks of digits 0, 1 and 2, induction from a finite window would give every odd
  number: this is the strict ternary scheme.
- The full-digit alternative allows any digit sigma on the line 3h + 2 sigma = N (klein, HYP-9454 (1)).

**State of the art, credited.**
- klein's root-response numerations attain every odd number in [69, 8,079,323] (THM-4865, HYP-9454 (6)).
- The full-digit sweep covers [307, 26551] with sole-exit marks of order <= 10 (THM-4862 (d)).
- The ready spectra R_1, R_2 have dyadic densities of about 0.65-0.77 (THM-4862 (f)).
- PREPARATION_REDUCTION: H-neutral preparation for digit s at h exists iff h is in R_s.
- Session191 (16:01-16:21) builds two-port ternary fans by fixed-core wrappers. The second generation fails: the
  profile (61, 63) has no core at any order.
- Session190 CP3 shows the S76 +-2 pair gives no ternary parent.

This session re-derived the any-mark one-step sweep from all classes of order <= 8 (section 3). It misses only
{1, 7, 11, 13, 21, 47, 67} up to 1,500. That is klein's height-two numeration, not new.

## 1. The end-count lemma (PROVED, elementary; not found elsewhere in the repository)

**Lemma.** Let T be a tournament with H = H(T) Hamiltonian paths.
- (i) If H >= 3, no vertex is the end of exactly H - 1 Hamiltonian paths.
- (ii) A vertex is the end of exactly H - 2 Hamiltonian paths iff T = TT_j => C3 (the last strong component is a
  3-cycle, and the earlier components are single vertices) and the vertex lies in that 3-cycle.
- If H = 1 (T transitive), every vertex ends 0 or 1 paths.

The same statements hold for starts, by taking the converse.

*Proof.*
- Let T - (last component) have strong components C_1 => ... => C_(k-1), and let C = C_k be the last one, of order m.
  Every Hamiltonian path ends in C, so e_v = 0 outside C.
- For v in C, e_v(T) = M e_v(C), with M = prod_{i<k} H(C_i) >= 1. Also H = M H(C).
- If m = 1, the single vertex of C ends all H paths. Then no vertex ends H - 1 or H - 2 (H is odd, so H - 2 >= 1 is
  impossible as an end count of any vertex).
- If m >= 3, C is strong. Each of its m vertices ends at least one of its Hamiltonian paths: open a Hamiltonian cycle
  of C before that vertex. So every vertex of C ends at least M paths of T.
- Hence a vertex a ending H - d paths leaves at least (m - 1) M <= d for the others.
  - d = 1 forces m <= 1 + 1/M < 3, impossible.
  - d = 2 forces (m, M) = (3, 1). Then C = C3, every earlier component is a single vertex, H = 3, and a is any vertex
    of the 3-cycle. Conversely, every vertex of the 3-cycle in TT_j => C3 ends exactly 1 = H - 2 paths.

QED.

## 2. Digit heredity under C2

**Native digits** (inherited: klein's tower lemmas e_x(T + x) = h - e_t and B_x(T + x) = B_t + h; re-verified here). In
the child T' = C2(B, a):
- sigma_x(T') = 0;
- sigma_y(T') = h(B) - e_a(B);
- e_x(T') = h(B), e_y(T') = B_a(B).

**Corollary (PROVED from the lemma).**
- The native mark y of C2(B, a) has digit 1 iff B is transitive and a is not its sink. It has digit 2 iff
  B = TT_j => C3 and a lies in the 3-cycle.
- In any grandchild C2(C2(B, a), u), the native mark y never has digit 1, since H(C2(B, a)) >= 3. It has digit 2
  exactly when C2(B, a) = TT => C3, i.e. B is transitive with its sink as mark, and u lies in that 3-cycle: three
  grandchildren per order.

**All marks** (FINITE-EXACT; every class B of order 1-7, every mark a, all marks of the child and of every grandchild;
`c2_digit_heredity_20261104.out`; order 7 pre-registered as HYP-9466 and confirmed):

| \|B\| | children | with a digit-1 mark | with a digit-2 mark | grandchildren | with a digit-1 mark | with a digit-2 mark |
|---|---|---|---|---|---|---|
| 2 | 2 | 1 | 1 | 8 | 0 | 3 |
| 3 | 6 | 2 | 4 | 30 | 0 | 3 |
| 4 | 16 | 3 | 4 | 96 | 0 | 3 |
| 5 | 60 | 4 | 4 | 420 | 0 | 3 |
| 6 | 336 | 5 | 4 | 2,688 | 0 | 3 |
| 7 | 3,192 | 6 | 4 | 28,728 | 0 | 3 |

- The digit-1 children are exactly the |B| - 1 transitive cases of the corollary.
- The three digit-2 grandchildren at each order are exactly the corollary's native cases. So **no old mark of any
  grandchild ever has digit 1 or 2** in this range.
- Typical grandchildren have no small digit at all: at |B| = 7 the smallest positive digit exceeds 20 in 28,463 of
  28,728 cases.

**Consequence for the ternary parent.** The strict ternary scheme dies by the second generation. Apart from the
transitive seeds, no C2 descendant carries digit 1, and digit 2 survives only in the three C3-bottom cases. A ternary
parent of N must therefore come from a fresh base, which is what klein's numerations and the full-digit sweep use. It
cannot come from the previous step's own marks. This complements:
- the preparation reduction (needing h in R_1 and R_2, densities about 0.75);
- Session191's fixed-core wrapper obstruction ((61, 63) at the second generation).

All three known self-reproducing routes stop at generation two.

## 3. Supporting data

**Explicit digit-1 family** (`ternary_parent_digit_families_20261104`; FINITE-EXACT).
- Construction: the unit-start rooted tournaments with a single-vertex core and a free transitive tail (FUNNELS F1),
  conversed and given a sole exit. They are digit-1 ready with H = h(C) + F(C, r) - 1. Every one was checked at
  t <= 9.
- At t <= 14 their values cover only about 35% of the odd numbers in the top range, with many gaps below 400.
- So no single explicit digit-1 family is a ternary parent.

**The any-mark one-step sweep** (`anymark_c2_digits_20261104`; a re-derivation of klein's height-two numeration).
- From all 6,880 + ... classes of order <= 8, the odd v <= 1,500 that are not of the form h + 2B_a are
  {1, 7, 11, 13, 21, 47, 67}.
- Only 65 of the 320 realised odd h <= 1,500 carry a digit-1 mark at orders <= 8, and 79 a digit-2 mark.

## 4. What a proof of full coverage now needs

Any C2-based proof must use digits that grow with h (fresh bases with large digits). A finite stock of small-digit
self-reproducing parents will not do, by section 2. Session191 rules out the fixed-core wrapper version.

The remaining open statement is HYP-9454 (1): the ready pairs (h, sigma) over all rooted tournaments meet every line
3h + 2 sigma = N for N >= 69. It is FINITE-EXACT to 8,079,323 by klein. A proof would need a two-dimensional
density or interval statement about the realisable pairs (h, B_a), for which no route is known.

## 5. Typing

| claim | type |
|---|---|
| End-count lemma (H - 1 impossible for H >= 3; H - 2 iff TT => C3) | PROVED (elementary) |
| Native digits of C2 children | inherited (klein THM-4865 lemmas), re-verified |
| Native-mark corollary (digit 1 only from transitive B; never at a grandchild's native mark) | PROVED |
| All-mark heredity table, orders 1-7 | FINITE-EXACT; order 7 pre-registered (HYP-9466), confirmed |
| Explicit digit-1 family sparse; any-mark sweep exceptions | FINITE-EXACT (the sweep re-derives klein) |
| "Every known self-reproducing route stops at generation two" | synthesis of PROVED/FINITE-EXACT pieces |

## 6. Pre-registration record

HYP-9466 (commit d96e81539f) pre-registered P1-P3 at |B| = 7 after the exploratory run at |B| <= 6. All hold:
3,192 native-digit checks with 0 failures; 6 and 4 small-digit children; 0 and 3 small-digit grandchildren.

## 7. Audit

(Filled after the independent audit.)

## 8. Reproduction

    gcc -O3 -march=native -o her c2_digit_heredity_20261104.c && ./her 7
    gcc -O3 -march=native -o any anymark_c2_digits_20261104.c && ./any 8 6 1500
    gcc -O3 -march=native -fopenmp -o fam ternary_parent_digit_families_20261104.c && ./fam 6 14 30000
