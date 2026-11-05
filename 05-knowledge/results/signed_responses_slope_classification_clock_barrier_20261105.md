# Signed one-apex responses: Session 198's clock in the cut-state language, the parity-and-slope classification, and the positive-count barrier

**Session:** opus-2026-11-05-S80. Worktree `C:/tmp/wt-s80` (check_docs through a `subst` drive). Canon: THM-4882.
Hypothesis: HYP-9475. Pre-registration: `04-computation/experiments/signed_slope_prereg_20261105.md` (P1-P5 at commit
0e5ba011be before the probe, P6 at commit c6d7e35dc4 before the recipe run). Independent audit: section 8.

**Owner prompt (verbatim).** "merge in the below as inspiration and try to extend ideas beyond it to proofs:" followed
by Session 198's packet `session198/rational-response-coverage/PROOF.md` ("Exact odd coverage for rational signed
source responses").

**Target.** Chan-Pak Conjecture 6.2 (every odd number except 7 and 21 is a Hamiltonian-path count) remains OPEN. This
note merges Session 198's signed-response theorem into the graded spine grammars of THM-4878 (S79), proves what the
signed image is for every rational slope (not only Session 198's parity class), and proves why the same mechanism
cannot program a positive count.

## 0. Summary

- **Merge (PROVED, elementary).** Session 198's source calculus O(D,E,U,d) = (D,0,2U,d), Z(D,E,U,d) = (D-E+U,d,2U-d,d)
  is the one-apex grammar G(1,0) of THM-4878 read in profile coordinates: a chain word b (apex x beats v_i iff b_i = 1)
  has (D, E, U, d) = (H, 1 - b_m, 1 + bin(b), 1), appending a bit 1 is O and a bit 0 is Z, and Session 198's g(n) is
  the one-apex count of the word of binary value n - 1 (checked exactly for n <= 4096 against the (H, A, B)
  automaton; the recurrences coincide). The clock W of Session 198 is a left eigenvector of M_W for the eigenvalue 1
  with a Jordan drift along the constant coordinate d; the sharpness interval (1/2, 1) is the one-apex cone.
- **Classification (THM-4882 (3); (i) is Session 198, (ii) and (iii) are new, PROVED).** For coprime positive a, b
  with 1/2 < a/b < 1 let S(a, b) = { a n - b g(n) : n >= 3 odd } (the signed response aU - bD over strong one-apex
  sources). Then S(a, b) is: (i) all odd integers if a is even and b odd; (ii) **all even integers if a and b are both
  odd** (Session 198's clock with a cyclic correction: the period word then starts with Z, the drop per period is
  K_cyc = K_lin - b, and the first period after the O drops by K_lin); (iii) **never all integers of its parity if b is
  even**: for a/b <= 3/4 every positive element is congruent to -(a + b) modulo 4, for a/b > 3/4 every element below
  -a - ab/(4a - 3b) is congruent to 2 - a (4 | b) or -a (b = 2 mod 4) modulo 4. Outside the interval one sign is
  missing (Session 198). So Session 198's "a even, b odd" is exactly the odd-image case, the both-odd case gives the
  even integers, and an even denominator is obstructed by a modulo-4 class on one sign side.
- **Verification (FINITE-EXACT).** The probe (n < 2^24, targets in [-999, 999]) and the exact recipe (40 random
  targets in [-10000, 10000] for each of 11 slopes of classes (i) and (ii), all hit exactly, final points strong)
  confirm the theorem's constructions; P3 CONFIRMED, P4 REFUTED (above 3/4 the negative side is confined, which the
  proof of (iii) then explained), P6 CONFIRMED.
- **Barrier (THM-4882 (4), PROVED, elementary).** In an apex grammar G(k, 0), H(p W^(h+1) q) >= H(p W^h q) for all
  words p, W, q (deleting one period leaves an induced subtournament), so no fixed block can lower a Hamiltonian-path
  count: the exact repetition programming of a signed response has no analogue for counts, and coverage of positive
  targets must come from digit freedom (S79's counting lemma) rather than from clocks. This is the structural reason
  behind Session 198's own caveat that its theorem "does not prescribe D or U independently, nor a positive count".

## 1. Prior art read (and what is reproduced)

- Session 198 `rational-response-coverage/PROOF.md` (13:42 today): the theorem for a even, b odd; the clock word
  from the bits of A = a(B-1)/b; K = aR - bC > 0 even; the modular seed (STRONG_SOURCE_INTERFACE: for every even
  M and odd residues D_0, U_0 there are strong unit-interface sources with those residues and D arbitrarily large,
  built from the one-apex gadgets of JOINT_MODULAR_GADGET / Session 197 `modular-gadget-universality`); sharpness of
  (1/2, 1); the 3/4 control (positive values 1 mod 4); the C2 parent/child offsets. Reproduced here: the drops and
  the recipe for (2,3), (4,5), (6,7), (4,7), (8,11).
- Session 197 `retained-unit-clock-reduction` (the Z operation on (h, e, R, d)), `palindrome-mersenne-clocks` and
  `modular-gadget-universality` (the recurrence f(2x+1) = f(x), f(2x) = f(x) + x + (x mod 2), f(n) = g(n+1)), the
  positive/unit/ternary clock packets (H_m = A 3^m + beta: transitive cloning, a different clock), Session 194's
  periodic counters (H(u d^k) = A + Jk + Kk^2 + (D + Ek) 2^(Lk) + G 3^(Lk) for two-apex blocks), HYP-9468's rooted
  triple (g(v), e(v), 1 + v), Session 16's four-coordinate transfer, THM-4876 (6) (the one-apex closed formula),
  THM-4878 (S79: the grammars, Hankel rank 4 of G(1,0)).
- No packet on the branch treats the both-odd or even-denominator slopes (git log and a grep of all Session 198
  markdown for "both odd", "even integers", "even denominator", at 14:30 today).

## 2. The merge: O/Z is the one-apex grammar

Let T(b) be the one-apex tournament of the chain word b = b_1 ... b_m (chain v_1 < ... < v_m transitive, apex x with
x -> v_i iff b_i = 1), rooted at x. Its profile in Session 198's conventions (D = H, E = paths ending at the root,
U = ordered two-path covers (P_1, x P_2) with P_1 possibly empty, d = H(T - x)) is

    D = H(T(b)),   E = 1 - b_m,   U = 1 + bin(b),   d = 1,

where bin(b) = sum_i b_i 2^(m-i). Proof: the chain is transitive, so T - x has one Hamiltonian path and a path ends
at x iff the chain's unique Hamiltonian path v_1 ... v_m is followed by x, i.e. b_m = 0; a cover (P_1, x P_2) is x
alone or x -> v_i followed by any subset of the later chain vertices in order, giving 1 + sum_{b_i = 1} 2^(m-i).
Appending a chain vertex with bit 1 (beaten by x and by every chain vertex: a sink) multiplies U by 2 and sets E = 0
without changing H: this is O. Appending a bit 0 (beating x) gives U' = 2U - 1 and E' = 1 and, by the one-apex
automaton of THM-4878 (letter 0 adds 2(A + B) = bin(b) + b_m = (U - 1) + (1 - E) = U - E to H), D' = D - E + U: this
is Z. Hence Session 198's orbit of the seed t = (1, 0, 2, 1) (the arc x -> v_1) is exactly G(1,0) with
(D, E, U, d) = (g(n), n mod 2, n, 1), n = U, and g(n) is the one-apex count of the word of value n - 1 (HYP-9468's
g(v), v = n - 1; THM-4876 (6)). The probe `signed_slope_probe_20261105.cpp` checks g(n) = H(word of value n - 1)
for n <= 4096 by the automaton: 0 mismatches.

In the cut-state language of THM-4878, a word W acts on the row vector (D, E, U, d) by a nonnegative matrix M_W, and
Session 198's transfer (3) says that on the slice E = d (strong inputs) M_W has the eigenvalue B = 2^k on U and the
eigenvalue 1 on D and d with a Jordan drift: the functional phi = aU - bD satisfies phi M_W = phi - K psi with psi the
d-coordinate, psi M_W = psi. The general lemma is trivial: phi(v M_W^h) = phi(v) - h K psi(v). Exact coverage then
needs only that the reachable set meets every residue class modulo K in the two coordinates, which is the
STRONG_SOURCE_INTERFACE theorem.

## 3. Lemmas on g

**Lemma 1.** For all n >= 1, g(n) is odd and g(n) <= n; for odd n >= 3, g(n) >= (n + 3)/2.
*Proof.* Induction on n. g(1) = 1. g(2n) = g(n) <= n < 2n. g(2n - 1) = g(n) + n - (n mod 2) <= 2n - (n mod 2), and
for even n this is at most 2n - 1 because g(n) = g(n/2) <= n/2 <= n - 1. Parity: g(2n - 1) - g(n) = n - (n mod 2) is
even. Lower bound: for odd m = 2n - 1 with n >= 2, g(m) = g(n) + n - (n mod 2) >= n + 1 when n is even (g(n) >= 1)
and >= (n - 1) + 3 = n + 2 when n is odd and >= 3 (g(n) >= 3 for odd n >= 3, since g(3) = 3 and g(2k - 1) >= k + 1 >= 3
for k >= 2 by the same induction), so g(m) >= n + 1 = (m + 3)/2. QED

(Lemma 1 is the one-apex cone of Session 198, D <= U <= 2(D - E), for the strong points, together with E <= d.)

**Lemma 2 (branches).** For odd n = 2m - 1 with m >= 2,

    F(n) := a n - b g(n) = (2a - b) m - a - b g(m) + b (m mod 2),

and for even m = 2m', F(n) = 2(2a - b) m' - a - b g(m').
*Proof.* g(2m - 1) = g(m) + m - (m mod 2), and g(2m') = g(m'). QED

## 4. The classification theorem

**Theorem (THM-4882 (3)).** Let a, b be coprime positive integers and S(a, b) = { F(n) : n >= 3 odd }.
1. If a/b <= 1/2 then S(a, b) contains no positive integer; if a/b >= 1 it contains no negative integer (Session 198
   section 4; from Lemma 1: F <= (a - b/2) n - 3b/2 < 0, respectively F >= (a - b) n >= 0).
2. (Session 198) If a is even, b is odd and 1/2 < a/b < 1, then S(a, b) is the set of all odd integers.
3. (New) If a and b are both odd and 1/2 < a/b < 1, then S(a, b) is the set of all even integers.
4. (New) If a is odd, b is even and 1/2 < a/b < 1, then S(a, b) consists of odd integers but misses infinitely many
   of them: when 4a <= 3b every positive element is congruent to -(a + b) modulo 4; when 4a > 3b every element
   below -N_0, N_0 = a + ab/(4a - 3b), is congruent to 2 - a modulo 4 if 4 | b and to -a modulo 4 if b = 2 mod 4.

*Proof of 4.* Parity: n odd, g odd, a odd, b even give F odd. Write n = 2m - 1 and use Lemma 2.
Case 4a <= 3b. For odd m >= 3, Lemma 1 gives F <= (2a - b) m - a - b (m + 3)/2 + b = (2a - 3b/2) m - a - b/2 < 0.
So positive elements come from even m = 2m', where F = 2(2a - b) m' - a - b g(m'). Here 2(2a - b) = 4a - 2b is a
multiple of 4 (b even), and b g(m') is congruent to b modulo 4 (g odd, b even: b(g - 1) is a multiple of 4). Hence
F = -(a + b) modulo 4.
Case 4a > 3b. For even m = 2m', F >= 2(2a - b) m' - a - b m' = (4a - 3b) m' - a, which is nonnegative for
m' >= a/(4a - 3b); for smaller m', F >= -a - b m' > -N_0. So every element below -N_0 comes from odd m, where
F = (2a - b) m - a + b - b g(m). If 4 | b then (2a - b) m = 2am - bm = 2 modulo 4 (a, m odd) and b, b g(m) are
multiples of 4, so F = 2 - a modulo 4. If b = 2 modulo 4 then 2a - b = 0 modulo 4, b(1 - g(m)) is a multiple of 4,
and F = -a modulo 4. In both cases the elements below -N_0 occupy one residue class modulo 4 and S misses the
other class there. QED

*Proof of 3.* Follow Session 198 section 2 with the first letter Z. Since b is odd there is k >= 1 with B = 2^k = 1
modulo b; put A = a(B - 1)/b, an integer. Both (B - 1)/b and a are odd, so A is odd; 1/2 < a/b < 1 gives
(B - 1)/2 < A < B - 1, hence B/2 <= A < B - 1 and the bit word W of A (least significant bit first, O for 0, Z for 1)
starts with Z, ends with Z and is mixed. Iterating O and Z over one period from an input (D, E, U, 1) gives

    D' = D + A U - C_lin - eps_0 E,   U' = B U - R,   E' = 1,

with R, C_lin as in Session 198 (3) (the term -eps_0 E is the ending count seen by the first letter when it is Z,
which Session 198 did not need because its first letter is O). Since bA = a(B - 1), the response obeys

    F' = a U' - b D' = F - K_lin + b eps_0 E,   K_lin = aR - bC_lin.

On a strong input (E = 1) the drop is K_cyc := K_lin - b = aR - bC_cyc, where C_cyc adds the cyclic term
eps_(k-1) eps_0 = 1 to C_lin; K_cyc > 0 by Session 198's run argument, which is formulated cyclically (their (5)
with cyclic indexing of the ending counts), and K_cyc is even because F is even on strong sources and W preserves
strongness. Given an even target r, choose odd residues D_0, U_0 modulo K_cyc with 2aU_0 - bD_0 = r + b modulo
K_cyc (write D_0 = 1 + 2x, U_0 = 1 + 2y; the congruence halves to 2ay - bx = (r + 2b - 2a)/2 modulo K_cyc/2, solvable
because gcd(2a, b) = 1), take a strong source with these residues and (2a - b) D >= r + b + K_cyc
(STRONG_SOURCE_INTERFACE), append one O (response S = 2aU - bD >= (2a - b) D, E = 0), then W exactly h times: the
first period drops by K_lin = K_cyc + b (E = 0), the others by K_cyc, so the final response is S - b - h K_cyc = r for
the positive integer h = (S - b - r)/K_cyc. The final point ends with Z: strong, s = E = d = 1, transitive root
deletion. QED

The constants for the simplest slopes: 3/5: W = ZOOZ, K_lin = 7, K_cyc = 2; 5/7: W = ZOZ, 11 and 4; 5/9: W = ZZOOOZ,
11 and 2; 7/9: W = ZOOOZZ, 11 and 2; 7/11: W = ZZOZOOOZOZ, 29 and 18; 9/13: 31 and 18. For Session 198's class
K_cyc = K_lin (first letter O): 2/3: OZ, 2; 4/5: OOZZ, 2; 6/7: OZZ, 4; 4/7: OOZ, 4; 8/11: OOOZOZZZOZ, 18.

## 5. Computations (FINITE-EXACT)

- `signed_slope_probe_20261105.cpp` (n < 2^24, targets in [-999, 999]): class (i) slopes attain 1000, 986, 981, 931,
  958, 843, 875, 856, 862, 868 of the 1000 odd targets for (2,3), (4,5), (6,7), (8,9), (4,7), (6,11), (8,11), (8,13),
  (10,11), (10,13) (the misses are large targets whose witnesses need n beyond 2^24: the recipe's h repetitions
  multiply U by B each time); class (ii) slopes attain 984, 980, 926, 935, 867, 867, 848, 865 of the 999 even targets
  for (3,5), (5,7), (5,9), (7,9), (7,11), (9,11), (7,13), (9,13), with the same edge profile; class (iii) slopes with
  a/b <= 3/4 have all 2^22 positive values in the single class -(a + b) modulo 4 ((3,4): 1; (5,8): 3; (7,10): 3;
  (7,12): 1; (9,14): 1; (11,16): 1; (9,16): 3), and those with a/b > 3/4 have both classes among positive values but
  one class below -200 ((7,8), (5,6), (9,10), (11,12), (13,14), (15,16): 3; (13,16), (11,14): 1), exactly as the
  theorem predicts; outside the interval one sign is absent ((1,3), (2,7), (4,3), (5,4)).
- `signed_clock_recipe_20261105.py`: for each of the 11 slopes of classes (i) and (ii), the per-period drops K_cyc
  (strong input) and K_lin (input after an O) are verified on 200 random strong orbit points each, and 40 random
  targets of the right parity in [-10000, 10000] are hit exactly by seed -> O -> W^h with the final point strong
  (up to 4,966 repetitions; U then has about 10,000 binary digits, so the witness is recorded as (seed n, h), as
  Session 198 says).
- Pre-registered outcomes: P1 supported (small-witness coverage; the theorem is Session 198's), P2 supported and
  then PROVED (3 above), P3 CONFIRMED and PROVED (4 above, first clause), P4 REFUTED (above 3/4 the parity
  hypothesis is still necessary: the negative side is confined; the first draft guessed the opposite with low
  confidence), P5 CONFIRMED, P6 CONFIRMED.

## 6. The positive-count barrier

**Lemma (THM-4882 (4)).** Let p, W, q be words of an apex grammar G(k, 0) (k apexes, transitive chain) and h >= 0.
Then H(T(p W^(h+1) q)) >= H(T(p W^h q)).
*Proof.* Deleting the chain vertices of one period of W from T(p W^(h+1) q) leaves the induced subtournament on the
remaining vertices, whose arcs are the chain order and the apex incidences of the remaining letters, i.e. T(p W^h q);
H is monotone under induced subtournaments (Redei insertion, THM-1370). QED
The same holds in G(k, w) when every slot of W closes inside the same period. Consequently a fixed repeated block can
only raise a count, no analogue of the signed drop -K exists for H, and the exact programming "seed, O, W^h" of
Session 198 cannot be transported to positive targets: the signed functional is programmable precisely because the
eigenvalue B = 2^k of U is cancelled in aU - bD by the choice bA = a(B - 1), while every readout of a physical count
is a nonnegative combination of cut-state counts in which the dominant eigenvalue survives. Coverage of positive
targets must therefore vary the word with the target (S79: the arcs are digits, with at least log N / (2 log log N)
of them), which is what G(0,2) and G(1,1) do in the censuses of THM-4878.

## 7. Beyond: what the clock says about the two-radix grades

In G(2,0) a block d has eigenvalues 1 (with a Jordan block of size up to 3), 2^L and 3^L (Session 194), so the
eigenvalue-1 left eigenspace carries two independent drift functionals (degree-1 and degree-2 counters). The
analogue of Session 198's theorem would program a pair of signed functionals of the 17-dimensional profile jointly;
HYP-9475 (B) records this as the next test, together with the question whether the recipe's witness size (n of
about B^h digits) can be reduced by mixing two clocks with coprime drops.

## 8. Audit

PENDING.
