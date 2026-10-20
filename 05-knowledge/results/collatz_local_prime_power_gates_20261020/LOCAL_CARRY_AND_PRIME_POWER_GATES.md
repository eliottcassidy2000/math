# Local carry content, bounded edits, and simultaneous prime-power gates

Date: 20 October 2026. Status: **PROVED elementary identities and two-odd-step classification; FINITE-EXACT one-swap census; general one-swap exclusion OPEN.** This report does not prove Collatz convergence, exclude all nontrivial positive cycles, or assert historical priority.

## 1. Inheritance, scope, and the retained coordinates

The starting mechanisms are THM-4604, `01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md`, and THM-4716, `01-canon/theorems/THM-4716-balanced-collatz-cycles-are-three-smooth-ladders-and-are-free.md`. The immediate predecessor is `05-knowledge/results/collatz_rotation_defect_ideal_20261020/LOCAL_ROTATION_DEFECT_IDEAL.md`, with `ONE_SWAP_TRINOMIAL_GATE.md` and their audits. Its connected-rotation-graph ideal theorem is inherited, not rediscovered here.

The hostile `000011` has clock 55, carry 80, and content 5. The corrected near miss is one-swap freeness: `0000101` has clock 119 and content 7. The least-used sidecar is the full valuation vector of the denominator, rather than its prime support. The concept board is carry content / local edit support / chronological rotation / prime-power depth / simultaneous witnesses / repetition. The anchor is local denominator recognition; the niche is the displacement-periodic two-one family; the wildcard is transferring only the local-global arithmetic mechanism, not a Collatz implication, to branch integrality.

For a chronological binary word w of length q with p ones, 0<p<q, put

    B(w) = sum_(i:w_i=1) 2^i 3^(number of ones after i),
    D = 2^q - 3^p,    g(w) = gcd(|D|,B(w)).

The shortcut branch is (3^p x+B(w))/2^q. Its rational fixed point is B(w)/D; its reduced denominator is |D|/g(w). All rotations have the same content. An integer 3x+1 cycle with this word exists exactly when g(w)=|D|, and it is positive exactly when D>0. This is a recognition criterion for a supplied word, not a universal exclusion theorem.

## 2. A single chronological edge gives the exact content

**Theorem A.** Let w' be the left rotation of w by one position. Strip any common prefix P and suffix S, writing w=PuS and w'=PvS. Then

    g(w) = gcd(|D|, B(w')-B(w)) = gcd(|D|, B(v)-B(u)).                 (A)

In particular, a connected graph of all rotations is unnecessary once D is retained: any one chronological neighboring edge suffices.

**Proof.** If the first bit is epsilon, the affine cocycle gives

    2 B(w') = 3^epsilon B(w) + epsilon D.

For epsilon=0, twice the carry difference is -B(w); for epsilon=1 it is B(w)+D. Since D is odd, either identity proves the first equality. Equal-shape context removal gives

    B(PvS)-B(PuS)=2^|P| 3^(#S) (B(v)-B(u)).

As gcd(D,6)=1, this proves the second equality. The two rotations are distinct because a word fixed by one-place rotation is constant, excluded by 0<p<q. ∎

Consequently, if the disagreement core has length ell and h ones, then

    g(w) <= |B(v)-B(u)| <= (2^(ell-h)-1)(3^h-2^h).                   (B)

The right bound follows by moving all ones to the left or right: the minimum and maximum carries are 3^h-2^h and 2^(ell-h)(3^h-2^h). Thus a bounded-width chronological disagreement gives a bounded content, independent of the common outside context. It excludes an integer cycle when |D| exceeds that bound. Small edit *count* alone is not small edit *span*; (B) retains both ell and h.

**Exact arbitrary-rotation refinement.** Let u be the prefix of length r, a=#u, and w_r the left rotation by r, 1<=r<q. Then

    2^r B(w_r) = 3^a B(w) + D B(u),
    gcd(|D|,B(w_r)-B(w)) = g(w) gcd(|D|/g(w),3^a-2^r).              (C)

The first identity follows by applying r shortcut branches to B(w)/D. For the second, reduce the first modulo D, factor g, and use gcd(B(w)/g,|D|/g)=1. This explains exactly why a nonchronological pair can overestimate content: it retains a subclock factor. If w_r=w in a periodic rotation, the formula still holds with a zero difference. For the one-place step, the subclock is ±1, so no extra factor survives. This special simplification uses multiplier 3; it is not automatic for every odd multiplier.

## 3. Bounded swaps give a sparse gate, not uniform exclusion

Fix the zero-intercept lower mechanical word

    c_i=floor((i+1)p/q)-floor(ip/q),  gcd(p,q)=1.

Suppose w is obtained by a sequence of t swaps of unequal cyclically adjacent bits. The positions and order of the swaps are part of the data. Each swap has one exceptional rotation cut, where its two positions occur at opposite ends of the linear word. Let E be the set of these cuts; |E|<=t, even if swaps reuse positions.

For any cut r outside E, telescoping the successive swaps gives

    B(w_r)-B(c_r) = sum_(j=1)^t epsilon_j 2^(a_j(r)) 3^(b_j(r)),

with the exponents computed at the j-th intermediate word, and epsilon_j=±1 according to the swap direction. The suffix counts must not be copied from the initial word after subsequent edits.

Order the balanced rotations by their intercepts. Each neighboring carry difference in that ordering is one signed 2–3-smooth monomial by THM-4716. At least max(0,q-1-2t) of its q-1 path edges have both endpoints outside E. For each such edge the corresponding difference for w is a sum of at most 2t+1 signed 2–3-smooth monomials. The content g(w) divides all these differences, hence divides their gcd together with D. If that latter gcd is proper in |D|, the supplied word is excluded.

This is a **PROVED sparse reduction**, not a bound depending only on t: the monomial exponents grow with q. The surviving path may be disconnected, so its gcd is not asserted to equal g. Adding a chronological edge invokes (A); adding a connected rotation graph invokes the inherited exact ideal theorem. The exceptional cuts cannot be discarded silently.

## 4. Exact two-odd-step theorem, including even periods

**Theorem B.** Suppose w has exactly two ones. Let d be the smaller cyclic distance between their positions, so 1<=d<=floor(q/2). Then

    g(w)=gcd(2^q-9,3+2^d)
        =gcd(3+2^d,2^(q-2d)-1).                                  (D)

Use absolute values for gcd if q=3. When q=2d the second entry in the final gcd is zero, and g=3+2^d.

**Proof.** Rotate to place ones at 0,d, yielding carry 3+2^d. Modulo that carry, 2^d=-3, so

    2^q-9 = 2^(q-2d) (2^d)^2-9 = 9(2^(q-2d)-1).

The carry is coprime to 3; cancel 9 inside the gcd. Rotation preserves g. ∎

For odd q=2k+1 and d=k-s, 0<=s<k,

    g(w)=gcd(3+2^(k-s),2^(2s+1)-1).                               (E)

For each fixed displacement s, the content is bounded by 2^(2s+1)-1 and is periodic in k with period dividing 2s+1. This is equality of the full gcd, including prime-power multiplicities, not only of prime supports. Periodicity follows because 2^(2s+1)=1 modulo that fixed odd modulus.

- s=0 is the balanced branch: g=1.
- For the **nonbalanced one-swap branch**, s=1 and k>=2, g is 7 if k=0 mod3, and 1 otherwise. Equivalently g=7 exactly when q=1 mod6. The balanced outputs of a swap remain in the s=0 branch and always have g=1.
- Higher displacement does not preserve the squarefree, small-prime picture: (s,k)=(10,27) gives 49; (36,99) gives 439; (58,102) gives 6559. Their periods are q=55,199,205 respectively.

For every q>=5, 3+2^floor(q/2)<2^q-9, so every two-one word has g<D and cannot encode a positive integer cycle. At q=4, d=2 is the repetition of the ordinary 1–2 cycle and has g=D=7; d=1 has g=1. At q=3, D=-1 belongs to the negative cycle, not a positive one. This elementary two-one exclusion is sharply scoped and no priority claim is made for excluding cycles with so few odd steps.

## 5. Prime powers, mixed primes, and the quantifier boundary

Write |D|=product_l l^(a_l). The exact criterion is

    integer cycle <=> for every l|D, v_l(B(w)) >= a_l.              (F)

Testing only B(w)=0 modulo every prime dividing D is insufficient. The explicit primitive word `11101000100` has

    D=1805=5*19^2,    B=475=5^2*19,    g=95, reduced denominator=19.

Every prime in D divides B, but the second factor of 19 is missing. A prime-support quotient destroys the very information that decides integrality.

A second loss is failure to keep the same word across different prime tests. At shape (q,p)=(15,8), D=26207=73*359. The complete finite necklace census has contents

    g=1: 418 necklaces;   g=73: 10;   g=359: 1;   g=26207: 0.

Each prime separately has a witness, but no common witness exists at this shape. Thus “for every prime there exists a word” cannot replace “there exists one word satisfying every prime-power condition.” This is the precise arithmetic analogy with mixed-prime symmetry refinement: separately attainable one-prime constraints need not glue. It is not a map identifying Collatz words with symmetry classes.

At (q,p)=(11,5), the necklace contents are 1:34, 5:5, 19:2, 95:1. At (7,2), they are 1:2, 7:1. In all three shapes gcd(q,p)=1, so every rotation orbit has q elements; these are true necklace counts, unlike the swap-occurrence census below.

## 6. Repetition and valuation inflation

For a nonconstant word w and integer r>=1, define

    H_r=(2^(rq)-3^(rp))/(2^q-3^p)
       =sum_(j=0)^(r-1) 2^(q(r-1-j)) 3^(pj) > 0.

The affine cocycle gives B(w^r)=H_r B(w), D(w^r)=H_r D(w), and g(w^r)=H_r g(w). Therefore the reduced denominator is unchanged. Large content alone can merely reflect repeated presentation, not improved integrality. The primitive/necklace-period sidecar is indispensable. For example 01 has D=1,B=2; 0101 has D=7,B=14, but both describe the same rational cycle.

## 7. Finite evidence and reproduction

`verify_local_integrality.py` uses no repository mathematical-code imports. Its recurrence and positional carry formulas are independent implementations. The recovered receipt records:

- 89,958 chronological neighboring-edge checks on every nonconstant binary word of lengths 2 through 12.
- 2,470 word-power cases: every nonconstant word of lengths 2 through 8, repetitions 1 through 5.
- 20,099 odd two-one displacement cases, 2<=k<=200, 0<=s<k.
- Complete fixed-shape word and necklace censuses at (7,2), (11,5), (15,8).
- 6,492,980 directed unequal cyclic swaps of zero-intercept lower mechanical words with coprime q,p and 3<=q<=400. Duplicate outputs are retained. Of these, 4,723,382 have D>0 and 1,769,598 have D<0. Nonunit content occurs 494,393 times; the maximum observed content is 3,924,169 at (q,p,i)=(350,79,186).
- No swap occurrence with |D|>1 has g=|D|. This is **FINITE-EXACT**, not an all-length theorem. The q<=400 universe must not be described as all necklaces, all near-balanced words, or a trajectory convergence test.

Reproduce from this directory:

    python3 verify_local_integrality.py > receipt.replay.json
    python3 -O verify_local_integrality.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json
    python3 independent_controls.py > independent-controls.replay.json
    cmp independent-controls.json independent-controls.replay.json

The independent script passes 223,248 checks: full affine compositions, arbitrary rotation factors, all-period two-one formulas through q=500, and bounded two-swap gates through q=12. Its direct positional q<=120 swap census has 175,668 occurrences, 12,924 nonunit contents, and no integral candidate with |D|>1. These proofs and checks do not resolve the all-length positive one-swap conjecture from the predecessor packet.

## 8. Next exact targets

1. Use the prescribed 2t+1-term relations to force a proper divisor for a precisely stated positive-clock one-swap family with p>=3. The two-one result does not dispose of this case.
2. Determine the exact minimal k-period of (E) from multiplicative orders modulo each prime power dividing its Mersenne modulus; 2s+1 is a guaranteed period, not always minimal.
3. Classify feasible simultaneous valuation vectors at fixed (q,p), retaining a common necklace witness and its primitive period.
4. Find a combinatorial condition forcing a short chronological disagreement core. Mere bounded swap count or a small collection of separate modular witnesses does not suffice.
