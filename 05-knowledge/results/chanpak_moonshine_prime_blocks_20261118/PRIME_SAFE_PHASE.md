# Prime-safe phase programming and a two-base absolute decoder

This follows [the crossed cyclic seam](POLES_TO_PRIMES.md). It separates the part of a prime request that can now be programmed freely from the absolute-height condition. Neither theorem here assumes that an integer close to a prime will be prime.

## 1. A physical low-bit program with the 3 and 5 addresses intact

The two signed-reader toggles with ordinary counts 1+3·2^j (j≥2) and 1+5·2^j (j≥3) are literal substitutions into small quotients. Their exact heights B are:
```
1+3·2^j:   B=−12+16·2^(j−1)+28·3^(j−1)
1+5·2^j:   B=−88+112·2^(j−2)+192·3^(j−2).
```
The optional-slot coefficients come from expanding the q=2 run vectors (104,100,28), (712,688,192) into the 1^k,2^k,3^k basis. Each additional visit within a given path needs the ordinary singleton through which the spine is re-entered, so the cutoff q+1 is sufficient. The first bit uses the literal count-11 tournament with row masks (16,17,3,7,12), height 100.

**Phase theorem.** Given L≥1 and any prescribed unit p modulo 15·2^L, one can construct an actual source with H≡p modulo 15·2^L. Moreover a fixed paired wrapper produces actual counts p' and 7p' with p'≡p modulo 15·2^L. The combined construction uses O(L²) vertices and computes both H and B exactly from the word of choices.

Start from H=1 and iterate j=1,...,L-1. If bit j disagrees, multiply H by a physical factor congruent to 1+2^j modulo 2^(j+1). At j=1 use 11; at j≡3 modulo 4 use 1+5·2^j; otherwise use 1+3·2^j. Each is a unit modulo 15. For any odd H this multiplication toggles precisely the requested bit among the bits processed so far. Realize multiplication by ordinal joins; it multiplies every slot reading, not just H.

Then choose a≥L with a≡1 modulo 4 and a≥2. The factor 1+3·2^a is neutral modulo 2^L and 3 and is 2 modulo 5, a generator of the units there. At most three copies set the desired mod-5 address. An odd exponent a'≥max(L,3) gives 1+5·2^a', which is neutral modulo 2^L and 5, and is 2 modulo 3. Use it once if needed. This proves joint phase programming, not three separate surjections.

Use the two archived signed-height scaffolds of slopes 60 and 420, intercepts 1 and 7 (the literal cores and mark masks are in the replay). Pad the common input by L ordinal singletons. The resulting counts are exactly p'=60·2^L B+H and 7p'. Because 15 divides 60, this preserves the whole requested CRT phase. The naive use of 1+2^j for all toggles would not do this: odd j would insert a factor divisible by 3 into H, spoiling prime eligibility on a slope divisible by 3.

The phase result alone does not guarantee primes or bound the distance to the original prescribed integer. The information still missing is how to choose B at the same time as H, within an absolute budget.

## 2. A separate exact family using both exponential bases

Let P_m be the m+2 vertex transitive tournament with its longest arc reversed. Its labeled adjacency determinant is
```
1−x_first x_last (product_middle (1+x_i)−1).
```
Every cycle uses the single reverse arc and a nonempty increasing selection of middle vertices. Odd-cycle packing therefore gives
```
f(q)=q^(m+2)+q[(q+1)^m−(q−1)^m],
H=1+2^m,     B=4·2^m+2·3^m−2.
```
This is the previously used long-cycle polynomial; the new application pairs it with the slope-60 wrappers and decodes from an absolute target. For
```
p=60·2^k B+1+2^m,
k+v_2(60B)>m,
```
the target itself fixes m=v_2(p−1); the unique k is found by checking whether (p−1−2^m)/(60B) is an integral positive power of two and checking the displayed inequality. The separately proved reader transfers build p and 7p without approximation. Verified examples include primes 13,457; 53,777; 27,525,137; and 217,252,097. Their odd mantissas of p−1 vary with the 3^m term, not just with the six fixed mantissas. There is no inference about the number of primes in this family.

A motif-level parallel to the archival pure-A Leech charts is exact but limited: the cycle of length h has univariate determinant 1−u^h; its cofactor after removing 1−u gives h at u=1, which is also the A_(h−1) Cartan determinant. The six selected h have h−1 dividing 12 (hence 24), making them pure-A Niemeier labels. The construction above retains the actual chord terms, and does not use the motif identity as a claim about prime placement.

## 3. The diagonal lattice and the arithmetic of a real repair

For equal incoming blocks in the crossed seam, put u=B−H. The exact count is the positive binary quadratic form
```
p=2u²+uH+2H²,   or   8p=(4B−3H)²+15H².
```
Thus one target-first branch is to solve this norm equation, then insist on a physical block of that (H,B), instead of using a formal norm representation as a tournament. For primes above 5 that appear in this branch, the count is 2 modulo 3 by the slot reflection congruence, and 2 or 3 modulo 5 by the norm identity, hence 2 or 8 modulo 15. The complementary residue branches need genuinely unequal blocks or other seams. TT2 and TT3 in both positions give literal primes 23 and 107.

For unequal blocks with equal q=2 score, assume a prime above 3 is at either end. Then neither block's H is divisible by 3: otherwise so is its B, and the crossed count is divisible by 3. Modulo 3 the scores (f(3)−H)/B recover the parity of each block order, so the orders have the same parity. The neutral ordinary difference Δ=h B_R−B_L x is divisible by 3. It is also divisible by 5: score equality identifies the same determinant with that of f(3) and f(2), and f(3)≡f(−2)=(-1)^order f(2) modulo 5. Packings show that B is divisible by 4 at odd orders ≥3 and by 8 at even orders ≥4. The slot-splitting inequality gives score ≥(order+1)/2, so the singleton (score 1) and the two-vertex block (score 2) cannot collide with larger same-parity orders. Nonzero prime-eligible repair steps are therefore multiples of 60 at odd incoming order parity and 120 at even parity. This makes 2,3,5 arise from actual reflected evaluations, not just a resemblance between short number lists.

The constructive split is now precise: the phase program solves the joint 2/3/5 unit-residue problem, the cycle family solves some absolute addresses, and the norm branch provides a different absolute-target inverse. A full proof still needs a physical profile theorem for the target-dependent heights; residue programming and arithmetic factorability do not themselves supply it.

Replay: `python3 crt_replay.py`, alongside the earlier construction files in this packet.
