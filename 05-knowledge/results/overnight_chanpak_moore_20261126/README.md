# Overnight tandem: invert the mutable contact before searching its suffix

2026-11-26, first stage. **PROVED** exact inverse and anchored pruning; **FINITE-EXACT** checks in the attached replays. Neither Chan–Pak nor Moore57 is resolved. This resumes the [last exact band and low-spoke theorem](../chanpak_exact_band_and_low_spokes_20261126/README.md) and the [physical Valiant cut](../moore57_quadratic_tensors_20261125/VALIANT-MARKED-BOARD-RANK.md).

## 1. Replace the cube-root estimate by a square-root window

The inverse here asks for literal one-edge crossings of a specified even midpoint E, not merely a residue or separate witnesses. In F2 take a binary rooted block C(v) with path count F(v), and an actual right rooted block with profile (h,e,R): its total paths, paths ending at the root, and ordered two-path covers whose second path starts at the root. The [seam proof](../chanpak_exact_band_and_low_spokes_20261126/README.md) gives H(v)=hF(v)−e[v even]+R(v+1). Always 0≤e≤h≤R: break each full path just before its first occurrence of the root to inject it into the indicated covers.

For a physical spoke reversal at bit i, take ℓ=2^i and write the smaller address v=2ℓq+r with 0≤r<ℓ, so the bit being flipped is zero. Pad both tournaments to the same spine size if the bit is newly leading. Let Z be the reverse complement of the i-bit suffix r, ε=Z mod 2 for i>0 (zero when i=0), and β=q+[q odd]. The padded concatenation identity gives

    F(v) = F(q)+β + 2qZ + F(r)−1,
    F(v+ℓ)−F(v) = Z+ε−β.

For the left rooted block the extra high zero vertices are beaten by all old spine vertices and by the root, so padding leaves the ordinary count, root-start count, and two-slot count unchanged. Thus this is a literal edge flip at fixed order, not insertion.

Write the required signed increment as 2s with s=±1, put γ=e for i=0 and γ=0 otherwise, and put w=Rℓ+hZ. Then the only candidate prefixes satisfy

    hβ = w + γ − 2s + hε,
    β is a nonnegative even integer, and q = β or β−1 (q≥0).

Substitution into the midpoint gives the exact value

    M = 2wq + R(r+1) + h(F(q)+β+F(r)−1) − e[v even] + s.

The suffix r is not a free correction after w is chosen: it is decoded from Z=(w−Rℓ)/h by bit reversal and complement. This is the important difference from choosing a modulus independent of the two contacts.

Because |q−w/h|≤2+2/h, F(x)≤x+1, h≤R, and r<ℓ, the error after the leading term satisfies |M−2w²/h|≤24w. For clarity: the first term costs at most 8w of error; the remaining nonconstant terms cost at most 4w+6h+5 together; and h≤w. Dividing the difference of squares by w+sqrt(hE/2) shows that a crossing of E requires

    |w−sqrt(hE/2)| ≤ 12h.

For fixed h,e,R, the powers ℓ are confined by Rℓ≤w<(R+h)ℓ≤2Rℓ. For each ℓ the permissible Z are an interval of at most 26 consecutive integers, before the divisibility and parity tests. Thus an exact search across *all* positions is possible without scanning any exponentially long suffix. The inverse also tests the original counts after each hit.

### Verified searches

* Every one of the 2,732 rooted profiles in the [inherited order-seven bank](../../../session203/coverage-recertification/session203-independent-audit/aud_bank7_inherited.txt) was checked for all spoke positions at exponents **31, 61, 89, 107, 127**. None crosses exactly. This rules out this entire F2 lane, not other edits or the target itself. The literal inverse is in [fixed_bank_inverse.py](fixed_bank_inverse.py).
* In F1 the unflipped block has an even positive address b=2^a d, d odd, so h=F(d−1), e=0, R=b. The same proof gives the tighter error bound 18w: the prefix deviation term costs at most 6w and the other terms at most 4w+3h+5. Hence use radius 9h. Set u=bℓ=d2^(a+i); then u≤w<1.5u since h≤d and a≥1. First try only the powers u meeting this narrow interval, and test h | u−2s before splitting the power between a and i. Then Z is determined within a short window, so each surviving split takes only a handful of exact reads. The [larger replay](inverse_large.out) excludes every b with odd part at most **2²⁷** and all spoke indices at exponent 61. After divisibility and window checks only 16 leading lanes, 78 dyadic splits, and 2,948 exact midpoint reads survive. This includes arbitrarily large dyadic padding: the midpoint bounds eliminate the rest, not an arbitrary exponent cap.

### Analytic pure-dyadic corollary

For b=2^a with a≥1 (h=1), there is **no** such crossing at any p≥7, on any right spoke. Put U=2^(p−1). The radius-9 window forces |w−U|<10, and u≤w<1.5u with u a power of two forces u=U. Then q=U+Z+ε−2s−δ, δ∈{0,1}. For s=−1 we have q≥U+1, and the midpoint exceeds E_p=2U²−U. For s=+1 and Z+ε≥2, q≥U−1 and β≥U. The possible negative contribution 2U(q−U) to M−E_p is absorbed by U+β alone; b(r+1)>0 makes the remainder strictly positive. For Z+ε=0 the suffix is all ones: Z=0, r=ℓ−1, and q=U−2 or U−3. Substitution gives M−E_p=−U/2 or −11U/4, respectively, using F(U−2)=U/2+1 and F(U−3)=U/4+1. No equality survives. This strengthens the earlier lowest-spoke result, but says nothing about an arbitrary left odd part.

## 1b. Two independent unit exits cannot live on the same binary arm

**Theorem.** For an arbitrary actual rooted right block in F2, if two spokes of the binary left block give exits −2 and +2 from the same tournament, the central count is **3 or 17**. This is an all-order obstruction. It is stronger than the prior result about the first spokes of the two different F1 arms.

**Proof.** Fix the central binary address v. For bit i let q_i=⌊v/2^(i+1)⌋, and set W_i=R2^i+hZ_i, where Z_i is the reverse complement of the i lower bits. The forward increment (always read from bit 0 to bit 1) is W_i+hε_i+γ_i−hβ(q_i), with γ_i=e when i=0 and zero otherwise. The actual increment is the forward increment times 1−2x_i. In particular the unit condition at the lower of the two bits gives

    h q_i ≤ W_i+h+2,

since hε_i+γ_i≤h. If the higher bit is j=i+m, then W_j≥2^m W_i and q_j=⌊q_i/2^m⌋. Its forward increment is at least

    (2^m−2^(−m))W_i − (1+2^(−m))h − 2^(1−m).

When i≥1, W_i≥2h, so m≥2 is impossible even for h=1; m=1 forces W_i≤h+2, hence h=1 and i=1. When i=0, m≥3 is impossible because W_i≥h; m=2 forces h=R=1. For i=0, m=1 the inequality gives R≤h+2. If h≥5 is odd, the upper unit condition gives h | 2R∓2. Among R=h,h+1,h+2 this leaves only R=h+1. Write y=x₀. The low unit condition is now hβ(q₀)=R+e±2; its size and even multiplier force β(q₀)=2, so q₀=1 or 2. With q₁=⌊q₀/2⌋ the high forward increment is 2(h+1)+2h(1−y)−hβ(q₁). It can have magnitude 2 only when y=1 and q₀=2, giving central v=5. Bit 1 is zero there, so the high actual change is +2. The low actual change is h−1−e≥−1, not −2.

It remains to check h=1 or 3, R≤h+2, 0≤e≤h, bit pairs (0,1),(0,2),(1,2), and v<32. The last bound follows directly from hq_i≤W_i+h+2 and the bounds on W_i or R above (in the i=1 case it gives v<28). The complete [integer check](same_arm.py) finds five rows. One has (h,e,R)=(1,0,1) and is impossible: H=1 forces transitivity by the odd-cycle packing formula, and for a transitive rooted block R is 2 raised to the number of vertices after the root, so R=1 forces e=1. The genuine rows are (1,1,1) with v=1,8,15 and marked pairs (0,1),(0,1),(1,2), respectively, and (1,0,2) with v=7 and pair (0,1). Their central totals are 3,17,17,17. The replay also verifies their literal orientations and both edits by Hamilton-path subset DP. This leaves opposite-arm edits, edits inside the right core, and larger construction grammars open.

## 2. Return to Moore: attach the cost to actual occupied slots

The shared-block change in Section 1 worked because the signed interaction was tracked before replacing the suffix by an absolute estimate. On the marked Moore board the corresponding exact quantity is the list of hits to fixed slots, not just a binary-overlap dimension or the coordinate edit count.

First fix all 50 existing non-marker contacts. The two holes (13,33) have 32 and 33 choices. Exhausting the 1,056 actual pairs finds minimum **three squares**, attained at (16,19), (16,20), and (29,16). Their center score profiles (56,0,28) are (1,2,0), (1,1,1), and (1,1,1), respectively. Thus a best unedited pair can balance its defects across all three centers or already make one center exact, but cannot fix any of those defects without moving an old owner. See [hole_profiles.py](hole_profiles.py).

For cap two, there must be one exact center, since every collision costs one or more and the three center scores are nonnegative integers. There is a useful finite description in any exact center: send each changed branch to the old owner of its selected slot, or to a vacancy symbol if the slot was initially unused. The changed owners and the two vacant slots have indegree one. The two original hole branches have indegree zero, and the two vacancy symbols outdegree zero. Every other changed owner has outdegree one. Therefore the motion is two directed chains from the two branch holes to the two initial vacancies, plus directed cycles. This applies separately to any exact center; the hard step is selecting the same physical contacts for both chains and the other centers. It suggests search by chains with edit interventions rather than treating all low-score signatures as unconstrained.

An exact variant of the anchored search also charges immutable collisions sooner. At a node, fixed newly chosen contacts and locked old owners already cannot move. Each open branch must independently pay at least its cheapest collision count against those immutable contacts, among the choices passing the existing optimistic edit bound. These counts refer to distinct pairs involving different open branches, so their sum cannot exceed the remaining square allowance. The new prune is in [anchored_immutable.cpp](anchored_immutable.cpp), which also verifies directly that different-branch contacts share at most one center slot, the condition allowing an owner bitset to count their collisions. The ordinary rooted-exchange completeness argument is still needed for using the anchored search as a nonexistence test within an edit cap; this is not an unrestricted negative.

To exploit the exact-center observation, [anchored_centers.cpp](anchored_centers.cpp) splits cap two into the six centerwise caps given by the permutations of (2,0,0) and (1,1,0). In a saturated center, any active owner hit by the next contact must be displaced, unless it is locked, in which case the choice is illegal. Parent transitions and immutable collision charges are checked in each center separately. If some allowed solution has a component of moved owners unreachable from the two vacancy branches under old-slot hits, reverting that entire component introduces no collision with retained moves in **any** center. Thus the completeness argument survives centerwise caps. All six searches excluded cap two through **10 reference edits**, agreeing with both unpartitioned versions; see [the receipt](moore_search.out). There is no cap-two exclusion at a larger budget claimed here.

A separate deterministic-seed annealing search ran 500 restarts of up to 400,000 attempted moves without beating score three. That is only a heuristic negative. Replayer [low_square_search.cpp](low_square_search.cpp) uses the existing board output from 'nested_core_check.py --board'.

## Next linked tests

On the Moore side, use the three exact-center branches separately only where total score ≤2, keeping the exposed two-chain endpoints. On the Chan–Pak side, extending the inverse to unbounded odd part needs a second independent contact reading, and the same-arm obstruction now tells us that a scalable two-exit design must look outside that binary arm. In both directions the distinction is the same: an arithmetic match after projecting away the contact addresses is not a physical lift.

## Replay commands

From this directory (C++17 and Python 3):

```sh
bank=../../../session203/coverage-recertification/session203-independent-audit/aud_bank7_inherited.txt
python3 fixed_bank_inverse.py "$bank" 31 61 89 107 127
python3 exact_bit_inverse.py 61 134217728
python3 same_arm.py "$bank"
python3 hole_profiles.py
python3 ../moore57_quadratic_tensors_20261125/nested_core_check.py --board > /tmp/overnight-moore-board.txt
g++ -O3 -std=c++17 anchored_centers.cpp -o /tmp/overnight-anchored-centers
for spec in '2 0 0' '0 2 0' '0 0 2' '1 1 0' '1 0 1' '0 1 1'; do
  /tmp/overnight-anchored-centers $spec 10 10 1000000000 < /tmp/overnight-moore-board.txt
done
```
