# An exact change of core that keeps its Hamilton count

The modular Heisenberg moves keep selected *residue* coordinates. There is a separate, exact physical switch in the bounded-profile inverse: change Q while preserving h as an integer, without formal inverses or a larger vertex budget. This does not contradict the binary address obstruction, which keeps a different coordinate (the absolute cut start weight).

On the transitive spine use the strong words `ABBBABB` and `ABABBAB`, in the [gap-polynomial convention](GAP_HEISENBERG.md). The resulting order-eight tournaments have
```
(h,Q)_− = (61,1744) = (61,16·109),
(h,Q)_+ = (61,1864) = (61, 8·233).
```
The equality h=61 is also a check in the prime/7p lane: unlike the smaller count-33 collision, it does not force the final gluing count to be divisible by 3.

Take an ordinal sum of k of these blocks, with i plus switches, and s additional singletons. Both h and Q multiply, so the second core W has exactly
```
h = 61^k,   Q = 2^s·1744^(k−i)·1864^i,   0≤i≤k.
```
With k,s fixed the full order remains 8k+s. The family fills every 2-adic valuation between 3k+s and 4k+s while holding h and order fixed. At modulus 5^j the ratio 1864/1744 is a unit with ratio−1 of valuation exactly 1. Its powers have exact order 5^(j−1) (the elementary lifting-the-exponent formula or the binomial theorem). For k≥5^(j−1)−1 and fixed s, the switches therefore fill every Q residue in one 5-adic unit coset mod 5^j while keeping the **integer** h fixed. This is not arbitrary control of the odd part of Q.

For the absolute target, set t=61^k and D=N−(M−L)t in the anchored equation N=LQ+(M−L)t. Test that D/L is a positive integer Q. If Q belongs to this family then i=v_233(Q) and s=v_2(Q)−4k+i, so the witness is recovered uniquely or rejected by one factor check. For k=0 the construction is just s singletons and one requires s≥1. In particular this is a physical solver for a sparse subset of bounded progressions and an explicit way of moving *inside* a progression with exactly fixed t. It does not show that a target-induced progression meets this family. The alternative switch has no independent additive increment at fixed h, which is what a full exact return would need.

Replay: the two base profiles are checked by literal subset-endpoint DP in [replay_profile_switch.py](replay_profile_switch.py).
