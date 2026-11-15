# Cancellation of both large turnback channels leaves a constant 48

The [return law](TURNBACK_FOCUS.md) gives a narrower literal construction if two old cores, rather than just two masks on one core, can change. With the return at the tip (L=j≥2), it reads
```
 H_j = 2^{j-1}U_1+(H(R)-J_1),qquad
 Q_{2,j} = 2^{j}(Q_2(R)-J_2)+2· 3^{j-1}U_2+2V_2.
```
To preserve (H_j) and obtain a constant difference in (Q_{2,j}) on an unbounded common clock, this makes an exact recipe: match (U_1, H(R)-J_1, Q_2(R)-J_2, U_2), but change (V_2). It is also necessary for those two identities in this grammar: the constant, (2^j), and (3^j) terms are linearly independent on an unbounded sequence of (j)'s, the dominant exponential forcing its coefficient to agree first.

A physical match exists using just one backward contact in each arm:
```
R_1 rows (172,225,2,70,79,92,133,60),    M_1=4
R_2 rows (86,92,152,1,232,143,172,11),  M_2=8
```
Their common ((H,Q_2,U_1,U_2,J_1,J_2)) is ((389,6032,406,3196,0,0)). The signed fields (V_2) are (-1508,-1484), respectively. For all (j≥2) the two ((8+j))-vertex spines have identical (H=389+406·2^{j-1}) and the second has (48) more ordered two-slot covers. The two growing channels cancel **exactly**; no asymptotic approximation is involved.

Wrap each returned spine (W) by two singletons in the cyclic orientation (W → p → q → W). The [ordinary outlet readout](../chanpak_carry_admission_20261115/OUTLET_WIDTH.md) gives (H_{wrap}=H(W)+Q_2(W)). Thus the same exact constant +48 becomes an ordinary Hamilton path shift at all (j≥2). The physical returned spines are strong: each old (R_i) is strong, can reach every tail vertex, and every tail vertex can reach the return tip, which reaches (R_i). In any case the cyclic singleton wrap is strong. The result is an exact readout lock, not yet a full matched ((h,S,d)) continuation or a closure of the wider folded window.

There is a small arithmetic floor for a constant (Q_2) difference in any such paired family of equal base orders and equal (H_j) for all large (j). Every odd-cycle packing contracts to a number (c) of blocks with (c ≡  n mod 2), so (Q_2 ≡ (-1)^nH mod 3). At all sufficiently large even orders the slot polynomial has exponents at least two, and all nonempty-packing coefficients are even, giving (8 divides  Q_2) (the empty (2^n) term is also divisible by (8)). Hence an all-(j) constant shift must be divisible by 24. The found 48 is twice that floor; no claim is made that 24 itself is achievable or excluded.

Replay: `python replay_turnback_lock.py` validates both kernels and literal Hamilton and slot counts before and after the wrap; the infinite equality follows from the derived law.
