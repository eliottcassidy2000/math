# One-sided recovery gives a coprime positional switch bank

The tangent came from Neményi's earlier [stress transport](../session285_five_web_20261112/stress/PROOF.md): an observable can be fixed on one boundary without fixing its response at the other. The graph-side distinction is exact, not an analogy to elasticity. The existing [matched Bézout return](../chanpak_tail_return_20261115/RECONTACT_ENERGY.md) switches blocks at fixed incoming and outgoing interfaces. Here a block is *invisible at its outgoing port* but visible to its incoming port. When an outgoing-invariant block is repeated, later changes can therefore see the growing prefix without knowing which earlier switches were chosen.

## Physical coordinates

A painted spine is a tournament B and a subset beaten by a common apex. Its data are ordinary Hamilton count h, count n after its own apex is added, two optional slots Q, and endpoint response (S,d,T,e). Here S,d are the total and odd-minus-even complement-weighted paths starting at a marked vertex, while T,e are the corresponding totals ending at an unmarked vertex. For painted ordinal concatenation, appending B to a left spine of data (g,N,A,D) yields

```
g' = h g
N' = h N + (n−h)g + T A + e D
A' = Q A + S g
D' = d g.
```

This is the [inherited cut-state law](../chanpak_selective_target_20261114/CUT_STATE.md). Fix two choices of painting on B having the same S,d. Changing a choice never affects the outgoing h,A,D, and so never changes the contribution of any subsequent choice to the difference. But if T differs, its local effect depends on the live incoming A. This strictly extends the mechanism in which all endpoint responses are held fixed.

Use the literal 7-vertex spine with bit rows
`(44,121,10,0,13,92,29)`.
The two choices label only vertex 1 (mask 2) or only vertex 5 (mask 32). Their exact data are:

|mask|h|Q|n|S|d|T|e|
|---:|---:|---:|---:|---:|---:|---:|---:|
|2|37|704|155|118|0|630|32|
|32|37|704|155|118|0|624|30|

## Independent positional digits on an identical outgoing interface

Take a single marked singleton, followed by one fixed copy of B painted with mask 2. Then append m cell pairs of copies of B, and a last unmarked singleton. Each cell, independently, is painted either (2,32) or (32,2). It is the same unpainted spine at every corner of this cube. Every resulting tournament is strong: the common apex reaches the first spine vertex, the last reaches the apex, and the joins are forward throughout the ordinal spine.

After initialization, normalized outgoing coordinates are x=A/g=822/37 and D/g=0. They are independent of choices. Each individual copy of B updates x to

```
F(x)=(704/37)x + 118/37.
```

For a cell standing at coordinate x, changing its order from (2,32) to (32,2) contributes exactly

```
6 (P/37) (F(x)−x),
```

where P=37^(2m+1) is the ordinary spine count of the entire board. To see this, using mask 32 in place of 2 changes the local N by −6A since D=0; the total downstream propagation of each local change is the product of ordinary counts to the right. Contributions from the two positions have opposite signs. Since

```
F(F^i(x))−F^i(x) = (704/37)^i · (704·785/37²),
```

the complete absolute count family is

```
N = N0(m) + 3315840 · sum_(j=0)^(m−1) ε_j 37^(2(m−1−j)) 704^(2j),   ε_j∈{0,1}.
```

N0 uses all cells (2,32) and is determined by the displayed integer recurrence, with no fitted term. The last singleton adds A+D and resets D to zero; it leaves every gap unchanged. Thus the 2^m counts share the full *outgoing* continuation data (P,A,0) and come from strong tournaments of order 14m+9. A common subsequent transitive extension leaves the same gap digits.

The ratio of successive gap weights exceeds 2, so the sumset is superincreasing. Its exact finite inverse is greedy from the largest weight, with a final zero-remainder check. This is a physical positional register, not a claim that arbitrary integers are reachable.

## Fixed content and primitive-compatible progression

Because gcd(37,704)=1, the gcd of all count differences on a fixed m-cell board is **exactly 3315840**, independent of m. This is a crucial difference from an equal-radix board whose ordinary factor divides its slot radix: there the common divisor can grow exponentially in m.

In fact infinitely many boards avoid every prime divisor of both 37 and 3315840. For primes other than 37 the baseline, divided by P in the indicated field, satisfies

|p|N0/P mod p|
|---:|---:|
|3|1+2m|
|5|2|
|11|4m|
|157|(69+77m)/37|

For example modulo 157 the normalized x is already 1 after initialization and remains 1, since 704+118≡37; each cell adds 1490/37 to N/g and the guard adds 1. Modulo 11, where 704≡0, x stays 2 and the cell increment is 2744/37≡4; the initialized N/g plus the guard is zero. Modulo 5, x stays 1, the pair adds zero, and initialization plus the guard gives 2. Modulo 3 the pair returns x=0, each term involving T vanishes, and the two local N increments sum to 2. Finally modulo 37 one has Q≡1 and the last pre-block A≡822≡8, so the final N is (T_last+Q)·8≡16 or 5. All counts are odd. As 3315840=2^7·3·5·11·157, take m≡2 mod 5181=3·11·157 to avoid all these primes at once. This is a necessary-arithmetic gate passed, not a primality or completeness theorem.

The position geometry is very far from an interval: it deliberately has only two digits for a ratio (704/37)^2. With the specified initializer and guard, the final count is at least the incoming A before the guard, which is at least 704^(2m+1) by the positive recurrence. Hence the number of this bank's distinct counts through X, even allowing every m, is O(X^(log(2)/(2 log(704)))), with exponent about 0.0529. The remaining positive-realization question is how to supply more load choices at a *single* incoming scale while still preserving the exact outgoing response. This suggests a testable search using the ratio Q/h, the number of paintings sharing S,d, and the range of their (T,e,n) values, not only the number of paintings or congruence classes attained.

Run `python verify_onesided.py` for raw subset counts of the spine and both first physical switch boards, all 256 switch choices at small m by exact transfer, and the primitive residue check.
