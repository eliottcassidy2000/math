# An exact orbit collapse for transitive middle cores, and three gap-two crossings

This note concerns the physical binary-cut gluing from the [absolute host decoder](../chanpak_heisenberg_target_20261115/ANCHOR_DECODER.md); it does not classify the full Chan–Pak spectrum. It keeps the absolute high and low pieces, not just remainders.

## 1. The all-transitive chamber has precisely one checksum over each possible high coefficient

Take the binary host spine word `1 0^g v` with (g≥1), (v) of width (k≥0). Let (b) be its binary value, (ρ=rev_k(b)) its reversed value, and (A(v)) the Hamilton count of its one-apex tournament (with (A(∅)=1)). The left core (U=TT_m), (m≥1), contributes (x=2^m-1), and the middle core (W) contributes its exact ((h,Q)). The gap formula is
```
N=2^gD+E,qquad D=xQ·2^k+h(2^k-ρ),qquad E=xQ(b+1)+h(A(v)+1-x).
```
Here (Q≥2h), so (E>0); choosing (2^g>E) makes ((D,E)) the actual quotient and remainder. The formula follows directly by substituting (B=2^{g+k}+b+1) and (A=2^g(2^k-ρ)+A(v)) in the binary-host transfer.

Suppose both cores are transitive, so (h=1, Q=2^s, s≥1). Write (j=s+k), (L=m+j), (M_m=2^L-2^{L-m}). Then
```
D=M_m+ℓ,qquad ℓ=2^k-ρ∈[1,2^k].
```
For given (m,j), allowing (k=0,…,j-1) gives exactly (M_m<D≤ M_m+2^{j-1}=M_{m+1}). In particular the bit length of (D) is (L). The intervals for (m=1,…,L-1) tile **precisely the nonpowers of two** at length (L≥2); each (D) fixes (m,j,ℓ) uniquely, including the exact interval endpoints.

For fixed (D), changing (k) to (k+1) replaces (ρ=2^k-ℓ) by (2^k+ρ), so (b) changes to (2b+1). This appends a terminal `1` to (v). Appending it contributes no new one-apex inversions, hence does not change (A(v)); equivalently that vertex can be slid from the physical transitive middle block to the host suffix without changing any arc. The factor (C=xQ=x2^{j-k}) is halved whilst (b+1) is doubled. Consequently (E=C(b+1)+A(v)+1-x) is unchanged in every representation. There is a canonical one with (s=1, k=j-1).

This gives a literal no-search inverse and obstruction for this lane: a high coefficient (D≥3) occurs iff it is not a power of two, and at most one exact low remainder is possible with that high coefficient, even after varying **all three** source orders. At any sufficiently large (g) this unique pair is realized. Also (E(2D)=E(D)): in the canonical representation of (2D), (m) is unchanged, (j,k) increase by one, and (ρ) doubles. Its reversal gains one initial zero and the value (b), and hence (A(v)), are unchanged. Varying the gap cannot cure this absence of independent checksum choices.

## 2. Three physical adjacent-odd crossings, using a cyclic second core

A cyclic middle core can make a genuine second checksum at the **same** high coefficient, differing by the minimal nonzero odd-count gap. Bits in the row records are zero-based outgoing adjacency. All checksums in a row are already positive and less than (2^g) for the indicated (g).

| Common (D) | Lower arm (m,(h,Q),v) | Higher arm (m,(h,Q),v) | Exact low pair | Sufficient (g) |
|---|---|---|---|---:|
| 1281 | (1,(1,2),111111110) | (1,(81,1200),∅) | 1279, 1281 | (g≥11) |
| 6241 | (2,(1,16),1111100) | (2,(73,2056),∅) | 6093, 6095 | (g≥13) |
| 16393 | (1,(1,1024),1110) | (4,(73,1088),∅) | 15369, 15371 | (g≥14) |

Respective cyclic middle row records:

- ((82,108,9,1,78,29,44)), (h=81,Q=1200).
- ((54,188,128,5,12,28,191,57)), (h=73,Q=2056).
- ((98,28,121,1,9,90,26)), (h=73,Q=1088).

The first and third middle cores are strong; the second has a source, and the integer profiles are independently checked via literal subset/end-point DP. The identity for each row itself holds for every (g≥1), even before the displayed low/high separation takes effect.

The first row has a particularly simple target-first form: it hits the two neighboring odd targets
```
1281(2^g+1)-2,  1281(2^g+1).
```
For a target (N) test whether either (N/1281-1) or ((N+2)/1281-1) is a power of two at least (2^{11}); if so that gives the actual (g) and the correct arm. This is exact selection at an absolute count, not just a coincidence of residues.

The same factor-minus-two mechanism works at any (k≥2) **if** a physical middle profile exists with (Q+h=5· 2^{k-1}+1). The lower arm uses (m=1,Q=2,h=1,v=1^{k-1}0), for which (b=2^k-2,ρ=2^{k-1}-1,A(v)=2^{k-1}+1); its leading coefficient is (5· 2^{k-1}+1) and its checksum is two lower. The higher arm uses the promised profile and the empty suffix. The condition (Q+h=5·2^{k-1}+1) is a real profile requirement, not guaranteed for arbitrary (k). It is certified at (k=9) above; no infinite certification in (k) is claimed.

## 3. Why even the cleanest crossing does not re-center itself on arbitrary targets

For the (D=1281) branch pair, replace the specific prefix `1 0^g` by an arbitrary nonempty prefix (p) of binary value (a≥1). Set (P=A(p)-1). Neither suffix generates a boundary adjacency correction: the nonempty suffix starts with 1 and the other is empty. The exact seam formula now gives
```
N_- =1279+1280a+P,   N_+=1281+1200a+81P.
```
Thus (N_+-N_-=2+80(P-a)). The pair is adjacent iff (P=a). Writing (P=P(a)) from the host recurrence, (P(2u)=P(u)+u+(u mod 2)), (P(2u+1)=P(u)), so (d(a)=a-P(a)≥0) from the binary-word inversion bound (A(p)≤ a+1), and
```
d(2u+1)=u+1+d(u),   d(2u)=d(u)-(u mod 2).
```
A positive odd input has (d(a)≥1), with equality only at (a=1). Removing trailing zeros lowers the deficit by just one at the first doubling and does not lower it on subsequent doublings. Consequently (P(a)=a) iff (a=2^g, g≥1). The two exact branches therefore reproduce precisely the claimed dyadic strip of adjacencies; arbitrary prefixes do **not** turn it into an interval. This statement is about the two fixed identified branches and does not rule out other cyclic witnesses for general (N).

## Replay

Run `python replay_checksum_crossing.py`. It checks the three profiles and literal small-gap tournament counts using subset/end-point DP, audits the exact quotient/low parameters, and checks the transitive chamber and prefix-rigidity laws over bounded data. The all-order conclusions above come from the symbolic formulas, not from this finite sample. The second and third rows are independent examples, not claimed minimal.
