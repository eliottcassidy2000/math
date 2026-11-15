# Sharp dyadic escape rails for three missing one-apex counts

A one-apex tournament is painted on a transitive spine: bit 1 means the apex beats the spine site. The earlier [absolute suffix decoder](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md) used Collatz's lesson of retaining the low-bit itinerary and its full address, not just a residue. The earlier [local-to-global failure](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_least_missing_jacobian_20261114/MODULAR_ESCAPE.md) gave arbitrary precision modulo every M, but not a sharp escape rate around the one-apex hole 13. The following fixed control word gives a tight order-of-magnitude escape with no variable modular search.

## Three literal rays

For every integer r≥1, the three strong one-apex tournaments below have the stated numbers of Hamiltonian paths:

| Spine word | Count | Order (including apex) |
|---|---:|---:|
| 11 · 0^r · 1110 | 7 + 27·2^r | r+7 |
| 11 · 0^r · 10110 | 13 + 57·2^r | r+8 |
| 11 · 0^r · 11100 | 21 + 75·2^r | r+8 |

In each case the initial 1 and final 0 give a spanning directed cycle. Thus the construction is literally strong at every depth, although it never realizes the intercept itself. The least achievable absolute count in the residue of each missing one-apex count e mod 2^r is Θ(2^r), and the minimum order is r+O(1): these statements are sharp in scale, not exact equality to the exhibited counts or orders.

## Exact normal form and sharp coefficient

Let P be a fixed prefix with at least one mark 1; let W be an arbitrary fixed suffix. Denote the prefix state by N=H(P), S=ordinary binary integer of P, and d=last bit of P. Here H is the count with the apex attached. The singleton cut rule gives the exact state update when appending a bit b:
```
(N,S,d) -> (N+(1-b)(S+d), 2S+b, b).
```
Write D(P)=S+1-N-d. It starts at 0, stays unchanged under a zero, and grows by S+d under a one. In particular it is nonnegative and even. After appending r≥1 zeros the state is (N+(2^r−1)S+d, 2^r S, 0). In a suffix W of length m put
```
t = Σ_(j=0..m−1) (1-W_j)2^j,    p(t)=H(W)−1 .
```
The recurrence gives p(0)=0, p(2u+1)=p(u), and p(2u)=p(u)+u+(u mod 2), so 0≤p(t)≤t. Appending W adds t times the entering value of S, plus p(t); that follows immediately by iterating the rule with entering d=0. Therefore
```
H(P 0^r W) = 1+p(t)−D(P) + S(P)(t+1)2^r.          (∗)
```
In particular the actual constant term can be an unavailable count. A decomposition which silently substitutes the value at r=0 for the constant term would miss this distinction.

The prefix 11 has S=3,D=2. For the three listed suffixes, t=8,18,24 and p(t)=8,14,22, respectively, so (∗) gives all three rays.

More strongly, 27,57,75 are the **least possible odd coefficient** of the 2^r term at the respective intercepts in every fixed layout P·0^r·W. Strip trailing zeros from P; D stays fixed and the removed powers of two in S are absorbed into the exponent. Strip leading zeros from W; t changes from 2^j(t'+1)−1 to t', and p is unchanged, absorbing another 2^j. These removals normalize to odd S and odd t+1 (if W is effectively empty its t=0 cannot give intercept e>1).

The intercept satisfies e=1+p(t)−D(P)≤t+1. If S=1, D(P)=0 and e would be a direct one-apex count, impossible here. If A=S(t+1) were below the coefficient in the table, S≤A/e leaves only S=3 among odd integers ≥3 in all three cases. That means D(P)=2, so p(t)=e+1. Checking only the even t up to 24 via the two-term recurrence gives the first solutions 8,18,24. This finite part is local to the optimality claim; the rails and the exclusion of larger unsampled parameter choices below these coefficients are exact.

## What is sharp and what is not proved

For completeness, the one-apex holes 7,13,21 can be seen without using the full conjecture. Trim leading zeros and trailing ones (they contribute only singleton ordinal components); any remaining cyclic word starts 1 and ends 0. For active spine length m≥3, the directed odd cycles with those endpoints alone contribute 2^(m−2) to H−1, by choosing an even subset of the m−2 interior vertices. Thus these targets can occur only at active m≤4,5,6 respectively. The replay checks these tiny finite universes; the exclusions 7 and 21 in *all* tournaments are separate inherited theorems. The one-apex 13 hole is not a full-spectrum hole.

Conversely the sum of all possible weighted start-end inversion pairs on m spine sites is 2^m−2 for m≥2: weight 2 at distance one and 2^(d−1) at distance d≥2. This bounds H≤2^m−1. Accordingly, once 2^r>e, any congruent witness for e modulo 2^r has H≥e+2^r and order at least r+2. The displayed rays give the upper bounds.

The occupied dyadic shells explain exactly how these one-apex congruence witnesses escape. They achieve linear order in precision, even for two targets that no tournament can realize. This is not a circumvention of the absolute obstruction. A constructive downward step still has to transform the signed intercept in (∗) into an ordinary count, or use a host family outside this layout. Nor is optimality among *all* modulo-2^r witnesses or all possible growing word shapes claimed.

Run `python verify_dyadic_escape.py` for literal subset-path checks, the small forbidden spectra, and the coefficient controls.

## A compiler for every positive odd intercept

For an arbitrary positive odd e, choose b≥1 with c=1+2^b>e and take W=1·0^b, of count c. Write c−e as a sum of its distinct binary bits 2^j, all j≥1. For each bit take P_j=1·0^j·1, whose defect D(P_j) is exactly 2^j. Construct the literal word
```
P_(j1) · 0^r · P_(j2) · 0^r · … · P_(jq) · 0^r · W.
```
Each variable run with X=2^r has polynomial update
```
(N,S,d) -> (N+(X−1)S+d, XS,0).
```
Thus the output H is a polynomial in X with integer coefficients. At formal X=0, a run resets S and d to zero. If a fixed P_j begins with S=d=0, at X=0 the composite P_j and its run change N by −D(P_j). The fixed final W adds c−1. Starting from N=1, the constant term is c−Σ 2^j=e. Every actual word (r≥1) remains strong: it starts with 1 and ends with 0. Its length is qr+O_e(1), with q=popcount(c−e); the endpoint-cycle lower bound and total inversion weight upper bound above imply a polynomial of degree exactly q, with a positive leading coefficient.

The compiler prescribes not just H mod 2^r but an ordinary binary address S mod 2^r: it is always S(W) there, because the final reset multiplies the entering S by X and each suffix bit does not undo divisibility by X. When e≠c, the limiting pair (H,S)=(e,S(W)) has no finite literal realization at that *exact* address, since the ordinary binary address fixes the bit word up to neutral initial zeros. This is a precise failure of taking a p-adic limit inside the family, now for a joint count/address observable. The distinction is the same one needed in a FRACTRAN-style simulation: a formal zeroing of a register permits a transition that no physical positive execution attains. The proof gives an explicit congruence compiler, not a proof of absolute coverage.
