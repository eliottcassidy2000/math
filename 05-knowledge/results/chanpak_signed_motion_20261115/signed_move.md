# Bounded signed differences by moving, not freezing, physical head counts

All windows below are literal strong tournaments at every specified height, not modular windows or a target-coverage claim.

## Readout and conventions

Take an arbitrary marked tournament `D,V` between nonempty transitive guards of lengths a,b. Append the sentinel t and the control x as in the [two-register head construction](https://github.com/eliottcassidy2000/math/blob/6d38e8ceb0ccc8a061270e52507820bd4ecf53e6/05-knowledge/results/chanpak_computational_tapes_20261115/HEAD.md): all tape vertices beat t, t beats x, and x beats all tape vertices except V. Denote the result by W(D,V;a,b). If (h,q,E,L,J) are, respectively, the Hamiltonian path count, the optional two-slot cover count, first-slot end in V, second-slot start in V, and both conditions, the head formula gives
```text
H(W) = q·2^(a+b) + E·2^a − L·2^b + h − J.
```
The transitive guards have H=1 and Q₂=2^length. Both guards nonempty implies W is strong: x reaches the first guard and hence all the tape, all tape vertices reach x through t, and x reaches t via that guard.

The core code stores a 1 when smaller beats larger, in lexicographic edge order. The mask stores membership by vertex index starting at 0. In tables, guard offsets (u,v) mean actual lengths (s+u,s+v), except in the explicitly one-sided window.

## Three repairs on one core of unchanged order

Use the single five-vertex core of code 4, with h=5 and q=64. The two *different* high counts E and L move on each step, and only their signed aggregate stays fixed.

| Mask | u | v | E | L | J | E·2^u−L·2^v | h−J |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 | 0 | 2 | 40 | 22 | 12 | −48 | −7 |
| 13 | 0 | 2 | 36 | 21 | 10 | −48 | −5 |
| 20 | 1 | 1 | 17 | 41 | 8 | −48 | −3 |
| 24 | 1 | 1 | 15 | 39 | 6 | −48 | −1 |

Hence for every s≥1 the same order 2s+9 realizes four consecutive odds
```text
256·4^s − 48·2^s − 7, −5, −3, −1  (relative to the common high part).
```
At s=1 they are 921,923,925,927. The active tournament and both control vertices are retained throughout. The middle step can be implemented on the same vertex set: move the first right-guard vertex to the end of the left guard, crossing exactly the five vertices of D. Equivalently fix TT_s and TT_{s+1} as the outer guards throughout, and put an unmarked pivot p in the same six-vertex repair region as D. It is after D in the first two states and before D in the last two. All other transitive guard relations are unchanged. In ascending output order, the three transitions require respectively 2, 8 and 2 arc reversals (the middle uses five guard/core reversals and three mask/control reversals). This gives a uniformly bounded local move script, not merely separate witnesses. It does not claim that the intermediate one-edge reversals each keep the count in the window.

Even if both guard attachments are to be frozen, two successive repairs are possible. Use the four-vertex core code 4, h=5, q=32; the common guards are (s,s+1).

| Mask | E | L | J | E−2L | h−J |
|---:|---:|---:|---:|---:|---:|
| 9 | 12 | 12 | 4 | −12 | 1 |
| 10 | 14 | 13 | 2 | −12 | 3 |
| 4 | 6 | 9 | 0 | −12 | 5 |

This gives `64·4^s−12·2^s+{1,3,5}`, starting at 233,235,237. Both guard lengths and all core arcs are fixed, and edits are restricted to arcs incident to x.

## Let even h and q move

The following is another four-count window at every s≥1. The middle steps change h,q,E,L all at once; the two high powers cancel separately by adjusting one guard.

| Core order | Code | Mask | u | v | (h,q,E,L,J) | h−J |
|---:|---:|---:|---:|---:|---|---:|
| 6 | 32759 | 7 | 1 | 2 | (9,168,35,112,20) | −11 |
| 6 | 32757 | 7 | 1 | 2 | (9,168,31,110,18) | −9 |
| 6 | 32622 | 25 | 0 | 2 | (33,336,150,132,40) | −7 |
| 6 | 32750 | 6 | 1 | 2 | (9,168,29,109,14) | −5 |

Every row has `q·2^(u+v)=1344` and `E·2^u−L·2^v=−378`, giving `1344·4^s−378·2^s+{-11,-9,-7,-5}`. The middle core has h=33 rather than 9 and removes a guard vertex; this window is not asserted to stay at a fixed order or use a fixed core.

A second, distinct mechanism lets the leading q and E move *against one another* rather than demanding that the two powers cancel separately. Make actual guards (s+u,1):

| Code (six vertices) | Mask | u | (h,q,E,L,J) | (2q+E)·2^u | h−2L−J |
|---:|---:|---:|---|---:|---:|
| 32615 | 48 | 0 | (37,376,138,96,20) | 890 | −175 |
| 32694 | 10 | 1 | (13,208,29,88,10) | 890 | −173 |
| 32759 | 37 | 1 | (9,168,109,71,38) | 890 | −171 |

Thus all five coordinates move at **each** repair, yet the total increases by precisely two: `890·2^s+{-175,-173,-171}`. Between the second and third rows the guard roles and lengths are unchanged and `2Δq+ΔE=−80+80=0`; opposite changes in the two growing contributions exactly cancel.

## Why eventual small differences with two growing guards lie on rays

Fix a finite bank of nonempty marked core states, and fix a gap bound K. Write C=h−J. Consider any pair of their outputs with guards a,b and a',b'. If all four are at least N and the gap is at most K, then for N large enough the high components are exactly equal. Indeed the difference of the high parts is divisible by 2^N, but its absolute value is at most K+|C−C'|, uniformly bounded over the bank.

Put Q=max q over the bank. If the leading values q·2^(a+b),q'·2^(a'+b') differ, their absolute difference is at least their maximum divided by 2Q². To see this, either one is at most half the other, or the difference of the outer exponents has magnitude giving a power of 2 at most 2Q. In the latter case their nonzero difference is at least the smaller *power of two*, and the maximum is no more than 2Q² times this power. Meanwhile each of the four E,L terms is at most its corresponding leading value times 2^(−N), since E,L≤q. Exact equality of high components therefore forces leading equality for large N, and then the signed linear terms agree.

Hence the sum displacement k=(a'+b')−(a+b) is fixed by the two positive q values (and is available only when their ratio is an integral power of two). For d=a'−a the remaining equation is
```text
(E−E'·2^d)·2^(a−b) = L−L'·2^(k−d).
```
E,L,E',L' are positive for nonempty masks. As d→+∞ the signs of the two outside coefficients are opposite; likewise as d→−∞. Thus only finitely many d are admissible. For each, either both coefficients vanish, giving independent synchronized guard growth on that pair, or a−b is fixed, giving a single diagonal ray (or no integer guard solution). On every such pattern the signed difference is precisely C'−C.

This classification is only for a fixed bank with all four guards growing. In particular the growing q and E can cancel each other along a one-sided ray because one guard remains short, as the explicit 890 family shows. Neither result extends to a bank that grows with the target. Any fixed-width window on a finite union of geometric rays remains sparse; the next issue is to choose and compose cancellations across genuinely changing cores without assuming that coarser modular coverage gives the desired absolute count.

## Reproduction

`python signed_move_check.py` independently enumerates all slot profiles of the stated masks with subset-DP, then explicitly builds and counts every full tournament at heights 1,2,3. The algebra from the verified profiles proves the infinite extensions. This is an existence check for each displayed witness, not an exhaustive census.

## Two active nontransitive blocks: cancel the amplification or bypass it

For adjacent marked blocks D₁,D₂, with nonempty transitive guards of lengths a,b, put T_i and S_i for the fractions' unnormalized numerators: the number of Hamiltonian paths terminating or starting in the mark of D_i. Specializing the exact multihead formula gives
```text
H = q₁q₂·2^(a+b) +(h₂E₁+q₁E₂)·2^a −(q₂L₁+h₁L₂)·2^b
    +h₁h₂−h₂J₁−h₁J₂−T₁S₂−L₁E₂.
```
The last two terms are genuinely joint: respectively the predecessor at x lies in D₁ and its successor in D₂, or vice versa. An individual J-only edit still has the multiplier from the other block; two edit choices can balance those multipliers, or the joint term can directly move the small remainder even when the multipliers share a factor.

### Different cores: balance the amplified local steps

| Order | Code | Mask | h | q | E | L | J | T | S |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 1464080 | 25 | 89 | 1240 | 574 | 436 | 166 | 48 | 32 |
| 7 | 1464080 | 35 | 89 | 1240 | 574 | 436 | 172 | 48 | 32 |
| 8 | 45442625 | 57 | 237 | 4080 | 1950 | 1950 | 770 | 122 | 122 |
| 8 | 45442625 | 101 | 237 | 4080 | 1950 | 1950 | 754 | 122 | 122 |

The two changes separately give −1422 and +1424. Taking both gives +2, identically in both independent guard lengths. The defining identity is 89·16−237·6=2. Thus amplified corrections can still be useful by being paired as a signed difference; neither correction is intrinsically a unit. The pair's count is
```text
5059200·2^(a+b)+2554038·2^a−1952430·2^b+{-942835,-942833}.
```
Complete enumeration of masks on the displayed cores shows they each have exactly two masks in the respective fixed (E,L,T,S) shelf. Moreover, across **all** pairs of nonempty masks on these two cores, grouping by the two signed high coefficients gives no three-count step-two run.

### Same core twice: move the interaction through the common factor

Take two adjacent copies of the five-vertex core code 955, with h=13 and q=104. The transition of ordered masks is (21,21)→(7,28). Here the one-point counts are **not** frozen in either block:

| Mask | E | L | J | T | S |
|---:|---:|---:|---:|---:|---:|
| 21 | 58 | 58 | 26 | 8 | 8 |
| 7 | 50 | 59 | 22 | 7 | 8 |
| 28 | 59 | 50 | 22 | 8 | 7 |

Both signed high coefficients are 6786 before and after. The three changing constant contributions are
```text
Δ[-13J₁−13J₂] = +104,
Δ[-T₁S₂]      =  +15,
Δ[-L₁E₂]      = −117.
```
Their sum equals +2. Since h₁=h₂=13, no linear combination of purely independent amplified repairs could produce that difference at all. The jointly exposed terms are what break the common 13 divisor, rather than a choice of coprime local scales. The tape has the same order, both guards and all core arcs fixed, and only control contacts change. Its counts are
```text
10816·2^(a+b)+6786(2^a−2^b)+{-3935,−3933}.
```
At (a,b)=(1,1) they are 39329,39331. The full mask census has no three-count run in a fixed high pair here either. The remaining target is a *bank* of interacting moves whose constant interactions produce more than a single unit at one high pair.

The same-core transition can be encapsulated as a single ten-vertex marked compound core, the ordinal sum of the two five-vertex pieces. Its whole profile keeps `(h,q,E,L,T,S)=(169,10816,6786,6786,104,104)` and changes only J: 4104→4102. The endpoint counts on the full core are T=h₁T₂ and S=h₂S₁, so both are retained. This synthesizes a continuation-ready unit bank from two pieces whose individual flat changes would each be multiplied by 13. It does not let arbitrarily many H=13 pieces share one control: the one-head contact bound still takes a factor 13 from every unexposed piece, and each summand exposes at most two. Three or more pieces with 13 dividing both h and q would force all output differences to retain 13.

A concise design rule behind the shared-factor example: suppose D has an anti-isomorphism that sends a mask W to W*, swapping E with L and T with S and preserving J. Start both copies at a balanced mask with E=L=e, T=S=t and J=j. If W has profile
```text
(E,L,T,S,J)=(e−(q/h)δ,e+δ,t−ε,t,j+η),
```
then (V,V)→(W,W*) preserves the two complete high coefficients. The signed change at the bottom is
```text
ΔH = −2hη + 2tε − ε² − 2eδ − δ².
```
Indeed the T₁S₂ correction sees both t−ε endpoints, and L₁E₂ sees both e+δ contact counts. This is an integer synthesis condition; no extra divisibility is imposed beyond integrality of its physical inputs. The code 955 example realizes δ=ε=1 and η=−4 (the numerical profile of mask 28 is the converse of mask 7). The measured two-site change travels through the difference of squares rather than through separately scaled local units.
