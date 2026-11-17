# Unrestricted signed flanks and a target-dependent replacement core

With the definitions of the marked interval law, a fixed center D+v is described by
`(h,P,L,R,L⁻,R⁻,κ)`, where κ=H(D+v). For arbitrary u,w≥0 take a left transitive spine whose v-beats word is the binary expansion of u and a right one whose v-beats word is the reverse complement of the binary expansion of w. Zero denotes omission. The exact law, now including unbalanced adjacent parity channels, is

    H=κ+h(F(u)+F(w)−2)+Ru+Lw+Puw−R⁻[u odd]−L⁻[w odd].

Indeed the signed left-start count is −1 precisely when the final binary bit is 1, and the signed right-end count is −1 precisely when the first right bit is 0. All the long-span terms are unchanged from the marked interval formula. The guessed requirement that both words be even was only a sufficient way to erase the adjacent sidecars. If two cores match the six non-κ coordinates, their exact difference survives **all** pairs of binary flanks. All displayed centers are strong with nontrivial v-contacts, so all constructed graphs are strong.

The former core D=(56,33,35,6,14,24), masks 6 and 9, has common (h,P,L,R,L⁻,R⁻)=(45,432,129,258,−15,−30) and κ=137,139. Removing the parity restriction alone still misses 257,631,8191,65537 and all their corresponding listed sevenfold companions. Changing the core supplies the following selections (masks are outgoing, and the value after the semicolon is the alternate mask if needed).

| rows, apex mask | u | w | H |
|---|---:|---:|---:|
| (4,1,2),1 | 13 | 1 | 257 |
| (4,1,2),1 | 11 | 4 | 631 |
| (4,1,2),1 | 1219 | 0 | 8191 |
| (4,1,2),1 | 65 | 83 | 65537 |
| (4,1,2),1 | 0 | 94342 | 458759 |
| (4,1,2),1 | 1 | 3370 | 57337 |
| (4,1,2),1 | 230 | 0 | 1799 |
| (2,4,9,3),1 | 82 | 1 | 4417 |
| (18,12,17,5,10),3;18 | 0 | 144 | 8191;8189 |
| (18,12,17,5,10),5;24 | 0 | 8064 | 458759;458761 |

The C3 center has profile (3,12,3,6,−1,−2,5); the order-four center's profile is (5,32,6,18,−2,−4,9). The order-five center shares the entire signed boundary `(15,120,42,63,−6,−9)` across all four indicated masks, with κ=43,41,41,43. This exhibits two new exact +2 marked recontact pairs across all binary flanks, including on a previously missed sevenfold prime target. It does not claim that every target admits one of these centers.

The inverse script uses the prior exact low/high binary split. For nontrivial contacts all coefficients A in the fixed-small-side inverse are positive; the newly signed parity corrections are bounded local constants and the (O(sqrt(N)\log N)) bound on candidate evaluations remains valid for a fixed center.

## Two-flank factor-seven extension

When v beats all of D, κ=h, L=P−h,R=0,L⁻=−h,R⁻=0. This follows from the vanishing parity character for optional slots and subtracting the term for the empty marked slot. Thus the whole read is

    H=h(F(u)+F(w)−1)+(P−h)w+Puw+h[w odd].

Take the previously verified paired ordinal blocks
`(h,P)=(2197,558848)` and its exact sevenfold multiple `(15379,3911936)`.
The two-flank reader gives exact strong counts p and 7p for every u≥0 and w>0. Strength follows from the common apex reaching the first block and the last flank reaching the apex, regardless of whether the center is internally strong. This is a replacement of the center, not a claim to retain it or to solve an arbitrary prime. At (u,w)=(1,5),(1,7),(2,6) the primes are 5586283,7812887,10061461 (trial division); their exact paired products are 39103981,54690209,70430227. With two arbitrary free contact words this operation differs from the previous two-run cyclic reader.
