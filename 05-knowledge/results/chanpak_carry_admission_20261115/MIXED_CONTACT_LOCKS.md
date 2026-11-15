# Exact mixed-contact windows beyond the homogeneous tariff

The [homogeneous-reader tariff](RICH_CORE_ADMISSION.md) forbids a cheap switch when lower profile observations are fixed. Instead change both the core and the contacts, and keep the one **mixed** response that the final suffix actually reads. This gives genuine consecutive runs of odd values, rather than subtracting shifts from unrelated windows.

## The one response needed by a painted terminal suffix

Take a spine tournament B, not necessarily transitive, and an apex q. Paint a vertex 1 if q beats it. Let h=H(B), N=H(B+q), and let S be the sum over nonempty directed paths of B starting at a painted vertex of the Hamilton count of their complements. This is the S of the [four-coordinate cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md). Append s≥1 ordinal singleton vertices, all painted 1, then a last singleton painted 0. Each 1 is a sink in the full tournament at its insertion: N is unchanged, S→2S+h, and signed start bias resets to h. The final 0 increases N by the new total start weight plus its signed odd-minus-even bias: this is twice the weight of odd-size eligible starts, equivalently the crossing odd-cycle weight at the new apex seam. Therefore the literal extended tournament satisfies
```
N_s = N + 2^s(S+h).
```
All starting cores below are strong and have mixed cuts. The enlarged graph is strong: from the old strong part one reaches the final vertex, which reaches the apex, and the old strong part is then reachable; any appended singleton reaches the final vertex and is reached from the old part.

This tells us what to seek backwards from a prescribed arithmetic lane. If the target is (odd leading slope)·2^j+(fixed odd checksum), admit any physical rows and painted mask for which S+h has that odd part and N is the desired checksum. Normalize powers of two by shifting s. Equal H or Q_2 across two arms is not required; neither a homogeneous reader nor an exact neutral move on a *fixed* core can substitute for this change of state.

## A genuinely new-profile triple at the same absolute slope

Masks use outgoing apex bits, indexed from zero. All three S+h equal 75.

| Rows of B | Mask | h | Q_2 | S | N |
|---|---:|---:|---:|---:|---:|
| (24,1,3,6,14) | 28 | 11 | 100 | 64 | 41 |
| (18,8,3,5,14) | 21 | 13 | 104 | 62 | 43 |
| (60,45,8,48,6,20) | 12 | 23 | 272 | 52 | 45 |

Thus every s≥1 yields all three consecutive odd values
```
75·2^s + {41,43,45}.
```
In particular at s=5 all three cross the nominal high place with one actual positive carry: their true quotient is 76 and low remainders are 9,11,13. For s≥6 they are carry-free. At an absolute target, for each of the three checksums test whether the excess divided by 75 is a power of two of exponent ≥1; the successful checksum selects the actual mask and core.

The third spine is strong and has two disjoint directed triangles, so its feedback vertex number is at least two. More strongly (H,Q_2)=(23,272) is **not** a profile of any ordinal sum of one-apex cores. The same bounded SCC-factor proof as in [RICH_CORE_ADMISSION.md](RICH_CORE_ADMISSION.md) certifies this: to contribute a factor at Q_2≤272, each strong one-apex SCC has order≤7 from Q_2≥2·3^(n−3). Enumerate all its opening/closing words and use multiplicativity. This is not just a flag requiring new geometry for one possibly old profile.

## Eight and nine consecutive odds at every high place

The same response also has an eight-value consecutive lane. The third column is S+h; the checksum is the marked-base count N. If the slope factor is 2^a, use s=j−a after choosing the arm.

| Rows of B | Mask | S+h | N |
|---|---:|---:|---:|
| (16,1,3,7,14) | 27 | 68 | 23 |
| (16,1,3,7,14) | 18 | 34 | 25 |
| (16,1,3,7,14) | 5 | 34 | 27 |
| (2,40,43,49,39,1) | 1 | 34 | 29 |
| (16,9,3,5,14) | 2 | 34 | 31 |
| (16,1,3,7,14) | 8 | 34 | 33 |
| (2,40,43,49,39,1) | 39 | 136 | 35 |
| (16,9,3,5,14) | 7 | 68 | 37 |

Hence 17·2^j+{23,25,…,37} is realized for every j≥4. The fourth and seventh spines are strong with two disjoint triangles, but their (H,Q_2) profile is old; this is a narrower richness statement than the slope-75 example.

Another lane uses the following seven-vertex strong spines, eight of them containing two disjoint triangles:

| Rows of B | Mask | S+h | N |
|---|---:|---:|---:|
| (114,8,99,5,110,74,10) | 58 | 510 | 235 |
| (16,93,9,33,12,23,61) | 83 | 510 | 237 |
| (48,9,3,53,70,22,47) | 32 | 255 | 239 |
| (8,61,89,80,1,93,19) | 96 | 255 | 241 |
| (106,88,115,36,41,2,56) | 16 | 255 | 243 |
| (40,57,51,52,97,64,15) | 8 | 255 | 245 |
| (10,16,35,6,45,75,31) | 29 | 510 | 247 |
| (68,9,2,53,71,23,46) | 8 | 255 | 249 |
| (30,64,42,2,14,91,29) | 32 | 255 | 251 |

Thus 255·2^j+{235,237,…,251} is realized for every j≥2. In each case the target-first test is the corresponding odd slope and an exact dyadic quotient after subtracting one of the listed offsets. The finite offsets in these lanes do not yield a growing interval at each height or guarantee any low pattern. In particular these examples do not solve universal placement.

Replay: `python replay_mixed_contact.py` independently computes each weighted start sum, literal base and extended Hamilton paths, and old-alphabet exclusion at the Q_2 cutoff.
