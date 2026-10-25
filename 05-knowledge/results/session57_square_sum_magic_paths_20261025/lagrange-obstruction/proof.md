# The Lagrange square-sum family has no nonzero rational doubly balanced member

## Scope and source

This statement concerns precisely the two-parameter family in equation (1.1) of Ajai Choudhry, *Two sets of integers such that all elements of the sumset of the two sets are perfect squares*, arXiv:2508.07806v1 (11 August 2025), https://arxiv.org/html/2508.07806v1#S1.E1 . Choudhry attributes this previously published family to Jean Lagrange through Guy's book, p. 267. It is not a claim about all of Choudhry's families, all square-sum K4,4 configurations, or the existence of magic squares of squares. The modular exclusion below is an independently derived exact certificate, not a theorem attributed to the paper.

The four first-part labels are

A1=9(18u²+50uv+37v²)²,
A2=36(u²+13uv+16v²)²,
A3=36(9u²+53uv+64v²)²,
A4=4(3u²−25uv−48v²)².

The second part is

B1=0,
B2=20(u−4v)(11u+19v)(u+4v)(u+v),
B3=36(3u+11v)(7u+12v)(3u+4v)(u+v),
B4=756(u−v)(3u+2v)(u+3v)(u+2v).

## Theorem

For rational u,v, the first quadruple admits a partition into two pairs of equal sum if and only if v=0. For any (u,v) other than (0,0), the two quadruples cannot both admit such a partition. Consequently, any decorated locus requiring both quadruples to satisfy parallelogram balance is missed by every nonzero rational parameter value of this family, even allowing arbitrary permutations within each part, exchange of parts, translation gauge, and nonzero square homothety.

## Proof

The three pair partitions of the first quadruple give exactly

A1+A2−A3−A4=−v F12(u,v),
A1+A3−A2−A4=F13(u,v),
A1+A4−A2−A3=−v F14(u,v),

where

F12=16608u³+102220u²v+205548uv²+135135v³,
F13=5760u⁴+50208u³v+168500u²v²+252948uv³+141345v⁴,
F14=19680u³+113996u²v+216300uv²+135135v³.

Modulo 47 these are respectively

17u³+42u²v+17uv²+10v³,
26u⁴+12u³v+5u²v²+41uv³+16v⁴,
34u³+21u²v+6uv²+10v³.

The table below lists their evaluations at every [t:1] in P¹(F47). All entries are nonzero; at [1:0] the evaluations are 17,26,34, also nonzero. Thus none of the three homogeneous forms has a projective zero over F47.

For completeness, this excludes rational projective zeros, not merely integer parameters in a finite range: if a nonzero rational pair (u,v) annihilated one of these homogeneous integer forms, clear denominators and divide the resulting integer pair by its gcd, obtaining coprime (U,V). Homogeneity preserves its zero. Since U,V are coprime, they are not both divisible by 47. Their reductions define a projective zero over F47, contradicting the certificate. If V is nonzero modulo 47 normalize to [UV^(-1):1]; otherwise use [1:0]. This covers denominators divisible by 47 as well.

It follows that for nonzero (u,v) the second balance can never vanish, while either of the other balances can vanish only when v=0. Conversely v=0 gives A=(2916,36,2916,36)u⁴, for which the first and third balances vanish. The origin also has v=0 and all labels zero.

The nonzero boundary v=0 has B=(0,220,2268,2268)u⁴. Its three balances are (−4316,−220,−220)u⁴, none zero when u≠0. Therefore no nonzero rational member is doubly balanced, including the boundary where A has repeated labels. At (0,0), both parts consist of four zero labels; this is the trivial degenerate solution, not a nondegenerate square-sum K4,4 realization. A zero homothety likewise collapses everything and is excluded.

A permutation merely selects one of the three pair partitions, already exhausted. Adding one common constant to all labels of a part cancels from every balance. A nonzero square homothety multiplies every balance by the same nonzero factor. Exchange of the bipartitions cannot repair the failure of the requirement on both parts. Thus all claimed equivalences preserve the obstruction.

## Complete finite-field certificate

| t | F12(t,1) | F13(t,1) | F14(t,1) |
|---:|---:|---:|---:|
| 0 | 10 | 16 | 10 |
| 1 | 39 | 6 | 24 |
| 2 | 19 | 19 | 2 |
| 3 | 5 | 29 | 7 |
| 4 | 5 | 23 | 8 |
| 5 | 27 | 1 | 21 |
| 6 | 32 | 23 | 15 |
| 7 | 28 | 21 | 6 |
| 8 | 23 | 34 | 10 |
| 9 | 25 | 20 | 43 |
| 10 | 42 | 44 | 27 |
| 11 | 35 | 43 | 25 |
| 12 | 12 | 14 | 6 |
| 13 | 28 | 14 | 33 |
| 14 | 44 | 19 | 28 |
| 15 | 21 | 18 | 7 |
| 16 | 14 | 13 | 33 |
| 17 | 31 | 19 | 28 |
| 18 | 33 | 17 | 8 |
| 19 | 28 | 1 | 36 |
| 20 | 24 | 25 | 34 |
| 21 | 29 | 15 | 18 |
| 22 | 4 | 4 | 4 |
| 23 | 4 | 38 | 8 |
| 24 | 37 | 35 | 46 |
| 25 | 17 | 20 | 40 |
| 26 | 46 | 31 | 6 |
| 27 | 38 | 25 | 7 |
| 28 | 1 | 19 | 12 |
| 29 | 37 | 43 | 37 |
| 30 | 13 | 46 | 4 |
| 31 | 31 | 37 | 23 |
| 32 | 5 | 38 | 16 |
| 33 | 37 | 37 | 46 |
| 34 | 41 | 35 | 35 |
| 35 | 25 | 46 | 46 |
| 36 | 44 | 3 | 1 |
| 37 | 12 | 40 | 10 |
| 38 | 31 | 22 | 42 |
| 39 | 15 | 15 | 19 |
| 40 | 19 | 4 | 4 |
| 41 | 4 | 34 | 13 |
| 42 | 25 | 22 | 15 |
| 43 | 43 | 39 | 26 |
| 44 | 19 | 28 | 15 |
| 45 | 8 | 39 | 45 |
| 46 | 18 | 41 | 38 |
| infinity | 17 | 26 | 34 |

Reproduce with `python verify.py` or `python -O verify.py`. Every mathematical check raises an explicit exception on failure; no verification relies on removable Python assertions. `certificate.json` contains all 144 projective evaluations and additional small-prime controls. `source-pin.json` pins the retrieved source bytes by SHA-256.
