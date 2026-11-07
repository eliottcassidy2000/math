# Exact finite scope of additional searches

Status: finite computations, independent audit requested. None of the negative searches below is an all-tournament or all-length impossibility theorem.

## Bounded-start source inventory and short periodic macros

The complete deletion bank through order8 gives63546 root-cut presentations of source order<=9 with 1<=s<=26 and s not divisible by3. The separately audited complete191536-class deletion order9 bank gives190262 more source-order10 presentations. Thus ALL_SOURCES_TO10.txt has253808 presentations, not necessarily distinct rooted isomorphism classes.

`all_root_sources.cpp` and `all_sources10.cpp` compile exact root-incidence coefficients; the latter performs11232768 cuts on21939 cores capable of a bounded positive start. These are not restricted by universal-pair reciprocity, since the weaker periodic criterion does not require it.

`periodic_group.cpp` considers every nonempty O/Z/J macro of length<=5 (363 macros) on either side, grouping initial adjusted bases X across residues. It asks whether repeated application of one fixed macro on each side keeps X equal, with at least one macro containing J. Among7822 shared initial bases it checks91406667 source/macro cases. There are ZERO matching keys even after the first THREE exact repeats. Therefore no modular-hash verification is needed; the negative answer is exact. The stored nineteen-term modular signature is only a later filter and plays no role in this zero-three-step result. Intermediate integer sizes for the three repeats fit signed64-bit arithmetic: the nine-state maximum grows by at most4 per letter, and source-order10 cover counts are bounded by10!*binom(12,2); multiplying by4^15 is below2^63.

This excludes fixed one-loop macros of this length on this source bank. It does not exclude an arbitrary finite-state paired grammar, asynchronous branching, other moves, or later sources.

## Left-reset rays and binary tests

The reset inventory includes all root cuts on deletion banks through order8 with 1<=d<=26 and d not divisible by3, even if the old start is0 or a multiple of3. It yields329 paired (X,a) keys and203 eventual ray signatures. No equal five-state binary Q/Z pair occurs in that inventory. Extending every reset by Q/Z words through length4 checks542066 presentations and139949 keys without a full binary match. Mixed O/Z/P/Q padding of the particular X=69 ray through length10 also yields no full binary pair.

## Structured deletion-core counter-rotations

These searches allow arbitrary root incidences on singleton runs around specific strongly connected deletion cores. They are broader than uniform common words. All tests enforce s=e=d and compare (a,rho,b), then the exact adjusted base. Partial matches must not be reported as paired-source witnesses.

1. `counterrotation.cpp`: transitive deletion d=1 versus TT_l + S4 + TT_r with d=5 and B2(S4)=32. Through18 total singleton pads on the S4 side:150994960 presentations,960 partial matches, zero full pair. Matching b determines the transitive length. A smaller defect control through10 pads gives h_A-h_B only -12 or-10, not the required+2. This observed pattern is not proved for all lengths.
2. `core_counterrotation.cpp`: either of the two strong order6 deletion cores with (H,B2)=(19,256), with N pads, versus S4 with N+3 pads. Through N=12:8126336 presentations,88 partial matches, zero pair. Smaller controls through N=8 give h_A-h_B=40, not required-10. Again no all-length claim.
3. `double_core_counterrotation.cpp`: S4 + TT_m + S4, with total N singleton pads including the exterior runs, versus S4 with N+5 pads. Through N=8:2751744 presentations,1268 partial matches, zero pair. Required h_A-h_B=-14.

The rooted ordinal-sum composition used by test3 is independently checked on literal graphs in verify_structural.py. If p=(h,s,e,L,R,d,b), q=(H,S,E,LL,RR,D,B), merge their roots and orient all left nonroots toward all right nonroots. The new state is

(dH+(R-d)LL+(e-d)S+(h-e-R+d)D,
sD,dE,bLL+D(L-b),BR+d(RR-B),dD,bB).

This is the inherited common-root ordered amalgamation specialized to the seven coordinates. Literal validation covers100 products of all labeled sources of orders2..3 and all256 root-cut products of S4.

The deletion-count choices 1,5,19,25 arise naturally by tabulating the complete small-core inventory and matching odd parts of B2. Turning this into an all-order core reduction additionally requires the cited strong-tournament lower bound used in earlier sessions and an explicit exhaustive odd-part table. This packet does NOT claim that unprovided reduction as a standalone theorem. Even with that reduction, arbitrary singleton-run lengths remain unbounded.
