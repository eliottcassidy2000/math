# An exact crossing at exponent 18 and two different pairs at exponent 21

**FINITE-EXACT.** There are two strong tournaments on 32 vertices differing in a single arc with Hamilton path counts **34,359,607,295** and **34,359,607,297**. Their midpoint is \(E_{18}=2^{17}(2^{18}-1)\). At \(E_{21}\) there are at least four different strong lower isomorphism classes on 40 vertices with a single-arc exit to the upper count. Their four upper classes also differ. Both exponents are composite. At exponent 61 the corresponding existence question remains open.

Use a core \(D\) with outgoing row masks in vertex order starting at zero. The apex beats precisely the vertices in the tag mask. The left transitive flank takes the bits of \(x\) from high to low and the right flank the complemented bits of \(y\) from low to high. A one means the apex beats the flank vertex. The left flank, core and right flank form an ordinal sum after deleting the apex. Each pair uses identical flanks and contacts and reverses only the indicated core arc.

| Exponent | Lower core \(D\) | Arc | Tags | \(x\) | \(y\) |
|---:|---|---|---:|---:|---:|
| 18 | 80, 21, 73, 131, 44, 15, 58, 119 | (4,6) | 89 | 7 | 749445 |
| 21 | 4, 17, 34, 7, 13, 27 | (2,3) | 22 | 7 | 923014801 |

The exact [two-flank transfer](../chanpak_gap_star_ports_20261122/README.md) uses the core path count \(h\), ordered two-slot count \(q\), and even/odd masses for tagged starts \((s_0,s_1)\) and untagged ends \((e_0,e_1)\). With \(k\) the apex-only count, its linear and constant terms are
\[
 (e_0+e_1)x+(s_0+s_1)y-(e_0-e_1)[x\ {\rm odd}]-(s_0-s_1)[y\ {\rm odd}]+k .
\]
The other terms \(h(f(x)+f(y))+qxy\) agree across each of these flips.

At exponent 18 the two profiles share \(h=413,q=6200,(e_0,e_1)=(1266,1527)\). The start pair changes from \((1042,1253)\) to \((1043,1252)\) and \(k\) rises from 1575 to 1579. Thus the gap is \(4-2[y\ {\rm odd}]\) for **every** choice of both flanks. The indicated addresses give the stated absolute heights. At exponent 21 the gap is \(-2+2[x\ {\rm odd}]+2[y\ {\rm odd}]\). The indicated pair has counts 2,199,022,206,975 and 2,199,022,206,977.

Each endpoint of the new exponent-21 pair has exactly one vertex of outdegree one. Each endpoint of the earlier 40-vertex pair has two. Thus the count pair and even its common order forget which marked tournament gave the crossing. Reversing every arc preserves both Hamilton counts and the specified edit. For the two constructions and their converses the numbers of vertices of outdegrees 1, 2 and 3 are respectively (1,2,0), (1,1,1), (2,1,0) and (0,0,5), reading each original directly followed by its converse. These four signatures distinguish both upper and lower classes at order 40. At exponent 18 taking the converse also gives a second lower and a second upper class, distinguished by two versus one vertices of outdegree one. At exponent 15 the corresponding distinction is one versus two. This proves lower bounds of two, two and four on the numbers of marked edges at the three exhibited orders. The first possible order and the total number of edges at order 40 are separate questions. [Earlier exponent-21 certificate](../chanpak_middle_out_square_sum_20261122/EXACT-SEAMS.md)

Run the linked Python script from this directory. It uses subset path DP on the two small cores to derive both profiles, evaluates the exact transfer, and constructs all full orientations to check strong connectivity, the unique edit and the degree distinction. [Certificate replay](other_seams.py)


## Rigidity throughout the consecutive band

The exponent-15 Farey construction closes the other gap. Taken together, there are literal single-arc crossings at every integer exponent from 3 through 21. The exponents giving Mersenne primes in this range are 3, 5, 7, 13, 17 and 19, meaning those exponents give even perfect centers. [Earlier exact exponent bank](../chanpak_mirimanoff_unit_lifts_20261122/README.md) · [Exponents 11 and 21](../chanpak_middle_out_square_sum_20261122/EXACT-SEAMS.md) · [Exponent 15](../chanpak_decimal_kernel_20261122/E15-FAREY.md)

**PROVED using finite integer input.** For exponents **3, 4, 6, 9, 10, 12, 13, 15, 16, 18, 19 and 21**, every tournament at the upper count, of any order, has trivial automorphism group. In particular the four upper classes at exponent 21 contribute at least four times 40! labelled realizations at that order. For the other exponents in the band the following are upper bounds on the set of possible automorphism primes.

| Exponents | Primes still allowed |
|---|---|
| 5, 8, 11, 20 | 7 |
| 7 | 11 |
| 14, 17 | 7, 11 |

The automorphism group acts freely on paths, so each of its prime divisors divides the height. The increasing odd-prime orbit lower bound from triangle packings has value 32,973,313,212,691 at 37, already exceeding the largest upper height in the band. It suffices to try the odd primes through 31 for divisibility and the orbit bound at each prime. This table also explains how the quotient by relabellings differs from the phase modulo 24. At odd exponents at least five, the endpoints always reduce to 15 and 17 modulo 24, and yet the orbit restrictions above depend on the exponent. The order-nine switching class already exhibits unit gaps in every even-center phase. Its relationship to marking a Monster involution and to the rank-24 Leech lattice is developed in the earlier switching study. [Orbit bound and Fermat/Platonic constraints](../chanpak_fermat_solids_seams_20261122/README.md) · [Switching phases, Monster and Baby Monster](../chanpak_mod24_moonshine_ports_20261122/README.md) · [Six Leech labels and their glue](../chanpak_leech_perfect_phases_20261122/README.md)
