# Nine consecutive exact flip centers at the perfect number with exponent 17

**PROVED by explicit constructions.** Write $E_p=2^{p-1}(2^p-1)$. For every even $j$ from $0$ through $16$, there are two strong tournaments differing in one arc whose Hamiltonian path counts are $E_{17}+j-1$ and $E_{17}+j+1$. Thus ten consecutive odd counts from **8,589,869,055 through 8,589,869,073** are represented, with a physical gap-two flip across every intermediate even center. The tournament pair may depend on $j$. Their verified orders in increasing center order are $35,33,31,36,33,34,35,34,33$.

The preselected windows were $[E_p-16,E_p+16]$ for $p=17,19$. These are the certified center offsets in those windows.

| Perfect number | Center offsets with strong one-arc witnesses |
|---|---|
| $E_{17}=8,589,869,056$ | $-8,0,2,4,6,8,10,12,14,16$ |
| $E_{19}=137,438,691,328$ | $-14,-12,-6,-4,0,2,6,8,12,16$ |

The second row gives thirteen consecutive represented odd values from $E_{19}-7$ through $E_{19}+17$. It certifies an actual flip only across the listed centers. An odd height shared by two constructions can have different realizers.

## How the middle of the band was filled

Fix two deletion cores related by one arc flip and use the same apex tags. With a left binary address $x$ fixed and the parity of the right address $y$ prescribed, the exact count is
$$
 K(x,y\bmod2)+h p(y)+(s+qx)y.
$$
Here $h$ counts paths in the deletion, $s$ is its tagged start mass, $q$ is the two-slot mass, and $K$ is obtained from the full tagged endpoint profile. Thus matching $h$ and $s+qx$ on both sides and making the constants differ by two preserves the physical flip at every right address of the prescribed parity. A fixed opposite arm lets the target select a different slope. At center $E_{17}+4$ the addresses $(x,y)=(17,79709)$ fill the gap between the preexisting one-arm witnesses. At $E_{17}+16$ the search uses a longer left address. [Exact two-arm count formula](README.md) · [Earlier perfect-center witnesses](PERFECT-TARGETS.md)

The [certificate script](check_perfect_band.py) lists the complete outgoing masks of each small lower core, its marked arc and tags, and both addresses. On each side it computes the small profile by subset dynamic programming and obtains the stated count from the exact identity. It also builds the two labelled hosts. A common source or sink is forced to the same end in every Hamiltonian path; removing it preserves both counts. The script performs these removals with common labels, then checks that each remaining tournament is strong and their symmetric adjacency difference is precisely one unordered pair. These are positive certificates with no dependency on the discovery bank.

Run `python3 05-knowledge/results/chanpak_gap_star_ports_20261122/check_perfect_band.py` from the repository root to print all twenty centers, the common order of each pair and both exact counts.
