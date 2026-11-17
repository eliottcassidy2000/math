# The normalized reader fails at the Fermat prime 257, but both witnesses exist

In the [normalized recontact construction](../chanpak_eight_split_20261117/NORMALIZED_SEED.md), a prime target p must have a strong F and tags A giving

    g=H(F),     p=3W+2g,     U+E=g.

Here W counts ordered two-slot covers with second slot starting at a tag; U also requires the first slot to end at a tag; E counts full paths ending at a tag. With a singleton F the only p is 5. With |F|>1, A cannot be a singleton: then U=0 would demand E=g, meaning every path ends at that vertex, impossible in a strong tournament (a Hamilton cycle can be cut at any vertex). So W≥2g, by cutting each full path before each tag. For p=257 this forces g≤32 and g≡1 mod 6. The only remaining possible strong counts are 13,19,25,31: the case g=7 is the known all-order hole. The exact minimum for strong order 8 is 45, and the floor increases thereafter, so F has at most 7 vertices. See [Busch's strong minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).

The attached standalone replay fixes a directed Hamilton frame and varies all other arcs. Any strong isomorphism class can be presented with such a frame. For every eligible tournament it tests **all** tag subsets of size at least two, computing subset path endpoints independently from the construction. Scalar matches and their remaining defects are:

| order | g | W | U+E−g | marked labelled frame presentations |
|---:|---:|---:|---:|---:|
| 5 | 13 | 77 | +40 | 13 |
| 6 | 25 | 69 | −6 | 25 |
| 7 | 25 | 69 | −4 | 25 |

No normalized match survives. These are frame presentations, not isomorphism-class multiplicities. Run `g++ -O3 -std=c++20 -Wall -Wextra -Werror target257.cpp -o target257 && ./target257`.

This is a gap in the construction, not a gap in the spectrum. The strong C₃[TT₈,1,1] realizes 1+2⁸=257 directly. The [moving Fermat reader](../chanpak_changed_contacts_20261117/MOVING_FERMAT_READER.md) gives a literal construction at 7·257=1799. Thus the defect really distinguishes different ways to realize the same arithmetic target.
