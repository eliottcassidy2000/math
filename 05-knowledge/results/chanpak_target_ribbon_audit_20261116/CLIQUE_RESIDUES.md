# No fixed-modulus hole in the one-apex conflict-clique family

**Proved, all moduli.** For every M≥1 and every class modulo M containing an odd integer, there is a strong tournament with a transitive vertex deletion and Hamiltonian path count in that class. Its odd-cycle conflict graph is a clique. Therefore a statement that this family misses half the *integers* cannot be explained by missing classes modulo a fixed modulus. This does not address exact coverage or density.

## Provenance and scope

[Session 25](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md) already proves all fixed moduli for broader strong one-apex tournaments; [Session 26](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md) separates this from the exact density limitations. The new claim is for this much narrower reflected contact family, with one common tape size for all classes, using an explicit two-packet modular compensation.

## Exact conflict-graph dictionary

Use the expansion of Grinberg–Stanley, [The Rédei–Berge symmetric function of a directed graph](https://arxiv.org/abs/2307.05569): H(T) is the sum of 2^j over vertex-disjoint packs of j nontrivial odd directed cycles. Build a simple graph whose vertices are the individual directed odd cycles of length at least three, taken up to cyclic rotation (distinct cyclic orders on the same support are distinct vertices). Two vertices are adjacent iff their cycles share a vertex. Its independence polynomial at 2 is exactly H. This dictionary does **not** allow an arbitrary interval or overlap graph, does not identify cycles by support, and does not discard long cycles. If deletion of one apex leaves a transitive tape, every directed cycle contains that apex, so the conflict graph is complete.

## Program every residue, not every absolute count

Use the [mirrored transitive seam](https://github.com/eliottcassidy2000/math/blob/a872455acffc668b3f577c1fb3941bc70e1cfd50/05-knowledge/results/chanpak_signed_seam_bank_20261116/DYADIC.md). Two TT_n halves carry a reflected binary mark word x. The left and right guards have length one, and the sentinel belongs to the transitive tape. The new head is the only apex. The head count is

    H = 2^(2n+2)+1 − M(x),
    M(x) = x²+[x odd]+2J(x),
    J(x) = Σ_{i>j, both bits set} (2^(i−j−1)+[i=j+1]).

To avoid notational collision write the desired modulus as 2^s m, m odd. It suffices to take s≥1 by adjoining a factor of two if necessary. Pick fixed n large enough for all the following slots, padding with zeros. Leave a gap of at least s+2 between distinct bit packets, including the separate low bit at index 0. All packet indices other than 0 are greater than s+1. Then cross-packet terms vanish modulo 2^(s−1) in J, and x² modulo 2^s is the low bit. Hence M(x)/2 modulo 2^(s−1) is the low bit plus the sum of the within-packet J's.

1. **Dyadic steering.** The low bit is independent. For each j=1,…,s−2 include or omit a two-mark packet at internal distance j+1. Its J contribution is exactly 2^j. These binary controls reach all residues of M(x)/2 modulo 2^(s−1), whatever the subsequent odd stage does. If s=1 this stage has no obligation; if s=2 the low bit alone suffices.
2. **Odd steering.** If m=1 stop. Otherwise choose D≥s+2 with 2^D≡1 mod m. An odd switch uses *two* packets. In one packet, the second mark moves from distance 1 ahead of the first to distance 1+D. In the other, it moves from distance 2+D back to distance 2. Neither crosses any other packet. Modulo m, moving by D preserves x and all pure power-of-two contributions, including cross-packet terms. Exactly one adjacency term disappears in the first packet, with no adjacency change in the second. Thus J drops by 1 modulo m. Dyadically the two changes are −2 and +2 modulo 2^(s−1), since c(1)=c(2)=2 and the long versions vanish. The two packets together preserve the dyadic residue and decrease M(x) by **2 modulo m**, an invertible shift.
3. Place m−1 switches initially off. Choose how many to turn on after steering the dyadic residue. The tape length n, guards and sentinel are fixed throughout. By CRT this programs every permissible residue of H.

The statement is about a fixed chosen modulus and one common (possibly very tall) tape for all its classes. The count size grows exponentially with tape height; nothing in the proof makes the representative a particular prescribed integer. For a power of two it needs only O(s²) tape vertices; the general bound is O(s² + m(s+D)). It illustrates exactly why congruence coverage for conflict cliques and eventual integer coverage are distinct questions. A two-apex family has to offer a better *absolute* selector, not merely fill a presumed dyadic gap.

## Reproduction

The verifier builds the packets, checks the full residue set in the exact arithmetic formula for several coprime moduli, and independently counts literal small tournaments by subset paths. Word shifts are specified directly, without an LP or statistical assumption.
