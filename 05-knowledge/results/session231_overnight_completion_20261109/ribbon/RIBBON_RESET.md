# A finite-menu digit obstruction for arbitrary ribbon prefixes

Status: PROVED and independently audited ([audit](../ribbon_audit/AUDIT.md)). The positive strong-state reset problem remains OPEN.

THM-4880(H), `01-canon/theorems/THM-4880-ribbon-tournaments-bounded-span-back-arcs-give-a-finite-state-transfer-for-hamiltonian-paths-strong-connectivity-is-regular-the-span-2-series-has-a-rank-14-integer-representation-irreducible-on-nonempty-words-47-needs-span-3.md`, records finite searches for constant/bounded residual digit blocks. This lemma does not rely on their state ranks or finite lengths.

## General spectator lemma
Suppose a construction F on tournaments has a domain closed under prepending ordinal copies of C3, and satisfies

    F(C3 => T) = C3 => F(T).

Then, for any scalar a, the set {H(F(T))-a H(T)} is bounded only if all its members are zero. Indeed a nonzero residual r generates residuals 3^t r.

For a finite menu F_1,...,F_m with the same property, if a is even, no bounded set D can satisfy: for every T in the domain there is j with H(F_j(T))-aH(T) in D. Fix any T. Each initial residual is an odd nonzero integer. On C3^t=>T all residual magnitudes are at least3^t. Thus none lies in D for large t. This remains true if the choice j depends on the entire prefix.

## Ribbon consequence at every width
Fix w>=1. A ribbon word has backwards arcs of span at most w. Let B_1,...,B_m be any fixed finite list of suffix words, of arbitrary finite lengths. On prefixes u of length at least w, appending B_j never changes edges before the last w prefix vertices. Consequently

    T((C3)^t => u B_j) = (C3)^t => T(u B_j).

The C3 components are legitimate span-one ribbon words: reverse both adjacent arcs in the three-vertex transitive order. Every ordinal join of ribbon words is another ribbon word of the same maximum span. When forming this ordinal concatenation, clear the formerly out-of-range bits in the first w letters of the right operand; they had no effect in its standalone encoding, and clearing them prevents unintended cross-boundary back arcs. The spectator lemma applies.

Thus NO finite bounded-length suffix menu implements an even-radix expansion with bounded additive digits on every ribbon prefix. This rules out doubling-plus-digits on the full state space, for all widths and all suffix lengths, not just the finite universe in THM-4880(H). For one block and any radix, a bounded residual image must be exactly {0}.

## Precise positive escape
The lemma does not exclude a digit construction on a proper invariant family of strong prefixes; prepending C3 leaves that domain. It also does not exclude operations reaching arbitrarily far back, such as an apex whose neighborhood depends on the existing prefix. Both escape routes must preserve a full port state, not merely H.

A viable next lemma is therefore: exhibit a strong-prefix state family S and two legal finite words B_0,B_1, possibly state-selected, with H(uB_i)=2H(u)+d_i, bounded d_i of opposite residue classes mod4, and uB_i still in S. Seed coverage and the exact inverse residues must then be checked. Merely finding some affine transitions on finite strong samples does not establish S-invariance or all-target coverage. In particular the global exclusions7 and21 mean no naive all-seed two-digit orbit can be admitted without a precise seed/exception rule.
