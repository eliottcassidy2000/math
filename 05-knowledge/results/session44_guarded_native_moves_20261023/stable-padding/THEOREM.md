# Endpoint-fixed partner transports survive arbitrary transitive padding

Status: PROVED all-size closure lemma and a native family at every order n>=8. This is a theorem about supplied actual observer packets, not a complete wiring classification.

## 1. The closure lemma

Let T be a tournament with two actual Hamiltonian paths P,Q. Suppose an isomorphism g:K_P(T)->K_Q(T) fixes their common first owner and their common last owner. Let T' be the ordinal sum

    R -> T -> S,

where R and S are transitive tournaments of arbitrary orders r,s>=0, and their orientations point from earlier blocks to later blocks. Prepend the unique directed order of R and append that of S to both paths, giving P',Q'. Then the extension of g fixing all new owners is an isomorphism

    K_(P')(T') -> K_(Q')(T').

Proof: within the old core this is exactly g's given property. Within each new transitive block, both partners keep the same consecutive path arcs and reverse the same other arcs. Between two distinct blocks all arcs are reversed except the actual Hamiltonian junctions. If R is nonempty, its last owner points to the first core owner on the retained path; that core owner is fixed by g. The corresponding last-core-to-first-S junction is likewise fixed. All remaining cross-block pairs are uniformly reversed, hence preserved. If a padding block is empty, its junction condition is absent. This is a direct literal-arc proof.

For one-sided prefix padding, fixing the common first core owner suffices; for suffix padding only, fixing the last suffices. The two-sided version retains both boundaries.

Every Hamiltonian path of T' consists of the unique R order, a Hamiltonian path of T, and the unique S order: a path cannot return to an earlier block. Therefore H(T')=H(T). The added source owners have distinct scores above every core score; the added sink owners have distinct scores below every core score. They are individually fixed by every automorphism, and restriction gives Aut(T')=Aut(T).

This closure concerns a partner isomorphism and its endpoint stabilizer. Equality of unpointed target classes alone would not justify extending a particular g across the new junctions.

## 2. A genuine uniformly guarded coupled move at every order n>=8

Use the audited original source1238, represented by canonical ordinary code29558785, with the two slots

    A0=(3,0,6,7,1,5,4,2), A1=(3,1,6,7,0,5,4,2),
    B0=(2,0,6,7,1,5,4,3), B1=(2,1,6,7,0,5,4,3).

The actual map g=(4,0,2,3,1,5,6,7) sends each State0 complement partner to its State1 partner, and fixes both owners2 and3. Thus the closure lemma applies simultaneously to both slots for ANY r,s>=0. Choose r+s=n-8.

Retain the actual new prefix and suffix, both old core endpoint owners in each slot, and equality of the six core interior owner variables between the slots. The two states remain compatible. They preserve literal shared arcs(6,7) and(5,4), now at clocks r+2 and r+5. The only changing shared owner variables are the old owners0 and1 at clocks r+1 and r+4.

Changing one slot alone is impossible: the unchanged other slot fixes the entire shared interior, while the prefix, suffix and core endpoints are pinned. These coordinates specify every position of the changed slot. This argument holds even if a padded target-class fibre acquires additional candidates. The known simultaneous two-slot change succeeds. Its sharp changed-slot budget remains two; its changed slot-clock count is four and its changed shared-variable count is two, independently of n.

Each slot preserves its actual target isomorphism class by the extended map g. Hence it preserves EVERY owner-class potential difference at the new order, in particular every actual mod8 harmonic character, with the source fixed. This does not assume or compute a numeric electrical atlas for n>8. The original numeric phase values5 and4 are not claimed to extend under padding. The statement is equality between the two states for whatever authenticated potential exists at that new order.

The original source has trivial automorphism group, so every padded source does too. The move therefore continues to cross distinct actual source-automorphism orbits. It is not explained by an invisible kernel rotation.

## 3. Converse, ports and global scope

Converse with path reversal exchanges prefix and suffix blocks and transports each actual clock j to n-1-j. The extended partner isomorphism and fixed arcs transport with them. Nothing is repinned to an unchanged printed clock. The selected packet still needs unused new literal ports or an explicitly injective global reassignment before it can extend into a fixed global remainder.

This supplies an all-size positive family under the strongest per-slot target-dependent guards. It does not show that every guarded fibre is connected, that these moves generate all native wiring, or that a sparse observer determines all electrical guards.

## 4. Exact tests and provenance

The original code and four paths come from the audited coupled packet and original owner atlas, not a synthetic tournament. The checker independently reconstructs its arcs and enumerates all40320 core permutations to verify H=409 and Aut=1. It checks direct target transports and actual pinned path states in every pair0<=r,s<=8 (81 hosts, orders8 through24), plus the one-sided32-owner paddings of order40. It also checks converse squares and that the two crossed slot combinations fail the shared-interior equality.

The proof establishes all r,s, not the finite range. Run check.py normally and with -O. The source pin document records the exact prior artifacts; no optimization-sensitive assertions or numeric target-class guesses are used.
