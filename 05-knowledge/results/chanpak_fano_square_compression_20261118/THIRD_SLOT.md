# The singleton separator is essential to this exact Paley rail

For the same two cores, let q_k count ordered decompositions into k nonempty directed paths. Their first four values are

    X:   27,  234,  900,  1800
    P:  189, 1638, 6678, 15120.

An independent enumerator can count q_k by concatenating the k slots into one permutation: every backwards adjacency in the permutation must be a slot cut; if there are b backwards adjacencies, exactly binom(n−1−b,k−1−b) cut sets work. Thus the first unscaled coordinate is q_3: 6678−7·900=378.

Replace the singleton separator b in C→A→b→C by a transitive block of width two, and also take A to have width two. C can now have three runs. Those need two separating runs in each of the other blocks, so the only word with three C slots is C A B C A B C. Neither A nor B has a vertex left for another run. Consequently at a=2, exactly one such word occurs, with coefficient (2²−2)²=4, and all words with one or two C slots still scale by seven. The exact excess at a=2, b=2 is therefore 1512. This is a literal continuation obstruction, not a statement against other width-two replacement cores.

The three-slot mismatch explains why proving that the primitive core profiles scale in two slots is sufficient for the singleton wrapper but not for arbitrary compositional reuse.
