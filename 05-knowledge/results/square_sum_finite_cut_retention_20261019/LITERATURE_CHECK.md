# Targeted primary-source literature check

Checked2026-10-19. This was a focused search, not an exhaustive novelty review.

- Robert Gerbicz's original2018 Mersenne Forum post is the primary construction source for ordinary interval square-sum Hamiltonicity: https://mersenneforum.org/showthread.php?p=477787 . Direct retrieval was blocked with403, so it was not re-read in full during this check.
- Gerbicz's author post#47 dated2023-10-15 in the same thread explicitly recommends nondecreasing2-opt moves for square-sum sequence pairs and retaining positional parity in the second member. Its text was available in indexed search results at https://www.mersenneforum.org/node/17331/page2 . This is clear prior use of the underlying heuristic; no novelty claim should be made for2-opt discovery.
- Elwyn Berlekamp and Richard Guy, *Fibonacci Plays Billiards*, arXiv:2002.03705, includes an appendix recording Gerbicz's ordinary square-necklace existence result and linking his construction. The full author manuscript was inspected: https://arxiv.org/html/2002.03705 . It does not supply a prescribed-path theorem in the material examined.
- The author's repository https://github.com/devisnik/square-paths describes cycle extension by cutting and reinserting, with a reported computational range up to1000000. This further confirms that incremental cutting/reinsertion is prior practice, not a new algorithmic idea here.
- Lee Trent's2022 senior thesis *Structure of Number Theoretic Graphs* is listed by its institutional repository as covering square-sum degrees, connectedness, Hamiltonicity and planarity: https://scholar.rose-hulman.edu/math_mstr/179/ . Only the repository record was reviewed in this pass, so no claim about absence of prescribed-path results in the full thesis is made.

Targeted queries for prescribed paths, prescribed edges and fixed edges in square-sum interval graphs did not locate a directly applicable general extension theorem. Many results called “square sum labeling” concern a different labeling problem, and results about the graph power G² concern a different graph. Neither should be used to assert an extension theorem here.

The present claims remain mathematical consequences of the explicit certificate construction. Historical priority of the exact prescribed-seed theorem or matching corollary is not established by this limited search.
