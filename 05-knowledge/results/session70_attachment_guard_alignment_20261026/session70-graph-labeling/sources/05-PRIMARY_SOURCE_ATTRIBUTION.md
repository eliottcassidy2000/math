# Focused primary-source attribution check

Checked 2026-10-19. Scope: the condition rho(A_i) >= s-i+1, independent rainbow selections joined by single replacements and closure-preserving double replacements, and its upper-tail graceful-tree corollary.

## Conclusion

**No exact prior statement or proof was identified in the sources checked. This is not a novelty certification.** The closest direct rainbow-spanning-tree announcement studies a stronger color-connectivity hypothesis with a stricter one-edge move model. The accessible announcement does not expose its theorem or induction. No checked primary source explicitly states the arbitrary fixed upper-tail integer-difference corollary. The safe publication wording remains: “We give a self-contained proof; no historical-priority claim is made.”

## Maezawa–Yamaguchi: directly relevant, but an abstract is not the full theorem

Shun-ichi Maezawa and Yutaro Yamaguchi, *On Reconfiguration Graphs Induced by Rainbow Spanning Trees*, ASIACOMB 2026 abstract book, printed page 160 (PDF page index 159), August 28, 2026:

https://asiacomb.com/documents/asiacomb2026-abstract-book-long.pdf#page=160

This official abstract defines adjacency by deleting one tree edge and inserting one edge. It explicitly restricts attention to edge-colored multigraphs in which **each individual color induces a connected spanning subgraph**. It announces investigation of connectivity, diameter, and circumference, without stating numerical bounds or a formal connectivity theorem. Thus its advertised hypotheses are not the descending rank condition: for n vertices, every color has rank n-1. Its moves are single-edge moves rather than simultaneous changes in up to two colors. The abstract contains no closure-preserving swap induction and no graceful-labeling application. This establishes related prior research, not equivalence or non-equivalence of the unpublished full results.

Earlier official lecture announcement, October 9, 2024:

https://www.en.sdu.edu.cn/info/1158/8226.htm

That abstract announces a sufficient connectivity condition but does not state it. Its linked Chinese announcement failed to load. The 2026 abstract above is more informative and should be preferred for comparison.

The authors' accessible publication lists were checked; no matching full manuscript was found there:

- https://smaezawa.org/papers
- https://ygussany.github.io/index.html

Absence from these lists does not establish that no manuscript exists.

A direct caution, independent of the cited sources: start with a path on r+1 vertices, r>=2, and put one parallel edge of every one of r colors on each path link. Every color class is a spanning path. Rainbow spanning trees are permutations assigning the r colors to the r links. No nontrivial single-edge flip is possible: a new edge must keep the deleted edge's color, and it must occupy the deleted link to remain a tree. Therefore every state is isolated under single-edge flips. The abstract's mention of investigating connectivity must not be rewritten as a universal single-edge connectivity theorem under its advertised hypothesis. Our simultaneous two-color moves can exchange the colors occupying two links.

## Wdowinski: precise rank theorem, different adjacency and excess hypothesis

Ronen Wdowinski, *Hall's theorem for reconfigurations and higher dimensional topological connectedness*, arXiv:2511.04863, Section 5.1, Theorem 5.3:

https://arxiv.org/html/2511.04863v1#S5.SS1

The theorem assumes r_M(X)+r_N(V-X) >= k+1 whenever V-X is an N-flat of rank at most k-1. Its vertices are common independent k-sets; adjacency requires their union to be M-independent of size k+1. For a partition matroid N with s nonempty color classes and k=s, this becomes the sufficient condition r_M(union_(i in I) A_i) >= |I|+1 for every nonempty I. This specialization is my direct derivation from the stated theorem. It does not follow from our staircase: rank of all classes may be exactly s. In particular, its union-independent adjacency has no edges between distinct bases of a rank-s matroid. Corollary 5.4 concerns sets smaller than a maximum common independent set. These statements therefore do not directly give the claimed full-size, two-color theorem. No exact staircase or graceful upper-tail statement was located in this paper.

## Lubiw–Pathak: labels move with the token

Anna Lubiw and Vinayak Pathak, *Reconfiguring Ordered Bases of a Matroid*, arXiv:1612.00958, introduction and Theorems 3–4:

https://arxiv.org/html/1612.00958

Their labels belong to selected basis positions and are transferred to the newly inserted element in an exchange. The replacement element is not restricted to a permanently fixed color class. They characterize reachability by matroid components and give O(r^1.5) exchanges, with an O(r log r) graphic bound. Consequently this paper cannot simply be cited as proving connectivity of rainbow bases in an already edge-colored graph. No staircase restriction or upper-tail graceful corollary is stated in the inspected text.

## Jamison–Kinnersley: group-sum colors and a spare color

Robert E. Jamison and William B. Kinnersley, *Rainbow spanning trees in abelian groups*, Journal of Algebraic Combinatorics 56 (2022), 5–21, DOI 10.1007/s10801-022-01120-7:

https://link.springer.com/article/10.1007/s10801-022-01120-7

The authors' university-hosted manuscript is indexed at:

https://math.uri.edu/~billk/pubs/abelian_rainbow.pdf

Search retrieval exposed its Section 4, including Theorem 4.1: for non-Boolean finite abelian groups, rainbow spanning trees can be changed into the star at zero by pivots. Edges carry group-sum colors, and a pivot inserts an edge of the currently unused color and removes an edge from the resulting cycle. These are not fixed integer absolute-difference colors with an unchanged upper tail. Direct PDF retrieval timed out, so the comparison is limited to the indexed primary manuscript excerpt and accessible publisher record; the full proof was not inspected.

## Additional scope limits

Hojin Chu's author-hosted talks page lists a forthcoming October 29–31, 2026 KMS talk with the same rainbow-spanning-tree reconfiguration title. Neither its abstract nor a manuscript was accessible through the checked links, so it supports no mathematical attribution claim:

https://www.hojinchu.com/talks/

Searches also covered rainbow bases, common-independent-set reconfiguration, two-color/two-colour exchanges, closure-preserving exchanges, graceful-tree reconfiguration and edge exchanges, and the exact upper-tail/staircase terminology. Search absence is weak evidence, especially for unpublished work or different terminology. A historical-priority claim would require a fuller bibliographic investigation and access to the Maezawa–Yamaguchi manuscript or detailed talk material. No author was contacted.
