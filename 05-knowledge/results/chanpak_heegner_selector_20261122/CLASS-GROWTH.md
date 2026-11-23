# One arc can change the degree of the class-counting polynomial

**PROVED.** For each of the five CM integers \(J_d\), consider the number \(a_n(h)\) of ordinary tournament isomorphism classes on \(n\) vertices with the fixed path count \(h\). It is eventually a polynomial in \(n\). Its degrees at the two endpoints of our physical flip are different in **every row**.

| \(d\) | Prime factors of \(J_d-1\) | Degree | Prime factors of \(J_d+1\) | Degree |
|---:|---|---:|---|---:|
|11|\(23\cdot31\cdot47\)|3|\(3\cdot11171\)|2|
|19|\(7^2\cdot17\cdot1063\)|3|\(941^2\)|2|
|43|\(11\cdot787\cdot102199\)|3|\(5\cdot23\cdot31\cdot41\cdot6053\)|5|
|67|\(147197952743\)|1|\(5\cdot47\cdot59\cdot10616513\)|4|
|163|\(3677\cdot4787\cdot91433\cdot163129\)|4|\(5\cdot103\cdot509781383768483\)|3|

At the separate rounded midpoint 33,506 the degrees are 2 and 5, from \(5\cdot6701\) and \(3^3\cdot17\cdot73\). Thus the 67 flip moves between an eventually linear and an eventually quartic isomorphism count even though the two marked tournaments on 31 vertices differ only at one arc.

**Proof.** The full size series at fixed \(h>1\) is the sum over ordered factorizations \(h=d_1\cdots d_s\) of the product of the finite strong-class polynomials \(A_{d_i}(z)\) divided by \((1-z)^{s+1}\). This follows by decomposing into strong components and inserting independently many singleton components in the \(s+1\) gaps. Positivity at \(z=1\) makes the degree equal to the largest attainable number \(s\) of nontrivial strong components. Every attained prime count occurs on a strong core, by the same factorization. All prime factors below \(2^{32}\) here except 7 occur by the inherited spectrum. The value 49 also occurs, and stripping its exterior singletons gives a strong core because a nontrivial split would require count 7. Thus in the only row containing sevens the attainable maximum is three by grouping 49. Four would require a separate strong component of count 7. [General class series](../chanpak_arc_fiber_grammar_20261122/README.md) · [Finite spectrum through the stated bound](../../../session203/coverage-recertification/README.md)

Two larger prime factors need their own witnesses. The new Fibonacci crossing supplies 147197952743 directly. For 509781383768483 choose the two one-apex words with bridge parameters \(b=29, e=12179923056399\). The literal two-apex bridge has order 51 and height
\[
  be+(1+f(b-1))(1+f(e-1))-[b,e\text{ both odd}]=509781383768483.
\]
This closes the outstanding factor in the concurrent ordinal-product construction for the Ramanujan upper neighbor. Trial division and the literal bridge formula are reproduced in the script. Together with the general prime-factor upper bound, this proves the table. [Bridge orientation and count](../chanpak_prime_frontier_selector_20261118/targeted/PROOF.md) · [Concurrent scalar Heegner constructions](../chanpak_heegner_absolute_targets_20261123/HEEGNER-TARGETS.md) · [Replay](class_growth.py)

**A coefficient bound for the Ramanujan lower endpoint.** The four strong factor cores from the concurrent construction have orders 10, 10, 12 and 12. Every one differs from its converse: the multiset of pairs consisting of a vertex score and the sorted scores of its outneighbors distinguishes all four pairs. The four prime heights distinguish the factors. Their \(4!\) orders and \(2^4\) choices of orientation therefore give 384 distinct classes already on 44 vertices. Insert exterior and intermediate singleton strings and add the separate strong 58-vertex witness. Coefficientwise,
\[
  \sum_{n\ge0}a_n(J_{163}-1)z^n
  \ \ge\ \frac{384z^{44}}{(1-z)^5}+\frac{z^{58}}{(1-z)^2}.
\]
The four-component and one-component words have different strong decompositions at every order. The physical flip at order 58 proves adjacency on top of the 44-vertex scalar constructions. Run `python3 -B class_growth.py` in this directory.
