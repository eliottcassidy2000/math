# The expanded primitive bank has the same generated language

Status: PROVED scope conversion, independently checked in `audit/AUDIT.md`. This note makes the audited factorization explicit; it introduces no new target computation.

Let G be the rooted-tournament language generated from the rooted singleton by arbitrary O/Z operations and the four rooted-triangle factors K_j, with the root retained. Then (T,r) belongs to G if and only if every strongly connected component of T-r is a singleton or a directed triangle.

In the forward direction, each factor appends its nonroot deletion component, with every earlier nonroot vertex pointing to every later nonroot vertex. Root deletion therefore leaves precisely the ordered singleton/triangle components.

Conversely, order the SCCs of T-r. For each component take its induced graph together with the same root and the original root incidences. The ordered root-amalgamation of these pieces reconstructs every original arc. A singleton component is O or Z. A triangle component is K_j up to an independent cyclic relabelling, where j is the number of triangle vertices beaten by the root. Thus the factorization is physical and exact, not just equality of scalar counts.

In particular, any primitive core whose deletion is an ordered join of ONE triangle and arbitrarily many singleton SCCs factors into O/Z letters and exactly one K_j. Associativity expands every word in that infinite primitive bank into G. If the bank includes all four K_j and permits arbitrary O/Z letters, the generated languages are equal. A subbank generates a subset, so the same negative target conclusion still holds.

The audited target inverse excludes H=55 and H=67 from all of G, at every total word length. Consequently every tournament with H=55 or H=67, after deletion of ANY chosen vertex, must have a strongly connected component of order at least four. This is not a global nonexistence claim for either target.

The prime palette {3,5,31,41,61} belongs only to the narrower NO-SPACER fixed-triangle family. It must not be imposed on G: the literal O1 graph already supplies H=13. For G and the expanded primitive bank, the audited new conclusion is the core-depth bound and the exact coprime-to-three omissions {7,55,67} through 67. Other deletion SCC types and other gluing operations remain outside the obstruction.
