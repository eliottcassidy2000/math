# The last forbidden pair protects a sink in the high block

**PROVED elementary grammar lemma.** This is a structural guard for the separate antisymmetric-star proof, not by itself an observer theorem.

Take the **last** occurrence of `be` or `ee` in a legal word. At its birth, use the strong bridge cut constructed in §5 of `FORBIDDEN_WORDS_ARE_STRONG_BRIDGE_CUTS.md`. Name the high triangle Q={r,u,c}, with

    r→u→c→r.

Here c is the newest vertex (created by the final e), u is the vertex created by the preceding b or e, and r is that preceding insertion's parent. The maximum pair at this moment is u→c. In particular c is the universal sink of Q−{r}={u,c}.

**Lemma.** Throughout the suffix after this last forbidden occurrence:

1. every insertion parent lies in the evolving high block Q;
2. c is never selected as parent;
3. c remains the universal sink of Q−{r}.

**Proof.** The persistent nontrivial strong bridge cut gives (1), by the indegree-one maximum/predecessor argument proved in the structural equivalence note.

If the suffix is empty, (2)–(3) are immediate. Otherwise its first letter is c or d: the state is B and e would create a later `ee`, contrary to the choice of the last forbidden occurrence. Neither c nor d selects the vertex c as parent: their parents are r and u, respectively. Thus c gains a second incoming edge and can never again have indegree one.

After this first safe insertion, every maximum is either the retained old vertex u (possible after d) or a later inserted vertex. The old vertex c loses to u. Suppose for contradiction that some later operation is the first to select c. Before this operation, c loses to every later inserted vertex, since each new vertex beats all old vertices except its chosen parent and c has not been chosen. Hence c loses to every current maximum. It cannot itself be a maximum (its indegree is already at least two); nor can it be the cycle predecessor of a maximum, since being a predecessor would require an arc from c to that maximum. These are the only legal parent forms. This contradiction proves (2).

Every newly inserted vertex lies in Q and, because its parent is not c, beats c. Initially u→c. Thus c remains a universal sink of Q−{r}, proving (3). □

The argument actually permits arbitrary later legal letters once the first c or d after the birth has occurred. Choosing the last forbidden pair is the clean way to ensure that first safe letter exists when the suffix is nonempty.

**Boundary.** A direct e immediately after the birth does select c, breaks this particular protected sink, and creates a new forbidden pair with a new protected sink. Therefore this lemma does not claim that one frozen star survives every suffix from an arbitrarily chosen earlier forbidden occurrence.
