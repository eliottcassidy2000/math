# THM-095: Arc-Flip Contraction Reduction

**Status:** PROVED
**Proved by:** kind-pasteur-2026-03-07-S39b
**Verified computationally:** n=3-6 exhaustive (all tournaments, all edges), n=7 sampled (2000 tournaments, 3 edges each), 0 failures
**Scope:** All tournaments

---

## Statement

Let T be a tournament on n vertices with arc e = (u -> v). Let T' be the tournament obtained by flipping e (replacing u -> v with v -> u). Define:

- **T/e**: contraction of e (merge u,v into w; w inherits IN-edges from u, OUT-edges from v)
- **T'/e'**: contraction of the flipped arc e' = (v -> u) (merge v,u into w; w inherits IN-edges from v, OUT-edges from u)

Then:

**H(T) - H(T') = H(T/e) - H(T'/e')**

This reduces the arc-flip H-difference to a problem on (n-1)-vertex digraphs.

---

## Proof

By THM-082 (deletion-contraction for Hamiltonian paths):

- H(T) = H(T \ e) + H(T/e)
- H(T') = H(T' \ e') + H(T'/e')

where T \ e (resp. T' \ e') is the digraph obtained by removing the arc.

**Key observation:** T \ e = T' \ e' as digraphs.

Both T and T' agree on all arcs except the pair {u,v}. Deleting the one arc between u and v from either yields the same digraph (all arcs of T except the u-v arc, with no arc between u and v).

Therefore:

H(T) - H(T') = [H(T \ e) + H(T/e)] - [H(T' \ e') + H(T'/e')]
             = H(T/e) - H(T'/e')

since H(T \ e) = H(T' \ e'). □

---

## Auxiliary Results

### THM-095a: Contraction Tournament Criterion

T/e (contraction of arc u -> v in tournament T) is itself a tournament if and only if u and v have identical out-neighborhoods on V \ {u,v}:

**T/e is a tournament iff T[u][x] = T[v][x] for all x in V \ {u,v}**

**Proof:** In T/e, the merged vertex w has: (w,x) iff T[v][x], and (x,w) iff T[x][u] = 1 - T[u][x]. For T/e to be a tournament, we need exactly one of {(w,x), (x,w)} for each x. That is:

T[v][x] + (1 - T[u][x]) = 1, i.e., T[v][x] = T[u][x].

Verified exhaustively at n=4,5 (0 errors). □

### THM-095b: Contraction Tournament Probability

The probability that T/e is a tournament (over uniform random tournament T and uniform random arc e) is (1/2)^{n-2}, achieved when all n-2 other vertices are beaten by both u and v, or both beat them.

**Computed values:**
- n=3: 50.0% = (1/2)^1
- n=4: 25.0% = (1/2)^2
- n=5: 12.5% = (1/2)^3
- n=6: 6.25% = (1/2)^4
- n=7: ~3.3% ≈ (1/2)^5

---

## Significance

1. **Inductive structure:** The arc-flip difference H(T) - H(T') is determined entirely by the contraction level. This could enable inductive proofs about H-value differences.

2. **Even H-difference:** Since H(T) - H(T') = H(T/e) - H(T'/e'), and H is always odd (Rédei), the difference is always even. This is consistent with the known fact that arc-flips change H by an even amount.

3. **Connection to OCF:** By OCF, H(T) = I(Omega(T), 2). The arc-flip contraction reduction relates I(Omega(T), 2) - I(Omega(T'), 2) to a (n-1)-vertex problem. This could give a structural explanation of how cycle structure changes under arc flips.

4. **Non-iterability:** The reduction cannot always be iterated because T/e is not always a tournament. It is a tournament only when u,v have identical out-neighborhoods, which happens with probability (1/2)^{n-2}.

---

## Scripts

- `04-computation/arc_flip_contraction.py`: Exhaustive verification n=3-6, structure analysis
- `04-computation/dc_tree_analysis.py`: Full DC tree analysis
