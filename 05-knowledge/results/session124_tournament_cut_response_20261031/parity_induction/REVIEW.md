# Independent review: parity induction from the raw cut law

**Verdict: valid and noncircular.** The proposed addendum supplies a self-contained proof of classical Rédei oddness and the integral cut weights. It makes no novelty claim and does not prove unrestricted Hamiltonian-path spectrum completeness.

## Dependency and quantifier audit

For a nonempty tournament T, deleting the added vertex from a Hamiltonian path gives the endpoint/two-path-cover expansion. For each original arc u→v, concatenation identifies W_uv with Hamiltonian paths using that arc. Summing outgoing minus incoming uses at each vertex gives the endpoint-divergence identity. These are integer-counting bijections and equalities; neither parity nor positivity of half-weights is used. They therefore establish the **raw** identity

`F_T(S)=H(T)+sum_(u in S,v outside S,u→v)(W_uv+W_vu)`

before defining or asserting integrality of b.

Fix n and assume oddness for **every** n-vertex tournament. For any pair u,v of any such T, reversing their arc produces another n-vertex tournament, so the induction hypothesis applies to both. In each cover counted by either W_uv or W_vu, the distinguished vertices belong to different component paths; their mutual edge is never used internally. Both cover counts are therefore unchanged by the reversal. The only Hamiltonian paths lost or gained are those using that arc, and cutting/concatenating gives

`H(T with uv reversed)−H(T)=W_vu−W_uv`.

The left side is even. Consequently W_uv+W_vu is even as well, and b_uv is an integer. The raw cut identity now shows that **every** extension of **every** n-vertex tournament has odd H. Every(n+1)-vertex tournament is such an extension, by deleting any chosen vertex and retaining its actual out-neighbor cut. This closes the induction; it is not restricted to a particular constructor family.

## Boundary cases and positivity

Use H(empty)=1 and H(one vertex)=1. The empty-parent cut law is checked directly, rather than by its empty endpoint expansion. At n=1 there are no pairs or edge reversals to consider; the raw identity directly gives H=1 for both orientations on two vertices. At n=2 the ordinary edge-flip argument applies.

After oddness has been established for all orders, the already-audited free-orientation block-insertion injection gives `b_uv>=H(T−{u,v})`. The right side is at least1 because it is a nonnegative odd path count; when n=2 the empty-remainder convention gives1 directly. Thus strict positivity is downstream of, not an input to, the parity induction.

No existing code or receipt was changed for this review. This is a proof audit of the stated addendum, not a replacement of previously frozen computational evidence.
