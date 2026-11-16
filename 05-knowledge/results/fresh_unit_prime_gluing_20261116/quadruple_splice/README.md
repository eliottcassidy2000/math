# A synchronized fresh +4 and the obstruction at deletion count 3

**Construction.** There is an unbounded family with strong induced extensions of counts $3H$ and $3H+4$, with both extensions tripling all three rooted continuation coordinates $(E,O,g)$. This is complementary to the [previous +2 splice](../triple_splice/README.md), not a third setting for its interface.

**Obstruction.** For any fixed continuation $(E,O,3)$ with $E\ne O$, no three consecutive odd counts can share that continuation, even allowing all orders, unrelated graphs and arbitrary choices of root. For the new +4 family the entire count fiber is explicit. See [DECODER.md](DECODER.md).

## Exact +4 splice

Masks below are outgoing neighbors with vertex 0 the low bit. Keep root 5 and inflate vertex 4. Both children preserve the source on the first six vertices and add one new vertex.

| | Outgoing row masks | $(H,E,O,g)$ at base size |
|---|---|---|
| $C$ | 30,60,24,16,32,13 | (15,11,12,1) |
| $A$ | 94,124,24,80,32,77,20 | (45,33,36,3) |
| $B$ | 94,60,88,16,32,13,58 | (49,33,36,3) |

$E,O$ count ordered path covers, with the first path starting at the root and its number of vertices even or odd. The second path may be empty. $g$ is the path count after deleting the root. Substituting any nonempty tournament $M$ into site 4 gives, with $u=H(M)$ and $R=P_M(2)$ the two optional ordered path-slot count,
$$
P(C)=(11R-7u,11R/2,11R/2+u,u),\qquad
P(A)=3P(C),\quad P(B)=3P(C)+(4u,0,0,0).
$$
A transitive block of length $s\ge1$ gives $u=1,R=2^s$. Thus the exact source count is $q=11\cdot2^s-7$. For $s=2$ the counts are $37\to111,115$, with common continuation $(66,69,3)$.

### Finite symbolic certificate

A nonterminal run of the marked module has to exit through vertex 5, so there are at most three runs across two paths. Exhaust the legal words made of the fixed vertices and marked runs. For Hamilton paths the coefficient pairs for one and two runs are $(15,11),(45,33),(49,33)$. Root deletion has only one run, with coefficients 1, 3 and 3. For the rooted ordered two-path covers, only patterns $(0,1),(1,0),(1,1)$ survive. In the source their unsigned coefficients are $(12,11,11)$, and their signed coefficients, using parity of the fixed vertices in the first path, are $(0,1,1)$. Both children have triple these tables.

With $\epsilon=(-1)^{|M|}$, the three unsigned filling counts are $u,u,R-2u$ and the corresponding signed filling counts are $u,\epsilon u,-(1+\epsilon)u$. For the third entry, include the empty-side terms in the signed sum and assign every cycle or singleton of a positive odd-cycle packing to a side. Moving one such part reverses its sign and preserves its weight. The full signed sum vanishes, leaving $-(1+\epsilon)u$ when the two empty terms are removed. This is the same [alternating cancellation as in the +2 splice](../triple_splice/README.md#arbitrary-modules-exact-digit-cost). Consequently $E+O=11R+u$ and $E-O=-u$ for the source, and both quantities triple in the children.

Run **python3 check.py** for the run-table certificate and independent rooted subset DP, including nontransitive modules. Run **python3 find_g3.py** for independent controls on the signed-offset decoder.

## The prior second-step exclusion is also insensitive to ordinal padding

The terminal-component lemma in [DECODER.md](DECODER.md) closes the padding and root-choice caveat in the [earlier carry note](../triple_splice/CARRY.md). The second-step nonzero targets 155, 157 and 161 are coprime to the required $g=27$ and have required $E-O=-27$. Any matching terminal strong component therefore has the entire rooted profile and retains the old strong source. The earlier literal test on strong hosts rules out those exact continuation matches at all orders, including arbitrary ordinal padding and any choice of root. The counts themselves are attainable, just not with the required continuation.
