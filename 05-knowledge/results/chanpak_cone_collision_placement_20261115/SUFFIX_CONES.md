# A finite suffix cone for all thirteen floor units

Let the units be numbered in the pinned order of [THM-5212(D)](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-5212-the-bottleneck-free-growth-bound-floor-units-compose-by-exact-hankel-laws-and-every-bottleneck-free-unit-word-grows-at-least-1.7699-per-vertex.md):
0.2, 1.2, 0.123, 2.13, 0.1.3, 0.3.1, 0.13.2, 0.1.123, 0.3.12, 0.12.3, 1.2.3, 1.12.3, 0.1.23. A concrete unit word is an initial vertex (ribbon letter 0) followed by those units. A unit word is *bottleneck-free* if its product has rank at least two in the minimal representation of the full thirteen-unit language.

- **All thirteen:** for every bottleneck-free word, $H\ge C\cdot 1.7349^n$.
- **Any mixture of the twelve other than 0.1.123:** the rate improves to $1.7699$.

Here $n$ counts vertices, not units. These theorems also hold in the larger accepted language where every four-unit window has reduced rank at least two (and shorter words are harmless). They do not cover unrestricted span-three ribbons, nor identify primes merely by checking rank. The all-thirteen rate $1.7699$ remains open in this packet.

## Replace the moving image plane by a fixed finite cone for each suffix

The ordinary forward transfer has 119 nonnegative path-fragment states. On adding a vertex one can start a singleton, append or prepend to a fragment, or join two distinct fragments through the new vertex when the necessary arcs exist. Each fragment keeps its head and tail ages 0, 1, 2, or old. An old tail always beats the new vertex; no new vertex beats an old head. A state with two old heads cannot contribute, since it would have at least two heads remaining forever. A final Hamiltonian path induces exactly this sequence of prefix fragments, and every accepted construction is one path. These semantics give an integer matrix $M_j$ for each unit, initial column $v_0$ after the initial vertex, and readout row $h$ selecting one fragment.

Use the Hankel projection $Z$ with rows $h$ and the seven rows $hM_j$ at indices $(0,1,2,3,5,6,7)$. Let $P$ have columns $v_0$ and $M_jv_0$ at indices $(1,0,2,3,4,7,9)$. The verifier reconstructs the whole representation from the fragment rule. It checks $ZP$ invertible and computes $R_j=ZM_jP(ZP)^{-1}$. The exact identities $R_jZ=ZM_j$ hold on every nonempty raw state; the empty state is gone after the initial vertex and never reappears. Thus the eight coordinates work on all states reached from the initial vertex. Reachability follows from the prefix basis, and observability from the suffix basis. This certifies minimality rather than extrapolating from sampled words.

For a three-unit suffix $w=(a,b,c)$ set
$$
K_w=\{ ZM_cM_bM_a x:\ x\ge0\}.
$$
This is one fixed finitely generated cone, using the 119 projected raw coordinate columns. It contains every actual state with suffix $w$. Appending $j$ slides the tag to $t=(b,c,j)$ and maps $K_w$ into $K_t$, since the discarded oldest factor $M_a$ acted nonnegatively. This is why the method works even though tracking *exact* image planes does not terminate.

Use every rank-at-least-two three-unit tag; keep the edge $w\to t$ iff $R_jR_cR_bR_a$ has rank at least two. Any rank-at-least-two whole product follows these edges because the rank of a product cannot exceed the rank of any contiguous factor. No converse is assumed.

## Exact comparisons and growth

The certificate assigns a linear row $\phi_w$ to each tag. On each of the 119 generators check
$$
H/2000\ \le \phi_w(ZM_cM_bM_ax)\ \le 2H,\qquad H=hM_cM_bM_ax
$$
and on each edge check
$$
\phi_t(ZM_jM_cM_bM_ax)\ge \rho^{|j|}\phi_w(ZM_cM_bM_ax)
$$
where $|j|$ is the number of vertices in the unit. Linearity extends all inequalities to $x\ge0$. The candidate functions can be found using an LP, but the files consist of integer rows divided by a common denominator; the verifier checks the scaled integer inequalities and reduced rank conditions without floating-point comparisons. Starting after the first three units, use $H\ge1$, multiply the inequalities along the word, and compare the final readout. The resulting explicit bound is $H\ge \rho^{n-10}/4000$. For shorter words the same bound follows directly from $H\ge1$.

For the twelve-unit rate the replay checks 1,193 tags, 11,134 edges, and 1,324,946 generator inequalities. The all-unit rate checks 1,471 tags, 14,725 edges, and 1,752,275 generator inequalities. The optional weaker all-unit certificate at 1.7101 uses only scalar multiples of $H$ as functionals, and already beats the floor $5^{1/3}$. Run, with Python, NumPy and SymPy, and assertions enabled:
```
python verify_suffix_cone.py CONE12.json
python verify_suffix_cone.py CONE13.json
python verify_suffix_cone.py CONE13_WEAK.json
```

## The omitted unit can be sparse without losing the stronger asymptotic rate

Suppose a bottleneck-free thirteen-unit word uses $b$ copies of 0.1.123. Split around those copies into at most $b+1$ runs of the twelve good units. On an arbitrary reachable prefix, the same twelve-unit cone bound applies from its first three good units: their initial readout is at least the readout before the run by the usual injective one-vertex extensions. Thus a run of $m$ vertices multiplies the running count by at least $1.7699^{m-9}/4000$ (also true for short runs by monotonicity). Each exceptional unit costs three vertices and never lowers the count.

Multiplying the run estimates gives, for $\rho=1.7699$,
$$
H\ge \frac{\rho^{n-1}}{4000\rho^9}\left(4000\rho^{12}\right)^{-b}.
$$
In particular $b=o(n)$ preserves the full asymptotic lower rate. This does not control words with frequent exceptional letters at the same rate; a stronger cone needs to handle their transition costs, rather than assuming they are infrequent.
