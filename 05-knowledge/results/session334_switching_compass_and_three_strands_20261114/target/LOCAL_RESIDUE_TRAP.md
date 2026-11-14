# An all-scale two-hub trap that a ternary product escapes

Let $W_{k,\ell}=11\,0^k\,2^\ell\,311$ have hubs $x\to y$, with letters $2[x\to v]+[y\to v]$, spine ordered left to right. The preceding session's smooth rail gives $H(W_{k,\ell})=1+36\cdot 2^k3^\ell$ for $k,\ell\ge 1$.

**Theorem.** For every $k,\ell\ge 10$, $W_{k,\ell}$ is a local minimum of the exact target score $w\mapsto  |H(w)-N|$ for $N=H(W_{k,\ell})+2$, among single-spine-letter insertions, deletions and substitutions, adjacent transpositions, and reversal of the arc between the two hubs. None of these moves even gives $H\equiv  N\pmod{288}$. Nevertheless $N$ is realized exactly by the ordinal sum of a directed triangle with $W_{k,\ell-1}$.

The word at the local minimum starts only two paths away from a target of unbounded magnitude. This is not a theorem about all arc reversals in the ambient tournament graph, compound edits, plateaux, or missing values. The exact witness uses an additional hub (one vertex of the triangle suffices) and is not strong.

## The finite check is uniform over both exponents

For a fixed order of the two hubs, assign successive spine vertices to the three increasing runs before, between, or after the hubs. A valid assignment gives exactly one Hamiltonian path. First and middle run states are $0,1,2$ for empty, ending with a vertex that beats the next hub, or ending with a vertex that loses to the next hub. The last-run state is just its nonemptiness. For the next spine vertex with $(a,b)$ the indicators for beating it by the first and second hubs, the first-run assignment always works (status becomes (1+a)); the middle-run assignment works if already nonempty or $a=1$ (status becomes (1+b)); the last-run assignment works if already nonempty or $b=1$. At the end, each of the first two runs must be empty or good, and an empty middle run requires the correctly oriented hub arc. Across both hub orders this gives exactly 36 states.

The powers $T_d^m$ of any one letter's transfer, entry by entry, are integer linear combinations of $1,2^m,3^m$ for $m\ge 1$. Indeed, fix a starting state and the set of runs the new batch of $m$ homogeneous letters actually enters. The eligibility and final state depend only on that set, and if enabled, the choices are the surjections onto that set. There are $1$, $2^m-2$, or $3^m-3\cdot 2^m+3$ such surjections. Thus modulo $288=32\cdot 9$ the powers are periodic with period dividing 24 once $m\ge 5$: $2^m$ and $3^m$ have vanished in their respective nonunit components, and the remaining orders divide 6 and 8. The zero-length case is kept explicitly.

There are three ways to cover all mutations without a bound on $k,\ell$:

1. Change a shell letter, insert at a shell or run boundary, or flip the hubs. Both high powers can be replaced by independently selected exponents from $6,\ldots ,29$.
2. Adjacent swaps of unlike letters lie at block boundaries. Unwind one letter from an altered run if present and use the same exponent representatives before unwinding, so the residual exponent still lies in the periodic regime. Swapping equal letters leaves the original count.
3. An edit within a run splits it into powers $a,b$ separated by one substituted or inserted letter (or nothing for a deletion). The untouched run has exponent at least 10, and $a+b\ge 9$ forces at least one split exponent to be at least 5. Thus independently admitting $a,b\in {0,\ldots ,28}$ with at least one $\ge 5$, and the untouched exponent in ${5,\ldots ,28}$, is a safe superset. We need no compatibility assertion about the overcounted representatives.

The replay evaluates these templates modulo 288. The seven groups have respectively $6,4,49,61,6,98,45$ residues (hub flip, shell deletion, shell replacement, boundary insertion, swap, split of the 0-run, split of the 2-run). The union has 125 residues, all odd and none equal to 3. The unchanged count is 1 mod 288, and $N$ is 3 mod 288. Since the score starts at 2 and all neighbors have odd counts, strict improvement could only mean attaining $N$ itself.

Using just one lens would not certify it. At $k=\ell=10$, with indices starting at 0, replacing the letter at index 14 by 1 gives count change $-349116574\equiv (2,5)\pmod{(32,9)}$, whereas replacing the letter at index 0 by 2 gives $40536336062\equiv (30,2)\pmod{(32,9)}$. Each modulus admits an apparent improvement, but they correspond to different physical words.

## The exact escape changes the operation, not the tolerance

The product law for ordinal sums follows because no path can return from the second block into the first; the connector from the final vertex of the first path to the initial vertex of the second always exists. Hence
$$
 H(C_3\to W_{k,\ell-1})=3H(W_{k,\ell-1})
 =3+36\cdot2^k3^\ell=H(W_{k,\ell})+2.
$$
The ternary strand of internal path assignments scales the nonconstant term without scaling the baseline. The triangle strand scales both. Those operations almost agree, and their exact difference corrects the target error. The statement is unbounded and exact, but does not show how to locate such an asymmetry for arbitrary $N$.

Run `python cert.py` from this directory (NumPy is the only nonstandard dependency).

## Even a whole nonlocal switching compass can miss the repair

If $k,\ell$ are *positive multiples of 24*, the same target $N=H(W)+2$ is absent from the full fixed-spine switching compass of $W$, including unequal hub switch bits. This still does not assert anything about switches with newly chosen hub pairs or a sequence mixing switches and edits.

Retain the named spine vertices. To keep them transitive, a switch must have at most one membership change in their order (any alternating triple becomes a directed triangle), so use a prefix cut modulo global complementation. With both hubs unchanged the result is the unturned suffix followed by the letterwise complement of the cut prefix. Switching both hubs further complements the whole word and leaves the hub arc unchanged. Switching just $y$ takes either case to the digitwise XOR-1 word and **also** reverses the hub arc. (An isolated hub-arc flip is not this switch.) These four choices at each cut exhaust the switching slice.

The seven run boundaries use representative exponents $24,24$. For an interior cut in one long block, represent its split exponents $a,b$ in $1,\ldots ,28$ with at least one $\ge 5$ and $a+b\equiv 0\pmod{24}$; the untouched block exponent is represented by 24. The homogeneous transfer periodicity applies to the complemented and XOR-transformed letters too. The union of all four phase readings over these 7 boundary charts and 2 interior charts has 49 residues mod 288, none equal to 3. Run `python compass_escape.py` to replay. The escape must therefore go outside even this nonlocal family, and the directed-triangle product still repairs it by exactly two.
