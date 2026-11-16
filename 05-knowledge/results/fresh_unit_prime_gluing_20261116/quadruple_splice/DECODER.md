# The signed inverse at root-deletion count 3

The [fresh +4 splice](README.md) gives fixed continuations $(33\cdot2^{s-1},33\cdot2^{s-1}+3,3)$ with $s\ge1$. The exact fiber of Hamilton counts is
$$
\begin{cases}
\{3q-4,3q,3q+4\}, & s=1,\\
\{3q-16,3q-4,3q,3q+4\}, & s\ge2,
\end{cases}
\qquad q=11\cdot2^s-7.
$$
A stronger statement holds for **every** continuation $(E,O,3)$ with $E\ne O$: it cannot carry three consecutive odd counts. This includes hosts which share no induced source and which choose unrelated roots.

## Reduction to the terminal component

Let $D=E-O$. If the root lies before a nonempty ordinal suffix, the suffix goes into two optional ordered slots, and its even and odd choices for the first slot cancel in the signed sum. This follows by assigning each odd cycle or singleton of a [positive odd-cycle packing](https://arxiv.org/abs/2307.05569) to a side. Hence $D\ne0$ forces a terminal root. The Hamilton count $a$ of the ordinal prefix multiplies all four profile coordinates of the terminal rooted component, because the first path stays in that terminal component and the prefix is wholly in the second slot. Ordinal singleton components have multiplier one. A terminal singleton root yields $H=g$ and $E=0,O=g$.

Thus for $g=3$ the only remaining cases are a strong rooted component with deletion count 3, or a prefix of count 3 multiplying a terminal component of deletion count 1. In the latter case the terminal component minus the root is transitive. If its contact word is the binary number $J$, with a 1 when the root beats the corresponding vertex in forward order, then $E+O=1+J$ before multiplication and $D=-1+(J\bmod2)$. Group by the first vertex after the root: its weight in the unsigned sum is a power of two given by the number of later vertices. All but the last such vertex have a nonempty potential suffix, cancelling their signed weight. Strongness forces the word to begin in 1 and end in 0, so for the prefix-multiplied case $D=-3$, and the count is precisely $3H_1(J)$, where $H_1(J)$ is the Hamilton count of the binary word over a transitive root deletion. Starting the word with extra zeros leaves $H_1(J)$ unchanged: such zeros are ordinal singleton predecessors of the whole graph.

## Exact reconstruction for a strong terminal component with deletion count 3

The root deletion is an ordinal sum of singleton components, a directed triangle, and more singletons. The nonsingleton component has count 3 and must be a triangle: a strong tournament on $n$ vertices has at least $n-2$ directed triangles, each contributing two over the empty packing. Let $U$ and $V$ be the binary outcontact words on the prefix and suffix, $b$ be the suffix length, and $k\in\{0,1,2,3\}$ be the number of outcontacts in the triangle. At each $k$ all contact choices are equivalent by triangle rotations. Grouping by the first vertex after the root gives
$$
J=(E+O)/3-1=2^b(4U+k)+V.
$$
A first vertex in the prefix is weighted by 12 times a later-singleton power of two (a triangle has $P(2)=12$); in the triangle it has weight $3\cdot2^b$, and in the suffix each weight is three times its later-singleton power of two. The singleton root path contributes 3. The same alternating cancellation for nonempty tails gives $D=-3+3(V\bmod2)$ if $b>0$. If $b=0$, a triangle contact contributes +1 each (the triangle's weighted rooted sign is $-1$, and the added root changes the sign), giving $D=-3+k$.

Strongness requires the root to reach the first old component and be reachable from the last one. This is also sufficient because the old components are ordinal. Thus a nonempty prefix starts in 1, or if empty then $k>0$; a nonempty suffix ends in 0, or if empty then $k<3$. Once $b$ is fixed, divide $J$ by $2^b$, and the quotient by 4, to recover the suffix of width $b$, the triangle contact count and the prefix with no leading zero. For strongness $4U+k>0$, so there are at most $1+\lfloor\log_2 J\rfloor$ possible rooted isomorphism classes. In particular $D=-1$ or $-2$ allows only $b=0$, and all strong cases have $D<0$.

## Closed signed offsets

Write $k=2c+d$, where $c,d$ are binary digits, so $J$ has word $U|cd|V$ (allow one initial zero when needed). In effect the triangle replaces the two bits $cd$. Let $e$ be the first bit of $V$, or 1 when $V$ is empty. Define $L(U)=U+\lfloor U/2\rfloor$. The exact count is
$$
H=3H_1(J)+2\Delta,\qquad
\Delta=
\begin{cases}
0 & cd=00\text{ or }11,\\
L(U)+e & cd=01,\\
-L(U)-1-e & cd=10 .
\end{cases}
$$

Here is an integer transfer proof. In the rooted odd-cycle insertion sum use $A_i$ for weighted nonempty paths of vertex parity $i$ that start at an outcontact and are available for extending through another old component. Use $B_i$ for paths stopped at an incontact. On each side of the selected path all unused vertices contribute their Hamilton counts, factored over the old SCCs. Ending a path is recorded only when its last vertex is first appended; subsequent unused components multiply the already stopped terms. The full half-increase is the even stopped count. On a prefix of singleton bits, $A_0+A_1=U$ by the update $A'_i=A_0+A_1+bi$ for a bit $b$; also $A_1-A_0=U\bmod2$. Thus $2A_0+A_1=L(U)$.

For a triangle tile the even/odd weighted menus for arbitrary paths, paths starting at an outcontact, and paths ending at an incontact are respectively $(3,6),(k,2k),(3-k,2(3-k))$. The menu imposing both endpoints is $(1,1)$ for $k=1,2$, else $(0,0)$. These are directly obtained from the eight triangle subsets. Replace the two singletons $c,d$ by this tile and compare to three times the singleton reading. The difference in the two active coordinates just after the tile is $(d-c,c-d)$. The difference in the even stopped coordinate is
$$
(d-c)(2A_0+A_1)-3(1-d)c+[c\ne d].
$$
After the first suffix bit $e$, the active difference vanishes (its coordinates sum to zero), and the even stopped coordinate changes by $(1-e)(c-d)$. With no suffix, take $e=1$ for zero extra change. Substitution gives the displayed $\Delta$. This proves the formula without assuming independence of contact choices.

## Why three consecutive odds are impossible

For the only multivalued signed case $D=-3$, the allowed half-offsets are from the formula above, with possibly 0 from a prefix of count 3. Distinct positive offsets are never adjacent. Indeed, two $01$ windows with $U>0$ cannot overlap by a bit. The later prefix $U'$ is at least $4U$, so the increase of $L$ exceeds any variation in the one-bit correction. A leading $01$ with $U=0$ has positive offset only when it gives 1, requiring that the significant word starts $11$. A later $01$ then cannot have $U=1$, so its offset is at least $L(2)=3$. Negative offsets are likewise never adjacent: successive $10$ windows do not overlap, and if one has $U=0$, the later prefix contains the initial $10$ and is at least 2. The early magnitude is at most 2 and the later at least 4.

The only remaining potential triple is $(-1,0,1)$. A $-1$ demands $U=0$, a $10$ window and $e=0$, so the significant word begins $100$. A $+1$ would require a leading $01$ with $U=0,e=1$ (significant word starts $11$), or a $01$ at $U=1,e=0$ (starts $101$). Neither is possible. Since the other nonzero signed cases have at most one count, this proves the global obstruction at $g=3$.

## The entire new +4 fiber

For the splice, $D=-3$ and $J=11\cdot2^s$ has word $1011$ followed by $s$ zeros. The baseline $H_1(J)$ is precisely $q$ (its rooted deletion is transitive, and the +4 source has the same unique strong $g=1$ word). With $b\le s-2$, the two bits are $00$ and the offset is zero. At the remaining five values the reconstruction is:

| $b$ | Prefix | $k$ | Suffix | $2\Delta$ |
|---|---|---:|---|---:|
| $s-1$ | 101 | 2 | $0^{s-1}$ | $-16$ |
| $s$ | 10 | 3 | $0^s$ | 0 |
| $s+1$ | 1 | 1 | $1\,0^s$ | 4 |
| $s+2$ | empty | 2 | $11\,0^s$ | $-4$ |
| $s+3$ | empty | 1 | $011\,0^s$ | 0 |

For $s=1,b=s-1$ the signed value is $-1$, so that row is inadmissible. Larger $b$ fails strongness. Thus every listed count has a strong witness, and ordinal padding cannot add a new count to the fiber. For $s\ge2$ the strong rooted isomorphism-class multiplicities are $s+1$ at $3q$, and one at each of $3q-16,3q-4,3q+4$. At $s=1$ the multiplicity is two at $3q$ and one at each of $3q\pm4$.

For a full three-digit search at nonzero signed response, this eliminates the entire $g=3$ layer, not just a particular contact, splice, or bounded host order. It does not eliminate other root-deletion counts or zero signed response.

The obstruction is sharp in window length: at $J=4,D=-3$, the continuation $(6,9,3)$ occurs at counts 13 and 15. By contrast, deletion count 1 with nonzero signed response determines a single Hamilton count, by the transitive-deletion decoder above. Thus a three-choice interface with nonzero signed response must leave both of the first two deletion layers.
