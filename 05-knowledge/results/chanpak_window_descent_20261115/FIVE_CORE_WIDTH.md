# The count-five core has at most one independent unit switch, at every height

**Theorem.** In the marked ordinal-apex switch framework of
[the carried-window note](CRT_WINDOWS.md), fix ordinary base count
$h=5$, two-slot count $Q$, and the unsigned marked start and unmarked
end totals $S,T$. Then $H(B+\mathrm{apex})$ takes at most two values
over all literal base tournaments $B$ and cut marks with these
parameters. If there are two, their difference is exactly two.
A fortiori this holds after matching the signed data as well.
The [literal repaired pair](CRT_WINDOWS.md) shows sharpness for the
full response. There is no residue-complete independent palette
modulo five, at any padding length.

## 1. Reading the core from both ends

By SCC product and the [Moon–Busch floor](https://doi.org/10.37236/1141),
every $H(B)=5$ base is $TT_a\to C\to TT_b$, where $C$ is the unique
strong four-vertex tournament. Use the outgoing rows $(6,12,8,1)$.
The uncolored two-slot count of $C$ is 32. For each marked subset $U$
of $C$, let $s,r$ be its total weighted starts and its total weighted
**marked** ends, respectively. The four separate weights are
$(7,9,6,5)$ for starts and $(5,6,9,7)$ for ends. These are computed
directly from path enumeration on subsets (replay). Also the marked
start sum $s$ is distinct for all 16 subsets.

Write $x,y$ for the binary integers formed by the $a$ left and $b$
right mark bits, respectively, in tournament order. Keep their
widths including leading zeroes. Put $T'=Q-5-T$ (the total weighted
marked end count, since the total nonempty end weight is $Q-5$).
Then
$$
Q=32\cdot2^{a+b},\quad
S=2^b(32x+s)+5y,\quad
T'=2^a(32\operatorname{rev}_b(y)+r)+5\operatorname{rev}_a(x).
$$
In the start calculation, earlier modules contribute their one-path
counts and later modules their two-slot counts; in the end calculation
the roles reverse. For a same-width split $(a,b)$, $S$ determines
$x$ by integer division by $32\cdot2^b$, because $s\le27$ and
$5y<5\cdot2^b$. Then $y$ is determined modulo $2^b$ since 5 is odd,
and $s$ determines the core mark. Hence any distinct pair in a fiber
must put its cores at different positions.

Compare two candidates with $a'=a+d,\ b'=b-d$, $d\ge1$. Congruences
of $S$ modulo $2^{b-d}$ and $T'$ modulo $2^a$ force literal common
outer bit words. The two candidates thus have local forms
$C_U\,z$ and $w\,C_{U'}$, respectively, with intervening words
$z,w$ of width $d$. After cancelling the common outer words, their
readouts agree iff
$$
2^d s+5z=32w+s',\qquad
32\operatorname{rev}_d(z)+r
=2^d r'+5\operatorname{rev}_d(w).
$$
Read the bits $z_i,w_i$ from left to right. These identities are
equivalent to a path of states
$$
 s_i=2s_{i-1}+5z_i-32w_i,\qquad
 r_i=(r_{i-1}+32z_i-5w_i)/2
$$
whose endpoints are the physical pairs $(s,r),(s',r')$ and whose
coordinates are integral. Integrality of intermediate $r_i$ follows
from solving the second displayed readout backwards from $r'$
with integer coefficients; integrality of $s_i$ is forward.
The $s_i$ lie in $[-5,32]$ by backward geometric summation from
$s_d\in[0,27]$, and the $r_i$ lie in the same interval by forward
geometric summation. Thus arbitrary displacement reduces to a
fixed finite graph, pruned to the states both reachable from and
able to reach a physical pair.

## 2. An exact accumulated count on that graph

Let $N_U$ be the marked apex count on $C$ alone. Let $d_s(U)$ be the
number of full paths of $C$ ending at a mark, and $d_e(U)$ the number
starting at an unmarked vertex. These are precisely the signed start
and end weights (the [cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md),
or signed odd-packing cancellation).

With an appended transitive bit $b$, an apex count receives no
increment if $b=1$, and receives $S+d_s$ if $b=0$; meanwhile
$S$ becomes $2S+5b$ and $d_s$ resets to $5b$. This is the singleton
case of the stated cut transfer. To count $w C_{U'}$, pass to the
converse, where the suffix has the reverse-complement of $w$,
and apply the same calculation. The exact local difference is
$$
\begin{aligned}
\Delta N={}&N_U-N_{U'}+d_s(U)(1-z_1)-d_e(U')w_d\\
&+\sum_{i=1}^d\big[(1-z_i)s_{i-1}-w_i(27-r_i)\big]
+5\big(\#10(z)-\#10(w)\big).
\end{aligned}
$$
Here each adjacent pattern $10$ is counted separately; this count
is invariant under reverse complement. To verify the unsigned
sum, the true forward marked start after the prefix of length
$i-1$ equals $s_{i-1}+32w_{<i}$. The true converse marked start
after the reversed complemented suffix $w_{>i}$ equals
$27-r_i+32\operatorname{rc}(z_{>i})$. These follow by telescoping
the two state recurrences. The terms with a factor 32 cancel in
the difference, since the pair contribution from $w_j(1-z_i)$
is $32\cdot2^{i-j-1}$ in both calculations for $j<i$.

Arbitrary common outer words only add
$$
f\,[d_e(U)-5(1-w_1)] + e\,[5z_d-d_s(U')],
$$
where $f$ is the last bit of the common left word (zero if empty)
and $e$ is one minus the first bit of the common right word
(zero if empty). Indeed $h,Q,S,T$ of the central modules agree;
crossing-cycle terms using unsigned channels then agree, and the
only signed changes are the interactions with those exterior
endpoint biases. Crossing from the common left word to the common
right one across the entire central module uses the same two-slot
count and no signed transport through a nonempty central module.
Allowing both exterior bits independently in $\{0,1\}$ is safe,
even for a fixed signed response requirement.

## 3. Finite closure certificate and consequences

[Replay](replay_five_core.py) computes the 16 core profiles by
literal permutations, prunes the bounded state graph, and explores
its **full closure** from each physical first edge. Its additional
coordinates are the previous two bits and the exact accumulated
integer in the displayed difference. The finite check processes 34
live unweighted states and 48 visited records in the weighted
closures, with all possible total differences exactly $-2,0,2$.
The queue is exhausted; there is no bound on displacement length.
Normal and optimized Python give the [same receipt](replay_five_core.txt).

Every two outputs in a fiber therefore differ by at most two;
their odd parity makes the width at most two in cardinality. This
is sharper than seeing the small initial two-valued example and
trying larger pads. Going from multiplier 3 to 5 does open a true
switch, but simply moving the count-five core further never
supplies an entire digit range. It is necessary to change the
strong core or the interface for this carry architecture.
