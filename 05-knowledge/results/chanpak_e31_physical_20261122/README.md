# A single arc straddles the perfect number at exponent 31

**PROVED.** Two strong tournaments on 60 vertices differing in exactly one arc have Hamiltonian-path counts
\[
        2{,}305{,}843{,}008{,}139{,}952{,}127
        \quad\hbox{and}\quad
        2{,}305{,}843{,}008{,}139{,}952{,}129 .
\]
Their midpoint is \(E_{31}=2^{30}(2^{31}-1)\). A second construction on 62 vertices uses two independent binary addresses. The shorter address selects a four-vertex nontransitive tournament on which exactly one full Hamiltonian path starts at a tag. This lets the longer address set the absolute count while the arc continues to change it by two. The same six-vertex core also gives a physical edge at \(E_{23}\), an integer of Euclid's form with composite Mersenne factor.

The earlier perfect-number band used unequal flanks at exponents 17 and 19. At 31 both odd neighbors were already realized separately by ordinal products. The earlier finite one-sided search left the single-arc question open. The new operation is to insert a nontransitive marked tournament between the flipped core and the growing flank. [Previous target and tested boundary](../chanpak_gap_star_ports_20261122/PERFECT-TARGETS.md) · [Previous scalar factorization](../chanpak_gap_star_ports_20261122/STAR-ARMS.md) · [Unequal-flank band](../chanpak_gap_star_ports_20261122/PERFECT-BAND.md)

The two-word construction admits a complete inverse for a prescribed odd lower endpoint, and its two endpoints together cover less than 13.3% of odd integers in upper density. Repeating fixed blocks gives a further obstruction across exponents. For flips incident to the apex, a different inverse gives at most four address pairs to test. The six Leech-related levels admit simultaneous congruence corrections of arbitrary precision, including an explicit lift from exponent 31 toward 61. [Two-word inverse and density](ALL-WORDS.md) · [Across exponents](../chanpak_exponent31_seams_20261122/README.md) · [Apex-arc inverse](../chanpak_perfect_spoke_inverse_20261122/README.md) · [Leech levels](../chanpak_leech_perfect_phases_20261122/README.md) · [Exponent 61](../chanpak_leech_perfect_phases_20261122/EXPONENT-61.md)

## The integer that controls the surviving arc

A deletion tournament \(D\) carries a set of tags. Write \(h_D=H(D)\), let \(k_D\) count full paths after attaching an apex beating exactly its tags, and let \(Q_D=\sum_{S\subseteq D}H(S)H(D-S)\), with \(H(\varnothing)=1\). For any nonempty subset, count its spanning paths starting at a tag, weighted by the complementary path count. Write \(c_D\) for their total and \(z_D\) for the same sum signed by \((-1)^{|S|}\). Write \(t_D,u_D\) for the unsigned and signed analogous sums over paths ending at an **untagged** vertex.

Suppose \(D,D'\) differ on one internal arc and have equal \(h,c\). Append a nonempty tagged tournament \(V\) in ordinal order. Use all the old tags and the selected new ones. Put \(\Delta k=k_{D'}-k_D\), \(\Delta z=z_{D'}-z_D\). In the glued tournament the two new \(h,c,z\) reads agree; the apex-count difference is exactly
\[
                  h_V\Delta k-R_V\Delta z,
\]
where \(R_V\) counts **full Hamiltonian paths of \(V\) starting at an untagged vertex**. Thus this is also the exact signed gap after *every* common transitive right flank. In particular the equality \(h_V\Delta k-R_V\Delta z=\pm2\) is necessary and sufficient for a unit gap in this construction.

**Proof.** The shared-apex ordinal law gives the glued reads
\[
 h'=h_Dh_V,\quad c'=Q_Vc_D+h_D c_V,\quad z'=h_Dz_V,\qquad
 k'=h_V k_D+h_D(k_V-h_V)+c_Dt_V+z_Du_V.
\]
The signed-end identity gives \(u_V=-R_V\). It can be seen by fixing the terminal vertex of the selected path, removing it and applying the resulting signed recursion: its signed-end read is minus the number of full paths starting there. Summing over untagged terminals gives the statement. For a suffix address \(y>0\) the absolute output is \(k'+c'y+h'p(y)-z'[y\text{ odd}]\). Its binary word at the common apex is the reverse complement of the standard bits of \(y\); address zero means the empty suffix. The two outputs therefore differ by the displayed determinant for all \(y\). The edge itself stays inside \(D\). [Ordinal arm identities](../chanpak_moonshine_prime_blocks_20261118/ALIGNED_ARMS.md) · [Signed end identity](../chanpak_branching_inverse_20261122/TWO-FLANKS.md) · [Suffix count and binary recurrence](../chanpak_gap_star_ports_20261122/README.md)

The binary recurrence here is \(p(0)=0\), \(p(2v+1)=p(v)\), \(p(2v)=p(v)+v+[v\text{ odd}]\). When \(\Delta k=\Delta z=2\), the output half-gap is simply the **number of full paths starting tagged in \(V\)**. This nonnegative integer counts the precise freedom retained by the mark. With one such path one may change the other path counts and the suffix coefficient without changing the gap. For a desired center \(M\), each qualifying \(V\) now gives one exact test on \(y\): put the lower adjusted intercept \(K=k'_{\rm low}-z'[y\text{ odd}]\), and solve \(M-1=K+c'y+h'p(y)\) by the two-half integer inverse in the preceding packet. Evaluating the two exact heights finishes a positive certificate.

## Three certificates for the prescribed center

Use the ordinal deletion \(D\to V\), the common apex and the suffix just defined. Outgoing binary masks use vertex zero for the low bit. The table gives the **lower** deletion, its flip, both tag masks and the right address. The other deletion is obtained by reversing the indicated arc of \(D\).

| \(D\) | Tags and arc | \(V\) | Tags | \((h',c',K)\) | \(y\) | Order |
|---|---|---|---:|---|---:|---:|
| 4,5,24,19,35,15 | 1, (1,3) | 2,8,3,5 | 8 | (125,1071,1021) | 1973896132951836 | 62 |
| 4,5,24,19,35,15 | 1, (1,3) | 0,5,17,7,11 | 17 | (125,2267,1021) | 986948066475918 | 62 |
| 96,13,25,49,67,22,46 | 4, (0,4) | 0,17,11,19,37,15 | 16 | (1749,53334,51861) | 42773862359093 | 60 |

For example the first deleted core has \((h,c,k,z)=(25,28,43,-10)\); the flip sends \((k,z)\) to \((45,-8)\). Its right factor has \((h,Q,c,z,t,u,k)=(5,32,7,-1,22,-4,11)\). Hence the surviving gap is \(2(5-4)=2\). In all rows the code constructs the labelled tournaments and checks strong connectivity, the single differing unordered pair and both exact counts. The count proof uses the exact transfer on each small deletion, with its data computed by integer subset dynamic programming.

Additional exact certificates give strong physical unit flips at center offsets \(-20,-16,-12,-4,-2,4,8,10,12,16,20\) from \(E_{31}\). The offsets \(-16,-12,+8,+10\) use the ordinal-graft bank, and the remaining seven use one-sided marked cores. Together with offset zero these realize every odd value from \(E_{31}-5\) through \(E_{31}+21\), a band of fourteen consecutive odd integers. Each listed center has its own actual arc; shared heights may have different representatives. The outgoing masks, suffix addresses and strong reduced orders (58–64) are in the [certificate table](e31_certificates.py). The [replay](check.py) checks each pair.

## Two freely chosen words and an explicit formula

Keep the first six-vertex \(D\). Build \(V(v)\) by adding a vertex \(a\) to a transitive spine. Read the ordinary binary expansion of any \(v\ge1\) from left to right, with a one where \(a\) beats a spine vertex. Give \(a\) the **only** tag. A full path starting at \(a\) must traverse the entire spine in increasing order, and its initial vertex is beaten by \(a\) because the leading bit is one. Thus **every** \(v\ge1\) and **every** right address \(y\ge0\) gives a single-arc pair with the two counts \(L(v,y),L(v,y)+2\).

Let \(b\) be the bit length of \(v\), let \(m_v=2^b-\mathrm{rev}_b(v)\), and put \(Q_v=Q(V(v))\). Then
\[
 L(v,y)=25\left[(1+p(v))(1+p(y))+(v+1)(y+1)-[v,y\text{ both even}]\right]
       +28\left[Q_v(y+1)-m_v\right]-10 .
\]
Indeed \(h_V=1+p(v)\) by the one-apex recurrence and \(c_V=v+1\): a path starting at \(a\) is either its singleton or chooses a first beaten spine vertex and any subsequent subset. It finishes a full path at \(a\) precisely when the last spine vertex beats \(a\), so \(z_V=-[v\text{ even}]\). The unsigned dual count of paths ending at \(a\) is \(m_v\) by reversal and complementation. Consequently \(t_V=Q_v-h_V-m_v\), \(u_V=1-h_V\). The even-path apex identity for one tag gives \(k_V=h_V+c_V+z_V\); a nonempty even path cannot start and end at its unique tag. Substitution into the ordinal formula proves the identity. [Even-path apex identity](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md)

The coefficient \(Q_v\) itself takes linear time in the number of bits in the arithmetic-operation model. Choose which of the two ordered path slots contains \(a\). Each spine vertex belongs before \(a\), after \(a\), or in the other slot. If the before part is nonempty, its last bit must be zero. If the after part is nonempty, its first bit must be one. The four states tracking validity of the current before part and whether the after part has started suffice; double their accepted total. Thus the formula applies to binary words of unbounded and independent lengths. With \((v,y)=(6,1973896132951836)\) it yields the 62-vertex certificate above (the spine labeling of \(V\) may be reordered). With \((250,184567259)\) it gives the two full counts \(35{,}184{,}367{,}894{,}527\) and \(35{,}184{,}367{,}894{,}529\) at \(E_{23}\), on 44 vertices. This is an exponent-independent search rule; existence at every exponent stays open.

The earlier true perfect centers at exponents 17 and 19 also fit this same transmission rule. Both earlier cores have \(\Delta k=\Delta z=2\), and both of their right addresses are even. In the actual reversed-complement suffix, the leading vertex is therefore tagged. Taking that vertex as the singleton receiver \(V\) leaves the address \(y/2\) for the remaining suffix. At exponent 31 the four-vertex nontransitive receiver changes the slope. Thus one criterion applies to all three perfect-number hits with the core and reader supplied by the target. [Earlier 17 and 19 data](../chanpak_gap_star_ports_20261122/PERFECT-TARGETS.md)

## Which marked tournaments transmit this unit change?

For the class \(\Delta k=\Delta z=2\) there is an exact recursive description of *all* admissible \(V\). Factor a nonempty tournament into its strongly connected components in ordinal order. Every path starts in the initial component \(C\). For a tag there, its number of starting paths in the whole tournament is its corresponding number in \(C\) times the product of all later component counts. Every vertex of a strongly connected tournament starts a full path by a Hamiltonian cycle; any nontrivial strongly connected tournament has more than one full path by rotating that cycle. Consequently exactly one tagged full path occurs **if and only if** every later component is a singleton, precisely one vertex of \(C\) is tagged, and that vertex starts a unique full path in \(C\). Tags in the later singleton components may be arbitrary.

For a strongly connected \(C\) with a marked vertex \(a\) and more than one vertex, delete \(a\) and tag its outneighbors. Full paths of \(C\) starting at \(a\) correspond exactly to full paths of \(C-a\) starting at these new tags. Apply the preceding rule to the components of \(C-a\); order strictly decreases until the singleton case. This supplies a terminating recognition procedure and a recursive construction in reverse when the required strong component is supplied. The transitive-spine reader is the case in which all the components after deleting the mark are singletons. The classification preserves a marked path rather than imposing a bound on the order of \(V\). The Hamiltonian-cycle fact used here follows from a longest cycle. Any outside vertex seeing both orientations on the cycle can be inserted at a change of orientation. If this never happens, strong connectivity gives an arc from a vertex dominated by the cycle to one dominating it, and the pair can be inserted consecutively. [Camion’s theorem, Proposition 4 in *Sur certains tournois*](https://www.numdam.org/item/M2AN_1972__6_1_27_0.pdf)

### Counting every unit-transmitting receiver

**Theorem.** Up to isomorphism preserving the tags, there are precisely \(F_{2n-1}\) right receivers of order \(n\ge1\) with one tagged full start, where \(F_0=0,F_1=1\). On a specified set of \(n\) labels their number is \(n!F_{2n-1}\). In particular the orders 2, 3 and 4 give exactly **2, 5 and 13** receiver types.

To see this directly, let \(a_k\) count strong rooted tournaments whose marked vertex starts exactly one path. Then \(a_1=1,a_2=0\). When rebuilding a strong rooted tournament of order \(n\ge3\) from the components of its root deletion, suppose the initial strong component has size \(k\). With no singleton tail, the new root beats exactly the unique starting mark in that component, and any of its other vertices beats the new root. With \(l=n-1-k>0\) tail vertices, the last one must beat the new root to maintain strongness. This suffices: the root reaches the initial component at its marked vertex and hence reaches every tail vertex, and all vertices can reach the last tail vertex. Its other \(l-1\) tail contacts are arbitrary. Component order and the selected child are recovered by deleting the root, so
\[
       a_n=a_{n-1}+\sum_{k=1}^{n-2}2^{n-k-2}a_k.
\]
Consequently \(a_3=1\) and \(a_n=3a_{n-1}-a_{n-2}=F_{2n-4}\) for \(n\ge4\). This recovers the earlier rooted Fibonacci count. Appending a sole-exit vertex at the root gives a tournament with exactly one Hamiltonian cycle, the connection to the previously recorded Garey count. [Earlier sole-exit correspondence, §2](../two_branch_rearming_sigma_sweep_20261104/PROOFS.md) · [Garey’s and Shelor’s theorems in the prior work](../uniquely_hamiltonian_one_arc_cut_shelor_20261104.md)

Now the new count of **all tagged receivers** is \(b_n=\sum_{k=1}^{n}a_k 2^{n-k}\). All \(n-k\) singleton vertices after the initial component may be tagged independently. Hence \(b_1=1,b_2=2\) and \(b_n=2b_{n-1}+a_n=F_{2n-1}\) by the Fibonacci recurrence. In every class the unique path starting tagged is fixed by each tag-preserving automorphism. Its ordered vertex list forces that automorphism to be the identity, proving the labelled formula. The coincidence with three levels in the Leech list counts concrete tournament receivers; the lattice construction has a separate predicate. Different receivers may still give the same integer output.

Run `python3 -B 05-knowledge/results/chanpak_e31_physical_20261122/check.py` from the repository root. Use `python3 -B 05-knowledge/results/chanpak_e31_physical_20261122/check.py --exponent 31 --reader 6` to recover all matching right addresses directly from the exponent and the chosen reader. The complete decoder at reader 6 returns exactly the one positive right address listed in the first row. The replay includes independent literal subset counts on both sides of four small two-word examples. Chan–Pak and the all-exponent perfect-center edge remain **OPEN**.
