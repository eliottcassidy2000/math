# A single arc straddles the perfect number at exponent 31

**PROVED.** Two strong tournaments on 60 vertices differing in exactly one arc have Hamiltonian-path counts
\[
        2{,}305{,}843{,}008{,}139{,}952{,}127
        \quad\hbox{and}\quad
        2{,}305{,}843{,}008{,}139{,}952{,}129 .
\]
Their midpoint is \(E_{31}=2^{30}(2^{31}-1)\). A second construction on 62 vertices uses two independent binary addresses. The shorter address selects a four-vertex nontransitive tournament on which exactly one full Hamiltonian path starts at a tag. This lets the longer address set the absolute count while the arc continues to change it by two. The same six-vertex core also gives a physical edge at \(E_{23}\), an integer of Euclid's form with composite Mersenne factor.

The earlier perfect-number band used unequal flanks at exponents 17 and 19. At 31 both odd neighbors were already realized separately by ordinal products. The earlier finite one-sided search left the single-arc question open. The new operation is to insert a nontransitive marked tournament between the flipped core and the growing flank. [Previous target and tested boundary](../chanpak_gap_star_ports_20261122/PERFECT-TARGETS.md) · [Previous scalar factorization](../chanpak_gap_star_ports_20261122/STAR-ARMS.md) · [Unequal-flank band](../chanpak_gap_star_ports_20261122/PERFECT-BAND.md)

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

The binary recurrence here is \(p(0)=0\), \(p(2v+1)=p(v)\), \(p(2v)=p(v)+v+[v\text{ odd}]\). When \(\Delta k=\Delta z=2\), the output half-gap is simply the **number of full paths starting tagged in \(V\)**. This positive integer counts the precise freedom retained by the mark. With one such path one may change the other path counts and the suffix coefficient without changing the gap. For a desired center \(M\), each qualifying \(V\) now gives one exact test on \(y\): put the lower adjusted intercept \(K=k'_{\rm low}-z'[y\text{ odd}]\), and solve \(M-1=K+c'y+h'p(y)\) by the two-half integer inverse in the preceding packet. Evaluating the two exact heights finishes a positive certificate.

## Three certificates for the prescribed center

Use the ordinal deletion \(D\to V\), the common apex and the suffix just defined. Outgoing binary masks use vertex zero for the low bit. The table gives the **lower** deletion, its flip, both tag masks and the right address. The other deletion is obtained by reversing the indicated arc of \(D\).

| \(D\) | Tags and arc | \(V\) | Tags | \((h',c',K)\) | \(y\) | Order |
|---|---|---|---:|---|---:|---:|
| 4,5,24,19,35,15 | 1, (1,3) | 2,8,3,5 | 8 | (125,1071,1021) | 1973896132951836 | 62 |
| 4,5,24,19,35,15 | 1, (1,3) | 0,5,17,7,11 | 17 | (125,2267,1021) | 986948066475918 | 62 |
| 96,13,25,49,67,22,46 | 4, (0,4) | 0,17,11,19,37,15 | 16 | (1749,53334,51861) | 42773862359093 | 60 |

For example the first deleted core has \((h,c,k,z)=(25,28,43,-10)\); the flip sends \((k,z)\) to \((45,-8)\). Its right factor has \((h,Q,c,z,t,u,k)=(5,32,7,-1,22,-4,11)\). Hence the surviving gap is \(2(5-4)=2\). In all rows the code constructs the labelled tournaments and checks strong connectivity, the single differing unordered pair and both exact counts. The count proof uses the exact transfer on each small deletion, with its data computed by integer subset dynamic programming.

## Two freely chosen words and an explicit formula

Keep the first six-vertex \(D\). Build \(V(v)\) by adding a vertex \(a\) to a transitive spine. Read the ordinary binary expansion of any \(v\ge1\) from left to right, with a one where \(a\) beats a spine vertex. Give \(a\) the **only** tag. A full path starting at \(a\) must traverse the entire spine in increasing order, and its initial vertex is beaten by \(a\) because the leading bit is one. Thus **every** \(v\ge1\) and **every** right address \(y\ge0\) gives a single-arc pair with the two counts \(L(v,y),L(v,y)+2\).

Let \(b\) be the bit length of \(v\), let \(m_v=2^b-\mathrm{rev}_b(v)\), and put \(Q_v=Q(V(v))\). Then
\[
 L(v,y)=25\left[(1+p(v))(1+p(y))+(v+1)(y+1)-[v,y\text{ both even}]\right]
       +28\left[Q_v(y+1)-m_v\right]-10 .
\]
Indeed \(h_V=1+p(v)\) by the one-apex recurrence and \(c_V=v+1\): a path starting at \(a\) is either its singleton or chooses a first beaten spine vertex and any subsequent subset. It finishes a full path at \(a\) precisely when the last spine vertex beats \(a\), so \(z_V=-[v\text{ even}]\). The unsigned dual count of paths ending at \(a\) is \(m_v\) by reversal and complementation. Consequently \(t_V=Q_v-h_V-m_v\), \(u_V=1-h_V\). The even-path apex identity for one tag gives \(k_V=h_V+c_V+z_V\); a nonempty even path cannot start and end at its unique tag. Substitution into the ordinal formula proves the identity. [Even-path apex identity](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md)

The coefficient \(Q_v\) itself takes linear time in the number of bits in the arithmetic-operation model. Choose which of the two ordered path slots contains \(a\). Each spine vertex belongs before \(a\), after \(a\), or in the other slot. If the before part is nonempty, its last bit must be zero. If the after part is nonempty, its first bit must be one. The four states tracking validity of the current before part and whether the after part has started suffice; double their accepted total. Thus the formula applies to binary words of unbounded and independent lengths. With \((v,y)=(6,1973896132951836)\) it yields the 62-vertex certificate above (the spine labeling of \(V\) may be reordered). With \((250,184567259)\) it gives the two full counts \(35{,}184{,}367{,}894{,}527\) and \(35{,}184{,}367{,}894{,}529\) at \(E_{23}\), on 44 vertices. This is an exponent-independent search rule; existence at every exponent stays open.

## Which marked tournaments transmit this unit change?

For the class \(\Delta k=\Delta z=2\) there is an exact recursive description of *all* admissible \(V\). Factor a nonempty tournament into its strongly connected components in ordinal order. Every path starts in the initial component \(C\). For a tag there, its number of starting paths in the whole tournament is its corresponding number in \(C\) times the product of all later component counts. Every vertex of a strongly connected tournament starts a full path by a Hamiltonian cycle; any nontrivial strongly connected tournament has more than one full path by rotating that cycle. Consequently exactly one tagged full path occurs **if and only if** every later component is a singleton, precisely one vertex of \(C\) is tagged, and that vertex starts a unique full path in \(C\). Tags in the later singleton components may be arbitrary.

For a strongly connected \(C\) with a marked vertex \(a\) and more than one vertex, delete \(a\) and tag its outneighbors. Full paths of \(C\) starting at \(a\) correspond exactly to full paths of \(C-a\) starting at these new tags. Apply the preceding rule to the components of \(C-a\); order strictly decreases until the singleton case. This supplies a terminating recognition procedure and a recursive construction in reverse when the required strong component is supplied. The transitive-spine reader is the case in which all the components after deleting the mark are singletons. The classification preserves a marked path rather than imposing a bound on the order of \(V\). [Camion and rooted paths in this workspace](../../../05-knowledge/reference/CORE-PAPERS.md)

Run `python3 -B 05-knowledge/results/chanpak_e31_physical_20261122/check.py` from the repository root. The replay includes independent literal subset counts on both sides of four small two-word examples. Chan–Pak and the all-exponent perfect-center edge remain **OPEN**.
