# A tooth diagram with sockets, a curried birth grammar, and unequal growth

**Proved:** the all-order pointed-incidence recursion, its reflection rule, an infinite obstruction to a reflection-equivariant one-parent tree, the exact rooted two-letter continuation observer, and a restriction on unequal growth with two homogeneous slots. **Open:** all-order Chan–Pak coverage and an unrestricted construction for the Mandelbrot \(3p\) quotient.

## 1. One tooth per class, one socket per pointed incidence

There is a direct construction of a tooth diagram at every order \(n\). It is an **incidence diagram**, not the unmodified [Applegate–Pol–Sloane cellular automaton](https://arxiv.org/abs/1004.3036): its teeth at level \(n\) are tournament isomorphism classes, and its wires are distinct pointed child classes \([T,v]\) of order \(n+1\), from \([T-v]\) to \([T]\). The distinction matters: different orbits of deleted vertices can have the same pair of endpoint classes. The [rendered first six levels](tooth_diagram.svg) use one segment per class, split each chiral pair across the mirror, place self-converse classes on the axis, and highlight the first class that cannot choose a parent compatibly with reflection.

The exact recursive rule is a bijection of finite sets:
\[
\bigsqcup_{[P],\,|P|=n} 2^{V(P)}/\operatorname{Aut}(P)
\ \longleftrightarrow\
\bigsqcup_{[T],\,|T|=n+1} V(T)/\operatorname{Aut}(T).
\]
On the left, a mask \(S\) specifies the old vertices that beat the newborn. On the right, remember the newborn as the distinguished vertex, then quotient by isomorphism. An isomorphism preserving the newborn restricts to an isomorphism of the old tournament carrying its mask, and conversely. So one can grow by masked births, glue isomorphic output teeth, but keep the distinct wires. Forgetting wires gives the \(A000568\) counts. This is a graded composition system, not a group: birth raises order and has no inverse within the grammar. Relabelings and converse are the reversible symmetries.

On wires, converse sends \([T,v]\) to \([T^{\mathrm{op}},v]\), and on a parent mask it complements the mask and replaces the parent by its converse. A fixed wire exists precisely when a reversing isomorphism of \(T\) can fix \(v\). Every automorphism group of a tournament has odd order (an involution exchanging a pair could not preserve that arc). The group including reversers therefore has Sylow 2 subgroups of order two whenever \(T\) is self-converse, and any reversing map has an odd power that is a reversing involution. Such an involution has at most one fixed vertex. Its unique fixed vertex in odd order, up to automorphism, follows from Sylow conjugacy: any conjugating reverser can be made an automorphism by multiplying by the source involution. Thus the number of fixed incoming wires is exactly **one** per odd-order self-converse tooth, and **zero** at even order. This is the precise reflected-socket version of the classical [self-converse structure](https://doi.org/10.4153/CMB-1979-004-6).

### Why just one predecessor breaks the mirror

The first counterexample is \(C_3\oplus C_3\) at order 6, where \(\oplus\) orders the strong components. Converse reverses this component word and fixes the class. Every deletion leaves either \(TT_2\oplus C_3\) or \(C_3\oplus TT_2\), two distinct chiral classes that converse swaps. A deterministic choice of one predecessor could not commute with converse. More generally \(C_3^{\oplus 2k}\) gives the same obstruction for every \(k\ge1\): the two singleton components left by a deletion would have to sit exactly at the center for the word to be a palindrome, and there are an odd number \(2k-1\) of triangles left. Deleting an entire orbit of the reversing involution instead preserves symmetry, at the cost of a two-vertex jump in even order.

The exact first-six-level class counts are \(1,1,2,4,12,56\); pointed incoming wire counts are \(1,2,4,12,48,296\); fixed wire counts are \(1,0,2,0,8,0\), counting the empty parent at level zero. The replay agrees independently from the parent-mask and pointed-child sides. This is an all-order definition with finite displayed levels.

## 2. First appearances are a filtration, not the same relation

For each attainable \(h\), retain **every** tooth with \(H=h\) at its minimal order. The [existing birth-fibre atlas](../birth_fibers_binary_growth_20261116/README.md) gives 320 distinct counts and 5,827 tied first classes through order eight, not arbitrary single representatives. It already records a different failure of one-step ancestry: both 63-path birth classes have no deletion that is itself at the first appearance of *its own* count. They can reach a 17-path birth by two deletions. So even before imposing reflection, it is not correct to replace the full socket graph by a tree of first births.

There are two independent operations here. Extension by a new directed neighborhood moves between levels. Ordered composition \(\oplus\) multiplies Hamiltonian counts, and converse reverses the word of strong components. Minimal order is the min-plus cost over factorizations, using the minimal strong order for each factor. The [known 7/21 exclusion](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) says those two birth fibres are empty; this grammar provides a way to enumerate classes without losing birth ties, not a proof that all other fibres are nonempty. In particular \(63\) is a literal strong birth, not a product derivation from either missing factor.

### Curry the futures, not the integer

Following the older [curried-packet discipline](../../../07-reflections/lrc14-curried-packet-functional-tower-codex-s170.md), associate to a marked source the function \(F_T(w)=H(w(T))\). For future contexts \(w\) consisting of ordinary sinks \(S\) and near-sinks \(N\) toward the retained root, the future of applying \(u\) is the residual function \(w\mapsto F_T(uw)\). This operation composes exactly with word concatenation. The unrooted tooth count is only \(F_T(\varnothing)\).

In this concrete subgrammar the four coordinates \((H,g,E,O)\) are a **complete observable state**: the [transfer maps](../prime_packing_inverse_20261116/README.md) act on them, and four scalar readings recover them. With extensions read left-to-right:
\[
\begin{aligned}
a&=F_T(\varnothing)=H,& b&=F_T(N)=H+2E,\\
c&=F_T(SN)=H+2(E+O),&
d&=F_T(NN)=H+2E+2(E+O-g).
\end{aligned}
\]
So \(E=(b-a)/2,\ O=(c-b)/2,\ g=(c-d)/2+E\). Equal profiles imply equal futures by composition, and equal futures imply equal profiles by these probes. This says nothing about completeness under arbitrary nonuniform contexts.

The two distinct 9-path first births \((2,4,9,19,7)\) and \((2,8,3,21,7)\), at roots 0 and 1 respectively (outgoing bit rows), read \((9,17,31,29)\) and \((9,17,27,33)\). Thus even the same birth count **and the same one-step read** do not determine the next generation. The differing socket/future data are measurable, not just aesthetic.

## 3. Unequal growth: two at-most-two-run homogeneous slots cannot use distinct rates

Let \(Q(p)=4^p+2^{p+1}+3\), the [Mandelbrot \(3p\) quotient](../crossover_rounding_mandelbrot_20261116/README.md). Mark two vertices in a fixed tournament pattern and replace them by homogeneous transitive slots of arbitrary lengths \(m,n\), each permitting at most two runs in any compressed Hamiltonian trace. If this realizes \(Q(p)\) infinitely often, then on all sufficiently large hits the schedules must be of the form
\[
m=p+u,\qquad n=p+v
\]
for one of finitely many fixed integer offset pairs \((u,v)\). This leaves the fixed-offset question at general skeleton size, and leaves higher-run slots and nonuniform contacts open.

**Physical input: both visible doubles can coexist.** A slot can be visited twice only if its vertex lies in a nontrivial strong component of the skeleton (the intervening walk gives a directed cycle). If the two marks lie in distinct nontrivial strong components, use a Hamiltonian cycle within each, repeat the marked endpoint, and concatenate the components. If they lie in the same strong component, say \(a\to b\), take a Hamiltonian cycle starting at \(b\) and ending again at \(b\), then prepend \(a\). The trace uses both marked vertices twice and every unmarked vertex once, with all adjacencies directed. Pad the other strong components ordinally. In particular, forbidding a double-double trace while keeping both individual double directions is impossible. This rules out the tempting additive architecture with lengths \(2p+u,p+v\) before searching for coefficients.

**Arithmetic input: distinct rates cannot hide in the crossed term.** Write the physical readout
\[
H(m,n)=A2^{m+n}+B2^m+C2^n+D.
\]
The trace formula gives \(A\) equal to the number of double-double traces, hence \(A>0\) if both slots can double. If only one can double, the readout is a single-slot affine exponential and cannot match \(Q\) infinitely often. If one width stays bounded on a subsequence, fix that width, absorb it into the finite core, and apply the same single-slot argument (or the [general single-slot obstruction](../homogeneous_mandelbrot_width_20261117/README.md)). Hence along any unbounded sequence of hits both widths go to infinity. The ratio \(A2^{m+n}/4^p\) tends to one; being a fixed integer times a power of two, it eventually equals one exactly. Thus \(A=2^r\) and \(m+n=2p-r\). Reducing the remaining identity modulo \(2^{\min(m,n,p+1)}\) forces \(D=3\). Put \(u=m-p,\ v=n-p\). Then
\[
u+v=-r,\qquad B2^u+C2^v=2.
\]
If both coefficients are nonzero, neither offset may escape to infinity: one term would diverge while the other tends to zero. If one coefficient vanishes, the other equation fixes its surviving offset immediately. If both vanish there is no hit. This proves the claimed finite offset reduction without any asymptotic assumption on the chosen schedules.

The constructive direction left by this analysis is narrower than “try unequal widths”: to get distinct growth rates one must admit higher-run terms or genuinely nonuniform contacts, rather than hope two homogeneous two-run slots can act as independent exponential registers.
