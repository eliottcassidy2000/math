# Select the receiver and suffix from the lower target

**PROVED, within the full Fibonacci construction.** For the six-vertex core with one marked arc used at \(E_{31}\), a prescribed positive odd lower count \(N\) gives a finite search returning **all** normalized rooted receivers and all their right addresses. It takes neither a receiver nor a suffix length as input. A positive output constructs a literal unit flip. An empty list excludes this construction at that target. The large search for \(E_{61}-1\) has not been exhausted. [Underlying arc and normalization](README.md) · [Previous inverse restricted to a binary spine](../chanpak_e31_physical_20261122/ALL-WORDS.md)

At order 123 the congruence modulo three now removes the triangle at every width. The two remaining receiver types require at least 22 bit changes and 22 equal adjacent pairs. The unrestricted suffix equations in these two charts are open. [Two-ended suffix inverse](SUFFIX-INVERSE.md) · [Reflected-address inverse](REFLECTED-SUFFIX.md)

Normalize by absorbing the singleton tail into the suffix and deleting forced sinks. For the remaining rooted strong receiver \(V\) (including the singleton), use the notation \(g,Q,s,m,e\) from the main proof. Here \((a,b,d)=(25,28,10)\). A completion at address \(y\ge0\), with \(j=y\bmod2\), has lower count
\[
 L_j=b\bigl((y+1)Q-m\bigr)
          +a\bigl((y+1)s+g(1+f(y))-e[ j=0]\bigr)-d.
\]
The same identity holds for the singly tagged cores described in the main theorem. We have \(Q\ge2m\) by exchanging slots (the two nonempty end-at-root events are disjoint), and \(g\ge e\). Consequently every solution obeys
\[
                  (2y+1)Q\le 2(N+d)/b.
\]
This is a target budget on the internal receiver as well as the remaining suffix.

When one vertex is added to a tournament, its two-slot count at least doubles. Insert the new vertex into the first or the second path of any old ordered path cover. A tournament vertex has an insertion site in any directed path, including an empty one. Using the first suitable site is injective after remembering the new slot, because deleting the new vertex recovers the original paths. The two slot choices produce disjoint sets. In particular \(Q(V)\ge2^{|V|}\), and a descendant reached by adding \(s'\) vertices costs at least \(2^{s'}Q(V)\).

Start with the singleton. The recursive construction of the strong receiver with distinguished root adds a new root and a tail of \(s'-1\) singleton vertices to a smaller rooted strong component. Skip the whole extension length when \(2^{s'}Q(V)>2(N+d)/b\), before enumerating any of its contact words. At every retained node solve \(N=K_j+cy+h f(y)\) with the exact suffix inverse in each parity and check address zero separately. In particular \(\gcd(c,h)\mid N-K_j\) is a cheap rejection implemented by the inverse. Tree depth is finite, its branching is finite, and deleting the root recovers its unique child and tail word. Thus no normalized receiver is repeated or skipped within the budget. The output counts marked receiver presentations; an isomorphism of the final unmarked tournament can forget the designated decomposition.

**FINITE-EXACT.** The stronger grammar makes a measurable difference. Among all lower targets through 3,000 it hits 50, compared with 45 when the permitted normalized receivers have transitive root deletion. The first added target is 859. The first added target requiring a positive suffix as well is 1,955. For both, the unique normalized receiver has outgoing masks \((2,4,9,3)\), rooted at vertex zero. Its deletion there is the cyclic triangle. Its profile gives \(c=1046,h=125,K_0=859,K_1=909\). Thus the two targets choose right addresses 0 and 1, giving physical edges \(859\to861\) on 11 vertices and \(1955\to1957\) on 12. Both orientations in each example are strong. The restricted comparison includes the singleton and all suffixes; it contains the prior transitive-spine word receivers after normalization. This shows that choosing the nested receiver can decide an absolute target that changing either binary word in the earlier construction could never attain.

Run `python3 -B receiver_select.py --lower 1955` here to obtain the receiver and the address. With no options it also checks all lower outputs through 3,000. The optional `--max-order` truncates the normalized receiver and makes the negative conclusion correspondingly finite. At \(E_{61}\) use the Mersenne carry and terminal-run tests in the [main note](README.md) for the large families already excluded.

## The target also selects the shape at a fixed order

Write \(n=|V|\), \(k=\operatorname{bitlength}(y)\) with \(k=0\) at zero, and \(m=7+n+k\) for the total order before removing any exterior singletons. We have \(2y+1\ge2^k\). A strong receiver of order \(n\ge3\) has at least \(n-2\) cyclic triangles. From the odd-cycle formula \(H(S)\ge1+2c_3(S)\). Summing \(H(S)H(\bar S)\) over all subsets, each triangle lies in one color of exactly \(2^{n-2}\) two-color assignments. Consequently
\[
 Q\ge 2^n(1+c_3(V)/2)\ge n2^{n-1},\qquad
 n\le\left\lfloor\frac{4(N+10)}{28\,2^{m-7}}\right\rfloor .
\]
The triangle bound also follows by taking successive strong one-vertex deletions down to a triangle: in a strong tournament the removed vertex is on a cyclic triangle. [Strong deletions, Moon](https://doi.org/10.4153/CMB-1966-038-7)

At \(N=E_{61}-1\), the six-vertex construction has **at most 123 vertices**. Indeed at 124 the bound would give \(n\le2\), and for the singleton its slope 81 already forces \(k\le115\). Orders above 124 are excluded by \(Q\ge2^n\) and the target budget. At 123, 122 and 121 the respective maximum nonsingleton receiver orders are 4, 9 and 18. At 123 only **three** rooted graphs pass the absolute height bounds:

| Receiver rows | Width of \(y\) | Forced initial zeros after the leading one | Forced terminal ones |
|---|---:|---:|---:|
| (0) | 115 | 0 | 0 |
| (2,4,1) | 113 | 2 | 0 |
| (2,12,9,1) | 112 | 8 | 5 |

For the table set \(Y=2^{k-1}\). Necessarily \(Y\le y\le (N-\min K_j)/c\), so compare the two end words to force the leading digits. Also \(f(y)\le (N-\min K_j-cY)/h\). A zero in bit \(i\) contributes at least \(\lfloor Y/2^{i+1}\rfloor\) to \(f(y)\), which forces the low bits recorded above. In particular the four-vertex receiver that succeeded at exponent 31, with slope 1071, is too large at order 123 here. The upper-order filter is necessary; none of the three surviving cases is asserted to give equality. The unrestricted question still allows different orders and cores. The selector accepts `--vertices m` and absorbs any forced terminal ones into its affine inverse before searching.

## Choose the free receiver contacts by a minimum cut

Fix a normalized child, the length of its next transitive tail, and any right address. All completions differ only in the contacts of one new root. There are \(r=\max(0,\ell-1)\) free contacts for a tail of length \(\ell\) in an admissible extension. The core, suffix, and orientations toward every other vertex are fixed. Apply the [directed cut formula for adjoining one vertex](../chanpak_minimal_core_grammar_20261122/README.md) to the entire lower graph with that new root deleted. The lower path count on the free contacts has form
\[
 L(x)=L(0)+\sum_i a_i x_i-\sum_{i<j} b_{ij}x_i x_j,\qquad b_{ij}>0.
\]
Indeed \(b_{ij}\) includes the two-vertex directed path on contacts \(i,j\), weighted by twice the positive complementary path count. The whole function is determined by \(1+r+\binom r2\) receiver profiles, using \(a_i=L(e_i)-L(0)\) and \(b_{ij}=L(e_i)+L(e_j)-L(0)-L(e_i+e_j)\). Each profile uses the scalar suffix formula, even for a long fixed address.

Give the arc \(i\to j\) capacity \(b_{ij}\) for \(i<j\) and take \(x_i=1\) on the source side. Its cut contribution is \(b_{ij}x_i(1-x_j)\). Subtract the row sum of these capacities from \(a_i\), obtaining a unary coefficient \(u_i\). A nonnegative \(u_i\) supplies capacity \(u_i\) from \(i\) to the sink. A negative one supplies capacity \(-u_i\) from the source to \(i\) and adds \(u_i\) to the constant. Pinning either value of any contact is just a side constraint on the cut. Its minimum plus the adjusted constant is the **attained minimum lower height** among all completions. Branch on unassigned contacts and discard a whole branch when this minimum exceeds the target. At full assignment equality gives the requested receiver. Every output carries the same literal marked unit flip. The matching argument for coordinate faces also shows that every such shell has at least \(\lfloor r/2\rfloor+1\) distinct unit-flip lower heights. It gives no guarantee that these heights are consecutive.

Run `python3 -B shell_cut.py` for an exhaustive comparison on three receiver shells. This cut selects the next receiver at fixed address. The preceding budget recursion selects the receiver and address jointly.
