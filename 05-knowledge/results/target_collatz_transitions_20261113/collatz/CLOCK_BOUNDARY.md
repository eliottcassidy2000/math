# The last clone: exact boundary debt and a native recharge

**Status: PROVED from a complete fixed-skeleton run certificate.** The
two-clock family and its negative-Collatz arc transfer are inherited from
[clone-clock CLOCKS, §§3–4](../../clone_clocks_negative_collatz_20261113/CLOCKS.md).
The digit-17 clock shape belongs to that packet's already explored clock
language; no priority claim is made for the shape itself. The result here
is its role as the actual zero-\(v\) boundary, the exact lost term, and a
literal native recharge at that boundary. The radius-two repair probe is
**FINITE-EXACT** and concerns only one named eight-vertex skeleton.

## 1. The inherited clock and its finite resource

The skeleton \(\Pi\) has ordered vertices

\[
a,v,c,u,d_0,d_1,x_0,f_0,f_1,
\]

all arcs forward except \(d_0,d_1\to a\) and \(f_0,f_1\to c\).
For nonempty tournament modules \(M,N\), its proved law is

\[
H(\Pi[u\leftarrow M,v\leftarrow N])
 =60P_M(3)P_N(2)+H(M)H(N).                              \tag{1.1}
\]

Here \(P_M(t)\) counts ordered \(t\)-slot spanning path covers, allowing
empty slots. With \(M=TT_j,N=TT_k\), \(j,k\ge1\), this is

\[
H(j,k)=1+60\,3^j2^k.
\]

For \(k\ge2\), take the last clone \(w\) of the \(v\)-block and reverse
only \(w\to c\). It becomes the first clone of the \(u\)-block, giving
\((j,k)\mapsto(j+1,k-1)\) and

\[
H\mapsto(3H-1)/2.                                      \tag{1.2}
\]

This is the **minus** Collatz odd step, equivalently the ordinary plus
Collatz odd step followed by one halving on the negative integer \(-H\).
It is not the ordinary positive
\((3H+1)/2\) step.

The clock's available transfer count is \(k-1\). Its arithmetic coordinate
has \(v_2(H-1)=k+2\): the minus map has three more odd-input steps when the
clock first reaches \(k=1\). Of those three, the first two still have odd
outputs and the third has even output, outside tournament path-count
parity. The published nonempty-module formula does not justify the first
two missing physical steps either.

## 2. The complete zero-site signature

Let \(\Gamma=\Pi-v\), retaining \(u\). In the order
\(a,c,u,d_0,d_1,x_0,f_0,f_1\), its outgoing bit rows are

```text
230, 60, 248, 241, 225, 192, 130, 2.
```

The outside in-neighbours of \(u\) are exactly \(a,c\). A path can visit a
substituted \(u\)-module in at most three runs: its first run can start the
path, and every later run consumes a distinct outside in-neighbour.

Let \(a_i\) count legal words on the seven outside vertices and \(i\)
indistinguishable copies of \(u\), with no adjacent copies. The complete
signature is

\[
\boxed{(a_1,a_2,a_3)=(163,180,60),\qquad a_i=0\ (i\ge4).} \tag{2.1}
\]

This is counted directly by the finite recurrence in
[clock_boundary.py](clock_boundary.py): the state consists of the subset
of used outside vertices, the number of used \(u\)-copies, and the last
symbol. Append an unused outside vertex exactly when the requisite arc
exists, or append \(u\) when its visit budget remains and the last symbol
is an in-neighbour. The empty-prefix state permits any first symbol. The
terminal value is one exactly when all outside vertices and all copies
have been used. There are only seven outside bits and at most three
copies; the structural support bound ensures the recurrence is complete.
The saved receipt gives its three integer outputs, independently of
inflated-module path counts.

For nonempty \(M\), write \(h=H(M)\), \(p=P_M(3)\), and let \(K_i(M)\)
count ordered covers into \(i\) **nonempty** paths. Using
\(K_1=h\), \(K_2=P_M(2)-2h\), and
\(K_3=p-3P_M(2)+3h\), (2.1) gives

\[
\boxed{H(\Gamma[u\leftarrow M])=60p-17h.}                \tag{2.2}
\]

This is an all-module identity, not a fit to a few transitive lengths.
One cannot obtain (2.2) by putting an empty module into (1.1): its
derivation used nonempty run sectors at \(v\), and the zero-visit sector
has to be counted separately.

## 3. The boundary flip has debt exactly \(18h\)

Start with \(\Pi[M,TT_1]\), count \(H=120p+h\). The same arc flip
\(v\to c\) is still a legal graph edit. Its result is exactly

\[
\Gamma[TT_1\triangleright M],
\]

because the old \(v\) becomes the first vertex of the enlarged \(u\)-block.
All internal arcs and vertex names of \(M\) are retained. Since ordinal
composition multiplies slot counts,

\[
H(TT_1\triangleright M)=h,\qquad
P_{TT_1\triangleright M}(3)=3p.
\]

Thus the literal output is

\[
\boxed{H_{\rm boundary}=180p-17h
       =\frac{3H-37h}{2}
       =\frac{3H-h}{2}-18h.}                            \tag{3.1}
\]

The middle term \((3H-h)/2\) is the interior clock-transfer law for a
general fixed \(M\); at \(h=1\) it is the minus Collatz step. The boundary
debt is precisely \(18h\), uniformly in the module and its size.

For transitive modules the literal examples are

\[
361\longrightarrow523\quad\text{instead of }541,
\qquad1081\longrightarrow1603\quad\text{instead of }1621.
\]

The new anchor is \(-17\): the output family is \(60\,3^{j+1}-17\).
Its numerical relation to the known negative-Collatz cycle containing
\(-17\) does not turn this graph edit into a further ordinary Collatz
step. The equality (3.1) is the actual statement.

## 4. A literal recharge, with its exact target change

On \(\Gamma[L]\), add a fresh vertex \(v'\) with

\[
a\to v'\to V(\Gamma[L])\setminus\{a\}.
\]

This is an ordinary root-child birth at the **same named vertex \(a\)**.
The resulting graph is exactly \(\Pi[L,TT_1]\). If \(h=H(L)\), its count
changes by

\[
\boxed{H(\Gamma[L])\longmapsto2H(\Gamma[L])+35h.}         \tag{4.1}
\]

Indeed (2.2) gives the old value \(60P_L(3)-17h\), while (1.1) gives the
new value \(120P_L(3)+h\). Combining the boundary flip and this birth,

\[
\Pi[M,TT_1]\longrightarrow\Gamma[TT_1\triangleright M]
\longrightarrow\Pi[TT_1\triangleright M,TT_1],
\]

returns to the same physical chart and gives

\[
\boxed{H\longmapsto3H-2h.}                              \tag{4.2}
\]

The block \(M\), every internal arc of it, and the named recharge vertex
\(a\) persist. The sole \(v\)-vertex is consumed into the \(u\)-block and
then replaced. Repeating is lawful for arbitrarily many rounds; its
scalar law is the already known ternary clock, now realized by an exact
boundary/recharge circuit. At \(h=1\), the numerical route is

\[
361\to523\to1081\to1603\to3241\to\cdots.
\]

It does not trace consecutive minus Collatz steps: the recharge is
\(2H+35\), and the complete circuit is \(3H-2\). In particular it is not
a target-neutral refill of the dyadic resource.

For a fixed amplitude in \(1+60\,2^k3^j\), the parameters \((j,k)\) are
uniquely determined by the count. Therefore no same-count re-presentation
**within this chart** can increase \(k\). A target-neutral refill must
change the external core, the amplitude, or the family of modules and
retain enough marked information to verify the new entry conditions.
The companion constructor's same-count rebase into \(\Pi\) is such a
change of physical presentation; it does not remove (3.1) at the far end.

## 5. Cheap hostile probe against a scalar-fit repair

To repair the missing \(18h\) universally in the zero-site skeleton, one
would need the complete signature \((181,180,60)\), with all higher
sectors zero, giving \(60P_M(3)+H(M)\).

The exact tested universe is all 28 single arc reversals and all 378
unordered two-arc reversals of the literal eight-vertex \(\Gamma\),
keeping the named module site \(u=2\). No single reversal gives even
the singleton-module count 181. Exactly two double reversals give that
count, but their full signatures are

\[
(181,260,132,18),\qquad(181,174,48,0).
\]

Thus neither is a uniform repair. Their literal rows and flips are
preserved in [clock_boundary_receipt.json](clock_boundary_receipt.json).
The mismatch has a precise missing coordinate: a one-point scalar test
sees only \(a_1\), whereas future module use sees all run sectors.
This excludes only the stated radius-two fixed-skeleton repair. It says
nothing about larger edits, more vertices, a different module site, or
other source families.

Run from the worktree root:

```bash
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/collatz/clock_boundary.py
```

The controls use \(TT_1,TT_2,TT_3,C_3\), and the strong four-vertex
tournament. Literal subset path DP checks counts; a separate repeated-site
recurrence certifies (2.1); explicit adjacency comparisons check both the
boundary flip and the recharge. Neither the proof nor the controls claim
exact coverage of arbitrary tournament targets or ordinary Collatz
convergence.
