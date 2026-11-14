# A maximal-depth target line with no ready source

The claim that every maximal-depth line has a physical ready source
is false even if **any rooted tournament** is allowed, not just the
current source bank. The concrete prime cofactor is \(p=25733\).
This does not exclude a wholly different construction for \(7p\).

## 1. Readiness intrinsically forces a one-outneighbor root

For a marked tournament \((T,v)\), write \(H=H(T)\) and
\[
R_{\pm}=\sum_{S\ni v}(\pm1)^{|S|}s_v(S)H(T-S),
\]
where \(s_v(S)\) counts Hamilton paths in \(T[S]\) beginning at
\(v\). Thus \(R_+\) counts pairs of ordered paths, the second
starting at \(v\), which cover the vertices exactly once. Each
Hamilton path gives one such split, and the extras are exactly
the pairs whose first path ends at a vertex beaten by \(v\).
This uses empty first paths, which always join.

**Outdegree lemma.** If \(R_+(T,v)=H(T)\), then
\(\deg^+(v)\le1\).

To prove it, restrict to the vertices reachable from \(v\), which
are its strong component and all succeeding components in the
ordered condensation of the tournament. Call this reachable
suffix \(W\). Suppose \(v\) beats distinct \(x,y\). Among
\(W-\{x\}\) and \(W-\{y\}\), at least one still has every vertex
reachable from \(v\). Indeed let \(D_x,D_y\) be their respective
inaccessible sets. They are disjoint: if a vertex \(z\) belonged
to both, a simple \(v\)-to-\(z\) path in \(W\) would contain
both \(x,y\). Shortcutting from \(v\) directly to whichever occurs
second would avoid the first, a contradiction. Also
\(y\notin D_x\), \(x\notin D_y\). If both were nonempty, take
\(a\in D_x,b\in D_y\). Vertex \(b\) is reachable after deleting
\(x\), so necessarily \(a\to b\) (otherwise \(a\) becomes
reachable too). Similarly \(b\to a\), impossible.

Suppose then \(W-\{x\}\) remains fully reachable. A tournament
fully reachable from \(v\) has a Hamilton path starting there:
start a Camion cycle at \(v\) in its first strong component
(or use its singleton), then concatenate Hamilton paths of all
later components. Use this path as the second path of a split.
All earlier components of \(T\), followed by \(x\), form the first
path, because earlier components dominate \(x\). The join fails
at \(x\to v\), showing \(R_+>H\). This proves the lemma.

Now define a ready signature \((a,b)\) by
\[
(H,R_+,R_-)=(a+b,a+b,-a), \qquad a>0,\ b>0.
\]
A sink \(v\) would have \(R_-=-H\), contrary to \(b>0\).
So \(v\) beats just one old vertex \(u\), and
\(T\) is the width-zero one-contact bridge over \(Q=T-v\).
The exact zero-width transfer gives
\[
(H(T),R_+(T,v),R_-(T,v))
=(h+P+L,\ h+P,\ -h-L)
\]
where \((h,P,L)=(H(Q),R_+(Q,u),R_-(Q,u))\). Hence **every**
ready source must come from
\[
                    (h,P,L)=(a,b,0)
\]
by a literal zero-width bridge. Conversely this signature on
\(Q\) suffices. In particular this reduction is not a mere
restriction to the positive-width source templates.

## 2. The exact finite obstruction at prime cofactor 25733

At \(p=25733\), maximal depth is \(k=0\) since
\(159^2\le p<907^2\). On the line \(47a+56b=p\) with \(a>0\)
odd, \(b\ge2a\) even, necessarily
\(a\equiv59\pmod{112}\) and \(a\le161\). The sole point is
\((a,b)=(59,410)\). A ready source thus requires a rooted \(Q\)
with \(H(Q)=59\), \(R_+(Q,u)=410\) (in fact with \(R_-=0\), but
we will not need this third condition).

Across the SCCs of \(Q\), ordinary \(H\) multiplies. Since 59
is prime, there is exactly one nonsingleton SCC \(G\),
with \(H(G)=59\); all others are singletons. By
[Busch's exact lower bound for strong tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/),
\(G\) has order at most 8: the floor at order 9 is 75
and increases thereafter. If \(u\) lies in \(G\), every later
singleton doubles \(R_+\) and every earlier singleton leaves it
alone. Thus \(R_+(G,u)\in\{410,205\}\) because \(410=2\cdot205\).
Neither value occurs, by the exact small-core replay below. If
\(u\) follows \(G\), the value is a multiple of 59, impossible.
If \(u\) precedes \(G\), then \(R_+=2^t F_2(G)\) for some \(t\ge0\),
where \(F_2(G)=\sum_S H(G[S])H(G-S)\).
Splitting each full path at every one of the \(|G|+1\) cuts gives
\(F_2(G)\ge(|G|+1)H(G)\); all split pairs are distinct. The same
replay shows \(|G|\ge6\), so \(F_2(G)\ge413>410\).

The same obstruction in fact defeats all positive terminal width
powers \(y\) at the length-seven seam within the old source cone.
The output must satisfy
\[
 180131=(57+68y)a+(68+81y)b,\qquad y=2^r,\ r\ge1.
\]
The cone forces \(180131\ge193+230y\), so \(y\le 512\).
For these nine candidate width powers, the only integral
odd/even cone point is \(y=4,(a,b)=(59,410)\); the replay checks
this by exact division. No assertion about all possible
ready signatures outside the cone is needed.

This excludes the whole maximal-depth ready line, rather than a
selected point or a named core bank. The sole nontrivial finite
certificate is the rooted \(H=59\) profile computation. The
replay generates the isomorphism classes by adding one vertex
in every possible way; score-preserving canonical labeling,
exhausting permutations inside score classes, makes the
quotient exact. It then uses literal path subset DP for every
relevant root. There are eight strong \(H=59\) cores of order 7
and four of order 8, and no smaller ones. The strong-order
theorem handles all larger cores, and all singleton padding
has been handled algebraically.

## 3. The numbers 11, 13, and 26 give an actual alternate source

Take the five-vertex tournament with outgoing masks
\((2,20,9,3,13)\), marked at vertex 1, with zero-based indices
and each bit representing a beaten vertex. Directly
\((H,R_+)=(11,13)\). Attach its width-one positive bridge:
the actual outgoing masks become \((98,52,105,99,109,64,2)\).
The rooted signature is \((H,R_+,R_-)=(37,37,-11)\), so its
deletion state is **\((11,26)\)**. This state is impossible
on the terrace: \(b/a=26/11\) forces width exponent 1,
\(q=0\) gives ratio 2, \(q=1\) would need an even first
predecessor \(d=2\), \(q=2\) would give \(3\mid11\),
and \(q\ge3\) gives first coordinate at least 17.

The terminal seven seam gives \(H=7(47\cdot11+56\cdot26)
=7\cdot1973\) on 17 vertices, also checked by literal
Hamilton-path DP of the full tournament. Increasing the
initial bridge's width to 7 keeps the same rooted core and
gives \((11,1664)\), whose cofactor \(93701\) is prime
and lies in the maximal-depth-zero range. The usual line
inverse recovers this actual \(a=11\) point as its least
positive choice. Thus the requested numbers correspond to
two genuinely used channels: ordinary rooted-core count 11
and rooted split count 13, with 26 the first outgoing
deletion coordinate they generate.

Run python ready_line_replay.py for the finite proof certificate
and both literal bridge constructions.
