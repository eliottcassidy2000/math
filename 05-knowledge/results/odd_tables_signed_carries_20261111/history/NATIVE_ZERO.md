# A native strong zero-digit return at the actual overshoot

Status: **PROVED** for the universal two-birth operation and the all-k count and marked \((H,F,s)\) formulas. **FINITE-EXACT / exact coefficient certificate** for the complete seven-coordinate response and the declared local repair search. The desired nonzero target \(2673\,2^k+647\) remains unconstructed by this patch.

## 1. The actual baseline, and a two-edge repair

The incoming [target-return packet](../../berggren_target_return_20261111/reset/TARGET_RETURN.md) supplies Q_k by inflating vertex7 of
\[
(162,1020,1017,1009,993,960,897,768,513,1)
\]
to \(TT_k\), then appending a vertex beaten only by vertex0. Four further births, each beaten only by the previous newborn, give
\[
H=3057D+739,\qquad D=2^k,
\]
where the desired next target is \(2673D+647\).

At \(k=1\), the actual 15-vertex overshoot has outgoing rows
\[
(1186,1020,1017,1009,993,960,897,768,513,1,
3070,5119,10239,20479,8191).
\]
Its count is 6853, and vertex14 has sole in-neighbor13.

Reverse the comparisons of vertex13 with vertices11 and12. Equivalently, when adding vertex13, let it be beaten only by11 instead of12; vertex14 is still beaten only by13. All older vertices and arcs are retained. The new rows are
\[
(1186,1020,1017,1009,993,960,897,768,513,1,
3070,13311,2047,22527,8191).
\]
The resulting count is5991. The same operation at every k gives the exact count
\[
2673D+645=3(891D+215),
\]
two below the desired nonzero target. It is a native **zero-digit** operation at that actual baseline, and its preserved interface makes it repeatable.

## 2. A universal two-birth multiplication by three

For a rooted tournament Y let \(h=H(Y)\), \(F\) be the weighted sum of paths ending at its root, and \(s\) the number of Hamilton paths starting at the root. Assume the root b has sole in-neighbor a, and a itself has indegree one. Call this the two-link indegree-one condition.

**Theorem.** Append c beaten only by a, then append d beaten only by c, marking d. The two-link condition is retained, and
\[
(h,F,s)\longmapsto(3h,\ 3h+2s,\ h). \tag{1}
\]
If Y is strong, the new tournament is strong. Thus (1) can be repeated indefinitely, each time using only two new vertices and preserving all old arcs.

**Proof.** Delete a,b from Y and call the remaining tournament X, with the sole in-neighbor of a marked. The degree assumptions identify Y literally as two successive births from X. Write \(h_X,F_X,s_X\) for its marked data. The inherited birth formulas give
\[
H(Y)=h_X+2F_X,\qquad s(Y,b)=F_X.
\]
Every Hamilton path of Y starting at a must next visit b, since b's only in-neighbor is a. What remains is any Hamilton path of X. Hence
\[
s(Y,a)=h_X=h-2s. \tag{2}
\]

The weighted path-end sum at a has three disjoint parts. The singleton path a contributes \(h_X\). A nontrivial path ending at a must enter from its sole predecessor in X. If b lies in the complementary path, b is its first vertex and the contribution is \(F_X\). If b lies in the path ending at a, b must be its first vertex and the contribution is again \(F_X\). Therefore
\[
F(Y,a)=h_X+2F_X=h. \tag{3}
\]

After the birth at a, its new marked triple is
\((h+F(Y,a)-s(Y,a),h+F(Y,a),h-s(Y,a))\). One more birth gives
\[
H'=h+2F(Y,a)=3h,\quad
F'=2h+2F(Y,a)-s(Y,a)=3h+2s,\quad
s'=F(Y,a)=h,
\]
proving (1). The new c has sole in-neighbor a and d has sole in-neighbor c. If Y is strong, every old vertex reaches a and hence c,d, while c and d each reach old vertices; the extension is therefore strong. ∎

In \((h,b,s)\) coordinates with \(b=F-h\), the two physical continuations are
\[
\begin{array}{c|ccc}
&h'&b'&s'\\\hline
\text{usual two births at the current root}&3h+2b&h-s&h+b\\
\text{two births through its predecessor}&3h&2s&h.
\end{array}
\]
Both retain the two-link interface. The second operation supplies unlimited exact ternary zero digits. It does not supply a free +2 correction: its new defect is \(2s\), rather than 1.

## 3. The all-k marked response at the repaired source

For the incoming source after its first two births,
\[
h=891D+215,\qquad s=297D+72.
\]
The theorem gives, with \(F=V\) in the seven-coordinate notation,
\[
\begin{aligned}
H&=2673D+645,\\
F&=3267D+789,\\
s&=891D+215,\\
b=F-H&=594D+144.
\end{aligned}
\]
The complete exact response is
\[
\begin{aligned}
H&=2673\,2^k+645,& E&=1,\\
S&=891\,2^k+215,& U&=169130\,3^k+95820\,2^k-6453,\\
V&=3267\,2^k+789,& d&=1485\,2^k+359,\\
K&=237530\,3^k+134800\,2^k-8982.
\end{aligned}
\]
The finite ordered-run certificate for the homogeneous placeholder gives the following coefficients \(A_1,A_2,A_3\); every higher coefficient through the maximum possible run count is zero:

|Observation|A1|A2|A3|
|---|---:|---:|---:|
|H|5991|2673|0|
|E|1|0|0|
|S|1997|891|0|
|U|692577|603210|169130|
|V|7323|3267|0|
|d|3329|1485|0|
|K|973208|847390|237530|

The exact substitution is \(\sum_j A_j j!\,{k\brace j}\), where the braces denote Stirling numbers of the second kind. Thus the three basis terms are \(1\), \(2^k-2\), and \(3^k-3\,2^k+3\). No polynomial interpolation or asymptotic fit is used.

The coefficient computation enumerates Hamilton paths on subsets of the fourteen nonplaceholder vertices with edge weight
\[
\mathbf1_{u\to v}+x\mathbf1_{u\to *\to v}.
\]
Leading and trailing placeholder runs are included explicitly. Products of these path polynomials over complementary subsets give the marked two-path covers. The ordered-run bijection assigns the k transitive placeholder vertices to the nonempty ordered runs, giving the displayed Stirling factor. This supplies the complete finite coefficient certificate for all k.

## 4. Exact local search scope

The cheap search held every edge incident to the final root fixed and tested all91 internal one-edge reversals and4095 internal two-edge reversals. It also tested all13 relocations of the root's sole in-neighbor. None has count5993.

Exactly two two-edge patches lie at distance two from5993; both have count5991 and the same complete ordered-run response above. One reverses \(\{11,13\},\{12,13\}\), giving the universal operation just proved. The other reverses \(\{0,10\},\{0,11\}\). Testing one further internal reversal at each of these two graphs adds182 narrowly targeted controls and still finds no5993. This does not exhaust arbitrary three-edge patches.

The companion [one-birth closure result](../tournament/BIRTH_CLOSURE.md) checks the separate proposed repair by adding one vertex to the first5991 graph. The current note's positive result is the unlimited native zero-digit operation; the +2 target gate remains open.

## Reproduction

Compile `native_repair.cpp` with C++17 and optimization to a temporary binary, and run it to reproduce `native_repair_receipt.json`. Run `python3 native_controls.py` from this directory for independent literal inflation controls and the universal two-link operation on small sources. The executable is not retained in the repository.
