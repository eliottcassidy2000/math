# A target-preserving source rewrite and two arithmetic sidecars

Status: **PROVED in this note**, relative only to the inherited exact fixed-core attachment formula. The accompanying finite checks exercise actual tournaments; they do not replace the proofs. No claim of complete odd-spectrum coverage is made.

## 1. The source rewrite

Use the fixed-core constructor and binary address convention of [`compiler.py`](../../session251_creative_completion_20261111/ternary/compiler.py). Its integer function satisfies

\[
g(0)=g(1)=1,\quad g(2q+1)=g(q),\quad
g(2q)=g(q)+q+(q\bmod2).
\]

In particular, every value of \(g\) is positive and odd, and
\(g(q)=1\) precisely when \(q=0\) or \(q=2^t-1\), \(t\ge1\).
The latter follows immediately by induction: a positive even input contributes at least two, while an odd input deletes its last digit.

Let \(Q=(18,8,3,5,14)\), with each entry the bit mask of a vertex's outgoing neighbors. Its Hamiltonian-path count is 13. With no left flank and core out-set \(O=3\), the upper of the inherited unit-switch twins has count

\[
T(B)=13g(B)+32B+6\beta+20,\qquad \beta=B\bmod2.
\]

The lower twin, with \(O=6\), has count \(L(B)=T(B)-2\). The all-word switch and rooted-state proof are in [the earlier continuation note](../../five_core_web_20261111/history/continuation/PROOF.md).

Consider these three replacement charts, always with left address zero:

| Source | Rows | Out-set | Right address | HP count |
|---|---|---:|---:|---|
| \(S_{5a}\) | \((2,4,9,19,7)\) | 20 | \(B\) | \(9g(B)+34B+2\beta+28\) |
| \(S_{5b}\) | \((2,4,25,3,11)\) | 18 | \(B\) | same |
| \(S_4\) | \((2,4,9,3)\) | 8 | \(4B+3\) | \(5g(B)+36B+36\) |

The last entry uses \(g(4B+3)=g(B)\). All rows specify literal tournaments, not a formal count semigroup.

**Theorem 1 (exact ray rebase).** For every integer \(B\ge0\), either of the first two replacement charts has count \(T(B)\) if and only if

\[
B=2^k-2\quad\text{for an integer }k\ge2.
\]

Exactly the same iff holds for the third chart. On this ray all four counts are

\[
T(2^k-2)=77\,2^{k-1}-31.
\]

**Proof.** The difference from either HP-9 chart is
\[
4g(B)-2B+4\beta-8.
\]
For odd \(B\), this is 2 modulo 4 and is nonzero. For \(B=2q\), it equals
\[
4\bigl(g(q)+(q\bmod2)-2\bigr).
\]
Since \(g(q)\) is positive and odd, it vanishes exactly when \(q\) is odd and \(g(q)=1\), hence exactly when \(q=2^{k-1}-1\), \(k\ge2\).

The difference from the HP-5 chart is
\[
8g(B)-4B+6\beta-16.
\]
For odd \(B\), it is 2 modulo 8; for \(B=2q\), it is
\(8(g(q)+(q\bmod2)-2)\). This gives the same zero set. Finally
\(g(2^k-2)=2^{k-1}+1\), which gives the displayed common count. ∎

This is a target-preserving rewrite from a source with HP count 13 to either HP-9 source at the **same tournament order**, or to an HP-5 source at **one additional vertex**. It therefore removes a simultaneous five-core failure for every twin pair indexed by this infinite ray: its upper member already has a five-core witness. It does not assert that the lower member is covered by this rewrite. [GRAMMAR.md](GRAMMAR.md) extends this to eight exact equality languages, and [ADDITIVE_CARRY.md](ADDITIVE_CARRY.md) constructs simultaneous rebases containing every finite prefix.

For the root added by the constructor, write the inherited state as \((H,E,U,d)\), where \(d\) is the HP count after deleting the root, \(E\) counts full paths ending at it, and \(U\) is its weighted start-path response. At \(B=2^k-2\), the four states are

\[
\begin{array}{c|c}
Q&(h,13,77\,2^{k-1}+13,13)\\
S_{5a},S_{5b}&(h,9,77\,2^{k-1}+9,9)\\
S_4&(h,5,77\,2^{k-1}+5,5).
\end{array}
\]

Here \(h=77\,2^{k-1}-31\). Adjacent differences are \((0,4,4,4)\). The inherited continuation
\[
Z(H,E,U,d)=(H+U-E,d,2U-d,d)
\]
fixes every difference \((0,c,c,c)\). Thus this equality persists under every future \(Z\), corresponding to \(k\mapsto k+1\). In contrast,
\(O(H,E,U,d)=(H,0,2U,d)\) changes the difference to \((0,0,2c,c)\), and the next \(Z\) exposes a count difference \(2c\). This identifies exactly why the rewrite is not automatically an arbitrary-word replacement.

The two HP-9 charts agree in their full rooted states for every right address \(B\) when \(A=0\). Introducing a left address distinguishes them: the first minus the second count is \(3A+(A\bmod2)\). The restricted interchangeability does not merge the two cores.

## 2. Shell fibres as reflection pairs

The inherited shell identity is
\[
D_k(n)=g(2^k+n)-g(n)=2^k-2\left\lceil\frac{\operatorname{rev}_k(n)}2\right\rceil,
\quad 0\le n<2^k.
\]
Interior fibres pair reversed words \(2j-1\leftrightarrow2j\). Conjugating this elementary reflection by bit reversal gives the explicit prefix rewrite
\[
0^t1s\ \longleftrightarrow\ 1^t0s,\qquad t\ge1.
\]
If the suffix has length \(\ell\) and value \(s\), so \(k=t+\ell+1\), its two addresses are
\[
n_0=2^\ell+s,\quad n_1=2^k-2^{\ell+1}+s,\quad
\delta n=2^k-3\,2^\ell>0.
\]

For \(\ell\ge1\), set \(C=2^\ell-1-\operatorname{rev}_\ell(s)\), and let \(a\) be the leading bit of the \(\ell\)-bit suffix. The elementary concatenation formula
\[
g(2^\ell p+s)=g(p)+g(s)-1+pC+(p\bmod2)(1-a)
\]
follows by induction on \(\ell\). Substituting \(p=1\) and \(p=2^{t+1}-2\) gives
\[
\delta g=2^t+(2^{t+1}-3)C-(1-a)>0.
\]
For the empty suffix, \(\delta g=2^t\).

Let \(\Delta=L(n_1)-L(n_0)\), and
\(J=32\,2^k+13D_k(n_0)\). The three independent choices of core cut, fibre partner, and high-bit attachment produce an actual scalar affine cube
\[
\{L(n_0)+2\epsilon+\Delta\eta+J\zeta:\epsilon,\eta,\zeta\in\{0,1\}\}.
\]
This is a legal-tournament construction; it is not an interval or target-coverage assertion.

**Lemma 2 (prefix-independent bridge).** The diagonal gap in this cube is
\[
L(2^k+n_0)-L(n_1)=
\begin{cases}
96\,2^\ell+39C+13(1-a),&\ell\ge1,\\
102,&\ell=0.
\end{cases}
\]
In particular it is independent of the arbitrarily long reflected prefix \(t\).

**Proof.** Since \(\operatorname{rev}_k(n_0)=2^{t+1}\operatorname{rev}_\ell(s)+2^t\), its shell jump is \(D_k(n_0)=2^t(1+2C)\). For \(\ell\ge1\), substitute this and the formula for \(\delta g\) into
\(J-\Delta=32\,2^k+13D_k-13\delta g-32\delta n\). The coefficients of \(2^t\) cancel. For \(\ell=0\), the parity term of \(L\) changes by \(-6\), giving 102. ∎

For the all-one suffix, \(C=0,a=1\), the bridge simplifies to \(96\,2^\ell\). For the empty suffix the addresses are \(2^k+1\) and \(2^k-2\), and their \(g\)-values are equal. This is the operation that led to Theorem 1's source rebase.

## 3. A divisor-box quotient of the binary energy

Set \(f(n)=n+1-(n\bmod2)-g(n)\). For four bits
\(n=u+2q+4r+8v\), fold the outer pair by \(a=u+v\). Direct expansion of the inherited weighted binary autocorrelation gives
\[
f(n)=4uv+2(u+v)(q+r)+2qr
=2a(a-1)+2a(q+r)+2qr.
\]

Thus the energy really descends to the exponent box
\(\{0,1,2\}\times\{0,1\}^2\), the divisor lattice of \(p^2qr\). This includes its bottom and top; the familiar ten proper nontrivial divisors are a further deletion, which we do not make. The quotient loses the orientation of the two lifts when \(a=1\). Exactly one extra bit \(v\) restores the target count:
\[
L(n)=38a+322v+90q+180r+31-13f(n).
\]

For any upper prefix \(P\ge0\), concatenation gives
\[
L(16P+n)=L(n)+13g(P)-13
+P(707-13\operatorname{rev}_4(n))+13(P\bmod2)(1-v).
\]
Consequently, the physical suffix swaps have exact, all-prefix effects
\[
\begin{array}{c|c}
(u,v):(1,0)\to(0,1)&322+91P-13(P\bmod2)\\
(q,r):(1,0)\to(0,1)&90+26P.
\end{array}
\]
These are actual address operations with known target changes. They are not target-preserving individually. The useful output is a finite quotient plus the minimal orientation sidecar, not a numerological identification of two collections.

## 4. Provenance and transfer boundary

* **Square-sum reflection:** [THM-4510, sum graphs as reflection orbits](../../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md) supplies the operation \(x\mapsto t-x\). Here the source is a two-element reversed-address fibre; the map is conjugation by bit reversal; the preserved predicate is the exact shell jump; the lost datum is the baseline count. Lemma 2 restores that datum explicitly.
* **Paired square-sum compiler:** [THM-4850, ninja pairs and base 49](../../../../01-canon/theorems/THM-4850-ninja-pairs-of-the-square-sum-graph-start-at-42-and-the-base-49-numeration.md) supplies the stronger methodological target: two witnesses with compatible sidecars and an operation for every next digit. Theorem 1 provides one closed ray, not that complete paired grammar.
* **Graceful gap insertion:** [THM-2787, signed-path-sum Weyl orbit and gap-tail leaf insertion](../../../../01-canon/theorems/THM-2787-signed-path-sum-weyl-orbit-and-gap-tail-leaf-insertion.md) and [THM-2793, rooted graceful gap reconstruction](../../../../01-canon/theorems/THM-2793-oriented-ramp-reference-and-rooted-graceful-gap-reconstruction.md) emphasize that an increment set needs its marked baseline and first-moment data. Here the exact rooted state replaces an insufficient scalar equality; \(OZ\) is the cheapest hostile test for dropping that sidecar.
* **Divisor fold:** [THM-4733, the divisor plane of p2qr and marked path tilings](../../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md) and [the arithmetic reset](../../trijection_reset_20261020/arithmetic.md) supply the outer-pair fold and the types \(p,p^3,p^2qr\). The new preserved predicate is the actual four-bit energy \(f\), and the missing orientation bit is exhibited, not discarded.

The strongest new result is Theorem 1: a literal target-preserving reduction of source complexity on an infinite, exactly classified binary language. The unresolved step is to enlarge this language with compatible rewrites so that simultaneous failure of both H13 twins and the binary shell becomes impossible.

## Reproduction

Run `python3 05-knowledge/results/h13_paired_descent_20261111/square_graceful/probe.py` from the repository root. The receipt states the exact finite universes and controls. The symbolic arguments above have arbitrary integer quantifiers.
