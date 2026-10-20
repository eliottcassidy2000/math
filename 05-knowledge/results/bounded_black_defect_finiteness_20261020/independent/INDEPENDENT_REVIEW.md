# Independent verification of the bounded-black-defect cycle mechanism

Date: 2026-10-20. Verdict: all algebraic claims are correct, and self-converseness is not needed. This is a verification of the counting mechanism and its degree consequence, not an independent verification of any external Hamiltonian-cycle lower bound.

## Definitions and scope

Let \(T\) be a tournament on \(n\geq 3\) vertices with a directed Hamiltonian cycle. Write \(A=\operatorname{Aut}(T)\). A Hamiltonian path is an ordered vertex word, not a word modulo reversal. Its *positional reversal* is the permutation interchanging its first and last vertices, its second and penultimate vertices, and so on. A path is blue precisely when this permutation is an anti-automorphism of \(T\).

The action of \(A\) on path words is free: an automorphism fixing the entire ordered word fixes every vertex. Blue status is invariant under \(A\), since positional reversal is conjugated by an automorphism. Thus if \(H\) is the number of Hamiltonian path words, \(h=H/|A|\), and \(b\) is the number of blue path orbits, then

\[
\delta=h-b,\qquad
\#\{\text{black path words}\}=\delta |A|.
\]

No self-converse assumption enters these identities. When \(T\) is not self-converse, no anti-automorphism exists, so \(b=0\). The singleton is a separate harmless case and has no nontrivial directed Hamiltonian cycle.

## 1. Exact cycle-cut orbit count

Fix a directed Hamiltonian cycle \(C=(v_0,\ldots,v_{n-1})\), let \(R(v_j)=v_{j+1}\), and take indices modulo \(n\). Put

\[
D=A\cap\langle R\rangle,\qquad d=|D|,\qquad q=n/d.
\]

Every automorphism of a tournament has odd order: an even-order automorphism has a nontrivial involutory power, and an involution swapping two vertices reverses the arc between them. Hence \(d\) is odd. As a cyclic subgroup, \(D=\langle R^q\rangle\).

Let \(P_i=(v_i,v_{i+1},\ldots,v_{i-1})\) be the cut before \(v_i\). The only vertex permutation taking \(P_i\) coordinatewise to \(P_j\) is \(R^{j-i}\). Therefore \(P_i\) and \(P_j\) share an \(A\)-orbit exactly when \(i\equiv j\pmod q\). The cuts meet exactly \(q\) path orbits.

Their positional reversals satisfy

\[
r_i(v_j)=v_{2i-1-j},\qquad r_i r_a=R^{2(i-a)}.
\]

If one cut \(P_a\) is blue, then \(P_i\) is blue if and only if \(R^{2(i-a)}\in A\), equivalently \(q\mid 2(i-a)\). Consequently the number of blue cut-path orbits is exactly \(\gcd(2,q)\), or is zero if there is no blue cut. Each black cut orbit is a global black path orbit, so

\[
q-\#\{\text{blue cut orbits}\}\leq\delta,
\qquad q\leq\delta+2.
\]

This is valid without self-converseness. In the non-self-converse case it strengthens to \(q\leq\delta\).

## 2. Black-cut cycles

Count directed Hamiltonian cycles up to cyclic rotation, equivalently by their successor permutations. Each Hamiltonian path word whose last vertex beats its first closes into exactly one such cycle.

For each cycle with at least one black cut, select one of its black cut words. Distinct cycles cannot select the same word because the closing edge and all successor edges are uniquely determined. This injects the set of black-cut cycles into the set of black path words. Therefore

\[
HC_{\mathrm{black}}(T)\leq\delta |A|.
\]

No canonical or equivariant choice is needed: there are finitely many cycles, so selecting the lexicographically first black cut in any fixed vertex labelling suffices.

## 3. All-blue cycles

If every cut of \(C\) is blue, then \(r_i,r_{i+1}\) are anti-automorphisms. Hence

\[
R^2=r_{i+1}r_i\in A.
\]

This implication is all that is required for counting. Its converse is not asserted: in even order, \(R^2\in A\) does not by itself force the cuts to be blue.

## 4. Exact full-cycle square-root fibers

For fixed \(g\in S_n\), let \(F_g\) consist only of *full \(n\)-cycles* \(R\) satisfying \(R^2=g\). This restriction is essential.

- If \(n\) is odd, \(R^2\) is itself an \(n\)-cycle, and \(R=g^{(n+1)/2}\). Thus \(|F_g|=1\) if \(g\) has cycle type \((n)\), and zero otherwise.
- If \(n=2m\), the square of a full \(n\)-cycle has exactly two disjoint \(m\)-cycles. Conversely, suppose \(g\) has those two cycles. Fix a vertex \(a\) in the first cycle. For each of the \(m\) choices of \(y\) in the second cycle, define
  \[
  R(g^k a)=g^k y,\qquad R(g^k y)=g^{k+1}a\quad(0\leq k<m).
  \]
  This gives a full \(2m\)-cycle with square \(g\). These are all possible roots: \(R\) commutes with \(g\), must interchange its two cycles, and is determined by \(R(a)=y\). Hence \(|F_g|=m=n/2\), and it is zero for every other cycle type.

In particular \(|F_g|\leq n\). Combining Sections 2–4 gives

\[
HC_{\mathrm{blue}}(T)\leq
\begin{cases}
|A|,& n\text{ odd},\\
(n/2)|A|,& n\text{ even},
\end{cases}
\qquad
HC(T)\leq(\delta+n)|A|.
\]

The displayed bounds are deliberately weaker than some available parity refinements. If \(4\mid n\), then no all-blue cycle exists: \(R^2\) would have even order \(n/2\), impossible in a tournament automorphism group. For odd \(n\), an all-blue cycle has \(R\in A\), and inverse pairing strengthens the blue-cycle bound to \((|A|-1)/2\).

## 5. Minimum-outdegree consequence

The subgroup \(D=\langle R^q\rangle\) has exactly \(q\) vertex orbits, each of size \(d\). It acts transitively on each orbit and preserves all tournament arcs. The induced tournament on each orbit is therefore regular, with internal outdegree \((d-1)/2\).

Writing \(\delta^+(T)\) for minimum outdegree, this proves

\[
\delta^+(T)\geq\frac{d-1}{2}
=\frac{n/q-1}{2}
\geq\frac{n/(\delta+2)-1}{2}.
\]

This lower bound holds separately for every chosen Hamiltonian cycle; one may choose the cycle giving the largest \(d\). No global vertex transitivity or global regularity is claimed.

## 6. Reproducible exhaustive tests

`verify_cycle_bounds.py` uses only Python's standard library, exact integer arithmetic, and exhaustive enumeration. It does not import any earlier probe code.

The script independently computes automorphism and anti-automorphism groups, Hamiltonian path words, literal blue reversals, complete path orbits, directed Hamiltonian cycles, cut orbit counts, and minimum outdegree. It checks every assertion above, including the exact \(0\)-or-\(\gcd(2,q)\) blue cut-orbit count and the injectivity of selecting one black cut per black-cut cycle.

For square roots, it enumerates every permutation in \(S_n\) for each \(1\leq n\leq 8\), every full \(n\)-cycle among them, and every fiber of the squaring map restricted to those full cycles. Thus missing cycle types are checked to have zero roots, not merely omitted from the report.

Run the default requested cutoff:

```text
python3 verify_cycle_bounds.py
```

Run the stronger order-six census, which includes non-self-converse strong tournaments:

```text
python3 verify_cycle_bounds.py --max-order 6
```

The complete results, invariant profiles, and individual strong-tournament records are written to `verification_results.json`. Console output is retained in `verification_results.out`.

### Initial completed census through order five

- 1,099 labelled tournaments were scanned.
- 570 nontrivial strong tournaments were tested.
- 794 Hamiltonian cycles and 3,942 cycle cuts were checked.
- These cycles split into 744 black-cut cycles and 50 all-blue cycles.
- All assertions passed.
- Every strong tournament through order five is self-converse. This is a limitation of this small cutoff, not an assumption in the proof or code; the order-six run below exercises the wider scope.

The square-root checks cover 46,233 permutations in total for orders one through eight. Their nonempty fiber sizes are respectively \(1,1,1,2,1,3,1,4\), exactly as proved.

### Completed stronger census through order six

- 33,867 labelled tournaments were scanned.
- 22,890 nontrivial strong tournaments were tested.
- This includes all 16,800 non-self-converse strong tournaments of order six.
- 62,234 Hamiltonian cycles and 372,582 cuts were checked.
- These cycles split into 61,704 black-cut cycles and 530 all-blue cycles.
- All assertions passed, including every minimum-outdegree inequality and every restricted square-root count.

The retained `verification_results.json` and `verification_results.out` are from this stronger run.

### Hostile control: \(R^2\in A\) does not imply all cuts blue

The complete order-six census contains 240 labelled strong self-converse tournaments having \(H=45\), \(|A|=3\), \(b=3\), \(\delta=12\), and five Hamiltonian cycles, all of which have black cuts. Nevertheless, two of their cycles have \(q=2\), hence \(R^2\in A\).

One explicit example has arcs

```text
0->3, 0->5,
1->0, 1->2, 1->5,
2->0, 2->4,
3->1, 3->2,
4->0, 4->1, 4->3,
5->2, 5->3, 5->4.
```

For its Hamiltonian cycle \((0,5,2,4,3,1)\), the successor and its square in one-line notation are

```text
R  = (5,0,4,1,3,2)
R² = (2,5,3,0,1,4).
```

The square is an automorphism, but none of the six cuts is blue. The example is record `n=6, mask=1332` in the JSON. In the code, pairs \(i<j\) are ordered lexicographically and bit one indicates \(i\to j\), so the mask is fully reproducible.

Portable packaging note: the retained JSON is supplied as deterministic `verification_results.json.gz`. Decompress it to obtain the exact `verification_results.json` produced by the reproduction command. No census data were removed.
