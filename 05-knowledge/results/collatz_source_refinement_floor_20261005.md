# Collatz source floors that survive refinement

**Status: PROVED scoped, with a FINITE-EXACT kernel certificate. Universal
positivity and universal Collatz remain OPEN.** This note gives a sharp
nonlinear lower bound through the designated inverse section, a way to
transport a floor without losing it merely by subdividing a certificate,
and explicit source-only floors on a decidable infinite binary tree. It also
recasts residue refinement as an energy budget for evaluating one fixed
source. The last reformulation identifies a concrete missing estimate; it
does not supply that estimate for every source.

The main inequality, for every odd \(n>1\) coprime to three, is

\[
 \boxed{\lambda(\rho(n))^2\ \geq\ {160\over441}\,W(n)^3.}                 \tag{1}
\]

Here \(W\) is the beta mixture from
[adaptive mixture flow, P2](collatz_adaptive_mixture_flow_20261005.md),
\(\lambda\) is its deficit probability on odd multiples of three, and
\(\rho(n)\) is the least predecessor divisible by three. Both the coefficient
and exponent are optimal in the senses stated below. A known positive floor
at \(n\) consequently gives a positive floor at this *particular* leaf that
persists in every finer residue class containing it. The inequality also
holds when both sides vanish, so it does not start a positive floor from zero.

## Inheritance and the concept board

The anchor is source positivity under arbitrarily deep refinement. The niche
is the fixed-source evaluation norm from moment selectors and electrical
networks. The wildcard is the use of continuation profiles and intrinsic
parent maps to retain certificate context.

- Closest proved mechanism: adaptive P2–P7, especially the positive difference
  array, exact incoming fibres, and the *complete* cost-eight kernel.
- Canonical hostile: the deleted geometric atom in
  [Fourier atom positivity](collatz_fourier_atom_positivity_20261005.md).
  Every finite residue class has positive mass, while the designated atom is zero.
- Corrected near miss: the bounded inverse-size section in
  [inverse predecessor sections](inverse_predecessor_sections_20261005.md)
  has no uniform positive *linear* relative-weight bound.
- Least-used sidecar: the signed difference between the counters of two
  checked common-future words, retained before taking any positive part.

The live board has five objects: **a labeled source, its mixture profile,
signed counter differences, nested evaluation energy, and a proper parent
map**. Every transfer below either keeps the label or explicitly states the
information it loses.

Incoming commits 2b38c40841, 5e291040bb, b92c07a9bd and e13c4223c9 were read
as mathematical input. The Pascal classification survives the audit; claims
that its discounted array is an exact split array and that it excludes every
finite-dimensional representation do not. The repairs are recorded in the
[Pascal note](collatz_pascal_boundary_leaf_section_20261005.md) and
[mistakes ledger](../../01-canon/MISTAKES.md). During this session the
independent audit in commit 905023c4dd reached the same corrections; its
incoming repairs were retained instead of duplicating that edit.

## Definitions and the complete finite boundary

Let \(U(n)=\operatorname{oddpart}(3n+1)\), and \(S(n)=4n+1\). A positive odd
source is rooted if its orbit first reaches 1 after a finite strict word
of valuations \(a_1,\ldots,a_\tau\). For a rooted source \(n>1\), put

\[
 L(n)=\tau-1,\qquad
 K(n)=\sum_{i=1}^{\tau}\left\lfloor{a_i-1\over2}\right\rfloor.
\]

Every such source has \(K\geq1\). For \(0<r<1\), define

\[
 f_r(n)=r^{K(n)}(1-r)^{L(n)},\qquad
 W(n)=\int_0^1 2(1-r)f_r(n)\,dr
     =w(L,K)={2K!(L+1)!\over(L+K+2)!}.                                \tag{2}
\]

Set both functions to zero at unrooted sources, and set their value at 1 to
one. These are definitions without a convergence assumption. In particular,
\(0\leq W(n)\leq1/3\) for every \(n>1\).

The adaptive note proves that the incoming deficit is
\(\lambda(n)=W(n)\) on odd multiples of three and zero at units. Its total
mass is one. If \(U(z)=n>1\) and \(v_2(3z+1)=a\), then

\[
 (L(z),K(z))=(L(n)+1,K(n)+j),\qquad j=\lfloor(a-1)/2\rfloor,            \tag{3}
\]

whenever rooted. Rootedness is equivalent on the two ends of the edge.

The designated section has the following exact table:

| \(n\bmod9\) | 1 | 2 | 4 | 5 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| \(a_0(n)\) | 6 | 5 | 4 | 1 | 2 | 3 |

Thus \(\rho(n)=(2^{a_0(n)}n-1)/3\), \(3\mid\rho(n)\), and
\(\rho(n)<64n/3<22n\). Its sibling cost \(j\) is at most two.

The inherited kernel contains every rooted *base* of cost \(K\leq8\), with
its forward parent and edge cost: 3,591 records. A base has initial valuation
one or two. Independent forward replay and exhaustive inverse boundary
closure verify completeness, with no assumed depth cutoff. Removing base 9
is a hostile control: its remaining records can still be individually valid,
but boundary closure fails.

For each \(q\leq8\), the map \(b\mapsto U(b)\) is a bijection between
nonroot bases of cost at most \(q\) and rooted unit sources \(n>1\) of cost at
most \(q\). The inverse uses the unique valuation in \(\{1,2\}\) permitted by
\(n\bmod3\), so it does not change \(K\). An independent enumeration by outer
sibling depth gives the same sets.

## A sharp floor through the inverse section

**Theorem SF1 (PROVED, using the complete finite boundary).** Equation (1)
holds for every odd unit \(n>1\). The largest possible coefficient in (1)
is \(160/441\). In a bound
\(\lambda(\rho(n))\geq cW(n)^\alpha\) with fixed \(c>0\), no exponent
\(\alpha<3/2\) works for all these sources.

*Proof of the bound.* It suffices to consider rooted \(n\). Since \(j\leq2\)
and \(w\) decreases in either coordinate,

\[
 {W(n)^3\over W(\rho(n))^2}
 \leq R(L,K):={w(L,K)^3\over w(L+1,K+2)^2}.                            \tag{4}
\]

For \(K\geq3\), this formal upper bound decreases strictly as \(L\) increases.
Indeed, put \(x=L+2\), \(k=K+1\). The denominator minus numerator in the ratio
\(R(L+1,K)/R(L,K)\) is

\[
\begin{split}
 &(x+k)^3(x+1)^2-x^3(x+k+3)^2\\
 &=(k-4)x^4+2(k^2-4)x^3
   +(3k+6k^2+k^3)x^2+(3k^2+2k^3)x+k^3>0.
\end{split}                                                        \tag{5}
\]

At \(L=0\),

\[
 R(0,K)={(K+3)^2(K+4)^2(K+5)^2\over2(K+1)^3(K+2)^3}.
\]

This decreases in \(K\), because

\[
 (K+3)^5-(K+6)^2(K+1)^3
 =15K^3+125K^2+285K+207>0.                                          \tag{6}
\]

Consequently all \(K\geq7\) satisfy
\(R(L,K)\leq R(0,7)=3025/1296<441/160\).

The complete cost-six kernel leaves exactly 398 rooted unit sources to check.
Their actual section bills have maximum \(441/160\), attained at 85.
This is a finite proof obligation with an exhaustive, independently replayed
certificate, not a search up to an ordinary-size cutoff. At that source,

\[
 \rho(85)=453,\quad W(85)=1/10,\quad \lambda(453)=2/105,
 \quad {W(85)^3\over\lambda(453)^2}={441\over160}.                     \tag{7}
\]

This proves the inequality and optimal coefficient.

For the exponent, take \(k=9t+1\), \(n_k=S^k(1)\). Then \(n_k\equiv5\bmod9\)
and the actual section is \(z_k=(2n_k-1)/3\). Directly,

\[
 W(n_k)={2\over(k+1)(k+2)},\qquad
 \lambda(z_k)={4\over(k+1)(k+2)(k+3)}.                               \tag{8}
\]

Thus \(\lambda(z_k)/W(n_k)^\alpha\) tends to zero for \(\alpha<3/2\).
This same family refuted the proposed uniform linear comparison. It now
identifies the exact threshold exponent. QED.

In particular, a floor \(W(n)\geq\eta>0\) yields
\(\lambda(\rho(n))\geq(4\sqrt{10}/21)\eta^{3/2}\). Both the source and the
selected leaf stay fixed while their representations are refined.

## Bounded inverse moves and endpoint mass

The sharp constant above uses the particular section table. A broader
estimate handles arbitrary legal bounded valuations and explains the
mechanism without optimizing factorial ratios.

**Theorem SF2 (PROVED with the stated finite kernels).** For \(1\leq J\leq5\)
and any actual inverse edge \(U(z)=n>1\) with
\(\lfloor(v_2(3z+1)-1)/2\rfloor\leq J\),

\[
 W(z)^2\geq {W(n)^3\over D_J},
 \qquad D_J=8(4/J+4^J).                                             \tag{9}
\]

The target need not have a known ROOT word for the inequality to be valid.
As before, applying it to obtain positivity requires a positive starting floor.

*Analytic part.* Suppose \(K\geq2J-1\). Write \(w=W(n)\) and
\(\varphi_J(r)=r^J(1-r)\). For \(t\leq2^{-J-1}\), split at \(r=1/2\).
The mixture mass where \(\varphi_J<t\) is at most

\[
 2\int_0^{(2t)^{1/J}}r^{2J-1}\,dr+
 2\int_{1-2^Jt}^1(1-r)\,dr
 \leq (4/J+4^J)t^2.                                                \tag{10}
\]

The first integral uses the arithmetic restriction \(K\geq2J-1\); arbitrary
nonnegative profiles would not give this estimate. Put
\(A_J=4/J+4^J\) and \(t=\sqrt{w/(2A_J)}\). Since \(w\leq1/3\), the cutoff
condition holds. At least \(w/2\) remains outside the bad set. Equation (3)
therefore gives \(W(z)\geq tw/2\), proving (9).

*Finite boundary.* For \(K\leq2J-2\), use the completed kernel, with the
largest possible sibling increment \(J\). The exact largest bills
\(w(L,K)^3/w(L+1,K+J)^2\) are:

| \(J\) | Cost cutoff | Exceptional units | Largest formal bill | \(D_J\) |
|---:|---:|---:|---:|---:|
| 1 | 0 | 0 | none | 64 |
| 2 | 2 | 5 | \(100/3\) | 144 |
| 3 | 4 | 40 | \(1225/12\) | \(1568/3\) |
| 4 | 6 | 398 | \(784/3\) | 2056 |
| 5 | 8 | 3590 | 588 | \(40992/5\) |

For example the five cost-two units are \(5,7,11,13,17\). This explains why
the earlier elementary bound \(\lambda(\rho(n))\geq W(n)^{3/2}/12\) was
already valid before optimizing to SF1. No claim for arbitrary \(J\) is
inferred from these five completed kernels.

## Signed bills preserve a floor under certificate subdivision

Suppose two checked finite actual words from \(n,h>1\) end at the same
positive odd source, with neither word passing through 1 before its endpoint.
The endpoint may itself be 1. Let \(a\) be the difference of the two sibling
cost sums and \(b\) the difference of their lengths. In either case,

\[
 f_r(n)=r^a(1-r)^b f_r(h).                                          \tag{11}
\]

For an endpoint other than 1 this follows by canceling the shared suffix.
For endpoint 1 the subtraction of one in both definitions of \(L\) cancels.
If the common future is unrooted, both profiles are zero. Guards and actual
word verification remain part of this statement.

**Theorem SF3 (PROVED).** If \(W(h)\geq\eta>0\), set \(m=\max(0,a,b)\).
For \(m=0\), \(W(n)\geq\eta\). In general,

\[
 W(n)\geq {\eta^{\,1+m/2}\over2^{\,2m+1}}.                           \tag{12}
\]

*Proof.* Every nonroot rooted profile satisfies \(f_r(h)\leq r\).
The mixture mass on \(r<s\) or \(r>1-s\) is at most \(2s^2\).
On the remaining interval, negative bill coordinates can only increase
the multiplier, and

\[
 r^a(1-r)^b\geq[r(1-r)]^m\geq(s/2)^m.
\]

Take \(s=\sqrt{\eta}/2\). At least \(\eta/2\) remains, proving (12).
For exact rational certificates one can instead choose dyadic
\(0<s\leq1/2\) with \(4s^2\leq\eta\), and record

\[
 \bigl(\eta-2s^2\bigr)\,[s(1-s)]^m.                                \tag{13}
\]

This is a positive rational floor. QED.

The important composition rule is to **add the signed pairs before quoting
the floor**. Bills telescope through any chain of checked common-future
relations. Subdividing the same certificate or inserting a checked relation
followed by its inverse leaves the total bill and the resulting floor
unchanged. Taking positive parts at every intermediate step destroys this
cancellation.

For a minimal hostile, \(U(3)=U(13)=5\). The two directions between 3 and 13
have bills \((1,0)\) and \((-1,0)\). Their sum is zero, preserving the initial
floor \(W(3)=1/6\). Repeatedly quoting the coarse one-edge bound loses mass
even though the represented source has not changed.

The existing guarded controller \(n=155+2048t\), \(h=111+1458t\) has bill
\((1,6)\), so \(W(n)\geq W(h)^4/8192\). For \(s\) consecutive legal uses,
the aggregate bill is \((s,6s)\); one quote has exponent \(1+3s\) and
denominator \(2^{12s+1}\). Repeatedly applying the already weakened scalar
bound would instead multiply exponents to \(4^s\). The aggregate estimate
still needs a grounded floor at the final child. It neither removes guards
nor proves that every source enters the controller.

This gives two separate meanings of refinement: *certificate subdivision*
preserves the aggregate bill; *residue refinement* preserves a floor attached
to the same atom. They should not be conflated.

## A fixed source and its evaluation energy

Write \(p_m=\lambda(6m+3)\). Let \(q_i\mid q_{i+1}\), \(q_i\to\infty\), and
keep one index \(m\) fixed. Let \(C_i=\{k:k\equiv m\bmod q_i\}\) and
\(c_i=\sum_{k\in C_i}p_k\). When \(c_i>0\), the function
\(u_i=\mathbf1_{C_i}/c_i\) represents evaluation at \(m\) on the space of
functions constant on these residue classes, with the \(L^2(p)\) inner product.
In particular,

\[
 \|u_i\|^2={1\over c_i},\qquad
 \|u_{i+1}-u_i\|^2={1\over c_{i+1}}-{1\over c_i}.                    \tag{14}
\]

The latter identity follows by integrating separately over \(C_{i+1}\) and
\(C_i\setminus C_{i+1}\). The increments are orthogonal: projection to an
earlier partition sends every later \(u_j\) to the earlier \(u_i\).

**Proposition SF4 (PROVED criterion).** A finite bound
\(\sum_i\|u_{i+1}-u_i\|^2\leq E_m\) gives the floor

\[
 p_m\geq {1\over c_0^{-1}+E_m}>0.                                  \tag{15}
\]

Conversely \(p_m>0\) implies a finite total energy, exactly
\(p_m^{-1}-c_0^{-1}\). This is equivalent to atom positivity. Its usefulness
is to identify the *fixed evaluation direction* whose energy must be bounded,
rather than treating positive finite matrices or average spectral bounds as
pointwise information.

In a nonredundant finite basis with positive Gram matrix \(G_i\) and evaluation
vector \(v_m\), the same norm is \(v_m^*G_i^{-1}v_m\). In residue indicator
coordinates it is simply \(1/c_i\). Singular or redundant bases require
restriction to the quotient by null vectors, with evaluation well defined
there; positive semidefiniteness alone is insufficient.

The hostile is exact. Start from \(p_k=2^{-k-1}\), delete its atom at 1, and
renormalize the others by \(4/3\). For \(q\geq2\),

\[
 c_q(1)={1\over3(2^q-1)},\qquad \|u_q\|^2=3(2^q-1).
\]

Every finite value is positive, but the evaluation energy diverges. For the
undeleted geometric measure, the same norm is
\(4(1-2^{-q})\), bounded by four and converging to the inverse atom.
No bound \(E_m<\infty\) for every Collatz source is asserted here.

The concurrent [localized resolvent floor](collatz_localized_resolvent_floor_20261005.md)
uses \(h_m(j)=4t/(1+t)^2\), \(t=2^{m-j}\), and the signed readout
\(9H_{m,d+1}-8H_{m,d}\), where \(H_{m,d}=\sum_jp_jh_m(j)^d\).
Its positive readout is a valid finite atom certificate. It is important
not to replace SF4's evaluation representers by the normalized densities
\(h_m^d/H_{m,d}\): these can have bounded energy while concentrating on a
neighbor. For target \(m=1\) and probability \(1/2\) at each of 0 and 2,
the target is missing, \(H_{1,d}=(8/9)^d\), and every such normalized density
has squared norm one. The signed readout correctly stays zero.
The required sidecar in the energy route is exact evaluation at the target,
not merely a source-centered shape.

## Source-only floors on a proper binary inverse tree

A related construction gives genuine unbounded coverage with a terminating
recognizer, while keeping its scope explicit.

Start at 5. At a unit parent \(n>1\), consider inverse exponents
\(\{2,4,6\}\) if \(n\equiv1\bmod3\), and \(\{3,5,7\}\) if
\(n\equiv2\bmod3\). Exactly one resulting predecessor is divisible by three;
retain the other two. They are distinct units, both greater than their parent,
and have \(U\)-image \(n\). Unique forward parents and strict increase make a
genuine binary tree \(\mathcal T\).

**Proposition SF5 (PROVED).** Membership in \(\mathcal T\) is decidable by
following \(U\), checking that the source is one of the two allowed children
and that the parent is smaller, stopping at 5. Every failure rejects.
For a member \(n\), let \(b=\operatorname{bitlength}(n)\). Then

\[
 W(n)\geq w(3b,9b+1)>0.                                            \tag{16}
\]

*Proof.* Recognition strictly decreases positive integers and hence always
terminates. At tree depth \(t\), equations (2) and (3) give \(L=t\),
\(K\leq3t+1\), since 5 has counters \((0,1)\). Every child \(c\) satisfies
\(c-1\geq(4/3)(n-1)\), hence \(n-1\geq4(4/3)^t\). Since
\((4/3)^3>2\), this implies \(t\leq3b\). Monotonicity of \(w\) in each
coordinate proves (16). QED.

For the labeled leaf \(\rho(n)\), (1) supplies a floor from (16). A simpler
all-rational weakening is \(w(3b,9b+1)^2/12\). If only the leaf \(z\) is given,
recover \(n=U(z)<2z\), run the terminating recognizer, and use
\(b=\operatorname{bitlength}(z)+1\) as an upper bound. Thus the certificate
can be computed from the source itself throughout this family.

The first twelve levels contain 4,095 distinct sources; all are checked.
The rooted source 7 is rejected because its forward step initially increases.
This is a proper infinite subfamily, not universal entry. The elementary
inverse-tree construction is not claimed as a literature novelty.

## Transfers from older topics

These connections transfer operations or proof obligations, not numerical
coincidences or entire theorems from another problem.

| Prior source | Collatz target and map | Preserved predicate | Lost information and required sidecar | Decisive test |
|---|---|---|---|---|
| [THM-2815, finite Laguerre carrier and radial selector access](../../01-canon/theorems/THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary.md) | Send a labeled coefficient selector to the fixed-source evaluation vector of SF4 | The norm of evaluation is an inverse Gram quadratic form | Finite quadrature nodes need not persist under refinement; retain the same source, compatible spaces, and a norm bound | Deleted geometric atom has finite positive Gram matrices but divergent evaluation norm |
| [THM-4015, character-sensitive Foster transference](../../01-canon/theorems/THM-4015-first-kind-character-sensitive-foster-transference.md) | Send effective resistance in a fixed direction to the cost of evaluating a specified atom | Positive quadratic energy and a direction-specific inverse operator | A trace or averaged Foster identity may find some cheap edge, not the specified source; retain its evaluation vector | Compare the bounded undeleted selector with the divergent deleted selector |
| [THM-2176, Gordian continuation profile and interaction cocycle](../../01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md) | Send a full continuation profile to the checked price profile, and telescoping differences to the signed bill of SF3 | Composition with context and cancellation before scalarization | A scalar price erases sign and intermediate context; keep both actual words and their guards | Insert the checked \(3\leftrightarrow13\) detour |
| [THM-3357, Berggren three-branch parent circuit](../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md) | Transfer intrinsic parent recognition and proper height to the two-branch tree of SF5 | Unique parent, a decreasing recognition algorithm, bounded local cost | An abstract branching tree does not identify arithmetic chambers or all integers; keep inverse exponents and a root stop | Reject 7 despite its known rootedness |
| [THM-3452, unequal-depth Hensel Heisenberg orbit law](../../01-canon/theorems/THM-3452-unequal-depth-noncommuting-smooth-hensel-heisenberg-orbit-law.md) | Candidate: retain carry coordinates when transporting nested source selectors | A proposed exact correspondence would have to preserve the fixed source and total energy | Finite residue surjectivity and mixed-depth rank forget point mass; an actual compatible map is still missing | A full residue census alone fails on the deleted-atom control |
| [THM-4482, ranks are tensions strategy cube](../../01-canon/theorems/THM-4482-ranks-are-tensions-strategy-cube.md) | Use dual obstruction witnesses to test each proposed weight class before searching it further | Scope-specific infeasibility rather than a global impossibility claim | A no-go theorem for one observer or bounded correction does not exclude an unbounded labeled representation | Check hypotheses before applying the inherited lookahead obstruction |

The first four transfers produced the explicit objects SF1–SF5 or their
tests. The Hensel transfer remains a proposed direction: no Collatz energy
bound has been imported from a different orbit theorem. The rank work
primarily improves the scope audit.

## Incoming deadlines and the next interface

Commits 905023c4dd, 5cf951f802, dedefebd4f and bfece2fb8a arrived while the
exact checks were running. Their mathematical packages were read before
integration. The [floor deadline theorem](collatz_floor_transport_deadlines_20261005.md)
and [threshold compiler](collatz_weight_threshold_receipts_20261005.md)
strengthen the operational consequence of this note. A true rational floor
\(W(n)\geq\eta>0\) forces \(N=L+K\leq B(\eta)\), where

\[
 B(\eta)=\max\{b\geq0:\eta(b+1)(b+2)\leq2\}.
\]

The actual ROOT word has length at most \(B(\eta)\). At the selected
predecessor, retaining the original floor and the actual single inverse edge
gives the better deadline \(B(\eta)+1\). Rebuilding a deadline from the
smaller scalar SF1 floor discards this information. Likewise, SF3 should
retain its complete word relation in addition to its scalar floor.

The concurrent [polynomial atom dual](collatz_atom_polynomial_dual_20261005.md)
and localized resolvent package give finite signed tests that can provide
the initial floor if their moment premises are independently established.
They complement SF4: a uniform bound on an exact evaluation norm and a
positive signed minorant are different certificate formats for the same atom.
Neither a positive Gram matrix nor a bounded norm for a nearby density can
substitute for either format.

The [refinement and shadows note](collatz_refinement_floor_shadows_20261005.md)
supplies actual rising-word cones and the exact segment-exit mechanism.
The integration audit retains those mechanisms but scopes three stronger
claims. An ordinary-size increase need not have expanding affine slope:
\(165\) reaches \(167\) in 17 steps with total valuation 27, so
\(3^{17}<2^{27}\); the positive carry matters. This refutes a wordwise
converse, not the possible existence of a different expanding ancestor.
A lower bound on the *whole* size tail is not by itself a necessary cost
bound for a residue-minus-tail test, because some of that same tail is
already present in the residue mass and cancels. Finally, segment
occupation must count visits with multiplicity if cycles are permitted.
The corrected note separates these points from its valid finite censuses.

These are useful route changes. SF1 avoids paying an ordinary-size tail to
transport a floor, and SF3 preserves the affine history until the final
quote. The remaining target is an independent positive input to one of
these interfaces, rather than improved normalization of a zero lower bound.

## Exact checks and remaining work

The [script](../../04-computation/experiments/collatz_source_refinement_floor_20261005.py)
uses only rational and integer arithmetic; its
[JSON](collatz_source_refinement_floor_20261005.json) records counts and hashes.
Normal and optimized Python runs agree on **250,069 exact checks**. The
repository documentation check and all eighteen local link targets pass.
Reproduce from the repository root:

    python3 04-computation/experiments/collatz_source_refinement_floor_20261005.py --json 05-knowledge/results/collatz_source_refinement_floor_20261005.json
    python3 -O 04-computation/experiments/collatz_source_refinement_floor_20261005.py

The universe includes the complete cost-eight boundary, all 10,922 odd unit
targets \(1<n<2^{15}\) with their 65,532 legal predecessors of valuation
at most twelve, the cost-six sharp-constant exceptions, twelve root-ray
sharp-exponent controls, 32 guarded controller instances, six nested residue
levels, and twelve full binary-tree levels. The polynomial identities in
(5)–(6) are proved algebraically above and also checked on a finite grid.
Hostiles include a deleted kernel record, the canceled-bill detour, the
deleted atom, and the rooted source excluded from the binary family.
Additional incoming controls check the neighboring-atom energy trap, the
contracting-slope rise \(165\to167\), cancellation between a residue and its
tail, and the visit-versus-set distinction on the negative two-cycle.

The next productive obligations are now more specific.

1. Find a source-dependent estimate of the energy increments in (14) whose
   sum is finite *without already using that source's ROOT word*. A finite
   positive Gram matrix, average resistance, or total mass bound does not do this.
2. Extend the proper-parent family by checked common-future controllers.
   Keep the total signed bill through route switches and use SF3 once.
   The extension must establish entry, not merely give another conditional
   transport inequality.
3. Seek complete cost boundaries beyond eight to support deeper inverse
   palettes. Completion must be proved by closure; a search timeout does not
   certify a finite kernel. Such an extension improves the transport library,
   not automatically its positive support.

The obtained source floors survive both meanings of refinement identified
above. Their present starting points cover specified infinite families.
An independent positive starting floor, or finite evaluation energy, for
every positive odd source remains OPEN.
