# Block reserves for the moving positive-part ledger

October 21, 2026. **PROVED:** exact native debit and associative excursion
summaries. **FINITE-EXACT:** an actual N4096 one-step obstruction to local
payment. **OPEN:** uniform block payment, uniform153N, and arbitrary entry.

The closest mechanism is the source-owned
[native annulus ledger](../session18_native_annulus/README.md), with the
[synchronous pending blocks](../future_tail_blocks_20261021b/README.md).
The canonical hostile remains N7 at cut8. The corrected near miss here is
to expect a positive-part moving potential to decrease at every native
step. The least-used sidecar is the maximum *prefix* debit inside a block,
which a block's final net debit does not retain.

Live concepts: original source norm; synchronous node; linked old/new
charge; current terminal; block net debit; internal reserve. Original
sources, mod15/mod35, integer clocks and the fixed frontier are unchanged.
There is no new terminal and no assumed unresolved arrival.

## 1. An exact stronger potential and its debit

Use Q(x,y)=x²−xy+y² and
δ(w,h)=Q(w+h)−Q(w). At a finite native cut D, let c_t be the frozen
old/new four-vector at terminal t and b_u the coalesced four-vector at
each actual pending node u. Define

\[
 \Psi_D=\Delta_D+\sum_{u\text{ pending}}\max(0,\delta(b_u)).       \tag{1}
\]

Then Ψ_D≥Δ_D and Ψ_0=153N. The latter is an exact equality: initially
each nonzero source is at its own node, each old source has δ=0, and
each new source contributes its original norm Q(h). Frontier sources
contribute the same quantity to Δ_0. Zero charges do not change the sum.

At the next native step let F_z be the set of pending nodes whose
actual successor is z, and put v_z=Σ_{u∈F_z}σ_u b_u, where σ_u is the
parity of that edge's *actual integer clock*. The exact debit at z is

\[
 d_z=\begin{cases}
  \delta(c_z+v_z)-\delta(c_z)-\sum_{u\in F_z}\delta(b_u)_+,&z\in V,\\
  \delta(v_z)_+-\sum_{u\in F_z}\delta(b_u)_+,&z\notin V.
 \end{cases}                                                \tag{2}
\]

Here x_+=max(0,x). Partitioning the pending sources by their actual next
node proves

\[
 \Psi_{D+1}-\Psi_D=\sum_z d_z.                              \tag{3}
\]

Equation(3) uses each block exactly once. No current source norm is paid
again when a block coalesces. A dual-zero block may be erased; a block
with w+h=0 alone may not. This is a full moving potential, stronger than
the physical first-arrival energy. It is not the already-proved all-future
envelope U_D, and it is not itself an all-future certificate at a single
cut.

## 2. Pointwise payment already fails in the intended range

At N4096, cut17, the current terminal811 has

\[
 c_{811}=(-6,-10,7,-36).
\]

The two pending nodes1081 and2163 have respectively

\[
 b_{1081}=(-2,-2,-9,-25),\qquad b_{2163}=(0,0,-1,0).
\]

Their actual edges are1081→811 with clock2 and2163→811 with clock3.
Their incoming vector is therefore b1081−b2163. Substitution into(2)
gives d811=2004>0. This is a source-specific native example, not an
arbitrarily chosen four-vector hostile.

The complete direct original-source replay gives

| Cut | Δ | Positive moving part | Ψ |
|---:|---:|---:|---:|
|17|166652|268827|435479|
|18|184903|254791|439694|

Thus the *global* one-step change is also positive:4215. Both cuts remain
well below626688=153·4096. The failure is of local monotonicity; it is
not a counterexample to the uniform budget.

The original discovery run followed all cuts0..150 at this N, with maximum
Ψ equal to its initial626688. It found further positive steps at23 and27.
That range is a finite exploratory observation. The independent checker
authenticates exactly cuts17 and18 by routing each original nonzero source
separately, without discovery's coalescence algorithm or a route suffix.

The complete source universe is the393216 primitive odd s<256·4096,
excluding s≡5 mod8. The checker retains joint register hashes, full
integer clocks, representative literal clock words, and counts/ranges
for every source participating in the displayed local fibres. It checks
the old named-charge equality and the exact annulus budget. A second
grouped calculation of(2) agrees with the direct two-cut result.

## 3. Larger blocks need a reserve, not merely a nonpositive endpoint

For any finite consecutive list of exact debits e1,...,er, retain the
two-coordinate summary

\[
 T=\sum_i e_i,\qquad
 P=\max\{0,e_1,e_1+e_2,\ldots,e_1+\cdots+e_r\}.             \tag{4}
\]

T is the net debit and P the minimum starting reserve that pays every
intermediate cut. Concatenation has the exact rule

\[
 (T_1,P_1)\star(T_2,P_2)
 =(T_1+T_2,\max(P_1,T_1+P_2)).                             \tag{5}
\]

This rule is associative because it is the summary of literal concatenation;
it also follows by expanding both parenthesizations. The identity is(0,0).
It preserves precisely the threshold predicate “all prefix debits≤B.”
It loses the positions and attribution of individual debits, so the native
source/clock ledger remains the authentication sidecar.

For example,(4,−5) and(−5,4) both have net−1, but require reserves4 and0.
Consequently proving only that every chosen block has nonpositive net
debit would not control the intermediate horizons. Equation(5) permits
blocks to pay together while keeping that missing peak coordinate.

**Conditional uniform route.** Partition every native trajectory of the
ledger into finite consecutive depth blocks whose endpoints are unbounded.
At the beginning d of a block, let B=153N−Ψ_d. If its authenticated
summary satisfies P≤B, every cut in that block obeys Δ≤Ψ≤153N and its
next reserve is exactly B−T. Applying this successively proves the bound
at every finite horizon. It does not require individual source arrival.
The existence of such uniformly payable blocks for all N≥4096 is OPEN.
The summary algebra is a way to state and check that obligation, not a
proof that it holds. A stopped-tail certificate may terminate this block
process as soon as one of the authenticated all-future gates succeeds.

This route is stronger than necessary: Ψ might exceed153N while Δ does
not. The angular stopping gate remains available if this moving-potential
route encounters such a counterexample.

## 4. Reproduce and scope

```
python -B 05-knowledge/results/uniform_prefix_20261021c/probe_potential.py 4096 --depth 150
python -B 05-knowledge/results/uniform_prefix_20261021c/check.py
python -B -O 05-knowledge/results/uniform_prefix_20261021c/check.py
```

The saved consumer receipt has406480 body checks and LF SHA256
`75694fcf2e9d37b03464e94d6e014f29dd15541346d44562db3caed25d22878e`.
The finite algebra controls exhaust all debit words on{−2,0,3} of lengths
0..6 and every split. The displayed native obstruction is independently
replayed from all source addresses. Neither the constant153 nor the
physical energy has been altered. The next proof obligation is a bound
on authenticated block excursions, or a uniform stopping-gate hypothesis
such as the linked coalescence condition in the
[tail-geometry packet](../tail_geometry_20261021c/README.md).

Independent full proof/consumer review and optimized replay by collatz_memory,
October21,2026: PASS, with the same receipt and406481 checks including equality.
