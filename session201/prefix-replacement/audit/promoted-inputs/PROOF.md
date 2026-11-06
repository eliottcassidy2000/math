# Constructive signed-fiber repair by maximal prefix retention

Session 201. **PROVED, independently audited in the stated scope.** An exact target-dependent
algorithm replaces a terminal segment of a supplied physical source recipe,
retaining the longest possible prefix, whenever the target has a physical
C1 section witness. Two infinite families previously hostile to every
append-only repair now have explicit replacements, with optimal retained
prefixes and optimal replacement lengths.

**OPEN:** all-odd C1/C2/C3 coverage. The algorithm decides and repairs actual
C1 membership; it does not assume that every odd target belongs to C1.
The meaningful new operation is replacement of a supplied witness, rather
than extending it indefinitely or substituting arbitrary cone coordinates.

## 1. Physical recipes and the exact inherited criterion

Use the actual nontrivial source orbit S from t=(1,0,2,1), under the
chronological O/Z/J operations in the frozen dependencies. Put

    f=U-D+E-d, L=6U-7D,
    T={actual q: E_q=2(U_q-D_q), d_q=3(U_q-D_q)}.

For an actual source s let Desc(s) include all its word extensions and itself.
The audited append-only classification gives

    A(s):={L(q):q in Desc(s) intersect T}
       = {D_v:v in Desc(s)}
         union ({d_s} if f_s=0)
         union {L(s.J^j):0<=j<=2, s.J^j in T}.                 (1)

It also proves that if f_s>0, A(s) is exactly the first set. The physical
compiler W(v)=v.OJJJJ lies in T and has target D_v; if f_s=0, s.ZJJJ lies
in T and has target d_s. The three short J completions are literal moves.
All section targets are positive odd integers. The full section image is
{1} union the C1 spectrum, not the full conjectural C1/C2/C3 union.

The new algorithm operates on a supplied chronological presentation. It
optimizes the number of terminal letters discarded from that presentation,
not arbitrary vertex deletions, root changes, graph isomorphisms or edge
reorientations.

## 2. A target-dependent finite table, with constructive witnesses

Fix a positive odd target n. If n>1, enumerate the actual source DAG

    S_n={s in S: U_s<=2n-3 and D_s<=n}.                       (2)

Every predecessor of a retained node is retained because U grows strictly
and D does not decrease. Start at t; if t exceeds the cutoff, the inventory
is empty. Mark all nodes with D=n as successful. Process remaining nodes in
decreasing U, marking s successful if at least one of O(s),Z(s),J(s) is a
successful retained node. Keep an actual successful suffix, optionally the
shortest such suffix inside this finite table.

For any actual queried prefix s, membership of n in D(Desc(s)) is decided
exactly as follows:

* If n<D_s, it fails.
* If n=D_s, it succeeds with the empty suffix.
* If n>D_s and U_s>2n-3, it fails.
* Otherwise it succeeds exactly when s is marked in (2).

For n=1 the equality test is sufficient, since all D-values are positive
odd and nondecreasing.

**Proof of completeness.** Suppose n>D_s occurs at a descendant. Strip
only the terminal O letters from its supplied extension. At least one new
Z/J remains, because D changed. Its positive-endpoint core is still a
descendant, has D=n, and has U<=2n-3. For Z this is
2D-U>=3d>=3. For J it follows from

    2D_new-U_new=D_old-2E_old+U_old>=3,

using D_old>=2E_old and U_old>=3, with the seed U_old=2 checked separately.
Every node on the retained path has no greater U and no greater D. Conversely
every marked path is an actual physical suffix. This is a finite search of
physical states, not a relaxation to their coordinate inequalities. QED.

Combine this decision with the deletion and three short-completion tests
in (1). On success the table immediately produces a physical suffix eta:

* a path to D=n, followed by OJJJJ;
* ZJJJ when f_s=0 and d_s=n;
* a successful J^j, j<=2.

Thus the test is constructive. A shortest path to D=n is not claimed to
give the shortest possible full section suffix. When that additional
optimization is desired, breadth-first search by word length can stop at
the known constructed suffix length; the search is finite and complete for
all shorter suffixes, without imposing an unjustified cover-height cutoff.

## 3. Longest retained prefix and optimal terminal rollback

Let w=c_1...c_l be the supplied physical word, with prefixes
s_i=t.c_1...c_i, including s_0=t. The requested target can be L(s_l)>0,
as in signed-fiber repair, or any specified positive odd n.

**Theorem 1 (maximal-prefix repair).** The feasible indices

    I_n(w)={i: n belongs to A(s_i)}

form an initial interval, possibly empty. If nonempty, let i_* be its largest
index. The test in Section2 constructs eta such that

    q=t.c_1...c_(i_*).eta belongs to T, L(q)=n.               (3)

Among all repairs formed by deleting a terminal segment of this supplied
word and appending fresh physical letters, exactly l-i_* is the minimum
number of deleted letters. Feasibility can be queried by binary search in
the prefix list, using at most one initial query and a logarithmic number
of further queries after the single target-dependent table is constructed.

**Proof.** Desc(s_(i+1)) is contained in Desc(s_i). Therefore
A(s_(i+1)) is contained in A(s_i), which proves the interval statement.
Each query is exact by (1) and Section2, and its successful suffix gives
(3). No larger index can be retained by any append-only completion, by the
definition of i_*. A smaller number of deleted terminal letters would retain
a larger index, proving optimality. QED.

The root query succeeds exactly when n is1 or an actual C1 count. On success,
the audited section inverse recovers an actual unsigned source with D=n.
For n>1, the last Z/J birth creating D=n also supplies a strong C1 graph:
close its parent by J and omit any trailing O letters. On failure, the
algorithm reports no C1 section repair; it does not report that n has no
C2 or C3 realization.

The theorem is conditional on actual target membership and makes that
condition decidable with a complete finite table. It is not an occupancy
proof. It supplies a canonical optimal cut in the given recipe even though
the original signed representative may be arbitrarily large or append-hostile.

## 4. Two infinite hostile families repaired optimally

Write b=(3,1,3,1) and v_r=t.Z.(OJ)^r=3^r b. Set lambda=3^r.
All equalities below are actual physical words, and their dependence on
lambda follows from the linear source operators. This is not division of
an arbitrary physical profile by its coordinate gcd.

### Target 3lambda: discard one letter and append two

The old witness

    s_r=v_r.OZ=lambda*(9,1,11,1)

has L=3lambda but admits no target-preserving append-only section repair.
Retain the prefix v_r.O and replace its last Z by JJ:

    q_r=v_r.OJJ=lambda*(15,6,18,9), L(q_r)=3lambda.           (4)

The two section equations hold exactly. Zero deletions are impossible by
the audited append obstruction, so one deletion is optimal. From the retained
prefix v_r.O no zero- or one-letter suffix reaches T: O has E=0; Z has E=d;
and J is balanced with positive E. Each is incompatible with the section
at positive target. Thus the replacement length two is also optimal.

### Target 33lambda: discard three letters and append seven

The old witness

    u_r=v_r.OOZ=lambda*(15,1,23,1)

has L=33lambda. Retain v_r and replace OOZ by ZZZZOJJ:

    p_r=v_r.ZZZZOJJ=lambda*(165,66,198,99),
    L(p_r)=33lambda.                                        (5)

To prove that no longer prefix can be retained, it suffices to exclude target
33 from A(b.O). This source has f>0, so its section-target set is exactly
its descendant D-values. The complete descendant core search through U=63
finds, for D<=33, precisely

    {3,9,15,19,27,31}.                                      (6)

The U-only inventory has55 profiles, of which34 also have D<=33. Its
completeness follows from Section2, so the
absence of33 is an all-height conclusion, not a guessed depth cutoff.
Linearity gives the corresponding exclusion of33lambda after v_r.O.
Every later prefix is even more restricted. Thus all three terminal letters
OOZ must be discarded, and (5) is an optimal retained-prefix repair.

The appended length seven is optimal for this retained prefix. Exhausting
all1+3+...+3^6=1093 words of length at most six from b finds no section output
with target33; ZZZZOJJ is a length-seven witness. This finite word check
settles a genuinely finite minimum-length question. Linearity transfers it
to every v_r. It is not used to prove all-height absence.

The repairs also give ordinary unsigned C1 witnesses: the balanced sources
v_r and v_r.ZZZZ have respectively D=3lambda and D=33lambda. The point of
(4)-(5) is stronger than those previously known counts: they optimally repair
the supplied bad signed representatives while retaining as much of each
chronological construction as possible.

## 5. Verification and the remaining existence step

The implementation keeps actual words, builds the complete target table,
and computes the optimal retained prefix. Declared tests compare binary-search
cuts against every prefix, reconstruct each returned section witness, retain
exact no-C1 controls such as7, and check the source-to-C1 conversion. The two
infinite repair identities are checked on many r; the all-order statements
come from linearity and the proofs above. Literal small tournament counts
independently verify the repaired profiles and the unsigned counts3 and33.
Normal and optimized Python receipts must agree.

The remaining all-target lemma is still nonemptiness of the relevant physical
C1/C2/C3 target fibers. This theorem does not create a successful root query
when its target is absent from C1. It turns a successful target into an exact,
optimal-prefix replacement of any supplied signed witness, and demonstrates
that replacement can resolve the earlier infinite append-only obstructions.
