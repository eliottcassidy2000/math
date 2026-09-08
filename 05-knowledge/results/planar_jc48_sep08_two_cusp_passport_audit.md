# Independent audit of the uniform two-cusp actual-sheet passport

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE + REPLAY AUDIT PASS.**
The [primary note](planar_jc48_sep08_two_cusp_passport.md) is a necessary
passport for the whole irreducible nonproperness support of a complex
planar Keller map. The degree-six example is an abstract finite-group
passport, not a realization of the actual curve group or a Keller map.

## 1. The local subsets are actual, and the equality sidecar is essential

I recovered the marked re-access and pure-deleted-boundary argument in
[the proved cusp passport §§2–3](planar_jc48_sep06_cusp_passport.md).
The relation g sigma g^-1=tau for g=sigma tau transports the retained
subset together with its meridian. The intersection of the resulting
actual subsets counts the cusp fibre: singleton local-group orbits
extend over the normal finite envelope, and pure deleted boundary
prevents removing only the cusp point. This does not identify all
inertia-fixed labels with retained labels. The one-cusp Euler conclusion
n<=1 is not imported into the two-cusp argument.

For the uniform lemma, if x and tau x both lie in U=A\B, the braid
relation forces tau x=tau^2 x. Then both generators fix x, so x lies in
gA=B, a contradiction. Thus tau U is contained in the outer block O,
giving 2n>=3a-d in every degree. At equality the bijection exchanges U
and O. Applying the braid again makes sigma exchange O and V. Hence
every nontrivial sigma-cycle is even and its whole fixed set is A.
This proves exactly the strengthened sidecar claimed; it does not
force transpositions unless both exchanged blocks have size one.

At a node, global conjugacy gives the same full fixed count a. The two
actual retained sets must therefore be the full fixed sets. Since the
node meridians commute, their support intersection is invariant under
either and is a union of nontrivial cycles. Its size is even under the
equality sidecar. Without that sidecar, deleted inertia-fixed letters
would invalidate this identification; the proof keeps this distinction.

## 2. The all-degree Euler table, including N=2

Normalization A1 with two cusps and N nodes has smooth-stratum Euler
characteristic -1-2N; the ambient complement has characteristic N.
Integrating actual counts gives 1=-a+n_1+n_2+W, W=sum omega_p.
The two cusp inequalities and the node lower bound give
N max(q,0)<=W<=q+1 for q=d-2a. Thus q>=-1. For N>=3 the positive
case is impossible. For N=2 the remaining scalar case q=1 has W=2
and each overlap one. Both cusps attain equality, so node parity
contradicts those individual overlaps. The explicit scalar hostile
d=3,a=1,n_1=n_2=0,W=2 confirms that the extra step is necessary.

For q=-1, W=0 and both cusp inequalities are equalities, forcing odd a.
For q=0, odd a gives two equal cusp counts and W=0; even a initially
also permits W=1 with both local inequalities equalities. Node parity
excludes that last case. The three displayed rows and all zero node
overlaps follow, with no finite degree cutoff. The node fibre is exactly
2a-d, explaining the respective zero/one entries.

## 3. Generator bounds and exact abstract controls

The generic-degree two/three/four exclusions remain the explicitly
cited, correctly typed inputs in the inherited note. They are not
deduced from the local passport: the degree-four S4 control passes all
local, Euler and three-generator conditions. At degree five local
equality forces a single transposition. A transitive group generated
by three such permutations would require a connected graph on five
vertices with only three edges, which is impossible.

At degree six, the same table gives a=3 and both cusp counts two. A
nontrivial meridian fixes at least three labels, leaving only a single
transposition or a three-cycle. The transposition case again fails the
three-generator connectivity bound. The cited connected-support
three-cycle lemma then identifies the transitive image with A6.
The explicit control indeed generates A6, has both correct transported
intersections, and has disjoint commuting node supports. It retains the
positive-meridian conjugacy class in one common group. Its omitted data
are the actual simultaneous access relations; the proof does not treat
this abstract existence as a representation of a specified complement.

## 4. Independent reproduction and a broader node check

Both normal and optimized replays match the frozen 674-byte output and
pass 1,051 always-active gates. I read the full source, including the
literal left-action convention, permutation products, exhaustive subset
loops, the scalar table filter, and complete generated-group searches.
The local universes on two, three and four labels have 2,15,136 rows.
The scalar universes through degree sixteen at N=2,3,4,5 each retain
the expected 15 ordered rows. The six-label four-cycle hostile prevents
silently strengthening equality to transpositions.

The producer's small commuting-pair check is inside a cusp-braid loop;
that does not by itself exhaust all commuting node pairs. The general
node-parity argument is analytic. As an additional independent control,
I enumerated every commuting permutation pair on d=2,3,4,5 with the
first permutation having only even nontrivial cycles, with no cusp-braid
filter. All 4,12,96,480 respective pairs have even support intersection.
These extra audit controls are not added to the producer's gate total.

    python3 -B 04-computation/planar_jc48_sep08_two_cusp_passport.py
    python3 -B -O 04-computation/planar_jc48_sep08_two_cusp_passport.py

Source SHA256: `c2794aac802fef9e0bebe5d47adf63bf49610de6bd1ca11b250d5e238053cee1`.

Output SHA256: `0bd7553781b01ddc2031fb4f4e7050404555e779ff36135fb3c4d53dc7e97880`.

The primary may be promoted to PROVED for the declared whole-support
hypothesis and necessary conditions. JC(2) remains OPEN.
