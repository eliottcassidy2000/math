# Complete positive selected-set tail seeds inside the triangle-chain family

**PROVED IN NOTE, specializing an inherited structural theorem.
FINITE-EXACT independent controls.** This classifies the selected-set,
transitive-tail constructor inside the K1/C3-SCC family at every order.
It does not classify arbitrary doubling constructions or prove universal
odd-slope seeds. The incoming dependency is
[Session258 CORE_CLASSIFICATION.md](../../session258_endpoint_seed_20261111/seed_decision/CORE_CLASSIFICATION.md),
merged from main a18596e598. Its general strong-core theorem is inherited;
the explicit family classification and marked enumeration below are this
packet's specialization.

## 1. Definitions and complete normal forms

For a tournament B and selected subset O, let K_ab count ordered pairs
of nonempty disjoint directed paths covering B, the first ending at a
and the second starting at b. Put

    W_B(O)=sum_(a,b in O) K_ab,
    U_B(O)=sum_(b in O) [s_b+sum_a K_ab],

where s_b counts full Hamilton paths starting at b. Thus U permits the
first path to be empty. A positive unit-offset tail seed means

    H(B)-W_B(O)=1,       U_B(O)>0.

The complete list among tournaments whose SCCs are K1 or C3 is:

1. `B=TT_a -> K1 -> TT_b`, selecting exactly that K1 vertex. Here
   `H=1, W=0, U=2^b`, for arbitrary a,b>=0.
2. `B=TT_a -> C3 -> TT_b`, selecting exactly two vertices of the C3 and
   no exterior vertex. Here `H=3, W=2, U=6*2^b`, for arbitrary a,b>=0.

There are no such seeds with two or more C3 components. Empty O, including
the empty source, is excluded by U>0. Allowing U=0 would add unselected
transitive sources and change the stated enumeration.

**Proof.** The inherited theorem says every positive unit seed has its
selected set in exactly one SCC, with every exterior SCC a singleton;
the selected core itself satisfies H-W=1, and b trailing singletons multiply
its U by2^b. Within this family the selected core is therefore either K1
or C3, and no other C3 may occur.

For K1, nonempty selection is the unique vertex; its H,W,U are1,0,1.
For a directed C3, K_ab=2 when a->b and0 otherwise, with zero diagonal.
Each vertex has s_b=1 and column sum2, so U_C3(O)=3|O|. The four selected
set sizes0,1,2,3 give W values0,0,2,6. Hence H-W=1 occurs precisely at
size2, with U=6. This proves both normal forms and their converse. The
inherited cross-component obstruction and exterior factorization supply
the all-order exclusion, not the bounded controls below. QED.

## 2. Exact multipliers, odd slopes and actual tails

The complete positive exact multiplier set in this family is

    {2^b:b>=0} union {6*2^b:b>=0}.

The possible odd parts, or odd slopes, are therefore exactly1 and3.
This does **not** say exact U=3 exists: slope3 starts at U=6, with one
factor2 already present. The only positive odd exact multiplier is U=1.

For fixed q the ordinary marked source generating function is

    z^(b+1)/(1-z), if q=2^b;
    z^(b+3)/(1-z), if q=6*2^b;
    0, otherwise.

Leading padding is arbitrary; the trailing length b is fixed by q.
The actual inherited tail constructor appends an unselected TT_t and an
apex that beats precisely O, with t>=1, and gives

    H(tail)=1+2^t U.

Thus the triangle-core seeds give `1+3*2^s` for s>=2, starting13,25;
they do not manufacture the forbidden value7. The singleton-core seed
starts3,5. These statements concern one common selected-set tail seed,
not unrelated realizations of the individual ray values.

## 3. Ordinary marked classes and labelled mass

The selected set is unordered. A transitive tournament is rigid, so
its selected vertex's position gives one class for each a,b. In a C3,
the three two-element subsets form one orbit under rotation, and the
stabilizer of a selected two-element subset is trivial. Consequently
there is one marked class for each a,b in each of the two normal forms.

The ordinary marked generating function is

    S(z)=z/(1-z)^2+z^3/(1-z)^2.

At order n>=1 its coefficient is `n+max(n-2,0)`:1,2,4,6,8,10,... .
Every positive marked seed is rigid, including the triangle-core case:
marking two vertices destroys its rotational automorphisms. Each ordinary
marked class therefore has exactly n! labelings. Its labelled **exponential**
generating function is the same S(z), while the actual labelled count
is n! times the displayed ordinary coefficient. This equality of generating
functions does not identify ordinary counts with labelled counts.

The unmarked family had the automorphism divisor3^r. It must not be
carried into this marked enumeration: the selected C3 subset has stabilizer1.

## 4. Constructor compatibility is stronger than unit extension content

The source C3->C3 has one-birth extension content1, yet no selected set
makes it a positive unit-offset tail seed. Indeed every source with at
least two triangles is excluded above, even if its extension-content
filtration reaches1 immediately. Coprime extension counts do not imply
the simultaneous identities H-W=1 and a specified U.

This also explains why the newly classified family cannot on its own
provide the incoming slope5 or slope29 seeds: those use different strong
five-vertex cores in Session258. Their existence is not contradicted by
this family restriction. The distinction is between actual source and
constructor interfaces, not a failure of ordinary primality or congruence.

## 5. Independent exact controls

Run `python -B TRIANGLE_seeds.py` and `python -O -B TRIANGLE_seeds.py`.
The checker builds every K1/C3 component word through order8 and tests
every selected subset. It obtains H, all endpoint counts and every K_ab
directly from permutations and split positions, independently of the
incoming endpoint-DP code. It checks the complete classification, exact
U law, ordinary orbit counts, labelled masses, empty-source convention,
and six literal one-/two-tail graphs. The receipt includes the inherited
theorem hash. These finite controls verify the specialization; the
structural theorem proves the unbounded quantifiers.
