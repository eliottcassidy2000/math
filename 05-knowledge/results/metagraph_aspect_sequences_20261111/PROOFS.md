# Exact response, symmetry, and the minimum seed order

Status: PROVED elementary identities; FINITE-EXACT census consequences use the
explicitly pinned complete class banks. Arbitrary odd-c and Chan–Pak are OPEN.

## 1. One computation for every one-vertex attachment

For a nonempty base B, let s_j and e_j count its spanning paths starting and
ending at j. Let C_ij count ordered covers of all vertices by two nonempty
directed paths, with the first ending at i and the second starting at j.
Thus C_ii=0. These are counts of actual paths, not walks or quotient orbits.

For O subset V(B), add v with v->O and B\O->v. Define

    U(O)=sum_(j in O) (s_j + sum_i C_ij),
    W(O)=sum_(i,j in O) C_ij,
    E(O)=sum_(j in O) e_j.

Classifying a Hamiltonian path by whether v starts, ends, or lies internally
gives the exact identity

    H(B+v)=H(B)+U(O)-E(O)-W(O).                  (1)

The interior case contributes sum_(i outside O,j in O) C_ij. The final-v case
contributes H-E(O), and the initial-v case contributes sum_(j in O) s_j.
This proves (1) with no assumptions of regularity, rigidity, or strongness.

Let end[S,i], start[S,j] be the usual path counts on the induced subset S.
Then

    C_ij=sum_(nonempty proper S) end[S,i] start[V\S,j].

The independent audit uses the equivalent scalar convolution
`u_j=sum_(nonempty S) start[S,j] H(B[V\S])`, with H(empty)=1.
Computing both subset tables and the C matrix takes O(n^2 2^n) arithmetic
operations. With O=R union {j},

    U(O)=U(R)+u_j,
    E(O)=E(R)+e_j,
    W(O)=W(R)+sum_(i in R)(C_ij+C_ji).

All 2^n actual attachments can therefore share one computation. Arithmetic
bit complexity also depends on the size of the resulting integer counts.
The independent audit traverses the subset cube in Gray order instead.

Converse sends attachment O to the complementary attachment in the reversed
base. Hence the multiset of child H values, and its support, is converse
invariant. It is a finite observation here that the support separates every
converse class through nine; the identity supplies no all-order injectivity.

The parent H need not be stored separately: it is the minimum of the support.
Every Hamiltonian path of B admits an insertion position for the new vertex
(at an end, or at the transition from a vertex beating v to one beaten by v).
Choose its first legal position to get an injection into the child paths.
Thus each child has at least H(B) paths, and a source/sink attachment attains
equality. The n<=8 registry includes H in its signature only for readability.

## 2. The seed observer retains direction

The inherited THM-5126,
[`THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md`](../../../01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md),
adds B -> TT_t -> v, t>=1, with the same O at v. It gives

    H=H(B)-W(O)+2^t U(O).

A positive unit selection satisfies W=H(B)-1 and U>0. If U=c 2^j with c odd,
it generates 1+c2^s for every s>=j+1. The group Aut(B,O) acts freely on all
ordered covering paths, so its order divides H(B)-W(O)=1. Therefore its
stabilizer is trivial; the number of such selected subsets divided by
|Aut(B)| counts their orbits exactly. This restriction is inherited from
the preceding odd-c packet, not claimed as a new symmetry theorem here.

Under converse W is unchanged, but U becomes

    U_dual(O)=sum_(j in O) (e_j+sum_i C_ji).

The number of unit subsets is preserved; their coefficient/exponent labels
need not be. The n6 H15 converse pair in `hostiles.json` has coefficient sets
{29,37,43} and {1,17,55}. This is why the order-nine producer computes BOTH
U directions and the independent audit visits EVERY ordinary source.
Simply doubling a converse-representative seed count cannot recover the
coefficient support. This prepublication correction is preserved in the board
and mistakes ledger.

Counting convention matters independently of converse. The earlier S21
`unitsets.c` reported19,002 selections of sizes2–4 at base order nine; our
19,011 positive selections additionally include the nine singleton selections
of the transitive base. Dividing each source's selections by its automorphism
order gives18,849 rooted isomorphism types. S21's separate `seeds.c` list counts
individual vertices in chosen root-tournament representatives and also includes
zero-slope offset-one roots. Its sequence2,5,7,16,31,68,207,1333 is therefore
not our positive orbit sequence1,2,4,9,22,59,194,1300 (base orders1–8).
Do not substitute our sequence into HYP-9658(D) without changing its definition.
The5,529 coefficient support and the33 exclusions were already complete in
S21's audit; the present all-subset scans independently replay that boundary.

## 3. Source padding preserves a unit seed exactly

Add a new source a to B and leave O unchanged (a is not selected). Every
spanning path starts with a, so H(a->B)=H(B).

In a W-cover, a must be first in the first path: it cannot be internal, and
placing it first in the second path would violate the selected-start condition.
Deleting a is a bijection to the old W-covers; its first path remains nonempty
because its endpoint lies in O. Thus W is unchanged.

In a U-cover the same deletion leaves either a nonempty first path or the
allowed empty first path. Conversely, prepend a to that first path, including
the empty one. A cover with empty first path in the padded base would start
its second path at a outside O and is not counted. Thus U is also unchanged.

Consequently every positive unit seed on at most nine base vertices has a
positive unit seed on exactly nine base vertices with identical c and exponent.
Both complete scans find exactly the same 33 missing coefficients <=4097;
`audit9.json` lists them. The preceding packet constructs a ten-vertex base
for every one of them. Minimum base order is therefore ten, and minimum root
order (base plus v, before the nonempty tail) is eleven, for each coefficient.
This proves optimality only within this precisely specified seed family.

## 4. Source-coherence defect is an actual sequence

At fixed n and h, let E_T be the set of one-vertex child counts and let
T range over ordinary n-vertex classes with H(T)=h. Define

    D(n,h)=# { {a,b}: a!=b, a,b in union_T E_T,
                       and no T has {a,b} subset E_T }.

The registry records sum_h D(n,h) and the number of h with D(n,h)>0.
Each value is invariant under relabeling and under replacing a source by its
converse. This uses an honest compatibility graph on target values; ties are
not arbitrarily oriented. It tests exactly the inference that union membership
supplies a common source. It says nothing about sources of another order.

The first failure is n=5,h=3,targets17 and29. The two exhibited attachment
witnesses are recounted in the final validation. The complete source list
proves no common n5 source. The all-height 105 obstruction is a distinct
inherited result; the present finite sequence does not replace its proof.
