# The two cyclic-triangle endcaps repair the h=3 backward branch

2026-10-25. **PROVED:** for n>=12, I={3,n-3,n-2}, and 3->1, the marked
closed blue P has a literal three-black detour. This closes one more
endpoint-four branch. It does not prove unrestricted black-core connectivity.
The [main packet](README.md) and its checker remain unchanged.

## Inheritance and exact carrier

The [Klein/closed-blue dictionary](../../black_detours_closed_blue_lines_klein_dictionary_20261020.md)
distinguishes a literal replacement for a marked path from merely finding
another Hamiltonian presentation. The
[rotation-symmetric endcap repair](../../rotation_owners_spectrum_20261020e/black.md)
retains arbitrary middle blocks but starts at endcap order five. Its
[regular companion](../../rotation_owners_spectrum_20261020e/regular.md) explicitly
does not supply an order-three macro. The carrier here also has nonuniform
incidence from the end modules to the middle, so it is not an ordered sum
of two outer strong components and the middle. Those earlier theorems do
not directly prove the present result.

The closest mechanism is the new backward ledger (D) in [README.md](README.md).
For h=3 it is legal and its first two moves remain black at n>=12, but its
last word can be blue. This is the canonical hostile. The corrected near
miss is treating the two triangles as ordinary outer SCC endcaps. The
less-used sidecar is the full permutation action of the blue final word
on the literal middle vertices. The board is marked target / retained
middle interval / two triangle modules / first-chord split / exact XOR.

Keep a=0,b=1,x=2,y=n-3,c=n-2,z=n-1, P=(0,...,z), and B=K_P(T).
The source has I={3,y,c} and 3->b. Put

    D3=(x,a,b,3,4,...,z).

Its reversal product g=r_D3 r_P is

    g=(a,x,b)(z,c,y),                                  (1)

fixing every vertex 3,...,n-4. Since P is blue in B, D3 is blue exactly
when g is an automorphism of B. Under the specified endpoint row this
has an exact carrier: L={a,b,x} and R={y,c,z} are cyclic-triangle modules;
L points to R; L points to middle vertex 3 and is beaten by every other
middle vertex; R is beaten by middle vertex n-4 and points to every
other middle vertex. The interval (3,...,n-4) remains an arbitrary blue
Hamiltonian presentation with its actual internal arcs retained.

For necessity, the first 3-cycle in (1) forces the rows of L to agree
at every fixed middle vertex. The known B(a) row determines their common
values. B(a) points to all three vertices of R, so the simultaneous two
rotations force every L-to-R arc forward. Reflection by r_P gives the
stated R-to-middle pattern. The original path edges make both induced
triangles cyclic. Conversely, the displayed carrier is preserved by (1),
so it indeed makes D3 blue. No quotient class or path-count statistic
has replaced the actual middle.

In particular this carrier forces B(b,4)=B(a,4)=0, hence T(b,4)=1.
This is the small extra arc needed for the following repair.

## A complete chord split, without testing for blue at runtime

If 4->b, use the original ledger (D). Its first two moves are black by
the proved k=4 degree/antipodal argument in the main packet. Its final
word is also black: were D3 blue, (1) would make the rows of a and b
equal at the fixed vertex 4. But

    B(a,4)=0, B(b,4)=1.

Thus the existing three-move ledger completes this branch.

If b->4, write W=(4,...,n-4), Rword=(y,c,z), and use

    E1=(a,b,W,Rword,x,3),
    E2=(b,W,a,Rword,x,3),
    E3=(b,x,3,W,a,Rword).                             (2)

This replacement was found by a six-block permutation menu on
`(a),(b),(x),(3),(W),(Rword)` under the exact carrier equations. Only
the resulting three words are part of the proof; no all-Hamiltonian-path
census or bound on possible middles is used.

E1 uses the new chord b->4 and the retained endpoint port z->x. Its
other joins are path edges. E2 inserts a after n-4 and before y: the
needed joins are the complements of a->(n-4) and y->a. E3 restores
b->x and 3->4 after two complements and retains the remaining joins.
The nonpath pairs b4,zx,(n-4)a,ay occur twice; every other edge has
odd multiplicity precisely when it belongs to P. Therefore

    E(E1) xor E(E2) xor E(E3)=E(P).                    (3)

All three black witnesses are local:

* E1 at T would equate T(a,x)=1 with T(b,3)=0.
* E2 at S=K_E1(T) would equate S(b,x)=0 with S(4,3)=1.
  Both original path edges bx and 34 were omitted by E1 and complemented.
* E3 at B would equate B(a,z)=1 with B(b,4)=0. The first complements
  the original closing edge z->a; the second complements the branch
  chord b->4. Both pairs are nonpath pairs of P.

Their positional reversals give exactly these reflected ordered pairs.
Thus (2) is a legal three-black detour with literal target B. It repairs
every forced two-C3 carrier and, more strongly, every source in this chord
branch regardless of whether the old final word was blue. The two chord
directions exhaust the source, proving the stated theorem.

## Controls and surviving frontier

`endcap_check.py` quantifies every reflected completion in both chord cells
for n=12..60, with additional orders 80,100,140. It checks the legal joins,
actual blue predicates at each intermediate, exact edge target, and two
literal arc-set replays per cell. A separate conditioned source check
verifies the exact two-C3 carrier, its forced b->4 chord, and the new repair.
The opposite b4 chord is infeasible under the carrier, which is retained
as a hostile to forgetting the module constraint. Run the checker normally
and with `-O`.

After this addendum the large-order endpoint-four source cells not handled
by this packet are h=3 with 1->3, and h=4 with 3->1, subject to the
inherited frontier guards. The scope is a local detour theorem, not a
claim that either residual is disconnected or that the global conjecture
has been resolved.
