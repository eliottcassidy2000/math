# The sum-of-two-squares endpoint has only two physical prime values

Session 199. **PROVED:** the physical source fixed points of the swap in the
quadratic norm x^2+4yz are classified exactly. Their diagonal C2 values are
13 and 41h^2, where h ranges over the Hamilton-path counts of balanced sources,
including the singleton. This is a restriction on that norm-fixed-point route,
not on arbitrary diagonal C2 values, independent source pairs, or the full
Hamilton-path spectrum.

This addendum depends on the unchanged, independently audited
`../source-norm-capacity/PROOF.md`, SHA256
`2a75ad55f56377cd430d16a05e4ca3be3dceb20c0bd70f0a17b14c892fe90508`.
It uses the physical O/Z/J grammar and the norm coordinates established there.
A nontrivial source has positive coordinates v=(a,f,e,p) and

    D=a+2e+p, E=e, U=a+f+2e+2p, d=e+p.

A balanced source means D=U. In this cone that is equivalent to f=p=0.
The singleton is separately allowed as a balanced input in the theorem below.
Words are chronological: JO(B) in operator notation means first O, then J.

## 1. Three exact grammar observations

Every nontrivial source has e+p>=1.

1. A last O has f'=a+2f+e+2p>0 and e'=0.
2. A last Z has a'=2a+f+e+p>=1, e'=e+p>=1, and
   a'-f'=2a+e+p>0.
3. A last J has f'=0, e'=a+e+p>=1, and a'=f+p>=0.

The seed b=(0,0,0,1) has e=0. These facts imply that every source with e>0
satisfies f<=a, and every nonseed source with e=0 ends with O in any chosen
complete source-word presentation. They also imply that a source with a=f=0
is either b or has a presentation J(B) with B balanced: its final letter
cannot be O or Z, and J's new a=f_old+p_old vanishes only when both vanish.
These are graph-construction statements from a complete physical word, not
inferences that all numerical solutions are sources.

## 2. Complete classification of the norm swap's physical fixed points

The norm coordinates from the preceding packet are

    x=D-E+U, y=(U+E)/2, z=U-d,
    N(v,v)=x^2+4yz.

The norm-preserving swap is (x,y,z)->(x,z,y). Its fixed-point condition is

    y=z  iff U-E-2d=0  iff a+f=e.                       (1)

**Theorem.** A nontrivial physical source satisfying (1) is either the
transitive edge or has a physical presentation J(O(B)), where B is a balanced
source and may be the singleton. Conversely all these sources satisfy (1).

Proof. The seed b is the edge and satisfies (1). Otherwise choose a complete
O/Z/J presentation and inspect its last letter.

If it is O, (1) would require 0=a'+f'=2(a+f+e+p), impossible.

If it is Z, (1) becomes e+p=2a+2f+e+p, hence a=f=0 in the predecessor.
By Section 1 that predecessor is b or J(B) for a balanced nontrivial B.
The first case gives the directed triangle J(O(singleton)). The second gives
Z(J(B))=J(O(B)). The equality here is literal rooted graph isomorphism: both
constructions append a sink to B and then a new root beating only B's old
root, with the two new vertices exchanged in their birth order. It is also
the inherited rank-two source identity ZJ=JO.

If the last letter is J, (1) becomes a+e+p=f+p in its predecessor, so f=a+e.
A predecessor with e>0 has f<=a by Section 1, a contradiction. Thus e=0.
If the predecessor is b, the output is the triangle. Otherwise it has a
presentation O(B). Substituting the O formulas into f=a gives

    a_B+2f_B+e_B+2p_B=a_B+e_B,

so f_B=p_B=0. Hence B is balanced and the output is J(O(B)).
This exhausts every physical presentation. The singleton boundary has been
handled separately rather than inserted into the positive cone.

Conversely, if B has D=U=h (including the singleton), the exact response of
J(O(B)) is (3h,h,3h,h). It satisfies U=E+2d and hence (1). QED.

This is an exact construction-class characterization. Different balanced
inputs can have the same h, and no new graph-isomorphism uniqueness claim is
needed here. The previously proved canonical balanced inverse is compatible
with this description.

## 3. The complete norm-fixed scalar image

The edge has (x,y,z)=(3,1,1), giving N=13. For J(O(B)), the profile in the
converse part gives

    (x,y,z)=(5h,2h,2h),
    N=25h^2+16h^2=41h^2.                               (2)

Thus the full physical fixed-point scalar image is exactly

    {13} union {41h^2 : B balanced, h=H(B)},             (3)

where the singleton supplies h=1. Its only prime values are 13 and 41:
every h>1 makes 41h^2 composite. Conversely both primes have the explicit
edge and triangle source witnesses, and their C2 grafts have those counts.

A different physical diagonal source need not be a swap fixed point; for
example the source word JOZ from the edge has profile (9,1,11,1) and gives
N(v,v)=601. This does not contradict (3). Likewise a failed lift of the
classical swap or outer norm involution is not a proof that the integer itself
is absent from any broader physical class.

The bridge's endpoint is now exact: the formal sum-of-two-squares fixed point
would need actual source membership, and enforcing that membership collapses
its prime outputs to two. No idoneal, Gaussian factorization, or prime-density
claim is imported.

## 4. Exact checks

Run `python3 verify.py` and `python3 -O verify.py` in this directory. The
checker authenticates its frozen source dependency by SHA256, independently
enumerates physical sources, confirms the iff classification and existence of a
numerical balanced predecessor for every qualifying profile in the declared
bounded universe, checks the formula on additional repeated balanced inputs,
and retains the nonfixed N601 control. The graph-level constructive proof
above gives the all-order result; bounded checks are not extrapolated.
