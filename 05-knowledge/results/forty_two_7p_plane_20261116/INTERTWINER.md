# Two full physical factor-seven intertwiners, not just anchors

The modulus-30 right-hand witness in README has full outgoing masks

    (256,269,281,97,299,327,407,319,8).

Keep the induced order-seven subset of mask 379, namely
{0,1,3,4,5,6,8}, with the same root 8. There are **no old arc
reversals**. Direct subset DP gives the full profiles

| witness | H | g | E | O |
|---|---:|---:|---:|---:|
| induced parent | 43 | 13 | 15 | 28 |
| two-vertex extension | 301 | 91 | 105 | 196 |

The entire four-coordinate vector, not just H, is scaled by
seven. This upgrades the right half of the 41/43 calculation:
the two **different E105 child cores** have a joint near-sink
ray; within that pair the actual H43->H301 induced lift is a
full intertwiner. This does not contradict the uniform near-sink
no-go, since the two added vertices at the seed do not have
the uniform block contact used in that obstruction.

Any common ordinary-sink or near-sink step acts on the four
coordinates by a homogeneous linear map. For a singleton
near-sink, the formula is (H+2E,g,E+O−g,E+O). For an ordinary
sink it is (H,g,E+O,E+O): the new sink can be placed at the end
of either part in the weighted path split. So every common
word in these moves retains exact sevenfold scaling. The
induced relationship survives too, because the new vertices
see every old nonroot in the prescribed uniform way; deleting
the same two exceptional seed vertices leaves the small
continued tournament. The same argument works for arbitrary
uniform near-sink blocks, whose transfer is homogeneous
linear in the source coordinates for each fixed right block.

In particular, every physically generated parent count
30·2^k+13 and every parent 30v+13u from the right block family
comes with a literal factor-seven lift, not merely a scalar
identity. This is an infinite positive intertwining theorem.
It applies to a specific rooted source and its continuations,
not to all prime witnesses.

There is a second exact full intertwiner for 41 at a different
slope. Take strong H287 rows
(256,5,329,307,295,455,411,287,2), with root 8. The induced parent
mask 367 includes the root. Its profile (41,13,14,27) is scaled
exactly to (287,91,98,189). Hence every actual parent
28·2^k+13 also comes with its induced 7-times lift.
This is not the left E105 witness in README; that witness is
an exact scalar anchor and has a different continuation.

## Why the two exact lifts do not already give a twin selector

The 41 and 43 parent instances above have the same g=13 but
even responses 14 and 15. For a nonempty binary word ending
in a near-sink, of length L and encoded integer x, the exact
response from a stable root is E'=2^L E+13x/2,
where 0<=x/2<2^(L-1). For distinct lengths the possible
response intervals on the two sides are disjoint:
one extra step multiplies the lower limit by two, while
the largest base upper limit is only 15+6.5. For equal
lengths the responses differ modulo 13, since 2^L is a
unit. The empty word fits its separate endpoints and
cannot coincide with a positive length either. Thus no two
such words synchronize these particular exact lifts, even
if different words are chosen. This is a precise reason
new source families or a third contact operation are needed.

Replay this fragment in the same folder:

    from check_rail import profile, induced, rails30
    from itertools import product
    def move(a,r,b):
        n=len(a)
        return tuple(row|((b or j!=r)<<n) for j,row in enumerate(a))+(0 if b else 1<<r,)
    for a,m in [(rails30[1],379),((256,5,329,307,295,455,411,287,2),367)]:
        s=[i for i in range(9) if m>>i&1]
        b=induced(a,s)
        for k in range(4):
            for seq in product((0,1),repeat=k):
                t,u=a,b
                for digit in seq:
                    t,u=move(t,8,digit),move(u,6,digit)
                assert induced(t,s+list(range(9,9+k)))==u
                assert tuple(7*x for x in profile(u,6))==profile(t,8)

## Modular targets: every odd modulus, with literal induced sevenfold lifts

The two exact lifts also give an all-order *local* selector.
On normalized parent coordinates I=H−2E and Z=2E, g=13, use
A=N (one near-sink) and B=SN (one ordinary sink followed by
one near-sink). Both end normalized and obey

    A(I,Z)=(I,2Z),      B(I,Z)=(I−Z,4Z+2g).

Work modulo an odd positive modulus M. These are permutations of
the finite profile space, so each inverse can be implemented by a
positive power. Composition below is right-to-left. Put

    F=B A^(-2):              (I,Z) -> (I−Z/4, Z+2g),
    C=(A F A^(-1)) F^(-2):   (I,Z) -> (I+3Z/8−g, Z).

The start Z=28 or 30 is a unit mod 13. Multiplication by 2 is
transitive on nonzero residues modulo 13 (order 12). Therefore A
can reach any specified nonzero residue mod 13 when 13 divides M,
and F adjusts all the remaining digits modulo M, since F adds 26.
When 13 does not divide M, F alone reaches every Z. This gives the
exact necessary restriction on Z: mod 13 both generators preserve
nonzero Z.

To tune I, first use F to move to some Z where 3Z/8−13 is a unit
modulo M. Modulo 3 this is automatic, and modulo 13 it is automatic
at unit Z. For every other prime factor of M it excludes just one
residue, which can be avoided by F and CRT simultaneously. Powers of
C then translate I by any desired amount, at fixed Z. Undo the
temporary F transport and land at the prescribed Z with arbitrary
I. Thus every H residue mod every odd M can be prescribed, and by
the full-profile scaling the very same word supplies the exact
induced larger count 7H. This does not guarantee the prescribed
*absolute* H or that H is prime.

For completeness, a small independent modular replay:

    for m in (3,5,7,9,13,19,39,169):
        for z in (28,30):
            pending=[(13%m,z%m)]
            reached=set(pending)
            for i,z in pending:
                for u in ((i,2*z),(i-z,4*z+26)):
                    w=tuple(t%m for t in u)
                    if w not in reached: reached.add(w);pending.append(w)
            assert len(reached)==m*(m if m%13 else m-m//13)
            assert {(i+z)%m for i,z in reached}==set(range(m))
