# Independent left and right growth cannot be a closed ternary digit grammar

**PROVED, elementary all-order obstruction.** Fix two words u,v. On attachment
words define L(w)=uw and R(w)=wv. These operations commute literally:

    L(R(w)) = uwv = R(L(w)).                           (1)

Suppose a nonempty family of physical words is closed under both operations,
and a scalar invariant h satisfies

    h(L(w)) = a h(w)+b,
    h(R(w)) = c h(w)+d                                (2)

for every word in that family. Necessarily

    (a-1)d = (c-1)b.                                  (3)

Indeed the two orders in (1), evaluated by (2), yield ac h(w)+ad+b and
ac h(w)+cb+d. Their equality gives (3). Closure through the two mixed
second-generation words is already enough.

In particular no such nonempty closed family can have independent fixed
prepend/append operations realizing 3h+2 and 3h+4: their offsets violate
2*4=2*2. This conclusion is independent of Hamilton-path counting or the
number of response coordinates. It also applies to any two independent
commuting clock increments, even when their physical blocks have unequal
lengths. The preserved invariant may be arbitrary, but its two asserted
updates cannot be incompatible affine maps.

## A concrete successful first fan and its forced failure

At the strong two-apex word w=120 with H=13, the asymmetric operations

    A(w)=w11,           B(w)=3w

give strong children

    H(A(120))=41=3*13+2,
    H(B(120))=43=3*13+4.

Thus asymmetric end growth does repair the one-birth obstruction at this
single parent. But both mixed operations give exactly the same physical word

    B(A(120)) = A(B(120)) = 312011,

whose actual count is 179. A closed ternary rule would require this same
word to have both 3*41+4=127 and 3*43+2=131. The failure is forced by literal
commutation, not merely the observed number 179.

The direct word counts and strongness are verified in `verify_extension.py`.
The formula remains a one-parent fan, not a recursively closed grammar.

## Scope and surviving route

The obstruction does not apply to state-dependent prefixes or suffixes,
re-marking apexes, deleting or rearranging old names, interleaving new letters,
or a finite-state rule whose branch actions genuinely fail to commute.
It says that changing from one clock to two *independent commuting* clocks
cannot by itself encode the order-sensitive ternary digits. A constructive
repair must preserve a coordinate carrying branch order and let subsequent
physical operations depend on it, or use intrinsically noncommuting growth.
This is a structural requirement, not a proof that the full 17-dimensional
two-apex grammar cannot realize all permitted odd values.
