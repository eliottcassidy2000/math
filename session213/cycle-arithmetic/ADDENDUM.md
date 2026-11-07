# Exact arithmetic of the four binary-cycle branches after deletion

**Status: PROVED compiler-specific identity and equality classification;
FINITE-EXACT controls are separate. No new spectrum completeness or
literature-priority claim.** This addendum leaves the first frozen packet
unchanged. It uses the definitions and elementary one-apex formula in
`PROOF.md`, and the unified family in Session212
`session212/cycle-bridge/PROOF.md` (immutable source pinned in `SOURCES.json`).

## 1. All four branches

Use Session212 parameters `m>=1`, `0<=q<2^m`, `tau in {0,1}`. Its two
comb tournaments have chain vertices `x_0,...,x_(m+1)` and extra vertices
`a,b`. The vertex `a` beats `x_1` and loses to every other chain vertex.
The coding apex `b` beats `x_0`; `x_1 -> b` iff `tau=1`; and
`x_(i+2) -> b` iff the low-to-high bit `b_i` of `q` is one. The last
comparison is `b -> a` in `B-`, or `a -> b` in `B+`.

Deleting `a` from either orientation leaves the same one-apex core
`U(m,q,tau)`, whose bitword is `(0,tau,b_0,...,b_(m-1))`. Recall
`D(q)=H(A(b))=1+2N(b)` from the first packet. The inherited cycle counts are

    c_minus = HC(B-) = q+tau,
    c_plus  = HC(B+) = q+tau+1.

**Theorem A.**

    H(U) = D(q) + (3-tau)q + (1-tau)(q mod 2) + 2tau.      (A1)

**Proof.** For `tau=1`, this is Theorem 2 of `PROOF.md`:
`H(U)=D(q)+2q+2`. For `tau=0`, the two zero sentinels contribute,
for each bit `b_i=1`, the weights `w(i+2)` and `w(i+1)` to `N`.
The first sum is `q`; the second is `(q+(q mod 2))/2`. The two
sentinels do not form a contributing pair, and the residual pairs give
`N(b)`. Thus `H(U)=D(q)+3q+(q mod 2)`. This proves (A1). QED.

## 2. Exact `2c+1` outcomes of this deletion

Only positive input cycle counts are considered below.

**Theorem B.** For the specified deletion `B -> B-a`:

1. For `B-`, `H(U)=2c_minus+1` iff `tau=1` and `c_minus` is a power of two.
2. For `B+`, `H(U)=2c_plus+1` iff either
   - `tau=0,q=1`, giving `c_plus=2`; or
   - `tau=1,q=3*2^j-1`, for some `j>=0`, giving `c_plus=3*2^j+1`.
   In the second case the permitted bit length satisfies `m>=j+2`.

**Proof of 1.** Subtract `2(q+tau)` from (A1):

    H(U)-2c_minus = D(q)+(1-tau)(q+(q mod 2)).

Since `D>=1`, equality to one for `tau=0` forces `q=0`, whose cycle
count is zero and is excluded here. For `tau=1`, equality is `D(q)=1`,
which the first packet classified as `q=2^k-1`. Hence `c_minus=2^k`.

**Proof of 2.** For `tau=0`, (A1) and `c_plus=q+1` reduce the target to
`q+(q mod 2)+D(q)=3`. At `q=0` the left side is one; at `q=1` it is
three; at `q=2` it is five. For `q>=3`, it is at least `q+1>=4`, so
only `q=1` works.

For `tau=1`, `c_plus=q+2` and the target is exactly `D(q)=3`, or
`N(b)=1`. There must be exactly one zero-before-one pair and it must be
adjacent. Indeed any bit between such a pair produces a second pair:
a zero pairs with the final one; a one pairs with the initial zero.
All bits before this `01` must be ones, and all bits after it must be
zeros, to avoid further pairs. Therefore

    b = 1^j 0 1 0^d,
    q = (2^j-1)+2^(j+1) = 3*2^j-1,
    m = j+2+d.

Conversely each such word has precisely one contributing adjacent pair,
so `N=1,D=3`. QED.

**Corollary C (precise density scope).** The set of positive input cycle
counts for which at least one of these four branches followed by deletion
of `a` gives `2c+1` is exactly

    {2^k : k>=0} union {3*2^j+1 : j>=0}.                    (A2)

This set has `O(log X)` elements at most `X`, hence density zero among
positive integers. The corresponding path targets are
`2^(k+1)+1` and `3*2^(j+1)+3`. These are familiar realizable sparse
families, not new global coverage. Their first terms are respectively
`3,5,9,17,33,...` and `9,15,27,51,99,...`.

This density statement concerns only the exact target equation in this
specified cycle family and specified vertex deletion. It is not a claim
about the image density of all one-apex path counts, all two-apex path
counts, or unrestricted tournament Hamiltonian-path counts. Deleting a
different vertex or modifying old arcs is outside its scope.

## 3. Controls

Run `python3 addendum_verify.py` and `python3 -O addendum_verify.py`; both
must reproduce `ADDENDUM_RESULT.json`. The verifier imports only the
first packet's independently implemented basic graph constructors and
dynamic programs. It checks all `m=1,...,8`, all `q`, both `tau`, and
both apex orientations: 2,040 cycle-compiler cases. It checks (A1), both
exact-equality classifications, and the equivalent `D=3` bit pattern.
These bounded checks support the implementation, not the all-order proof.
