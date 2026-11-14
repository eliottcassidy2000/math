# All singleton contacts: a complete mod-three skeleton

This specializes the [exact cut-state transfer](CUT_STATE.md) to
contacts that are singleton strongly connected components of the
base. Every one of the r outcontacts of the apex has this property.
Everything between consecutive contacts is an arbitrary ordinal
tournament seen entirely on the losing side of the apex. An
initial all-losing prefix factors in the scalar; a trailing
all-losing block is handled like any gap. Empty gaps are allowed.

For a nonempty all-losing block X, write a=H(X), q=P_X(2).
Its ending total is q−a and its signed odd-minus-even total is a
by the signed cover identity. It applies
[
 (h,N,S,d)\mapsto(ah, aN+(q-a)S+ad, qS, 0).
]
An outcontact singleton applies
[
 (h,N,S,d)\mapsto(h, N, 2S+h, h).
]
Start with (1,1,0,0) for the empty base. This also explains why
one neutral contact at the right can matter to a *later* suffix,
even though it contributes nothing to the count at the moment.

For r=2, ignoring the initial prefix, the exact normal form is
the four-block chart p→A→q→D of [SPARSITY_AND_INVERSE](
SPARSITY_AND_INVERSE.md), not just one example of two marked
singletons. Every such pair lies in that form. In particular its
boundary formulas matter for adjacent contacts or a terminal
contact. For a prime output with a genuine cycling core, the
initial prefix must have H=1, since it contributes a product
factor and the apex core has count at least three.

## Replace every losing block by at most two vertices modulo three

Given a marked singleton-contact source B, form a transitive
painted word w by retaining every outcontact and replacing each
nonempty all-losing block by a losing singleton if its size is
odd and a losing pair if its size is even. Omit empty blocks.
Let G(w) be the path count after adding the apex to this word.
Then
[
 H(B+v)\equiv H(B)\,G(w)\pmod 3.
]
This is an exact statement for arbitrary many contacts and
arbitrary actual losing modules. No modular division is needed.
Indeed the ordered-cover odd-cycle expansion gives
q≡(-1)^|X| a mod 3, since all packing blocks are odd. Thus the
entire first displayed transition for X is a times the transition
for its one- or two-vertex surrogate modulo 3. For the even
surrogate the two all-losing singleton transitions send
(N,S,d) to (N+3S+d,4S,0), identical mod 3 to the required even
case. Scalars from previous blocks factor since all coordinates
are homogeneous, proving the claim even when 3 divides H(B).

Equivalently, use a three-register finite automaton over F3,
starting at (N,S,d)=(1,0,0):
* contact: (N,1−S,1);
* nonempty even losing gap: (N+d,S,0);
* odd losing gap: (N+S+d,−S,0).

This is a complete necessary **and sufficient** test for avoiding
divisibility by 3 within an already built chart: the source H
must be coprime to 3 and the first automaton output nonzero. It
is not an absolute target test, and contact vertices contained
in nonsingleton components can escape it. It does let a
target-driven search reject whole synchronized families before
ever evaluating their exact integer residual. For two contacts
and both gaps nonempty it reduces to the odd/even test in the
four-block inverse.
