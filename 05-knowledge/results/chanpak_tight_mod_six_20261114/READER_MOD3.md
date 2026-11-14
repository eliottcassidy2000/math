# Why parity controls modulus three in every homogeneous reader

A reader replaces one distinguished vertex of a fixed quotient by
an arbitrary nonempty tournament A; all edges between A and the
other vertices or modules are uniform as specified by that quotient.
Write h=H(A), n=|A|. Compress maximal visits to A as in the
[run decoder](../chanpak_family_selection_20261114/PRIME_QUOTIENT.md).
The output count has the form J(A)=Σ_j c_j G_j(A), where G_j
counts j ordered nonempty covering paths and c_j are nonnegative
integers independent of A.

**Universal mod-3 test.** With TT_m transitive of size m,

    J(A) ≡ h · J(TT₁) mod 3   if n is odd,
    J(A) ≡ h · J(TT₂) mod 3   if n is even.

Only two tiny test sources are needed for this particular residue
of the entire arbitrary reader. The analogous statement for parity
is J(A)≡1 mod 2 for any literal nonempty reader output.

**Proof.** Permuting the labelled runs freely gives j! | G_j.
So for mod 3 only G₁=h and G₂=F₂−2h can survive. By
[Grinberg–Stanley’s odd-cycle packing formula](https://arxiv.org/src/2307.05569v1/anc/redeiberge-long.pdf#page=14),
each power in F_q(A) has exponent e of the same parity as n:
collapsing one odd cycle changes the number of atoms by an even
number. Thus F₂≡(−1)^n h mod 3. This gives G₂≡0 for n odd and
G₂≡2h for n even. Evaluate the coefficients c₁ and c₁+2c₂ at
TT₁ and TT₂ respectively.

There is a sharp reason for 2 and 3 rather than any two convenient
small primes. For every integer q, all tournaments of order n have
F_q(A)≡q^n H(A) modulo p for prime p=2 or p=3, because the
nonzero elements have square one (and all powers are positive).
No larger prime has this universal property: at q=2 the directed
triangle gives the difference (2³+2·2)−2³·3=−12. Thus the collapse
of *all* weighted color powers to just parity is special.

For any three uniform cyclic modules A,B,C, apply the test
to each module in turn. If e among their orders are even, the
result is their product of path counts times the cyclic readout
on e copies of TT₂ and 3−e copies of TT₁. At e=0,1,2,3
those actual small counts are 3,5,11,45. Therefore a uniform
cyclic join is automatically divisible by three when the
module orders all have the same parity, and otherwise its
count is congruent to −H(A)H(B)H(C). In particular the factor
of three can be forced by parity alone, without equal modules
or a three-fold automorphism.

## A physical branch test, not division by three

The singleton triangular wrap W(A)=C₃[A,TT₁,TT₁] reads h+F₂,
hence has a mandatory factor of 3 when n is odd. The reader
C₃[A,TT₂,TT₂] reads 4F₃+5F₂−11h and has a mandatory factor
of 3 when n is even. Their two tests are

| Reader | TT₁ | TT₂ |
|---|---:|---:|
| W | 3 | 5 |
| Two pairs | 11 | 45 |

For a target not divisible by 3 and a source with 3∤h, one can
immediately reject the wrong geometry according to the parity
of its intended order. Neutral singleton padding changes n parity
and doubles F₂ without changing h, so that parity is a physical
choice, not part of the scalar count alone.

But mandatory divisibility does **not** construct a tournament
realizing the quotient by 3. For example the two-pair reader has
coefficient vector (11,17,4) in G₁,G₂,G₃. A fixed homogeneous
reader of its exact quotient for every even-order source would
need coefficients (11/3,17/3,4/3): test even-order transitive
sources, for which G_j(TT_n)=j!S(n,j) has dominant term j^n
and these functions are linearly independent on even n. The
coefficients of any literal fixed reader are integers. Thus no
such fixed readout performs the divide, even with an arbitrary
fixed quotient. This leaves target-dependent changes and
nonuniform contacts open.

## Related mechanisms, and their limits

* **Collatz.** The parity of an exponent selects which powers of
  2 undo “3x+1” modulo 3. Here that inverse move becomes a literal
  binary append with exact tournament recurrences; see
  [COLLATZ_RAYS.md](COLLATZ_RAYS.md). [Tao's work](https://arxiv.org/abs/1909.03562)
  uses a 3-adic distribution for a different purpose, control of
  typical orbit sizes. Neither calculation proves or reduces one
  conjecture to the other.
* **Prime sums and gaps.** Apart from 2 and 3, primes lie on the
  two signs ±1 mod 6. For binary Goldbach at an even target,
  opposite signs are necessary when the target is 0 mod 6;
  the same sign is necessary at 2 or 4, unless one summand is
  3. Twin primes above 3 must straddle a multiple of 6. These
  are local filters, not existence in a bounded interval. The
  distinction is the same one exposed by the two tournament
  Collatz rays: every joint residue can be present while small
  absolute witnesses are absent.
* **Primitive Pythagorean triples.** From squares mod 4 and mod
  3, exactly one leg is even and exactly one leg is divisible by
  3. This is a structural partition of every primitive witness,
  not a selector for every proposed hypotenuse. What transfers
  is the order of proof: classify the local modes first, then
  require genuine objects of the correct size.

The two-prime mechanism is a small unit group plus a physical
recursion, not a coincidence of recurring numbers. The positive
height-window constructions are in [TIGHT_COMPLEMENT.md](TIGHT_COMPLEMENT.md).
