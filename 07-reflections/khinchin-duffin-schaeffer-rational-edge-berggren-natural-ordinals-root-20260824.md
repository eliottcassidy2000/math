# Khinchin, Duffin--Schaeffer, rational edges, and Berggren natural ordinals

**Session:** root / 2026-08-24  
**Canonical outcomes:** THM-4056 and THM-4057  
**Headline status:** the discrete bridges below are proved; irrationality of
`e+pi` and of Khinchin's constant remain open; no prize problem is closed.

## 1. Portfolio and inheritance pass

### Anchor

Turn the proposed natural-number rational-edge picture into a correctly typed
object, then ask whether Khinchin/Duffin--Schaeffer input supplies arithmetic
or tournament structure.

### Niche

Use the odd-square/Berggren chart to make the natural-number ordinal explicit,
including inverse descent and parity content rather than only a scalar rank.

### Wildcard

Treat metric approximation mass as a periodic discrete clock and compare its
phase, lift, reciprocal, and hidden-coordinate losses with the repo's clock
and tournament work.

### Closest proved mechanism

- THM-3756: coprime odd roots `(q,d)`, exact Berggren affine descent, and the
  warning that outer odd-square rank is not a node address.
- THM-3509: all reduced fractions, the `J(x)=(1-x)/(1+x)` leg-swap fibre, and
  the distinction between a six-edge carrier and a tournament.
- THM-2294: antisymmetric data orient; symmetric data colour.
- THM-3744: a continued-fraction scalar loses the ordered recurrence/carry.
- HYP-2212's proved sharpening: for `S=e+pi,P=e*pi,D=e-pi`, no two are
  algebraic.

### Canonical hostile and corrected near miss

- Coprime pairs are missing natural edges such as `{2,4}`; all reduced
  rationals include both orientations. Neither object is a tournament.
- Metric full measure is not pointwise membership. The golden-ratio hostile
  in THM-4056 is explicit.
- Reciprocal rationals share a custom projective Euclidean coefficient word,
  but not the standard finite Khinchin word. The exact hostile is
  `3/5=[0;1,1,2]` versus `5/3=[1;1,2]`.
- The Stern--Brocot and Calkin--Wilf trees share reciprocal word reflection
  but do not assign the same rational to a mixed word.

### Least-used sidecar

Prime-power denominator depth. Reduced-residue phases and the factor
`phi(q)/q` see only `rad(q)`, whereas natural-edge lift multiplicity sees the
full denominator. The collision `q=2` versus `q=4` became the decisive
sidecar test.

## 2. Live concept board

The session kept six objects live and compared every new signal with all six:

1. scale-decorated rational edges `(p,q;g)`;
2. reciprocal/Stern mirror and direction bit;
3. Calkin--Wilf natural ordinal `k`;
4. Berggren branch word, parity, and Pythagorean content;
5. Duffin--Schaeffer phase/lift clocks;
6. Vieta trace/product/discriminant and named-point irrationality gates.

The main gain came from refusing to quotient them simultaneously. Every
useful theorem says exactly which coordinate survives one quotient and which
must be restored as a sidecar.

## 3. Triangular integers recover the odd coordinate, not the address

For the integer extension `T(z)=z(z+1)/2`,

```text
T(z+h)-T(z-h)=h(2z+1).
```

At `h=2`, this is the motivating identity

```text
T(z+2)-T(z-2)=2+4z.
```

Putting `z=r-1` and dividing by two extracts the `r`th positive odd number
`2r-1`. THM-3756 already turns the two odd square roots adjacent to a primitive
Pythagorean hypotenuse into a lossless two-coordinate chart. The new session
did not promote “the nth odd square is n” into a unique node label; instead it
found the missing address coordinate in the Calkin--Wilf heap.

For standard Euclid parameters `(m,n)`,

```text
(q,d)=(m+n,m-n),
(m,n)=((q+d)/2,(q-d)/2).
```

Thus the triangular odd-coordinate extractor, the reduced rational `d/q`,
and the Berggren parameter pair are exact coordinate changes. The ordinal
`k`, not `q` alone, supplies unique ancestry.

## 4. The rational-edge proposal has an exact three-way type split

### Reduced arcs

All `p/q` give a bidirected coprimality graph. Restricting to `p<q` gives an
acyclic orientation of an incomplete graph.

### Complete natural tournament

Orient every `a<b` upward and retain

```text
g=gcd(a,b),       (p,q)=(a/g,b/g).
```

The primitive type fibre is exactly `floor(N/q)`, giving

```text
C(N,2)=sum_(q=2)^N phi(q)floor(N/q).
```

### Arbitrary tournament

Every tournament on the positive naturals is an arbitrary bit

```text
epsilon(p,q;g)
```

on primitive reciprocal pairs and scale. A rational-only selector exists
exactly when the bit is dilation-invariant. Directed ratio products telescope,
so rational magnitude is a coboundary and cannot see cycle curl.

This is why the proposed edge picture is useful but does not automatically
inherit tournament theorems.

## 5. Natural-number breakthrough for the Berggren tree

Let `s(k)` be Stern's diatomic sequence and enumerate the Calkin--Wilf tree by

```text
k -> s(k)/s(k+1).
```

Reciprocal reflection within a dyadic layer is

```text
k*=3*2^floor(log2 k)-1-k.
```

The primitive-Pythagorean parameter subset is exactly

```text
k=3 or 5 mod 6.
```

On this set, the Berggren branches become

```text
A(k)=2k-1,
B(k)=4k*+3,
C(k)=4k+3.
```

Their images are recognized by terminal radix classes `1 mod 4`, `3 mod 8`,
and `7 mod 8`, giving a complete inverse descent. At Berggren depth `d`, the
number of nodes pulled back to Calkin--Wilf depth `1+d+j` is

```text
binom(d,j)2^j.
```

This is a stronger and cleaner natural-number map than outer odd-square rank:
it is injective, carries ancestry, and makes reciprocal reflection affine on
each dyadic layer. Its density `1/3` is encoding-dependent and cannot be
transferred to the ternary heap without the index map.

## 6. The tournament/Berggren collision is sharp at the root

Each parent-child Berggren incidence forms a natural three-vertex triangle
after adding the remaining closure. Its closure gcd is controlled by the
parity of one Euclid endpoint, so the depth-`d` incidences split exactly into

```text
primitive K3      (3^(d+1)+(-1)^d)/2,
scale-two closure (3^(d+1)-(-1)^d)/2.
```

The natural triangle is transitive. Flipping its unique least-to-greatest edge
makes it cyclic; that edge is the closure in `A,B` and the child in `C`. The
unique collision among its normalized unordered Pythagorean edge shadows
occurs at the root `A` branch:

```text
(1/2,2/3,1/3),       U(1/2)=U(1/3)={3,4,5}.
```

After flipping `1->3`, the triangle is `1->2->3->1`, while every gcd and the
entire unordered-Pythagorean multiset is unchanged. The missing data are
exactly direction, ordered legs, and scale.

## 7. Duffin--Schaeffer mass has two clocks, not one

For finite denominator support and `L=lcm(1,...,Q)`, the reduced-residue
phase clock is

```text
P(t)=2 sum_q psi(q)1_[gcd(t,q)=1].
```

The natural-edge lift clock is

```text
R(t)=2 sum_(q|t) psi(q)phi(q).
```

Both have complete-period mean

```text
Psi_Q=sum_q 2phi(q)psi(q)/q,
```

the finite raw Duffin--Schaeffer mass. They are not pointwise equal: for
`Q=2,psi(2)=1`, their periods are `(2,0)` and `(0,2)`.

The phase clock has a signed divisor expansion but remembers only sums over
equal radicals. Moving unit weight from `q=2` to `q=4` preserves phase and
mean, while the lift periods change from `(0,2,0,2)` to `(0,0,0,4)`. Thus
mean equality is a first-moment correspondence, not an overlap theorem.

This mirrors several prior repo lessons:

- LRC clock mass does not decide interval containment or owner arrival;
- Rule 30 phase coincidence does not recover a hidden response bit;
- sequence support does not recover indexed multiplicity;
- a symmetric Jacobian observer does not recover the hidden source jet.

The common meta-pattern is now evidence-rich, but this session did not add it
to META-PATTERNS because the existing carrier/sidecar patterns already cover
the mechanism.

## 8. Khinchin and named irrationality: exact positive and negative results

### Positive

If an irrational's infinite continued-fraction digit geometric mean has a
finite limit, then its irrationality exponent is exactly two. This is the
lawful Khinchin-to-Duffin--Schaeffer bridge.

Membership in `W'(psi)` with `psi(q)->0` certifies irrationality, and

```text
W'(q -> 1/q)=R\Q.
```

### Negative

The golden number `alpha=(sqrt(5)-1)/2` is outside
`W'(q->1/(4q))` even though the Duffin--Schaeffer sum diverges and that set
has full measure. Hence metric abundance cannot name a point.

Euler's digits occur in blocks `(1,2k,1)`, whose product through `3n` digits
is `2^n n!`; the digit mean for `e` diverges. A named familiar irrational need
not be Khinchin-typical.

### `e+pi` proof gate

Assuming `S=e+pi` algebraic forces `P=e*pi` transcendental. A symmetric
recurrence coefficient which remains nonconstant in `P` cannot be
denominator-cleared from `S` alone. Independent linear forms can combine only
if each decay rate beats the other form's coefficient growth and a
synchronized nonvanishing statement survives. No existing packet meets those
requirements.

## 9. Connection ledger

| Source | Target | Map | Preserved | Destroyed | Needed sidecar / cheapest test |
|---|---|---|---|---|---|
| natural edge `(a,b)` | reduced rational | divide by gcd | slope, order | scale | retain `g`; test fibre `floor(N/q)` |
| rational edge | Stern--Brocot node | `(p,q)->p/q` | reciprocal mirror | tournament completeness | selector bit on every scale |
| Calkin--Wilf ordinal | Berggren node | `k=3,5 mod 6` | full parameter pair, ancestry | ternary depth | radix inverse (25) |
| odd roots `(q,d)` | Pythagorean triple | quadratic light cone | primitive triple | orientation after leg fold | sign of `q-d`, ordered legs |
| DS intervals | phase clock | reduced residue activation | raw mean | overlaps, prime-power depth | radical/depth sidecar; `2` vs `4` hostile |
| digit mean | irrationality exponent | infinite coherent CF limit | exponent two | Markov constant, named membership | ordered tail; golden hostile |
| `(e,pi)` | `(S,P,D)` | Vieta | symmetric algebra | owner and sign | trace-pure recurrence or synchronized forms |

## 10. Frontier consequences

The two canonical theorems generated six concrete tasks now routed to the
backlog: trace-purity and paired-form audits for `e+pi`; a prime-power-resolved
metric phase clock; a Berggren/Gauss transfer operator; finite-state reciprocal
tournament selectors; completion-fibre ranges; and the parity-restricted odd
Berggren limsup problem.

No LRC(14), planar Jacobian, Rule 30, `abc`, or IUT theorem follows from the
shared clock/Vieta/tree grammar. The useful advance is a set of exact carriers,
inverse maps, and minimal loss witnesses which prevents those frontiers from
being joined by analogy alone.

## 11. Post-merge generated breakthrough: the depth packet is inverse arithmetic

The first task generated by the merge closed as THM-4059. For reduced `p/q`,
write `pr-qh=1` with `r` the least positive inverse modulo `q`. Then

```text
D(p/q) = pq+pr+rh mod 2.
```

This is not a numerical fit: the Farey-parent matrix, reduced modulo two,
acts on the three nonzero vectors of `F_2^2`; its permutation sign is the
parity of the Stern word. The formula turns the tree checkerboard into an
exact denominator-packet weight.

For odd `q`, the sign is `(-1)^(a+a_inverse)`. For even `q`, it is the
negative of the parity of the Bezout carry `(a*a_inverse-1)/q`. The latter is
exactly the sheet bit deciding whether the inverse lifts from modulus `q` to
`2q` without moving or after adding `q`. Thus the prime-power clock and Stern
depth meet at a genuine dyadic lift coordinate.

The exact packet sum `S(q)` has inversion and negation symmetry and satisfies

```text
S(q)=phi(q) mod 4,  q>2.
```

For odd `q`, its odd- and even-numerator signed half-shells are both `S(q)/2`.
This strengthens the earlier `phi(q)/2` Pythagorean cardinality split to a
signed Stern-depth split.

The natural-number tournament connection is now exact. The lower-star
imbalance at apex `N` is

```text
B(N)=sum_(q|N,q>1) S(q),
```

and the full cyclic compiler including its declared zero packet is
`A(N)=1*S(N)`, with Möbius inversion recovering `S`. This is the first direct
map from an exact denominator packet to a tournament score in this session.

Two seductive upgrades fail minimally:

- the depth weight is not a unit-group character (`q=9`, `2*2=4`);
- `S` is not multiplicative (`S(2)S(5)=0`, but `S(10)=-4`).

Accordingly the associated Dirichlet series factors by convolution but has
no Euler product. A depth-signed Duffin--Schaeffer shell is a checkerboard
observable, not positive mass; its cancellation cannot prove coverage or
name `e+pi`. These failures define the next two tasks: bound inverse-parity
packet cancellation, and relate the lower-star sequence to genuine global
cycle statistics without confusing an apex star with a full infinite degree.
