# Independent audit: all-period golden fibres and changing-modulus lifts

**Verdict: PASS (independently proved and exactly checked).** This is an audit
recommendation, not a canonical promotion. The root session owns promotion.
No Collatz convergence, ordinary-integer classification, or priority claim is
made. The all-period conclusion follows from the proofs below, not from the
bounded word census.

## Audited objects and provenance

- Producer: `/tmp/session74-clocks/RESULT.md`, Sections 2–4.
- Prior result: `01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`
  at commit `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`.
  The producer's copied source is byte-identical to `git show` at that pin;
  SHA256 `3a3a04a17d40d366a594889d0d2651eed867d91ad108dd3e631eec247260203e`.
- Its previous word-fibre assertions were FINITE-EXACT through length 22,
  with the historical independent audit reporting word checks through 29.
- Independent implementation: `audit.py`, using Python standard-library
  integers and `Fraction`, with no imports of producer code and no floats.
  Runtime checks survive `python -O`.

## 1. Complete carry bounds, including the cyclic seam

Let O = Z[φ], φ² = φ+1, ψ = 1−φ = −1/φ. Let W_n consist of rooted binary
words of length n≥1 such that w_i w_(i+1)=0, indices modulo n. Set
F(w)=Σ_(i=0)^(n−1) w_i φ^(n−1−i).

For a collision modulo φ^n−1, define z_0=(F(u)−F(v))/(φ^n−1) in O and
z_(i+1)=φz_i−(u_i−v_i). Direct iteration gives z_n=z_0. The recurrence
therefore defines a bi-infinite n-periodic sequence in O. The signs and
reading order agree with the claimed left-shift convention.

Solving this recurrence forward in the real embedding gives

    z_i = Σ_(j≥0) (u_(i+j)−v_(i+j)) φ^(−j−1).

Each of the two nonnegative expansions belongs to [0,1]. Indeed, group
positions (0,1),(2,3),…; each pair contributes at most φ^(−2j−1), since
its digits cannot both be one, including pairs across the cyclic seam.
The sum of these bounds is φ^−1/(1−φ^−2)=1. Hence |z_i|≤1.

If M=max_i |z_i'|, conjugating the recurrence gives M≤|ψ|M+1, so
|z_i'|≤1/(1−φ^−1)=φ². For z=a+bφ, we obtain

    |b|√5 = |z−z'| ≤ 1+φ² = φ√5.

As 1<φ<2, b is −1, 0, or 1. The real bound gives |a|≤1+φ<3.
Checking these integer possibilities against both bounds gives exactly

    S = {0, ±1, ±(φ−1), ±(φ−2)}.

Thus the seven-carry list is complete without any empirical search cutoff.
The producer's optional coarse bound |z|≤φ follows instead from |u_i−v_i|≤1.
It yields |b|<2 and |a|<4, and exactly the two extra carries ±φ. The
producer's nine-carry enumeration ranges therefore also exhaust its region.

## 2. Independent finite graph and why SCCs suffice for every n

Use states (a,b,p,q), representing z_i=a+bφ and previous digits
p=u_(i−1), q=v_(i−1). A transition with current digits (u,v) is

    (a,b,p,q) -> (b−u+v, a+b, u,v),
    pu=qv=0,

provided its new carry is in S. This produces **28 vertices and 41 edges**.
Independent Tarjan decomposition gives precisely the following SCCs capable
of containing a closed walk (all tuple entries are integers):

1. D = {(0,0,0,0), (0,0,1,1)};
2. P = {(1,0,0,0), (−1,1,1,0)};
3. N = {(−1,0,0,0), (1,−1,0,1)};
4. A = {(−2,1,1,0), (2,−1,0,1)}.

D has the three internal edges 00→00, 00→11, 11→00, and hence exactly
identical legal word pairs. Each of P,N,A has exactly its two directed
cross-edges and no internal self-loop. They represent respectively an
alternating word versus zero, zero versus an alternating word, and the
two opposite alternating words. Their rooted rotations include every
ordered nonidentical pairing of those three words.

A closed walk cannot combine SCCs: after contracting SCCs, every cross-edge
strictly decreases the integer rank supplied in `receipt.json`; a closed
walk cannot decrease that rank. Thus there is no gap from treating an SCC
union as a single recurrent possibility. Every non-diagonal closed walk
has even length and uses one of the three two-cycles. This proves the
collision classification for every n, not just sampled n.

The producer's seven-row hand transition table and elimination argument
were also checked and are valid. In particular, a visit to carry +1 forces
(1,0), then (0,0), indefinitely; the tempting second (1,0) is forbidden by
no-11. A periodic visit to ±(φ−2) forces alternation backwards because each
has only the other as a predecessor; periodicity makes that a full two-cycle.

For even n, the two alternating readings really do vanish:

    F((01)^(n/2)) = (φ^n−1)/(φ²−1) = (φ^n−1)/φ,
    F((10)^(n/2)) = φ^n−1.

Here φ is a unit in O, since φ(φ−1)=1. Therefore the unique non-singleton
fibre is {0^n,(01)^(n/2),(10)^(n/2)} over zero.

## 3. Surjectivity, n=1 and n=2, and primitive periods

The adjacency matrix [[1,1],[1,0]] counts rooted closed walks, giving
|W_n|=trace(A^n)=L_n. This includes n=1: only the word 0 is legal, since
the digit 1 is adjacent to itself across the cyclic seam.

Multiplication by a+bφ on basis (1,φ) has matrix [[a,b],[b,a+b]] and
norm a²+ab−b². Since φ^n−1>0 and ψ^n−1<0 for every n≥1,

    |O/(φ^n−1)| = −(φ^n−1)(ψ^n−1)
                 = L_n−1−(−1)^n.

The collision classification therefore gives exactly the full quotient's
cardinality: L_n for odd n, L_n−2 for even n. This proves surjectivity.
At n=1, both sets are singletons. At n=2, W_2={00,01,10} and R_2 is the
one-element ring, so the stated three-way collision is the entire map.

Left shift satisfies F(shift w)=φF(w)−w_0(φ^n−1), proving equivariance.
These are rooted words / periodic points of period **dividing** n, not
necklaces and not exclusively primitive words. Away from the exceptional
zero fibre, uniqueness and equivariance preserve least period. In
particular, exact-period n words and exact-period n quotient labels are in
bijection for n≥3; the least-period-two sources are the exceptional collapse.

## 4. Rational standard-Collatz source bijection, separately proved

Use C(x)=3x+1 for odd x and x/2 for even x on Q∩Z_2, not the shortcut map.
For a legal word, h=#zeros≥1 and k=#ones. Its dictated affine composition is

    F_w(x)=(3^k x+B_w)/2^h,
    x_w=B_w/(2^h−3^k).

The denominator is a nonzero odd integer. Every cyclic rotation likewise
has a unique fixed point in Q∩Z_2. Cyclic composition shows its fixed
points satisfy the prescribed affine transitions x_i -> x_(i+1).
Those transitions really have the intended parities: if w_i=0 then
x_i=2x_(i+1) is even; if w_i=1, the following digit is zero, hence
x_(i+1) is even and x_i=(x_(i+1)−1)/3 is odd. This proves every prescribed
parity, without assuming in advance that a formal branch fixed point is
an actual Collatz point. Conversely, an actual C^n-fixed point determines
its unique legal parity word and this same affine fixed-point equation.

This proves the source bijection. The three exceptional words correspond
to 0,−1,−2. It does not classify which other rational sources are integers,
nor imply convergence of positive integers.

## 5. Period extension and the elevenfold hostile

For m,r≥1 set g_m=φ^m−1, Q=Σ_(j=0)^(r−1)φ^(jm), so g_(rm)=Qg_m.
The map

    J: O/(g_m) -> O/(g_(rm)),  a -> Qa

is a well-defined injective O-module homomorphism. If Qa=Qg_m b in O,
cancellation of nonzero Q proves a=g_m b. It commutes with φ and repeated
word reading is exactly multiplied by Q. Thus it preserves the same
rational source, while for even m the original three-source zero fibre
still needs to be retained separately.

The ordinary quotient map ρ in the reverse direction satisfies ρJ=r,
not generally the identity. At m=5, R_5≅F_11 and φ=4; r=2 is the old
factor-two colour swap. At r=11, all eleven projected labels become zero
but the eleven lifted labels remain distinct. Independent exact norm
calculation confirms |R_55|=312119004989. No enumeration of this ambient
ring was performed or needed.

**Additional precise algebraic corollary (independently proved):**
φ^m Q≡Q modulo g_(rm), hence Q²≡rQ. Consequently

    J(a)J(b) = r J(ab).

J is therefore generally not a ring homomorphism, not merely a non-unital
one. Its image is the ideal (Q)/(Qg_m). For m=5,r=11, multiplication of
any two elements of this eleven-element ideal is zero, because 11=0 in
R_5. The lift retains all eleven additive labels in a **square-zero ideal**,
while the naive return projection loses them. This is an optional corollary,
not an extra prerequisite for the fibre theorem.

### Smith return-loss formula: PASS

For m≥1, φ^m−1=a+bφ with a=F_(m−1)−1 and b=F_m. The integral
multiplication matrix [[a,b],[b,a+b]] has entry-gcd d=gcd(a,b)>0 and
absolute determinant N=|R_m|. The two Smith invariant factors are d and
N/d. In particular d² divides N because all matrix entries are multiples
of d. Therefore the additive group is Z/d ⊕ Z/(N/d), and multiplication
by r has exactly

    K(m,r)=gcd(r,d) gcd(r,N/d)

kernel elements. Its nonempty fibres all have that size and its image has
N/K elements. The return is injective exactly when gcd(r,N)=1. This
covers the one-element rings at m=1,2 without exceptions. Direct quotient
enumeration in the independent script verifies 264 cases, independently
of the Smith calculation and producer's residue-key enumeration.

### Clock sidecar and the 2/3/11 operator palette: PASS with precise scope

Repeating a length-m presentation r times preserves its underlying rational
point and least period, but multiplies traversal length, odd/even counts,
and additive accumulated traversal costs by r. Commutation with one-step
C does not erase this integer-clock information.

On R_5=F_11 the three return operators for r=2,3,11 are respectively ×2,
×3, and zero. With H={1,3,4,5,9}, ×2 swaps H and −H and reverses the Paley
tournament x→y iff y−x is in H. Multiplication by 3 preserves this tournament,
has order five, fixes zero, and rotates **each** of the two nonzero five-cycles;
indeed 3=4^4 mod 11, so this is four standard-Collatz shifts on labels.
The r=11 ideal lift is still injective although the return is zero.
This is a valid palette of return operators at lengths 5→10,15,55, not an
explanation of why some independent bracket classification selected these
three primes. Other nonsquares, squares, and multiples of 11 share the same
coarse return types. Source weights are transported, not replaced by a new
uniform ambient distribution.

The root synthesis `/tmp/session74-root/THREE_RETURN_OPERATORS.md` was
read and its mathematics accepted. Its checker was independently rerun and
matches `palette.normal.json` byte-for-byte (220 directed-pair arc tests,
total source mass 66). One wording guard was requested: say r repetitions
of the original length-five **presentation**, not r traversals of the minimal
source orbit without qualification, since the zero source has least period
one. This does not affect any count, fibre, or operator theorem. The same
presentation-versus-minimal-orbit distinction applies to general m.

## 6. Exact checks and reproduction

Run from any current working directory:

    python3 -B PATH/audit.py > receipt.json
    python3 -B -O PATH/audit.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

Both runs passed **105,461 checks**, producing byte-identical receipts.
Receipt SHA256:
`0b191472543d422b8ad7a65e360b9c7365af468d7a81af798e6df346b0cd9a71`.
The receipt contains every edge, all SCCs with a decreasing condensation
rank, exact word/source/quotient checks for n=1,…,16, extension checks for
m=1,…,10 and r=1,…,8, and successful checks of all elevenfold projected labels and 121 products
inside the square-zero image. It also records 264 independently enumerated
Smith return-kernel cases (m=1,…,12, r=1,…,22). The source code includes illegal internal-11
and cyclic-seam collision controls.

The producer's verifier was separately rerun successfully; its generated
receipt matches the producer's retained receipt byte-for-byte, SHA256
`8be63c669fae6bf6c5dd8f5f13b52b80b109cb28bfb1b669cf16bc3470260958`.
That rerun is a reproduction check, not the independence basis for this audit.

The separate producer Smith-control verifier was rerun successfully; its
receipt matches the retained `return_loss.json` byte-for-byte (396 cases),
SHA256 `feaac99e670cec179d285f2bf208adc1ff1e7f6793d068fe0565c469962572dc`.

No mathematical blocker to producer Sections 2–4 was found, including its
added Smith paragraph and clock sidecar. The presentation-versus-minimal-orbit
wording guard above is the only requested editorial precision repair.
The independent proof and exact graph support upgrading the bounded fibre
assertion to an all-period theorem, subject to the root session's explicit
promotion decision and normal provenance checks.
