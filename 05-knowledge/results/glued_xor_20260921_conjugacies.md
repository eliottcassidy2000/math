# Which signed Collatz maps are actually conjugate?

**Status: PROVED elementary classifications and domain distinctions;
FINITE-EXACT controls.** No complete classification of odd-parameter
Collatz cycles or convergence is claimed. The even-parameter cycle theorem
below concerns the explicitly defined **full odd-part map**, a different
system. No novelty claim is made.

## 1. Inheritance and the map definitions

The closest proved mechanisms are odd dilation and signed conjugacy in
[signed parameter strata, §§1–3](arithmetic_braids2_20260917_signed_cycles.md)
and the exact two-sheet model in
[the Catalan note, §2](catalan_elliptic_20260921_catalan.md).
The canonical hostile is the fixed source `1` for `U_(-1)`: any proposed
conjugacy must carry it to a fixed target. The corrected near miss is
applying a valid raw affine identity after silently changing the halving
rule. The least-used sidecar is the retained two-adic shell.

Anchor: the user's proposed opposite-sign pairs. Niche: all integer-affine
embeddings for the specified maps. Wildcard: whether three matched pairs
of known cycles supply a commuting threefold symmetry.

| Live concept | Exact retained coordinate | Cheapest hostile |
|---|---|---|
| Raw affine map | Fixed center and arithmetic sublattice | A raw arrow need not preserve a valuation |
| Odd acceleration | Actual power-of-two division | Compare valuation cells one and two |
| Full odd-part map | Parameter parity | Even parameter becomes raw after one step |
| Shell normalization | Common `v2(n)=v2(b)` | Full odd-part erases even scaling |
| Shortcut map | Both parity branches | Translating only the `3n+b` branch fails |
| Matched known cycles | Least periods `1,2,7` | Three cycles are not one period-three cycle |

Use these distinct definitions, with valuations taken on absolute values:

```text
A_b(n)=3n+b                       raw, n in Z;
O_b(n)=(3n+b)/2^v2(3n+b)          full odd-part, 3n+b != 0;
U_b(n)=O_b(n)                     n,b odd; hence v2 >= 1;
S_b(n)= n/2       if n even,
         (3n+b)/2 if n odd        shortcut, b odd, n in Z.
```

For `O`, the exponent can be zero. For even `b`, inserting an odd `n`
in `(3n+b)/2` is not an integer operation. Thus there is no odd-accelerated
`U_b` with even parameter under the stated positive-exponent convention.

Throughout, `H(n)=an+c` has integer coefficients and `a!=0`. An identity
`H f=g H` is a conjugacy onto the invariant image of H; it is an ambient
bijection of Z only when `a=+/-1`. No arrow reversal is being asserted.

## 2. Raw maps do realize the requested sign pairs, on sublattices

**PROVED complete raw classification:**

```text
H A_b=A_d H  iff  d=ab-2c.                              (G1)
```

Compare constant coefficients in `3an+ab+c=3an+3c+d`.
Equivalently H sends the raw fixed center `-b/2` to `-d/2`.
Bijective integer-affine raw conjugacy is possible exactly when `b,d`
have the same parity. Rational translation has no such parity restriction,
but changes the integer lattice.

Here are actual embeddings for the user's pairs:

```text
A_(-1) on n>=1  -> A_(-2) on negative n= -1 mod4,
                         H(n)=3-4n;
A_(+1) on n<=-1 -> A_(+2) on positive n= 1 mod4,
                         H(n)=-3-4n.                   (G2)
```

These are invariant arithmetic progressions, not the whole opposite
half-line. Among integer-affine raw embeddings into the stated opposite
half-line, absolute slope four is minimal: in the first case write
`a=-2r`, so (G1) forces `c=r+1`; then `H(1)=1-r<0` iff `r>=2`.
The second case similarly has `c=-r-1` and `H(-1)=r-1>0` iff `r>=2`.

The exact accelerated failure is already at these boundary inputs:

```text
H(U_(-1)(1))=-1,   O_(-2)(H(1))=O_(-2)(-1)=-5;
H(U_(+1)(-1))=1,   O_(+2)(H(-1))=O_(+2)(1)=5.           (G3)
```

The raw source numerator has valuation one in each example; the target
numerator has valuation zero. The raw coordinate change preserves the
affine arrow but destroys the halving guard.

## 3. Accelerated affine rigidity: translations cannot survive

**PROVED.** For odd parameters `b,d`, an integer-affine H mapping odd
integers to odd integers satisfies `H U_b=U_d H` on their entire valid
source domain iff

```text
c=0, a odd, d=ab.                                      (G4)
```

The same conclusion follows if the identity is required on every odd
integer of an unbounded positive or negative ray. It is a statement about
all inputs in that ray, not about a single sparse trajectory. For the
full odd-part maps `O_b,O_d`, defined on all valid integer inputs with
arbitrary integer parameters, exactly the same classification holds.

Proof. Fix a source valuation `k` and put `m=O_b(n)`. Then
`3n+b=2^k m`, and the proposed identity gives, for some target valuation l,

```text
a 2^k m+E = 2^l(am+c),       E=3c+d-ab.                 (G5)
```

Every fixed valuation cell contains an infinite arithmetic progression,
so m is unbounded there. The quotient of the two sides before the power
`2^l` tends to `2^k`. Powers of two are isolated, hence eventually `l=k`
on that cell. Consequently `E=2^k c`. For U use cells `k=1,2`; for O use
cells `k=0,1`. In either case `c=0`, then `d=ab`. Since every image is
odd, the equality forces a odd. Conversely odd dilation preserves
the valuation exactly and proves (G4).

Thus the only ambient affine bijections are identity or negation with
the corresponding parameter change. In particular

```text
U_b(-n)=-U_(-b)(n),
O_b(-n)=-O_(-b)(n),
S_b(-n)=-S_(-b)(n).                                    (G6)
```

Each is forward conjugacy. It does not exchange outgoing and incoming
arrows; the latter operation would produce a multivalued inverse relation.

## 4. Even parameters for full odd-part have a complete simple dynamics

For even b, every O-image is odd and, on odd inputs,

```text
O_b(m)=3m+b,
O_b^t(m)=3^t(m+b/2)-b/2.                               (G7)
```

**PROVED complete cycle classification for this model:** the only possible
cycle is the fixed point `m=-b/2`; it is in the odd domain exactly when
`b=2 mod4`. After its first image, every valid trajectory either is at
that odd fixed point and stays there, or escapes to positive or negative
infinity. Indeed a periodic point in (G7) must
satisfy `(3^t-1)(m+b/2)=0`, and the remaining coefficient is nonzero for
every other point. No even point can lie on a cycle because all outputs
are odd.

The first-image qualification is essential: the nonfixed even start two
satisfies `O_(-2)(2)=1` and then stays fixed. When `b=2 beta` with beta
odd, the entire fixed basin consists exactly of the integer values
`n=-beta(2^k+2)/3`, `k>=0`. This is the equation `3n+b=2^k(-b/2)`;
no other odd state can enter the fixed point under the injective raw map.

In particular `O_(-2)` has its unique fixed point at positive one, while
`O_(+2)` has its unique fixed point at negative one. On the user's proposed
target sectors, `O_(-2)(n)<n` for negative odd n and `O_(+2)(n)>n` for
positive odd n. There are no cycles in either sector.

Therefore **no map at all intertwining the whole requested source and
target sectors** can send `U_(-1)` on positive odds to `O_(-2)` on negative
odds: the source fixed point one would need a target fixed point there.
The same obstruction applies to `U_(+1)` on negative odds and `O_(+2)`
on positive odds. This is stronger than merely excluding affine conjugacy.

## 5. Retaining the two-adic shell restores even scaling

Write `b=2^s beta`, `n=2^r u`, with odd signed beta,u and `r,s>=0`.
For nonzero numerators, the exact three-way rule is

```text
O_b(n)=3u+2^(s-r)beta       if r<s,
       U_beta(u)           if r=s,
       3*2^(r-s)u+beta     if r>s.                     (G8)
```

For unequal valuations the summand with smaller valuation determines the
numerator valuation. At equality there is additional cancellation and
the odd acceleration reappears. This explains why forgetting the shell
changes the map rather than merely changing its scale.

On the fixed shell `v2(n)=r`, with parameter `2^r beta`, define instead

```text
V_(2^r beta,r)(n)=(3n+2^r beta)/2^(v2(3n+2^r beta)-r).
```

The exponent is at least one and the output remains on that shell. Then

```text
V_(2^r beta,r)(2^r u)=2^r U_beta(u),
O_(2^r beta)(2^r u)=U_beta(u).                          (G9)
```

Thus even dilation conjugates U to V, while full odd-part O erases that
dilation. The correct opposite-sign shell correspondences are

```text
positive U_(-1) --n -> -2n--> negative V_(+2,1),
negative U_(+1) --n -> -2n--> positive V_(-2,1).         (G10)
```

Notice the parameter signs in (G10). The known three cycles transport
exactly under these maps, with the same ordered exponents; their being
the only cycles at arbitrary period is not claimed.

Nor is taking an input's odd part by itself a dynamical quotient: inputs
one and two have the same odd part, but `O_1(1)=1` and `O_1(2)=7`.
To recover a map one must retain the shell and parameter coordinate.

## 6. A parity-swapped shortcut is exact only when both branches change

Define, for even d and odd e,

```text
P_(d,e)(m)=(3m+d)/2 if m even,
          (m+e)/2  if m odd.
```

**PROVED complete integer-affine classification:**

```text
H S_b=P_(d,e) H iff a odd, c=e, d=ab-c.                 (G11)
H S_b=S_d H     iff a odd, c=0, d=ab.                   (G12)
```

Here b,d in (G12) are odd. For (G11), an even slope maps both source
parities into one target branch, which cannot match the two distinct
source slopes. An odd slope preserves both parity classes; swapping their
roles requires c odd. Coefficient comparison on each class gives `c=e`
and `d=ab-c`. The same comparison without swapping gives (G12).

In particular translation by one gives a genuine `3m-2` branch:

```text
H(n)=n+1:  S_(-1) -> P_(-2,+1),
H(n)=n-1:  S_(+1) -> P_(+2,-1).                         (G13)
```

These translations preserve the relevant sign half-lines: positive source
states in the first case go to positive states at least two; negative
source states in the second go to negative states at most minus two.
Their fixed points move from one to two and minus one to minus two.
They do not give the opposite-sign identifications requested in the
question. The odd branch is essential: specifying only `3m-2` or `3m+2`
does not specify a shortcut system.

## 7. Three matched cycles do not supply a threefold commuting action

Consider only the known positive `U_(-1)` cycles of lengths `1,2,7`, and
their signed `U_(+1)` partners. On these 20 tagged states, the dynamics
permutation has cycle type `1^2 2^2 7^2`. Any commuting permutation must
preserve least period, because a bijection commuting with iteration
preserves the least positive return time.

On a pair of length-j cycles, a commuting permutation can rotate each
cycle independently and can swap the cycles; there are `2j^2` choices.
Thus the full centralizer on this specified finite core is

```text
S2 x (C2 wr S2) x (C7 wr S2),
order 2*8*98=1568=2^5*7^2.                             (G14)
```

It contains no element of order three. Negation with parameter reflection
is one commuting involution matching the three pairs. Three known cycles
is a count of components, not an orbit of period three. This statement is
only about the indicated finite core; it imposes no classification on
unknown cycles elsewhere. A sixfold action on some other finite object
needs its own map and preserved predicate before it can be transferred.

## 8. Connection contract and reproduction

The raw connection (G2) preserves arrows on sublattices and loses valuation
guards. The accelerated connection (G4) preserves the guards and allows
only odd dilation. The restored even connection (G9) preserves an explicit
shell instead of discarding it. The shortcut connection (G13) preserves
both branches after shifting both. These are different exact statements;
none is a cycle-completeness or convergence theorem for odd parameters.

[The exact script](../../04-computation/experiments/glued_xor_20260921_conjugacies.py)
and [JSON](../../04-computation/experiments/glued_xor_20260921_conjugacies.json)
declare every finite box. Run from the repository root:

```text
python 04-computation/experiments/glued_xor_20260921_conjugacies.py
python -O 04-computation/experiments/glued_xor_20260921_conjugacies.py
```

Controls include independent finite affine classifiers, raw sign-pair
examples, both minimal fixed-point hostiles, valuation shells through
five, all three cases of (G8), exact even-parameter iterates, and both
parity branches of the shortcut maps. All checks use explicit exceptions.
The universal classifications rely on the proofs, not the finite search.

Normal and optimized runs agree. The finite classifiers test 5,888
accelerated candidates, 15,210 full odd-part candidates, 3,240 ordinary
shortcut candidates, and 10,800 swapped shortcut candidates. The script
also constructs all 1,568 commuting permutations of the known core;
only the identity has cube equal to the identity.

Independent peer audit: **PASS after correction.** The initial trajectory
wording omitted nonfixed even preimages of a fixed odd point; `O_(-2)(2)=1`
is the minimal witness, now retained with its entire basin formula.
A preliminary 22-state count was also corrected to `2(1+2+7)=20`;
the cycle type and centralizer formula were unchanged. The peer checked
the revised proofs and independently replayed `main()` against the JSON.
