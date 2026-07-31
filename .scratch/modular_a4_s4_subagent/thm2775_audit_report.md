# Independent hostile audit of THM-2775

## Verdict

**ACCEPT**, with one nonblocking but important wording hardening:
Section 3 should say that the half-Hadamard forms recover the four
**centered root states**, not the four absolute roots.  The exact identity is

```text
H D = 4 I_4-J_4,
h_i=4r_i-(r_1+r_2+r_3+r_4).
```

The opposite-edge map has the diagonal translation line as its kernel.
Recovering the absolute roots additionally requires the common sum `e_1`.
The candidate already displays the correct formula, so this is a heading and
scope repair rather than a flaw in the matrix/group theorem.

Suggested audit label:

```text
modular-a4-s4-independent-hostile-audit-2026-07-28
(independently derived all four D3 matrices from root substitutions;
verified triangle relations, W(D3)/A4 closures, centered Hadamard sheet
actions, diagonal-translation hostile, THM-2769 Newton/local 110 row,
face-square inertia equality, scope, and normal/-O/stored/hash replays:
ACCEPT, with centered-root wording guard)
```

## 1. Independent derivation from the root variables

Put

```text
D = [ 1  1 -1 -1
      1 -1  1 -1
      1 -1 -1  1 ].
```

Then `d=D r` and `D D^T=4I_3`.  If a sheet substitution has permutation
matrix `P`, its induced opposite-edge matrix is forced by

```text
D P=M D,
M=(D P D^T)/4.
```

Applying this formula, without importing the candidate companion, gives
exactly

```text
X_S(d)=(d1,-d3,-d2),       Y_S(d)=(d2,d3,d1),
X_A(d)=(d1,-d2,-d3),       Y_A(d)=(-d3,d1,-d2).
```

These are respectively the substitutions `(12),(234),(12)(34),(123)`.
Integer multiplication gives

```text
X_S^2=Y_S^3=I,     ord(X_SY_S)=4,
(X_SY_S)^2=diag(-1,-1,1),

X_A^2=Y_A^3=I,     ord(X_AY_A)=3.
```

Exhaustive closure gives all `24` even-sign signed permutations for the
`S` frame and the `12` determinant-`+1` elements for the `A` frame.
Thus the images are `W(D3)=S4` and its `A4` subgroup exactly.

## 2. Half-Hadamard reconstruction and its sharp loss

For

```text
H = [ 1  1  1
      1 -1 -1
     -1  1 -1
     -1 -1  1 ],
```

direct multiplication gives `HD=4I-J`.  Hence the four rows transform
under the derived matrices as

```text
X_S:(12),     Y_S:(234),
X_A:(12)(34), Y_A:(123).
```

This proves the literal sheet-permutation statement on centered roots.
The hostile `r -> r+c(1,1,1,1)` leaves every `d_i` and every `h_i`
unchanged.  No absolute-root, affine-centre, or common-translation claim may
be inferred without retaining `e_1`.

## 3. Face square versus the THM-2769 boundary

For the proved hostile

```text
R_t(U)=U^3-4U^2+16tU-64t^2,
```

the Newton points are `(0,2),(1,1),(2,0),(3,0)`, giving root valuations
`(1,1,0)`.  Substitution `U=tZ` leaves the residual quadratic

```text
Z^2-4Z+16,
```

whose discriminant is `-48`, while the valuation-zero root at `U=4` has
derivative `16`.  Thus no ramification rescaling is hidden and the normalized
parity row is genuinely `110`.

Over the completed complex local field, both valuation-one units are
squares, so the unique local quadratic inertia flips their two square roots
and fixes the third.  Its matrix is

```text
diag(-1,-1,1)=(X_SY_S)^2.
```

This confirms the candidate's face-square/inertia identification after the
allowed relabelling.  It does not say the modular relation forces that
inertia to occur; THM-2769 is the separate geometric realization.

## 4. Scope audit

The candidate correctly stops at a finite marked generator frame.  It does
not:

- choose such a marking canonically for every quartic extension;
- retain the common root translation;
- force the `V4` boundary row to be `000`;
- construct the THM-2655 unit/class-group torsor;
- identify a Keller/Jelonek divisor or prove `JC(2)`/`DC(2)`;
- construct an LRC carrier or prove `LRC(14)`.

The full even code remains `{000,110,101,011}`.  The finite group relations
only preserve membership in that code.

## 5. Replay audit

Canonical THM-2775 normal, optimized, and stored outputs match byte-for-byte:

```text
script 758edc501f1f0cd6939e9e634611043cb42670e92286ba12c5b657c98f085e84
output 4c62a2409bd5ecfec67200d820bb03840ab5f4d650214e784a075e898fae01b5
```

The THM-2769 dependency also matches in normal, optimized, and stored modes:

```text
script 01d6b992cf58e1900eec0f88a22bff0ef32682954a45a00314459029c02c0fe2
output 1f2a73699c26ea65386b81a18890a8e9fdc8a7e16bf20e1f5ca0c2e4280f922b
```

The independent companion has no truth-bearing assertions and reproduces
its stored transcript in normal and optimized modes.
