# Converse preserves the full labelled descent kernel

**Corrected and independently audited; see ../promotions.json for the reviewed byte bindings.**
The original B1 sentence in RESULT.md said that recovering a named Phi-kernel
vector required the orientation gauge. This was false. Its proof established
`Phi(M(T^op))=Phi(M(T))^t` and equality of nullities but failed to audit the
additional row-zero tournament identity making Phi(M) symmetric.

The root's incoming-commit recovery supplied
`session321_selective_returns_and_fano_wiring_20261113/BSD_PARITY.md`, introduced
at `012251a1e7` and read at `0eab57c432`. Its SHA256 is
`ff315c0645bc37b42170b725a0aca77eebcd608401c6c25ee8d35ffbe74d5e0c`.
This is a stronger current proved mechanism; the bridge is credited as its
consumer, not as a new independent discovery of full-matrix invariance.

Explicitly, `M^t=M+J+I`, `MJ=0`, `JM=(n−1)J`, `J²=nJ` imply
`Phi(M)^t=Phi(M)` over F2. Consequently the converse preserves the entire
labelled Phi matrix and its embedded kernel. A joint-path reconstruction
up to converse recovers those objects, not only their dimension. M itself
still uses orientation. The renewed independent review additionally clarified
the parity of the full2n-coordinate Monsky embedding: every Phi-kernel vector
has even coordinate sum, so `Mop x=Mx+(n mod2)x`. Thus the embedding changes
by `(0,(n mod2)x)` and is already invariant at even n. RESULT.md now records
this exact distinction rather than assigning a universal gauge loss to it.
The B2 same-Q separation and B3 F8 subgroup obstruction are unchanged.

| Artifact | Superseded SHA256 | Repaired SHA256 |
|---|---|---|
| RESULT.md | `79b53e748063b365c537e1ee22e5a489353f64fa7a08d8bc40e14041de7262ba` | `29fc9ca464171bd1f5ef30602786cac2bbe222bda2b5b7cab9deb34e7e9699b8` |
| check_bridges.py | `10dad0c9e6b809e0652e4aacb0bde6b04c112c3832343c7dcef650019c3950cf` | `89fae3be3a595d5de0b7427485092c24a5099e9cba3894718a08746775f3da95` |
| Each normal/-O receipt | `8e77556bfe17bf5d61f4f4bf4d086ccb18e41b45ee34771c4fae2b5ab391896b` | `1d98df0f4bbc99b7f4dc6969eb8f6291a22b6ac21a6998c40b4e5ca3b28a82f9` |

Intermediate corrected RESULT.md hash
`a749d5f3c1d159b736905d2d1b40e2db3a3f65b9990b08e638855012295fd371`
established the full kernel identity but preceded the explicit even/odd
embedding distinction. It is superseded by the final hash in the table.

The repaired checker adds equality of the literal Phi matrices and complete
kernel sets under converse for all1099 labelled tournaments of orders1..5,
and after the independent two-slot decoder for every labelled tournament of
orders3..5. It also tests the kernel's even coordinate sum and the exact
embedding parity formula throughout orders1..5. Both modes pass **19,447**
always-active gates, byte-identically.
No optimization-disabled assertion is used. Earlier review hashes refer to
their earlier scope and do not silently certify this amended result.

Failure mechanism: a valid quotient proof can still overstate information
loss when additional identities on the source object have not been applied.
Here transpose was not a genuine remaining ambiguity on the constrained Phi
matrices. The correction strengthens the consumer while retracting the false
gauge claim.
