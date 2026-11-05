# The sharp dyadic residue image of self-converse C2 tournaments

Session 197. This is a separate constructive corollary of the joint-profile theorem and the diagonal classification. **It is a local residue theorem, not an exact-target sufficiency theorem.** It is not promoted until both dependencies and this corollary have passed independent review.

## 1. Named dependencies

1. The adjacent PROOF.md, SHA256 49bfc941828bf680dcaca1d5e4ec448f3b689d0f72c71f164bcb8b6bf0088c64, gives actual strong start-one sources with independently arbitrary odd residues of (D,U)=(H,R) modulo 2^k and exact (E,d)=(1,1), using at most 4k^2+3k+2 vertices. One or two universal sinks give (E,d)=(0,1) and multiply U by 2 or 4.
2. The frozen `session197-hybrid-constructions/self-converse-c2/PROOF.md`, SHA256 3e6d6e7344413504f299329fec05304600deef9432911f7aa9fecbd8fd4cb908, proves the physical converse-swap and diagonal classification. For every actual nontrivial start-one source V, its type-N diagonal N(V,V) is a nonregular strong self-converse tournament with exactly two Hamilton cycles. If V has (D,E,U,d)=(D,0,U,1), its Hamilton count is

       H=N(V,V)=D^2+2DU+3U^2-2U
               =(D+U)^2+2U(U-1).                  (1)

   It also proves that every nonregular self-converse strong C2 tournament has H=1 modulo 4. The only exceptional graph is the regular five-tournament R5, with H=15.

The diagonal identity is a specialization of the independently reviewed full-C2 bilinear normal form, published at f01df1766885c807aa4acae17685c2f4777f3ead. Source equality here means using the same actual rooted graph twice, not merely equal response quadruples from unrelated sources.

## 2. Odd square roots modulo powers of two

For k>=3 an integer a is an odd square modulo 2^k if and only if a=1 modulo 8. Necessity follows by squaring an odd integer. For sufficiency begin with x=1 modulo 8. If x^2=a modulo 2^j, j>=3, replacing x by x+t2^(j-1), t in {0,1}, changes its square by t2^j modulo 2^(j+1): the remaining square term is divisible by 2^(j+1), and x is odd. Choose t to correct the next bit. Induction supplies an odd root modulo 2^k. This is an effective bit-lifting construction.

## 3. Exact dyadic image and physical witnesses

**Theorem.** For every k>=2, the exact residue image modulo 2^k of Hamilton counts of strong self-converse tournaments with exactly two Hamilton cycles is

    {r: r=1 modulo 4} union {15 modulo 2^k}.        (2)

Each residue in the first set occurs for infinitely many distinct counts of nonregular graphs. Modulo two the image is the single odd residue. Modulus one is automatic.

**Proof of sufficiency.** First suppose k>=3 and r=1 modulo 4. There are two cases.

- If r=5 modulo 8, prescribe U=2 modulo 2^k. Then r-2U(U-1)=r-4 is 1 modulo 8.
- If r=1 modulo 8, prescribe U=4 modulo 2^k. Then r-2U(U-1)=r-24 is 1 modulo 8.

Take an odd square root x of this difference, and prescribe D=x-U modulo 2^k. This D is odd. The joint-profile theorem supplies an actual strong start-one source with H=D and R=1 modulo 2^k, and with endpoint and deletion counts one. Append one universal sink in the first case, or two in the second. The resulting actual nontrivial source V has E=0, d=1, the desired D residue, and U=2 or4 respectively. It need not be strong, which is permitted by the full-C2 source theorem.

Use exactly this rooted graph on both sides of N(V,V). Formula (1) gives H=r modulo 2^k, and the diagonal theorem supplies strongness, self-converseness and exactly two Hamilton cycles. If the joint-profile source has order at most 4k^2+3k+2, the final N graph has order at most

    8k^2+6k+9.                                    (3)

For k=1 or2, lift the prescribed allowed nonregular residue to 1 modulo 8 and use the k=3 construction. More uniformly, (3) holds with k replaced by max(k,3).

Arbitrarily many neutral loops in the joint-profile construction produce unbounded L, hence unbounded source Hamilton counts D. The positive form in (1), with U positive and even, is at least D^2, so the resulting counts are unbounded as well. This proves infinitely many distinct counts in each permitted nonregular residue class. R5 realizes the exceptional residue.

**Necessity.** The diagonal classification proves H=1 modulo4 for every nonregular graph. R5 has H=15 and is the only exception. Thus there are no other residues in (2).

## 4. Consequences and boundaries

A positive integer N passes every dyadic congruence obstruction for this self-converse C2 class if and only if

    N=1 modulo4, or N=15.                         (4)

For N=3 modulo4, (2) would require N=15 modulo every power of two, forcing N=15. This does not say that every N=1 modulo4 is realized exactly. Nor does it assert prime infinitude in any constructed ray or progression. The infinitely many counts above are unrestricted composite or prime integers.

In particular prime candidates other than the exceptional composite15 must be 1 modulo4; seven-times-prime candidates must have their prime factor 3 modulo4. These are necessary residue statements only.

The proof constructs only type-N diagonals. Type-M diagonals are not needed for residue surjectivity. Universal sinks are used inside the source before applying the strong C2 construction; no sink is appended to the final graph. The explicit converse map exchanges the two copies of V and fixes the one central N vertex, so the output is physically self-converse, not merely numerically symmetric.

## 5. Verification

`verify_self_converse_c2.py` uses the adjacent standalone `verify.py` only for the exact one-apex word constructor and literal source operations. It independently implements odd-square lifting, the type-N adjacency construction and its converse permutation. It checks every dyadic target residue through the declared finite bound and large exact bit examples. Small literal final graphs have their Hamilton paths and cycles counted separately by subset dynamic programming. The arithmetic and physical controls corroborate the all-order proof; no finite scan is used to establish (2).

Run the verifier normally and with `-O`. Neither run writes files. Keep this addendum separate from the frozen joint-profile proof and the frozen diagonal classification.
