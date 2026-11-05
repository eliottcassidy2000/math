# Complete ready-to-ready classification for the rooted growth grammar

Session 195, 2026-11-05. **PROVED all-order structural theorem.** This strengthens the separate one-gadget binary-return barrier. It concerns a particular marked interface, not unrestricted scalar Hamilton-path coverage. The `{7,21}` spectrum remains open.

## 1. Statement

Use the physical state `(h,e,R,d)`, with `h=H(T)`, `e=e_r(T)`, rooted two-path count R, and `d=H(T-r)`. A state is binary-ready with sign σ if

    e=d, R=h+d+σ, σ in {-1,+1}.

The allowed elementary operations are:

- K: add a vertex that beats only the old root and retain that root;
- D: add a universal sink and retain the root;
- any root-switch gadget G: form T->G, reverse one arc from the old root r to the new gadget root b, and mark b.

Gadgets may have arbitrary orders and may be chosen adaptively with full knowledge of the input and target. They need not belong to a fixed or finite bank.

**Theorem.** Start at any physical binary-ready state with d>=3. A finite word in these operations ends binary-ready if and only if every operation in the word is K. In that case its sign and deletion count are unchanged and, after k operations,

    H=2^k(h+σ)-σ.                                 (1)

Thus even arbitrarily many root switches while outside the ready class cannot repair this binary-ready interface. A change of state, reconstruction of old arcs, or another operation outside this grammar is necessary to rearm it.

## 2. Dependencies and exact formulas

The companion [arbitrary-gadget binary-return barrier](../binary-return-barrier/PROOF.md), producer proof SHA-256 `a4f2e2517d410a29a4b07dca70ab57bfe698a97e5927689ff8bcdb9380150389`, proves these elementary facts and their structural prerequisites:

1. `0<=e<=d<=h<=R` and `R+e>=h+d`.
2. If e=d, then either the root is a universal sink, with `(h,e,R,d)=(d,d,d,d)`, or h>=3d. This uses the proved Session 194 weak-deletion ratio after SCC factorization.
3. From any physical endpoint-complete state e=d with d>=3 and offset `τ=R-h-d` not equal to ±1, no K/D word can end binary-ready.
4. From a binary-ready state with d>=3, one arbitrary root-switch gadget followed by any K/D word cannot end binary-ready.

Write

    c=R-h, f=R+e-h-d.                              (2)

Both are nonnegative. The exact elementary formulas give

    K: (e,c,f,d) -> (d,f,f,d),
    D: e'=0, c'=h+2c, f'=h-d+2c, d'=d.            (3)

Counts h never decrease, because every operation preserves the old tournament as an induced subtournament.

For a rooted gadget let its profile be `(g,s,L,E,Q,B,delta,d_G)` as in the companion proof and put

    U=Q+E-g-d_G, V=B-L, W=s-delta.

One has U,V,W>=0 and the exact identity

    V=sum_{u in N_G^+(b)} L_u(G-b).                 (4)

Consequently a non-sink gadget root has V>=d_G>=1. A sink gadget root has U=V=W=0 and g=E=Q=d_G. The singleton is a sink but has s=delta=1; every non-singleton sink has s=delta=0.

The root-switch outputs satisfy

    e'=E h-delta e,
    c'=(Q-g)h+s e+V R,
    f'=U h+W e+V R,
    d'=d_G h.                                    (5)

The remaining profile inequality needed below is `Q>=g`.

## 3. A small-endpoint gap with no graph-order cutoff

**Lemma 1.** For every rooted physical tournament with h>=9,

    h-e is either 0 or at least 3.                 (6)

**Proof.** If the root is not in the last strong component, e=0. If its last strong component is a singleton, the root is a universal sink and e=h. Otherwise let C be its nontrivial last strong component and let a be the product of the earlier component counts. Then h=a H(C) and e=a e_r(C).

If |C|>=4, every other vertex of C can end a Hamilton path, so `H(C)-e_r(C)>=|C|-1>=3`. If |C|=3, then C is a directed triangle, giving H(C)=3 and e_r(C)=1. Since h=3a>=9, one has a>=3 and h-e=2a>=6. QED.

Endpoint positivity in a strong tournament follows from endpoint-preserving path insertion. This argument does not use a fixed-order census or a minimum-order theorem.

## 4. The invariant gap class

Let the gap class consist of physical rooted states satisfying

    h>=9, d>=3,
    each of e,c,f belongs to {0} union [3,infinity). (7)

This class can contain binary-ready states. For example `(17,3,21,3)` has `(e,c,f)=(3,4,4)`. Invariance alone is therefore not the exclusion argument; the last-root-switch lemma in section 5 is essential.

**Lemma 2.** The gap class is invariant under K, D, and every rooted gadget switch of arbitrary order.

**Proof.** For K this is immediate from (3): e'=d>=3, c'=f, f'=f, and d is unchanged.

For D, c'=h+2c>=9. If c>=3, then `f'=h-d+2c>=6`. If c=0, the physical cone gives `f=e-d>=0` while e<=d, so e=d and f=0. The endpoint-complete structural dichotomy from section 2 gives either h=d or h>=3d. Thus `f'=h-d` is either zero or at least 2d>=6. This also treats all neutral sink histories, and e'=0.

For a non-sink gadget root, (4)-(5) give c',f'>=V R>=h>=9. The gadget is non-singleton, so e'=E h is either zero or at least nine. For a non-singleton sink gadget, c'=f'=0 and e'=g h>=9. For the singleton gadget, (5) gives c'=e, f'=0, and e'=h-e; Lemma 1 supplies the necessary gap for e'. In every case d'=d_G h>=9. QED.

**Lemma 3.** Starting from any binary-ready physical state with d>=3, K preserves binary readiness, while either D or any root-switch gadget enters the gap class.

**Proof.** K gives equation (1) for one step. A binary-ready state cannot be a universal sink, since its offset is ±1 rather than -d. Therefore h>=3d>=9. Its endpoint count e=d is at least three. Its c=f=d+σ is at least two.

For D, `c'=h+2(d+σ)>=13` and `f'=h+d+2σ>=10`, while e'=0. For a non-sink or non-singleton sink gadget, the preceding proof works unchanged. For the singleton, c'=e=d>=3, f'=0, and `e'=h-d>=2d>=6`. QED.

The exceptional starting possibility c=f=2 is real: `[2,12,9,17,7]`, rooted at vertex 1, has `(h,e,R,d)=(13,3,15,3)` and sign -1. Lemma 3 explicitly handles it rather than assuming every ready input was already in the gap class.

## 5. No last root switch can return from the gap class

**Lemma 4.** From an input in the gap class, one arbitrary rooted gadget followed by any finite K/D word cannot end binary-ready.

**Proof.** Let z be the state after the gadget. The offset after a following K is

    τ_K=(U+V-d_G)h+W e+V c.                       (8)

If the gadget root is non-sink, all coefficients are nonnegative integers, because V>=d_G. Every term is zero or at least three by (7), so τ_K is zero or at least three. If the gadget root is a sink, (8) gives `τ_K=-d_G h<=-9`.

For a following DK the offset is

    τ_DK=A h+s e+(B+V)c,
    A=(g-E)+(L-E)+2U+2V>=0.                     (9)

Again every coefficient is a nonnegative integer, so the value is zero or at least three. A direct binary-ready gadget output would stay ready after K and has therefore already been excluded.

The output deletion count is d_z=d_G h>=9. If the first K in a padding suffix follows at least two D operations, its offset is `2^m R_z-h_z-2d_z>=h_z>=d_z>=9`. If it follows zero or one D, equations (8)-(9) exclude ±1. The companion padding lemma then excludes every delayed repair after that first K. A suffix ending in D has e=0, and the empty suffix was handled above. QED.

The possibility τ=0 in these algebraic bounds causes no problem: it is excluded from readiness, and the padding lemma applies to every offset other than ±1. In actual tournament states parity may exclude some formal zero cases, but that fact is not needed here.

## 6. Proof of the complete classification

A word of only K operations remains binary-ready and obeys (1).

Conversely, consider a word containing an operation other than K. Its prefix before the first such operation consists of K's and is binary-ready by the K transition. At the first D or root switch it enters the gap class by Lemma 3, and remains there by Lemma 2.

If the word contains no root switch, it is a K/D word with at least one D. After the first D-block followed by K, the companion padding inequality gives offset at least d>=3. Every later repair is excluded by the padding lemma; a final D is not ready.

Otherwise take the **last** root switch. If it is also the first non-K operation, its input is binary-ready and the already-proved arbitrary-gadget one-switch theorem excludes its entire K/D suffix. If an earlier non-K operation exists, its input is in the gap class, and Lemma 4 excludes its entire K/D suffix. These cases exhaust every finite word. QED.

## 7. Consequences and precise limits

From a fixed mature ready input (d>=3), the exact ready-to-ready count set is a single dyadic ray (1). From any fixed finite collection of such inputs, it is a finite union of these rays and has O(1+log X) members up to X. Its reciprocal sum converges, so it misses infinitely many ordinary-prime atoms and infinitely many 7p atoms. This statement permits arbitrary, even target-dependent gadgets during the attempted derivation.

The d>=3 restriction matters. Starting from a rooted C3 with `(3,1,3,1)`, the physical word A,D,K reaches `(17,3,21,3)`, where A is the singleton root-switch gadget. Thus an initial small-deletion reset can succeed. The theorem says that after reaching this mature ready state, no further non-K excursion can return.

This does **not** bound the scalar spectrum of the whole four-coordinate grammar. A target may be constructed at an unready root, using a different interface, or by a derivation that never enters a mature ready state. It also does not restrict operations that change old arcs, delete or rebuild a core, or reroot arbitrarily at existing vertices. The full `{7,21}` spectrum and an all-order scalar coverage proof for the enlarged grammar remain open.

## 8. Verification

`verify.py` is a standalone hostile-control program, not an all-target proof by extrapolation. It recomputes physical profiles on every labelled tournament through order six; checks Lemma 1, gap invariance, mature-ready entry, and the last-root-switch formulas against every one of the 67 audited physical gadget profiles; and checks arbitrarily long reachable histories within a stated R-bound, recording whether any non-K operation occurred. Singleton gadgets, endpoint-zero states, and neutral K/D histories are included explicitly. The all-order quantifiers follow from sections 2-6.

Run

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The exact finite control parameters and operation counts are reported in the receipt. They are not a claimed extension of the known Hamilton-path spectrum census.
