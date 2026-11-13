# Exact target repairs with a universal chain entry

**PROVED scoped construction, using inherited endpoint laws; FINITE-EXACT
base certificates. Chan–Pak and arbitrary-target entry remain OPEN.**
The new contribution is a physical23-core connection and its general
prefix theorem. The R/N and retained-chain operations are inherited, as
are the small unrooted tournament classes themselves.

## 1. A single off-root arc controls a later obstruction

Use outgoing rows, root r=5, and the six-vertex strong tournaments

    K9  = (4,33,26,35,43,5),     z9  = (23,2, 9,60),
    K11 = (4,37,24,35,43,5),     z11 = (23,2,11,60).

Here z=(H,s,q,F), where s counts paths starting at r, q=H(T-r), and
F counts ordered spanning two-path covers with the first path nonempty
ending at r and the second possibly empty. Reversing the single arc
between physical vertices1 and2 exchanges these graphs. It keeps every
edge incident with r, the full count23, and the first birth target81.
Literal permutations, deletion permutations and cuts independently check
the two profiles. This is a physical point of the inherited current/next
target fibre Delta=(0,v,w,v), with v=0,w=2.

The direction useful after adding a prefix is **K11→K9**. An unprefixed
first-R gcd test favours the opposite direction3→1; that local preference
does not determine which source survives the later chart entry.

## 2. Arbitrary nonempty prefixes and a retained incoming port

Let L be any nonempty tournament and write

    h=H(L),       f=P_L(2),

where P_L(2) counts two ordered possibly empty spanning paths. Both h and
f are positive, h is odd, f is even, and f≥(|L|+1)h≥2h.
Put L before either six-core in an ordinal join. Retain r and any named
prefix vertex p. In both graphs p→r, and all prefix names and arcs survive
the one-arc change. The rooted states are

    z_q=(23h,0,qh,60f),              q in {9,11}.          (1)

Full counts and deletion counts multiply. A full path cannot start in
the terminal six-core. In a two-bin cover, both bins may traverse L,
giving f choices, followed by a rooted cover of the six-core. This proves
(1) directly, with its endpoint convention intact.

Add x with r→x→every other old vertex and retain r (operation R). The
two resulting graphs are strong and have the **same exact target**

    M=23h+60f,
    Rz_q=(M,qh,qh,120f-qh).                            (2)

Strongness follows because x reaches L, L reaches the strong six-core,
and that core reaches r→x. The same old root and the same newborn x are
available on both sides. The switch is still the one off-root arc.

## 3. Complete entry divisors and an infinite exact repair

Use the inherited physical L_k operation, k≥2, from
[CHAIN_REUSE](../h17_exact_repairs_20261112/constructors/CHAIN_REUSE.md):
k births at the old root, retain the newest chain child w, then birth y
at the old root and z at y; keep the live pair(y,w), activating w.
For input (H,s,q,F), set d=F-q, c=H-s+4q-F and v=2^(k-1). Its output is

    C(Q,S)=(Q+3S,S,Q,4Q+2S),
    S=vd,       Q=3S+c.                               (3)

Applying (3) to the actual sources(2) gives

    d_q=120f-2qh,       c_q=(23+4q)h-60f,
    gcd(Q,S)=gcd((23+3q)h,60f-qh).                     (4)

Indeed c_q is odd; remove powers of two from S, and add c_q to
60f-qh. This proves the exact divisor **for every k**, not a sampled
clock assertion. In particular

    D9  = gcd(50h,60f- 9h),
    D11 = gcd(56h,60f-11h).                            (5)

Assume the explicit prefix guard

    gcd(h,30f)=1.                                    (6)

Then D9=1. A common odd prime with h is excluded by(6); the only remaining
prime in50 is5, and60f-9h is nonzero modulo5. Likewise

    D11 = 7 if f≡h (mod7), and1 otherwise.             (7)

There are no hidden powers of7: once the h factor is removed,56 has only
one7. Thus K11→K9 preserves M and the root star, but repairs the entry
divisor7→1 exactly on the stated locus. At other admissible prefixes it
preserves primitiveness. This is a complete conditional prefix theorem.

After entry, the inherited chart map for any next ell≥2 is

    S'=2^(ell-1)(3Q+2S),      Q'=Q+3S',
    gcd(Q',S')=gcd(Q,S).                              (8)

The actual retained ports make every subsequent choice legal. Positivity,
odd Q, even S and strongness persist, so (6) supplies indefinitely reusable
primitive chain entry. This proves the operational guard, not unit readiness
F-H=1 or2. The chain changes the full target after the neutral repair.

## 4. A compensated third core removes the same obstruction

A third strong23-core, aligned at the same root, is

    K3=(14,32,26,2,43,13),       z3=(23,2,3,30).

It differs from K9 on exactly {0,1},{0,3},{3,5}. Insert a singleton tau
between L and K3. The prefix is now L→tau, whose pair is (h,2f). Hence
the completed source has profile(23h,0,3h,60f), and its first R still
has precisely target M in(2). Keep every old name, r, the prefix port p
and the common R newborn. The new singleton is explicit; individual
intermediate edits are not asserted H-neutral.

Its entry content is

    D3=gcd(32h,60f-3h)=1 under(6).                    (9)

This is a dyadic compensation of the **marked** response: one singleton
doubles f, while the six-core's F coefficient halves60→30. Equality of H
alone would not have proved equality of the next target. The factor32 in
(9) makes the primitive-entry test especially transparent.

## 5. Explicit unbounded families and exact target inversion

For L=TT_a, a≥1, h=1 and f=2^a. The common targets are

    M_a=60·2^a+23 =143,263,503,983,1943,3863,7703,…

K11 has content7 exactly when a≡0 mod3, since2 has order3 modulo7.
K9 and the compensated K3 enter primitively at every a and every k≥2.
At a1,k2, K9 yields the literal chart state(2603,444,1271,5972);
one further L_2 yields(57683,9402,29477,136712), again with content1.

The prefix can have unbounded nontrivial cores, rather than only an
ever-longer transitive block. Let L be any ordinal arrangement of b copies
of the nearly transitive R6 and a singletons, with a+b>0. The literal
R6 pair is(H,P2)=(17,224), so

    h=17^b, f=224^b 2^a,
    M_(a,b)=23·17^b+60·224^b 2^a.                    (10)

Because224=2^5·7 and17 is coprime to30·224, guard(6) holds at every
prefix step. This gives an unbounded strong-output grammar, a known
preparation, and entry into (8) at every selected(a,b). The b1,a0
target is13831. It is not a finite union of the original fixed-core
dyadic charts; no claim of density or of arbitrary-target coverage follows.
More precisely it is not even contained in finitely many rays A_i2^k+B_i
with fixed positive integer A_i and integer B_i. Fix b and let a increase.
If B_i differs from23·17^b, then the2-adic valuation of M_(a,b)-B_i is
eventually fixed, forcing k to be fixed and giving at most one sufficiently
large a. Thus some B_i must equal23·17^b for each b. Infinitely many
distinct offsets require infinitely many rays. This is the same valuation
mechanism as the inherited finite-atlas obstruction, applied to an actual
entered family.
For a requested N, test b≤log_17(N/23) and whether
(N-23·17^b)/(60·224^b) is a permitted power of two. This terminates and
returns an actual graph whenever N belongs to(10).

More broadly, [SELECTOR.md](SELECTOR.md) gives an exact finite search for
the entire guarded prefix language M=23H(L)+60P_L(2), with a complete
implemented small-target head. This selector is only for that language.
It does not prove every admissible odd integer has such a representation.

Reproduce: `python3 -B 05-knowledge/results/surplus_selector_20261113/repair23.py`.
