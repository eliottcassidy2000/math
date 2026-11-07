# Typed doubling recursions, and a physical lift of the even-response selector

Status: PROVED elementary deductions from inherited formulas; FINITE-EXACT independent subset-DP checks. No claim of unrestricted Chan–Pak coverage or literature priority. The inherited and newly deduced statements are distinguished below.

## 1. What actually doubles

Different operations have different types:

* Skew-Sylvester D sends vertex order n to 2n and sends each spectral quantity lambda²+1 to 2(lambda²+1). It does not double H.
* A single interior birth in the classical one-back-arc NT_s doubles delta=(H−1)/2. This has H' = 2H−1, but only on that family.
* The nested selector lift below doubles the single-arc response H(V)−H(U)=2N for every positive N, and preserves both endpoint graphs as induced subgraphs. It does not double either endpoint's H.
* Blue tiling Mode B sends order n to n+2; its dimension increases by n, its cardinality multiplies by 2^n, and its glue dimension increases by one, so the glue cardinality doubles.
* A numeric address map q→2q need not preserve adjacency, geometry, or Hamiltonian path count.

H(T) is always odd (classical tournament parity). Consequently H→2H is impossible. Also H→2H+1 cannot be universal since 3→7 is impossible, and H→2H−1 cannot be universal since 11→21 is impossible (repo THM-1370, exact slug `THM-1370-h-spectrum-omits-7-21-all-n.md`). Family-specific formulas avoid this obstruction but do not solve coverage.

## 2. A nested physical lift doubles every single-arc response

Use Session221 physical packet's convention: T(x) consists of an ordered chain v_0→...→v_(m−1) and apex a, with x_j=1 iff v_j→a. At a marked 10 at positions i,i+1, swapping to01 corresponds, after swapping vertex names v_i,v_(i+1), to reversing just the chain arc v_i→v_(i+1). The inherited response is

    H(T(x^swap))−H(T(x)) = 2N,
    N=1+L+R,
    L=sum_(j<=i−2) (1−x_j) 2^(i−j−2),
    R=sum_(j>=i+3) x_j 2^(j−i−3).

**Nested doubling theorem.** Suppose i>=1 and R=0. Insert a new 0 immediately BEFORE position i−1, the immediate left neighbor of the marked pair, in both endpoint words. Preserve every old chain vertex's relative order and apex incidence. Then the new marked pair is at i+1,i+2, and its response is exactly 4N.

Proof: old contributing left vertices have their distance to the marked pair increased by one, doubling every old summand. The new 0 sits two positions before the new marked pair and contributes one. The old immediate left neighbor remains immediate and contributes zero. Thus L'=2L+1. All right distances are unchanged, so R'=R=0. Therefore N'=1+(2L+1)=2N. Deleting the newly inserted chain vertex recovers the original endpoint tournament literally; no old arc is changed. The relabeling transposition that identifies the second word with the literal one-arc-flipped tournament acts only on the same marked two vertices, and commutes with the induced embedding. Hence this is a nested construction of pairs differing in one arc, not merely an equality between arithmetic evaluators. Iterating gives response 2^(k+1)N on pairs with exactly k extra vertices.

If both endpoint words begin0 and end1, both tournaments are strong. The insertion preserves these boundary values, even when it is before position0. Thus the nested doubling lift preserves strongness whenever the inputs are strong.

**All positive N have seeds.** Let d be the digitwise complement of the standard binary representation of N−1 (empty if N=1). For N>=2 use

    x=d1101,    y=d1011.

The marked10 in x is its penultimate pair before the last1. The last1 is response-null (immediate right neighbor), R=0, and L=N−1. Both words start0 and end1, so both endpoints are strong. Inserting0 immediately after d produces d01101 and d01011, the canonical seed for 2N because complementing the bits of 2N−1 appends0 to d. For N=1 use the separate strong seed0101/0011; insertion before its first bit gives00101/00011. Repeated insertion then gives a nested strong tower for N=1 too.

Examples (first pair, then its one-vertex extension):

    N=1: 0101 / 0011 has H=9 / 11;
         00101 / 00011 has H=19 / 23.
    N=2: 01101 / 01011 has H=15 / 19;
         001101 / 001011 has H=37 / 45.

The distinct N=2 seeds illustrate that N does not determine either absolute count. The theorem doubles a difference coordinate on a marked pair, not an absolute-count scalar. This directly connects Session221 exclusion tapes to a genuine recursive physical *2 operation; its extra one-bit input is the left binary charge.

Source: `05-knowledge/results/session221_computation_to_construction_20261107/physical/PROOF.md`, equations(1)–(3),(8),(11). Inherited response formula and selector are explicitly credited; the induced nested doubling lift is the deduction here.

## 3. Classical defect doubling in a single tournament

Let NT_s, s>=3, be the transitive tournament on0,...,s−1 with only the endpoint arc reversed: s−1→0. Then

    H(NT_s)=1+2^(s−2),  delta_s=(H−1)/2=2^(s−3).

Proof: one Hamilton path avoids the reversed arc, the increasing word (legal for s>=3). Every path using s−1→0 is uniquely an increasing run ending at s−1, followed by an increasing run beginning at0. Each of the s−2 interior vertices independently chooses the first or second run, giving2^(s−2) paths. Insert one interior chain vertex, retaining the endpoints and every old arc; the independent assignment gains one binary choice. Thus delta_(s+1)=2delta_s and H_(s+1)=2H_s−1. Two interior births give delta_(s+2)=4delta_s and H_(s+2)=4H_s−3.

This is inherited classical content, routed by THM-4748 `THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`, citing Moon1972/THM-250. Its role is an explicit realizable *2 witness, not a new theorem or a universal constructor. Examples H=3,5,9,17,33,65,... .

## 4. D(T): doubling order, and a necessary tiling-gauge correction

THM-447 `THM-447-skew-sylvester-doubling.md` defines

    M_D = [[M,M+I],[M−I,−M]],
    M_D² = diag(2M²−I,2M²−I).

Its canonical path is P followed by the reversed P in the second copy. Order doubles; lambda²+1 doubles; H has no such law. The singleton tower has orders1,2,4,8 and H=1,1,3,189.

THM-452 `THM-452-doubling-tiling-coset-law.md` has an indexing typo in section6: the five weights0,0,2,13,61 belong to orders1,2,4,8,16, not the four listed orders2,4,8,16.

There is also a convention issue. These weights and the stated law

    w_(2n)=2w_n+(n−1)+binom(n,2)

count FORWARD nonconsecutive arcs along the canonical path as1. The n−1 non-path twin arcs are all forward. Under THM-554's explicitly stated convention bit1=BACKWARD relative to the path, the weight is the complementary count and obeys instead

    wback_(2n)=2wback_n+binom(n,2),

with values0,0,1,8,44 at orders1,2,4,8,16. Global complementation preserves the equation (1+sigma)t=c_n, since sigma(1)=1. However, THM-452's description of its solution fiber as B+c_n is itself incorrect: c_n is sigma-fixed and therefore lies in B, so B+c_n=B. The correct fiber is

    {t:(1+sigma)t=c_n}=t_0+B,

where t_0 is any one actual D-image, not the defect c_n. This follows because subtracting t_0 turns the affine equation into the kernel equation. For n>=2, c_n is nonzero, so this fiber is disjoint from B; D-images are not blue. The transpose of a D-image cannot itself be a D-image: equality on the diagonal copy blocks forces the original input, whereas the nonzero cross defect changes the cross block. These two nonexistence statements require n>=2. For n=1 the tiling is empty, c_1=0, and the D-image is blue and is its own transpose. The report therefore distinguishes the invariant defect EQUATION from the incorrect inherited choice of coset representative.

The two diagonal tiling subtriangles carry the input and its mirror, while off-twin cross coordinates appear in complementary pairs. Thus for n>=2 the image lies in the affine fiber t_0+B disjoint from blue; the sigma defect is constant. The all-ones twin seam under the forward convention is not a Hamilton-count multiplier.

## 5. Mode B and where a true factor2 occurs

THM-477 `THM-477-the-blue-code.md` gives B_n=ker(1+sigma) and its dual im(1+sigma). Their quotient has dimension f_n=floor((n−1)/2), represented by fixed hypotenuse tiles. Under the centered embedding from n to n+2, the old tile(x,y) becomes(x+1,y+1). The old reflection is preserved. The new two-leg border has n−1 sigma-pairs and one fixed apex, giving

    B_(n+2) ≅ B_n ⊕ F_2^(n−1) ⊕ F_2,
    dim B_(n+2)−dim B_n=n,
    (B_(n+2)/B_(n+2)^perp) ≅ (B_n/B_n^perp) ⊕ F_2.

Thus the BLUE FAMILY cardinality grows by2^n, while its hypotenuse glue cardinality doubles. This is a direct, rigorous +2→*2 bridge, but on the quotient's cardinality, not H or the full number of tilings. New paired border bits lie in the dual, so only the new fixed apex survives in the glue.

## 6. Reproducible finite evidence

Run `python3 verify.py`. It writes `RESULT.json` and checks:

* All255 binary apex words of lengths0..7: subset Hamilton-path DP versus1+2C.
* All1100 base-path framed tournaments of orders1..6: defect equation, c_n being blue, exact n>=2 nonzero guard, representative coset, global-complement invariance, and image injectivity.
* All75 labelled tournaments of orders1..4: D validity, exact squared dominance law, odd H of both endpoints (D outputs through order8).
* The NT_s formula for all3<=s<=12 using subset DP.
* Both weak and strong selector lifts for all1<=N<=32,64 pairs in total, including literal induced containment.
* Blue-code ModeB border/reflection decomposition and both dimension increments for every2<=n<=100.
* Canonical tower paths and the FORWARD weight law through order16, and H through order8.

Finite tests do not establish Chan–Pak coverage. The symbolic proofs above establish the scoped all-order identities. The cheapest hostile controls are H parity, missing7/21, the two distinct response4 seeds, and the forward/backward weight mismatch.

## 7. Active-check and independent-audit receipts

All verifier conditions use explicit `require` calls that raise RuntimeError on failure; none rely on Python assert. Run `python3 verify.py` and `python3 -O verify.py` from this artifact directory. The normal and optimized runs produce byte-identical RESULT.json and console receipts; `REPLAY_RECEIPT.json` records the hashes and active-check count.

The independent session auditor additionally verified the nested operation at every eligible marked10 site in all words of lengths3..8 (354 cases, including31 strong). It also independently audited all1100 framed D-images through base order6 and tower weights through32, and caught the coset-representative error and n>=2 guard described above. Those audit results are attributed to the sibling audit packet rather than represented as checks implemented by this verifier.
