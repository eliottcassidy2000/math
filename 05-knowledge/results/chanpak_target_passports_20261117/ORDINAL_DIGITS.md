# Exact ordinal digits for one-exit profiles

All graphs are tournaments. Let R be a nonempty tournament with a designated tag set M. A new vertex q beats precisely M in R. Add a module B beating q and losing to all of R. For R, define:

* J=H(R) and P=P2(R), where the latter counts ordered path pairs covering R and permits empty paths.
* W counts ordered covers (X,Z) with Z nonempty and starting tagged, X optional. Z_ends (written Zeta below) is the dual count: X ends tagged and is nonempty, Z optional.
* E counts full Hamiltonian paths ending tagged, S counts full paths starting tagged.
* U counts two-path covers with X ending tagged and Z starting tagged. K=J−U.

Since B has only one external exit, there are at most two consecutive runs from B in any full Hamiltonian path. Its one- and two-run skeleton counts are 2W+K and W: with one run, it can be terminal after a path of R+q not ending at q, or it can immediately exit to q followed by the remainder (including the case q ends the whole path). These give L=W−U and W+J choices respectively. With two runs, the first exits to q and then through one tagged vertex, and the second ends the full path. These give W choices. Consequently the full read is
    H = W P2(B)+K H(B).
In particular B=TT_k gives H=W*2^k+K.

## Exact transfer under marked ordinal composition

Append a tagged tournament A ordinally to R. Let its stats be (h,p,w,z,e,s,u), in the order H, P2, weighted starts, weighted ends, unweighted ends, unweighted starts, both-tagged covers. The transfer is
    J' = hJ,
    W' = pW + Jw,
    E' = Je,
    U' = hU + Ju + Es + Wz.
For W, condition on whether the marked start lies in R or A. In the first case A is freely split across both paths; in the second all R goes to the first path. For U, locate the ending tag of the first path and starting tag of the second: old/old leaves all of A to the second path; new/new leaves all of R on the first path; old/new gives E*s; new/old uses a weighted split of each subgraph, giving W*z. These cases are disjoint, including when one of the optional *partial* paths is empty. Thus
    K' = hK − Ju − Es − Wz.

Now suppose each ordinal block offers two tag sets with the same (h,p,w,z,e,s), but local K values differing by delta_i. The stack has a common J,W,E and all mixed terms are exactly unchanged. Each chosen local replacement changes stack K by delta_i times the product of h_j over all other j, regardless of other choices. Therefore
    H_epsilon(TT_k) = W*2^k + K0 + sum_i epsilon_i delta_i (J/h_i).
This is independent exact programming, without assuming that arbitrary contact toggles have vanishing mixed defect. For arbitrary common B, replace 2^k by P2(B) and multiply the intercept by H(B).

When the h_i are pairwise coprime and delta_i is a unit mod h_i, a prescribed residual after subtracting W*2^k+K0 selects its unique physical word, if any: reduce modulo each h_i; only zero or delta_i*(J/h_i) is permitted, and that chooses the corresponding bit. Then **must check the integer equality**; congruences alone admit integer shifts by J. If the h_i come from a fixed finite set and modules are repeated, the lattice step is at least J/lcm(h_i); repeated modules cannot furnish a constant-step growing interval by this mechanism.

## Checked physical family

The seven local fixtures are literal in replay_ordinal_digits.py, with pairwise coprime prime full counts 137,89,151,523,659,389,521. Using zero-based bank order (0,1,3,2,5,6,4) gives
    J = 128606212764656219
    W = 8187244045210794521306510
    K0 = -2883491385431820369445837.
Exactly 20 of the 128 independent bit choices have K divisible by 7. All produce distinct strong sevenfold outputs at each fixed k≥1. For k≥4 their quotients meet every odd residue modulo 16 at the same height; for every k≥1 their quotients also meet every residue modulo 9, separately. The width of the band is fixed, so these modular statements are not claims about an interval of consecutive quotients. One ordinary block out of the seven has a sink, but q reaches the first strong block via either of its tag sets, all later blocks are downstream, and all R reaches q via B. Strength also follows by explicit reachability in the witness compiler.

A smaller exact twin fixture uses R=(56,97,51,54,98,64,13) and q tagged masks 70 and 69. Here J=101 and W=553 are common, K=−35 and −21. The two one-exit graphs differ only in interchanging q's arcs toward vertices 0 and 1. Their counts are 7(79*2^k−5) and 7(79*2^k−3) at all k≥1, with literal twin-prime quotients 311 and313 at k=2. When B is arbitrary, the absolute difference is 14 H(B). This pair does not need equal full stackable boundary; its equality of W is enough for the two outer reads.

Replay: python replay_ordinal_digits.py and python replay_joint_digits.py. Computes local literal subset counts, compares directly on an independently assembled 15-vertex ordinal sum, checks exact CRT for all 20 addresses, and reconstructs the residue sets. No finite check is used in place of the transfer proof.

## Joint residue control with repeated banks

In the displayed bank order the seven absolute increments have gcd exactly 2, checked directly in the replay. Let J,W,K₀ denote the baseline full-stack values above, and let `D={Σ εᵢδᵢJ/hᵢ}` be its palette of even changes. Joining r copies of the entire bank gives a baseline `W_r,K_r`; choosing a word in each copy leaves W_r fixed and changes K_r by `J^(r−1)(d₁+⋯+d_r)`, with d_i∈D. Also 7 divides W_r since 7 divides W and W composes by addition of integer multiples. The literal graph stays strong: either tag choice reaches a vertex of the first strong block, which reaches every downstream block, and all exterior vertices return to q through the final module.

Suppose a prescribed modulus is 2m with gcd(J,14m)=1. In the group of even residues modulo 14m, the palette D contains 0 and its seven individual increments, which generate the whole group (their gcd is 2). The directed walks adding one palette entry can reach every residue in at most 7m−1 steps: the finite directed Cayley graph is strongly connected, as inverses of any generator are nonnegative multiples of it, and a shortest walk repeats no vertex. Pad with zero increments. Multiplication by J^(r−1) is invertible mod 14m, so for r=7m−1 and **any** fixed final height k≥1, choose contact words with

    K_r + J^(r−1)Σd_i ≡ 7c − W_r2^k (mod 14m)

for any odd c. This constructs a literal strong output divisible by 7 whose quotient is c mod 2m, irrespective of primality. Neither the core bank nor the final height depends on the desired residue once m is fixed; only tag bits change.

For a concrete efficient case, three repetitions already suffice for **all 72 odd classes mod 144 simultaneously** at each fixed k≥1. The one-, two- and three-copy palette sumsets modulo 1008=7·144 have 108, 500 and 504 residues respectively, and the last is exactly the even subgroup. The new replay implements residue decoding, verifies positivity and exact divisibility, and reconstructs a strong physical example. This improves the separate mod-16 and mod-9 statements above; it does not give an absolute window of consecutive counts.

The relative scale is itself constrained. For a stack of l blocks from a fixed finite bank with first W positive, W/J grows at least as a positive constant times 2^(l−1), by `W'/J'=(p/h)W/J+w/h` and p≥2h. Total switch width divided by J is at most l·max_i(|δ_i|/h_i). Since U≤W, the actual output is at least W(2^k−1)+J for k≥1. Thus the relative switch width shrinks exponentially with l and k. Joint residue control is a physical achievement here, but the remaining obstacle is choosing the absolute baseline.
