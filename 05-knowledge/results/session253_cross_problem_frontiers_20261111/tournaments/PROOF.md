# Session253: modular substitution germs and a constructive 3-adic cover

**Status: PROVED HERE**, using the scoped inherited path-cover substitution and reciprocity theorems. The accompanying bounded tests concern these new consequences, not an inherited census. No claim of external literature priority, universal integer realization, or Chan–Pak closure.

## 1. Recovery and live board

The closest unused mechanism is THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: maximal block runs make substitution polynomial in ordered path-cover counts. THM-3166, `THM-3166-falling-factorial-order-join-path-colour-transform.md`, diagonalizes ordinal join and supplies negative-colour reciprocity. THM-3181, `THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`, strengthens this to polynomial parity. These are inherited as scoped proved mechanisms, not audited again.

A corrected near miss is scalar multiplicativity outside transitive quotients (THM-2195): exact integer H loses repeated block visits. A native hostile is C3[C3,C3,C3], whose H=3159 is not a scalar product rule. The least-used sidecar is the factorial in the full run profile. It kills deep runs modulo a fixed prime power, even though it cannot kill them integrally.

Board: (i) full substitution path jets; (ii) parity reciprocity; (iii) finite residue states and costed witnesses; (iv) five-core Chan–Pak interfaces; (v) non-P-recursive quotient blowup growth. The successful operation is reduction modulo a prime power before attempting scalar compression. This leaves a true compositional invariant rather than a scalar approximation.

Two recovered close targets remain separate: THM-3226/3248's Q4 radius/non-P-recursiveness wall is not solved by eventual modular periodicity; THM-1950's determinant discrepancy inequality is already closed and is not a current proof target. The live new target is to optimize integer cost inside a complete modular constructor, not reprove discrepancy bounds.

## 2. Factorial jet theorem

For a nonempty tournament S let p_S(k) count unordered spanning covers by k nonempty directed paths. Set F_S(k)=k!p_S(k), H(S)=F_S(1), and Q_S(t)=sum_k p_S(k)(t)_k. All quantities here have integer coefficients/values. A singleton path is allowed.

Fix a prime p and a>=1, and set m=p^a and

    kappa = min {k>=1 : v_p(k!) >= a}.

The modular ordered jet is J_m(S)=(F_S(1),...,F_S(kappa-1)) modulo m.

**Theorem 1.** J_m is a congruence for arbitrary tournament substitution. In particular, replacing any factor S by S' with the same J_m in any finite iterated substitution context preserves H modulo m and the entire resulting J_m.

**Proof.** For a quotient tournament R on q vertices let

    W_R(x)=1^T X(I-A_R X)^(-1)1.

The inherited maximal-run bijection, with the d resulting paths ordered, is

    F_{R[S_1,...,S_q]}(d)
      = sum_{c_1,...,c_q>=1} [x_1^c_1 ... x_q^c_q] W_R(x)^d
          product_i F_{S_i}(c_i).

Coefficients are integers. Each F_S(c) is divisible by c!. Consequently a summand vanishes modulo p^a whenever any c_i>=kappa. More sharply it vanishes if sum_i v_p(c_i!)>=a. The surviving finite polynomial depends only on J_m of each factor. For output d>=kappa, F(d) itself vanishes modulo m. Applying the same rule inductively proves the context assertion. QED.

The naive state bound is p raised to sum_{1<=k<kappa}(a-v_p(k!)); actual tournament jets form a smaller subset. This is an upper bound on a finite interface, not an efficient algorithm for computing H of an arbitrary input tournament. For fixed quotient, prime, and exponent, however, the transition table is finite and independent of factor order.

The factorial lattice diag(1!,2!,...)Z^k is a genuine integral source of this truncation. It should not be identified with the Session252 seven-coordinate parity lattice without an explicit map between the observers.

## 3. Half-grid state over an odd prime

For odd p the entire modular jet is determined by order parity and

    Q_S(1),...,Q_S((p-1)/2) modulo p.

Indeed polynomial parity Q_S(-t)=(-1)^|S| Q_S(t) gives all evaluations on F_p from this half-grid and Q_S(0)=0. For k<p, Newton inversion gives

    F_S(k)=Delta^k Q_S(0)
          =sum_{j=0}^k (-1)^(k-j) binom(k,j) Q_S(j).

For k>=p the factorial makes F_S(k)=0 modulo p. Thus there are at most 2 p^((p-1)/2) states with the parity sidecar. The sidecar is necessary for this reconstruction in general.

This is a finite-field statement about evaluations and finite differences. It does not claim that (t)_p vanishes as a polynomial in F_p[t]; it vanishes on F_p.

## 4. Exceptional scalar-plus-parity compression at 3

**Theorem 2.** Let R be any tournament quotient and S_i any nonempty tournaments. Put epsilon_i=1 if |S_i| is odd, and epsilon_i=2 if it is even. With T_j denoting the transitive tournament of order j,

    H(R[S_1,...,S_q])
      = H(R[T_epsilon_1,...,T_epsilon_q]) product_i H(S_i) (mod 3).

**Proof.** Since Q_S(-1)=(-1)^|S| H(S) and (−1)_k=(-1)^k k!, modulo 3 only k=1,2 remain. Therefore

    F_S(2)=(1+(-1)^|S|)H(S) (mod 3).

The complete mod-3 jet is consequently H(S)(1,0) for odd order, and H(S)(1,2) for even order. These vectors are precisely H(S) times the jets of T_1 and T_2. Each factor occurs once in every product in the substitution formula, so the scalars H(S_i) factor out. QED.

Consequences:

1. A factor with 3 dividing H cannot be rescued by any substitution context, including contexts that are themselves iterated substitutions. Its zero jet is absorbing.
2. If all factors have odd order, the expression reduces to H(R) product H(S_i) modulo 3.
3. The complete quotient obstruction is the finite Boolean germ epsilon -> H(R[T_epsilon_i]) mod 3. It has 2^q evaluations, not an unbounded search over factor tournaments.
4. For R=C3, the germ is 0 on 111 and 222, and 2 on all six mixed parity words. This follows immediately from the balanced-content kernel: for one even block H=5, for two even blocks H=11, for zero even blocks H=3, and for three even blocks H=45.
5. Ordinally adjoining a singleton to any factor preserves its H and flips its parity. Thus a zero C3 germ with three 3-unit factors is neutralized by changing the parity of one block. This creates actual new integer outputs, not just a renamed residue state.

**Sharp hostile to a tempting generalization.** The four-vertex strong Q4 with arcs 1->2,2->3,3->1,3->4,4->1,4->2 has F=(5,22,36,24). Thus 5|H(Q4), but H(C3[Q4,T_1,T_1])=3*5+22=37. Divisibility by 5 is not substitution-absorbing. At p=5 the second half-grid coordinate Q(2) survives, explaining exactly the missing information.

## 5. An exact integer constructor and all 3-power residue classes

Here is the constructive continuation beyond a modular filter. Let S be any nonempty tournament, r>=0, and S=>T_r denote ordinal padding on the right by r singleton vertices (r=0 means S itself). Define

    B_r(S)=C3[S=>T_r,T_1,T_1].

**Theorem 3.**

    H(B_r(S))=H(S)+2^r Q_S(2),
    |B_r(S)|=|S|+r+2.

**Proof.** Only block-run contents (1,1,1) and (2,1,1) enter C3[U,T_1,T_1], giving H=3H(U)+F_U(2). Ordinal join gives Q_U(t)=t^r Q_S(t). Hence H(U)=H(S) and F_U(2)=Q_U(2)-2H(U)=2^r Q_S(2)-2H(S). Substitution gives the formula. QED.

Taking S=T_1 gives H=1+2^(r+1), on r+3 vertices. Taking S=Q4 gives Q_S(2)=2H+F(2)=32, so H=5+2^(r+5), on r+6 vertices. All these are actual strong tournaments because their quotient is strong and all blocks are nonempty.

**Theorem 4.** For every a>=1, these two explicit families together attain every residue class modulo 3^a. Each residue has a representative with r<2*3^(a-1), hence with at most 2*3^(a-1)+5 vertices.

**Proof.** The order of 2 modulo 3^a is 2*3^(a-1). To see this without an external input, 4=1+3 and an elementary induction yields v_3(4^(3^j)-1)=j+1; the order of 2 is even and this forces the stated order. Thus 2^r runs through every unit modulo 3^a in one period. The first family has residues 1+u for u a unit, exactly those not congruent to 1 modulo 3. The second has residues 5+u for u a unit, exactly those not congruent to 2 modulo 3. The two sets cover every class. QED.

One may select the branch deterministically: if N is not 1 modulo 3 use seed T_1 and solve 2^(r+1)=N-1 modulo 3^a; otherwise use Q4 and solve 2^(r+5)=N-5. The right side is then a unit and a solution r exists within one period. A finite discrete-log calculation supplies the graph address.

This is distinct from changing the odd-cycle fugacity to 3. THM-4950's ternary H_3 is another invariant; here H is always the ordinary Hamiltonian-path count, reduced modulo powers of 3.

**Remaining lift obstruction.** Counts grow exponentially in r. The theorem does not give 0<=H<3^a, or H=N for an arbitrary integer N. In fact these families have at most O(log X) actual values up to X, so residue completeness cannot be mistaken for integer-density or Chan–Pak coverage. This exhibits the same precise bounded-lift issue as Session252's Fermat modulus construction, now for every 3-power via a native quotient operation.

## 6. Uniform eventual periods for all fixed quotient constructors

Let lambda(p^a) be the exponent of the units modulo p^a. For every integer b and r>=a,

    b^(r+lambda(p^a))=b^r (mod p^a).

For a unit this is the defining exponent property. For a nonunit both sides vanish, since r v_p(b)>=a.

**Theorem 5.** Fix a tournament S. The full J_{p^a} of its r-fold ordinal join is periodic with period dividing lambda(p^a) for r>=a. Hence in ANY fixed iterated substitution context whose leaves are ordinal powers of fixed tournaments, H and J_{p^a} are separately eventually periodic in each exponent, with that same universal threshold and period. In particular this applies to arbitrary fixed quotient blowups by transitive blocks T_r.

**Proof.** The path-colour polynomial of the r-fold join is Q_S(t)^r. For k<kappa,

    F_{S^{=>r}}(k)
      =sum_{j=0}^k (-1)^(k-j) binom(k,j) Q_S(j)^r.

Apply the elementary power congruence termwise, and then Theorem 1 to each context. QED.

This makes modular search bounded independently of the large factor sizes. It coexists with the inherited non-P-recursive growth results: periodic residues do not bound complex singularities or imply an integer recurrence.

## 7. Tests and exact scope

Run `python3 probe.py` and `python3 -O probe.py` from this directory; both outputs agree. The script uses no assertion-dependent checks. An endpoint subset DP computes H, and an unordered set-partition DP computes path profiles.

- 870 direct vertex-level substitutions: quotients C3 and transitive T3, factors every labelled tournament on 1,2,3 vertices, restricted to total order <=7. Compare Theorem 2 against the parity skeleton. This is a bounded new-theorem control, not an isomorphism census.
- 162 low-jet period checks: seeds T1,C3,Q4; prime powers 2,4,8,3,9,5; nine exponents beginning at the proved threshold.
- All 8 C3 and 16 Q4 parity germs; the C3 statement is also proved analytically.
- A small search through labelled orders <=4 gives the prime-5 hostile, directly checked at six vertices.
- 12 direct padded-constructor checks, seeds T1,C3,Q4 and r=0,1,2,3.
- All residue classes for the two-seed construction for a=1,...,7. The all-a claim is arithmetic proof, not extrapolation.

No large inherited census, source re-audit, or Git mutation was performed. Read-only startup coordination could not run because the supplied checkout has no .machine-id; no registration was fabricated.

## 8. Concrete next work

1. Compile modulo 9 and 27 the finite path jets of alternative small core constructors. Modulo 3 the absorbing ideal is completely classified by H; modulo 9 higher runs can distinguish which zero residues are repairable outside a fixed factor context.
2. Pair each modular state with minimum graph order and actual integer output in a costed constructor search. The exact two-seed cover gives a baseline upper bound of 2*3^(a-1)+5 vertices; a useful advance must reduce output cost, not merely regain residues.
3. Compare the full path-jet transitions with the seven-state flank compiler only after explicitly deriving the map. The shared lattice language alone is not a transfer.
4. A conditional integer transfer is immediate: if for every target N a modulus 3^a>N and a constructed witness with H=N mod3^a and 0<=H<3^a can be supplied, then H=N. This note supplies the modular witness but not that output bound.
