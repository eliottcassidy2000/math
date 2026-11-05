# A sharp residue-bank length, its carry ledger, and a nonlinear escape route

**PROVED elementary results; FINITE-EXACT nonlinear census; independently audited.** Session 192 follow-up, 2026-11-05. These are specified one-apex tournament presentations. The global Hamiltonian-path spectrum conjecture with exceptions `{7,21}` remains open. No exact-density conclusion is claimed.

## 1. Inheritance, objects, and the question

For a binary word `w=b_1...b_n`, let `T(w)` have chain vertices `1,...,n`, oriented increasingly, and apex `x` with `x->i` iff `b_i=1`. Write `H(w)` for its number of directed Hamiltonian vertex words. Here `n` is the **chain length**, not tournament order. For n>=2, strongness is equivalent to first bit 1 and last bit 0.

The inherited S25/S191 recurrence starts at `(H,A,p)=(1,0,0)` and is

    append 1: (H,A,p) -> (H,       2A+1,1),
    append 0: (H,A,p) -> (H+A+p,   2A,  0).

The adjacent-switch identity from S191 is

    H(u10v)-H(u01v)=A(u)+2-p(u)+C(v)-epsilon(v),        (1)
    C(v)=sum_(j:v_j=0) 2^(j-1),
    epsilon(v)=[v nonempty and v_1=0].

Empty prefixes have `A=p=0`. In particular every such switch strictly increases `H` by a positive even integer.

Closest proved mechanism: [the dyadic selector](../dyadic/PROOF.md), sections 2, 4, 5, uses independently insulated adjacent switches and a two-ended dyadic selector. Its bank has length `2k^2-2`. The canonical hostile is the all-order one-apex exclusion of 13, together with S26's positive missing-density theorem. The corrected near miss is substituting growing `k` into a fixed-`k` square-root estimate without its bank-length cost.

Concept board: independent switch packing; information budget; exact mixed derivatives; quotient versus carry; nonlinear modular automata. The anchor is simultaneous size/residue control. The niche is the exact carry ledger. The wildcard is whether nonlinear residue decoding can beat the independent-switch packing cost.

Two earlier unrelated repository ideas guide actual tests here:

- THM-4604, `01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md`, distinguishes a coarse invariant from the carry needed for integer incidence. The transfer here is methodological: retain the *actual tournament carry* after reduction modulo `2^k`, rather than import any Collatz theorem.
- THM-3512, `01-canon/theorems/THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary.md`, section 6, counts the information omitted by an emitted prefix. Here plain counting bounds the number of programmable residue bits, and the finite automaton retains `(H mod 2^k,A mod 2^k)`. No Rule 30 state-complexity theorem is asserted for tournaments.

## 2. The exact switch interaction and sharp packing theorem

Let `P_0=01`, `P_1=10`. For fixed words `u,v,z`, set `F(e,f)=H(u P_e v P_f z)`. If `g=|v|`, then

    F(1,1)-F(1,0)-F(0,1)+F(0,0)
      = 2       if g=0,
        2^g     if g>=1.                              (2)

Proof: apply (1) to the later switch. Changing the earlier `01` to `10` increases its prefix integer by `2^g`. If `g>=1`, the last bit of this prefix is fixed, giving `2^g`. If `g=0`, the prefix's last bit simultaneously changes from 1 to 0, adding one more. This touching-pair boundary is essential: the mixed difference of `0101` is **2**, not 1.

**Architecture.** A bank consists of fixed binary positions and `m` disjoint labelled adjacent pairs, each independently chosen as `01` or `10`. Every vector `e in {0,1}^m` is allowed. Require an additive residue law

    H(u B(e) v)-H(u B(0) v) = sum_i c_i e_i mod 2^k,  (3)

with coefficients independent of both exterior words `u,v`, and require all odd residues to occur. Throughout this theorem `k>=2`. This excludes overlapping switches, variable-length coding, nonlinear dependence among choices, and target-dependent exterior corrections.

**Theorem.** Every such bank has

    |B| >= k^2+2k-2.                                  (4)

This bound is attained for every `k>=2`.

Proof of lower bound:

1. There are `2^(k-1)` odd residues, so `m>=k-1`.
2. An additive law has zero mixed second differences modulo `2^k`. If two switch starting positions are separated by distance `d`, then `g=d-2`. Equation (2) excludes touching pairs and every `d<k+2`. Thus consecutive switch starts are at least `k+2` apart.
3. Let `a` be the number of fixed positions before the first switch. For `a>=1`, changing the exterior prefix integer by 1 changes that switch increment by `2^a`, so exterior independence requires `a>=k`. If `a=0`, its exterior-dependent term is `A-p`; exterior integers 0 and 2 change it by 2, which is nonzero modulo `2^k`. The suffix argument is the same: a nonempty fixed suffix of length `b` multiplies the exterior contribution `C` by `2^b`, while an empty suffix leaves `C-epsilon`, which can change by 2. Hence `b>=k` after the last switch.
4. Consequently

       |B| >= 2k+2+(m-1)(k+2)
            >= 2k+2+(k-2)(k+2)=k^2+2k-2.

The theorem is a lower bound for this precisely stated architecture, not for all residue-programming constructions.

## 3. A bank attaining the bound

Take `D_0=0^k`. Recursively for `j=1,...,k-1`, put

    q_j = (2^(j-1)-1-floor(A(D_(j-1))/2)) mod 2^(k-1),

and write the `k-1` low-to-high binary digits of `q_j` as `q_(j,0),...,q_(j,k-2)`. Define the next shared spacer by

    D_j = 1 (1-q_(j,0)) ... (1-q_(j,k-2)).             (5)

Thus `C(D_j)=2q_j`, and its first bit is 1. Form

    B_k(e)=D_0 P_(e_1) D_1 ... P_(e_(k-1)) D_(k-1).  (6)

Its length is exactly

    L_k = k+(k-1)(k+2)=k^2+2k-2.                      (7)

Each switch sees exactly the preceding and following `k`-bit spacer modulo `2^k`; all more distant exterior bits vanish. Its increment is

    A(D_(j-1))+2-last(D_(j-1))+C(D_j)
      = 2(1+floor(A(D_(j-1))/2)+q_j)
      = 2^j mod 2^k.                                 (8)

Therefore, for **every** exterior and every choice vector,

    H(u B_k(e) v)=H(u B_k(0) v)+sum_(j=1)^(k-1)2^j e_j
                     mod 2^k.                        (9)

This is a bijection onto the odd residues. For `k=1`, take the empty bank and `L_1=0`.

The shared-spacer bank improves the previous length by `k(k-2)` bits. For example, at `k=8`, it uses 78 bits rather than 126. Its quadratic leading coefficient is optimal inside the architecture in section 2.

## 4. The omitted integer carries are explicit

Fix any exterior `u,v`, let `H_0=H(uB_k(0)v)`, and let

    delta_j=H(uB_k(unit_j)v)-H_0,
    c_j=(delta_j-2^j)/2^k >=0.

The count `H` is a degree-at-most-two polynomial in the attachment bits, by the direct path-slot formula. Each disjoint pair is affine in its choice bit. Hence its exact expansion in the choice bits stops at degree two. The distance of switch starts `i<j` is `(j-i)(k+2)`, so (2) proves

    H(uB_k(e)v)=H_0+sum_j 2^j e_j+2^k K(e),           (10)

    K(e)=sum_j c_j e_j
         +sum_(i<j) 2^((j-i-1)(k+2)) e_i e_j.         (11)

All coefficients are nonnegative integers. In particular adjacent selected switches create an unavoidable unit carry in `K`; selecting two nonadjacent switches creates a larger carry. The residue observer deletes exactly this ledger. Thus reducing the independent bank's length does not make its physical increments equal their small residue increments.

This is the promised transfer from the earlier carry viewpoint: source and target both have a coarse arithmetic projection; the map here is `H -> H mod 2^k`; the lost coordinate is the integer `K(e)`; equations (10)-(11) restore it. There is no asserted conjugacy to a Collatz map.

## 5. Simultaneous control with growing k, and the certificate barrier

Import only S192's already-proved two-ended selector: for any middle length `L`, any `a,b>=1`, `n=a+b+L`, and any real `h in [2^(n-2),2^n]`, it supplies strong prefix/suffix choices with

    |H-h| <= 2^(n-a)+2^(n-b)+floor(n/2)+1             (12)

uniformly over **every** middle word of that length.

Combining it with (6)-(9), for every integer `k>=1`, every `n>=L_k+2`, every odd residue `r mod 2^k`, and every such real target `h`, there is a strong word of length `n` with

    H=r mod 2^k,
    |H-h| <= E_(n,k),
    E_(n,k)=2^ceil((n+L_k)/2)+2^floor((n+L_k)/2)
             +floor(n/2)+1.                          (13)

These quantifiers permit growing `k` with its explicit cost. For `k>=2`, the admissible range is `n>=k^2+2k`, and

    E_(n,k)=O(2^((n+k^2+2k)/2)+n).

For fixed `k`, this is the residue-specific square-root net inherited from S192, with a smaller constant. For example, if `k=alpha*sqrt(n)+O(1)`, `0<alpha<1`, the relative error is at most `2^(-(1-alpha^2)n/2+O(sqrt(n)))`, while the modulus is only `2^(alpha*sqrt(n)+O(1))`. This is still far from forcing exact equality.

There is a general obstruction to using **this unchanged error certificate**, even after replacing the bank by an arbitrary nonlinear binary middle. Programming all odd residues through a fixed exterior needs at least `2^(k-1)` possible middle words, so necessarily `L>=k-1`. But for `a,b>=1`,

    E=2^L(2^a+2^b)+floor(n/2)+1
      >=2^(L+2)>=2^(k+1).                            (14)

Thus this particular independent size bound can never certify `|H-h|<2^k`, the sufficient exactness threshold when `H=h mod 2^k`. This does **not** rule out sharper target-dependent estimates, correlated ends/middle, actual cancellations, or entirely different tournament families. At `k=1` the certificate is also insufficient directly.

## 6. Nonlinear modular automata escape the switch lower bound, finitely

The exact modular continuation state is `(H mod 2^k,A mod 2^k)`; `p=A mod2`. Keeping its nonlinear dependence on all middle bits permits a shorter architecture.

**Insulation lemma.** For every word `t` and every exterior `u,v`,

    H(u 0^k t 1^k v)
      = H(u 0^k 1^k v)+H(t)-1 mod 2^k.                (15)

Indeed, after `0^k`, the integer state is zero modulo `2^k` and the previous bit is zero. Reading `t` therefore increments `H` by `H(t)-1` modulo `2^k`. Reading `1^k` leaves `H` unchanged and sets the integer state to `-1 mod 2^k` and its final bit to 1. The suffix response is now independent of `t`.

So **any** fixed-length residue-complete word family of length `m` yields a fixed-exterior bank of length `2k+m`. There is no independent additive-switch requirement. Leading zeros and trailing ones leave `H` unchanged exactly, so all words of length at most `m` can be padded to length exactly `m`.

The exhaustive finite search, its minimal-length table, and its scope are recorded separately in `FINITE_SEARCH.md`. It certifies complete nonlinear banks through `k=14`, including a 47-bit bank at that modulus, compared with 222 bits for the sharp independent-switch architecture. No all-k linear bound is claimed.

## 7. Verification and scope

Run:

    python3 verify.py
    python3 -O verify.py

Both runs must give identical stdout and receipt bytes. No third-party dependency and no `assert` gates. The independent controls include the direct path-slot formula for every word of length at most 12, vertex-subset Hamiltonian-path DP through chain length 7, all mixed switch differences for three short exterior/gap words, the touching `0101` hostile, full bank cubes through `k=12` in four exteriors, exact carry ledgers, simultaneous selectors through `k=16`, and the nonlinear census through `k=14`.

The nonlinear census is a search over binary words, not unlabelled tournament classes. It concerns finite residues, not integer support density. The pre-existing all-order one-apex exclusion of 13 and S26 positive missing density remain intact. The exact `{7,21}` spectrum is not settled.
