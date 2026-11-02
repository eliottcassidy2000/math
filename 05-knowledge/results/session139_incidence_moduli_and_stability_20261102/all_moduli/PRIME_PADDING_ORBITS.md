# Exact source-padding phase cycles at odd primes

**Status: PROVED for the formal parity-function space at every odd prime; actual tournament counts are established at3,5,7,11 by the audited complete-state theorem.** No new search for actual states is performed. The existing main proof and source pins are unchanged.

## 1. The carrier and all-prime fixed-point formula

Fix an odd prime p. The exact prime contextual state of a nonempty tournament T is its path-cover polynomial function `f=F_T:F_p->F_p`. It belongs to

    X_p=V_+ union V_-,
    V_e={f:f(0)=0, f(-x)=e f(x)},       e in {+1,-1}.

The two spaces intersect only at zero, so `|X_p|=2p^((p-1)/2)-1`. Source padding `P(T)=K1->T` acts by

    (Pf)(x)=x f(x).

This is a permutation of X_p: P exchanges the two parity spaces and has inverse P^(p-2). Its order is exactly p-1, witnessed by the actual transitive seed f(x)=x.

**Fixed-point theorem.** For every positive integer d,

    Fix_p(d)=#{f in X_p:P^d f=f}
      =1                              if d is odd,
      =2p^(gcd(d,p-1)/2)-1            if d is even.             (P1)

**Proof.** An odd power switches parity, so a fixed function must lie in both spaces and hence vanish. For even d, a fixed function may be nonzero exactly where `x^d=1`. There are gcd(d,p-1) such nonzero field elements, paired freely by x↔-x. For each parity, values at one representative of each pair are arbitrary. Count the two spaces and subtract the zero function counted twice. QED.

The parity restriction is indispensable: in the full space of all functions vanishing at0, the nonzero indicator of{1} is fixed by P. It is not parity-pure and is not an admissible tournament state. Likewise, retaining two separate zero parity tags would overcount the fixed state.

## 2. Exact periods and actual prime cycle tables

Let `E_p(ell)` count states of least padding period ell. It vanishes unless ell divides p-1. For each divisor ell,

    E_p(ell)=sum_(d|ell) mu(ell/d) Fix_p(d).                   (P2)

This is ordinary Möbius inversion of `Fix_p(d)=sum_(ell|d) E_p(ell)`. Thus E_p(1)=1, all odd periods greater than1 are absent, and the number of cycles of length ell is E_p(ell)/ell.

There is also a direct support proof. Choose representatives `1,...,(p-1)/2` of the opposite pairs. A nonempty support S contributes `2(p-1)^|S|` functions, all with least period

    lcm_(x in S) [2, ord_p(x)].                              (P3)

This independently recovers(P2) without enumerating function values.

The audited complete contextual-state theorem in `PROOF.md`, Section5, makes the following **actual all-tournament contextual** cycle tables exact:

| p | least period | states | cycles |
|---|---:|---:|---:|
|3|1|1|1|
|3|2|4|2|
|5|1|1|1|
|5|2|8|4|
|5|4|40|10|
|7|1|1|1|
|7|2|12|6|
|7|6|672|112|
|11|1|1|1|
|11|2|20|10|
|11|10|322080|32208|

All these classes have actual width-three representatives. At other odd primes, (P1)–(P3) classify the entire **formal** parity-function space; the actual tournament subset is a union of these padding cycles, but equality with the full formal space is not proved here.

The10 length-two cycles at p=11 and10 length-four cycles at p=5 count cycles of **source-padding on contextual residue classes across all orders**. They cannot be identified with the ten merged five-vertex metagraph nodes through F alone: the independently verified order-five census has only **eight distinct full path-cover profiles**, so any map factoring through F has at most eight values on those ten nodes. An additional sidecar would be necessary. Sources: Session138 `independent_zero_framework_later/AUDIT.md` and `zero_framework_later/minimality_literal.out` under `05-knowledge/results/session138_stability_and_actual_contexts_20261102/`; the ten merged-node count is also recorded in `01-canon/theorems/THM-4857-the-gf2-automaton-of-a-tournament-classifies-the-merged-metagraph-exactly-through-order-five-circulants-are-additive-cellular-automata-and-paley-tournaments-give-quadratic-residue-codes.md`. This is a cardinality obstruction, not a numerical analogy.

## 3. The explicit all-modulus padding matrix

Write `d_c=c!(c-1)!`, `kappa(M)=min{c:M|d_c}`, and

    r_c=p_T(c) mod m_c,       m_c=M/gcd(M,d_c),
    1<=c<kappa(M).

The identity `x(x)_c=(x)_(c+1)+c(x)_c` gives the exact source-padding update

    r'_c=c r_c+r_(c-1) mod m_c,       r_0=0.                 (P4)

It is well-defined: m_c divides m_(c-1), so the predecessor coordinate has enough precision. In weighted coordinates `q_c=d_c p_T(c) mod M`, the same lower-bidiagonal matrix is

    q'_1=q_1,
    q'_c=c q_c+c(c-1)q_(c-1) mod M,       c>=2.             (P5)

No division by a nonunit occurs. The sharp universal law from `PROOF.md`, Section3, is therefore an explicit finite matrix identity on actual contextual states:

    P^(A(M)-1+L(M))=P^(A(M)-1),                             (P6)

with the previously proved minimal index A(M)-1 and period L(M). At a prime, A=1 and L=p-1, agreeing with the permutation cycles above. Formula(P6) also holds on the full formal reduced-coordinate product: lift a tuple to an integer falling-factorial polynomial, expand into integer monomials, and apply the transitive clock termwise. This algebraic extension does not assert that every formal tuple has a tournament representative.

## 4. Compact analytic verification

`verify_padding_orbits.py` compares Möbius fixed-point counts with the independent support/order formula(P3), checks integrality and cycle totals, and checks(P4)/(P6) on formal coordinate-basis vectors. It includes the parity-loss and doubled-zero hostiles. It does not enumerate additional tournaments or promote formal state counts at new primes.

    python3 -B verify_padding_orbits.py > padding.normal.json
    python3 -B -O verify_padding_orbits.py > padding.optimized.json
    cmp padding.normal.json padding.optimized.json

Dependency: the frozen `PROOF.md` (SHA256 `74aa0348207b196ddc64712d34adebbb5552b01c198e562d5e02e0ffd4885a3c`) supplies exact contextual coordinates, the sharp all-modulus clock, and actual completeness only at3,5,7,11. This corollary changes neither integer-spectrum status nor real stability bounds.
