# Binary LFSR graphs, primary-height tomography, and supernatural clocks

2026-10-20. **PROVED within the stated elementary scopes; FINITE-EXACT controls.** Historical novelty is not claimed. The classical cycle-structure ingredients are attributed below. The repo-specific contributions are the exact information-loss classification and its sharp degree-18 hostile, the typed inverse-limit construction, and the source-adjunction application in the companion note.

## 0. Inheritance and connection contract

Closest proved carriers:

- `01-canon/theorems/THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier.md`: a full Hasse jet faithfully retains the finite sheet module, and translation is triangular unipotent.
- `01-canon/theorems/THM-2043-period14-parity-hasse-jet-completeness.md`: complete local coordinates cannot recover information already discarded before reduction.
- `01-canon/theorems/THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar.md`: an integer recurrence modulo a modulus is eventually periodic, and is purely periodic when its constant coefficient is a unit.
- `01-canon/theorems/THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase.md`: full Fourier phase retains a periodic Boolean observable, while its density does not.
- `05-knowledge/results/tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md`: source adjunction, source-free cores, and Hasse digit clocks. The companion `TOURNAMENT_SOURCE_ADJUNCTION_ODOMETER.md` supplies an actual operation-preserving bridge.

The user-supplied plane-coloring blueprint, §8, constructs a finite observable quotient before using a characteristic polynomial. The same order of operations is essential here: state dynamics, an observed output sequence, and the full module may have different minimal polynomials. The separate observable/Hankel packet handles that quotient. The supplied countable-coloring blueprint does not identify its “linear graph” hypothesis with an F2-linear state transition; no such identification is used.

**Map:** cyclic F2[x]-module → finite functional graph → period counts → supernatural clock after a specified compatible limit.

**Preserved:** the first arrow keeps the successor relation; period counts completely classify invertible unlabelled functional graphs; the limit clock keeps the closure of powers of the specified operator.

**Lost:** addition, polynomial labels, primary exponents beyond dyadic observations, a chosen scalar output, and the multiplicity of orbit closures. A clock never identifies the whole phase space.

**Restoration:** retain the F2[x]-module or full labelled primary Hasse jets; retain the output functional for sequence questions; retain the transition maps and topology for limits.

Concept board: cyclic modules; odd root orders; Hasse thickness; observed versus hidden coordinates; source adjunction; direct versus inverse limits. Hostiles: a zero constant term, equal characteristic polynomials for noncyclic matrices, conjugate-looking cycle profiles with different primary heights, and direct limits with only finite-order points.

### Classical attribution

Chang–Ezerman–Ling–Wang, [*The Cycle Structure of LFSR with Arbitrary Characteristic Polynomial over Finite Fields*](https://arxiv.org/abs/1612.07928), published in *Cryptography and Communications* 10 (2018), 1183–1202, already determines arbitrary-polynomial LFSR cycle structures. Xu–Zou, [*Linear Dynamical Systems over Finite Rings*](https://arxiv.org/abs/0810.3164), explicitly records the elementary-divisor classification over finite fields. The proofs below are self-contained; no priority claim is made for those ingredients. Limit-periodic/procyclic terminology is also classical; see [Gan, *An Exposition of the Connection between Limit-Periodic Potentials and Profinite Groups*](https://arxiv.org/abs/1001.1748).

## 1. The precise 2^n-state object

Let f be monic of degree n≥1 over F2. Set

    A_f = F2[x]/(f),       T_f(v)=xv.

Using the coefficient basis 1,x,…,x^(n−1), T_f is a companion transition. Its characteristic and minimal polynomial are both f. Changing shift-register convention may transpose or reciprocate the displayed matrix; the module convention here fixes that ambiguity.

The state graph has **2^n vertices**, each with exactly one outgoing edge v→T_f(v). The operator's matrix has size n, not 2^n.

    T_f is invertible  ⇔  f(0)=1.

Indeed x is invertible modulo f exactly when gcd(x,f)=1. When f(0)=0 there are trees feeding cycles; it is incorrect to describe every LFSR graph as a disjoint union of cycles.

Write

    f = x^h g,    g(0)=1,    gcd(x,g)=1.

Allow g=1. CRT gives

    (A_f,T_f) ≅ (F2[x]/x^h, multiplication by x)
                  × (F2[x]/g, multiplication by x).               (1)

The first factor is the length-h nilpotent shift. Every periodic point has nilpotent coordinate zero. Thus there are 2^(n−h) periodic states, and the maximum preperiod is h. The number of states of exact preperiod t is

    t=0:       2^(n−h),
    1≤t≤h:     2^(n−h+t−1).                                    (2)

This counts all basins together. It does not say every individual rooted tree is a complete binary tree of height h.

**Proof of (2).** The nilpotent factor has 2^t states killed by x^t for 0≤t≤h, hence 2^t−2^(t−1) of exact killing time t≥1. Multiply by the 2^(n−h) choices in the invertible factor. ∎

## 2. Fixed-point and cycle formulas

For every m≥1,

    F_f(m) := #Fix(T_f^m) = 2^deg(gcd(f,x^m−1)).                 (3)

To see this, put d=gcd(f,x^m−1), f=df1, and x^m−1=db. The kernel condition is f1|bv; since gcd(f1,b)=1, it is f1|v. There are 2^deg d such residues modulo f. This proof works without invertibility.

Let P_f(m) count states of **exact** period m and C_f(m) count cycles of length m. Integer Möbius inversion gives

    P_f(m) = Σ_(d|m) μ(m/d) F_f(d),
    C_f(m) = P_f(m)/m.                                        (4)

Transient states contribute no fixed points to (3). Dividing by m counts each cycle once.

Factor the invertible part as

    g = ∏_i p_i^(e_i),

where p_i are distinct monic irreducibles other than x. Write d_i=deg p_i, and let r_i be the multiplicative order of a root of p_i in F_(2^d_i)^×. The integer r_i is odd and divides 2^d_i−1.

Write m=2^a u with u odd. In F2[x],

    x^m−1 = (x^u−1)^(2^a),

and x^u−1 is squarefree. Therefore

    log_2 F_f(2^a u)
       = Σ_(i: r_i|u) d_i min(e_i,2^a).                       (5)

The exact operator order on the invertible factor is

    L(g) = lcm_i [r_i · 2^ceil(log_2 e_i)],                    (6)

with L(1)=1. Each irreducible-power factor divides x^m−1 precisely when r_i|u and 2^a≥e_i. A state attaining L(g) is 1 modulo g, so (6) is a realized orbit length, not only an upper bound.

For an individual v in A_g, its annihilator is g/gcd(g,v), and its period is the order (6) of that annihilator. This distinguishes one chosen seed from the clock of the full register.

### Full Hasse jets retain the lost thickness

If α is one root of an irreducible p of degree d, then

    F2[x]/p^e → F_(2^d)[ε]/ε^e,
    F ↦ F(α+ε) = Σ_(j<e) D^[j]F(α) ε^j                         (7)

is an F2-algebra isomorphism. Vanishing of the displayed jet means (x−α)^e|F over the splitting field; Galois conjugation gives p^e|F. Equal dimensions give surjectivity.

In these coordinates T is multiplication by α+ε=α(1+α^(−1)ε). Hence

    J_j(TF) = α J_j(F) + J_(j−1)(F),    J_(−1)=0.              (8)

The odd clock comes from α and the dyadic clock from nilpotent jet depth. This is a precise characteristic-two version of the translation carrier in THM-2201, not an identification of its characteristic-thirteen geometry with this graph.

## 3. Exact graph classification by dyadic primary-height profiles

For every odd r and integer a≥0 define

    H_(g,r)(a) = Σ_(i:r_i=r) d_i min(e_i,2^a).                 (9)

Only finitely many r occur, and each profile is eventually constant.

### Theorem A: graph tomography and its complete loss ledger

For companion systems f=x^h g and f'=x^h' g', the following are equivalent:

1. their directed unlabelled state graphs are isomorphic;
2. h=h' and their periodic cycle counts agree at every length;
3. h=h' and H_(g,r)(a)=H_(g',r)(a) for every odd r and a≥0.

**Proof.** (1)⇒(2): graph isomorphism preserves cycles and maximum preperiod. For (2)⇒(1), pair cycles of the invertible factors preserving their cyclic orientations, then take the product with the identical nilpotent shift in (1). Formula (5) gives (3)⇒(2) by (4). Conversely, cycle counts give F(2^a u) for every odd u, and integer Möbius inversion in the odd divisor poset gives

    H_(g,r)(a) = Σ_(d|r) μ(r/d) log_2 F_f(2^a d).            (10)

Thus (2)⇒(3). ∎

All roots of a fixed odd order r have the same degree ord_r(2), with the convention degree 1 for r=1. Consequently the graph remembers, for each r, only the dyadically sampled sums Σ min(e_i,2^a). It forgets which irreducible labels carried those exponents, and can even forget the exponent multiset itself.

The graph is still a complete dynamical invariant **as an unlabelled finite set with one successor**. The claim is not that equal cycle counts fail to give graph isomorphism in the invertible case; the failure is their inability to give F2-linear conjugacy or recover the primary module.

### Theorem B: the sharp degree-18 primary-partition hostile

Put

    p=x^3+x+1,       q=x^3+x^2+1,
    f=p^2 q^4,      g=p^3 q^3.                               (11)

Both irreducibles have root order 7. The companion graphs of f and g are isomorphic, with

    C(1)=1,   C(7)=9,   C(14)=288,   C(28)=9216,              (12)

and no other cycles. Their characteristic polynomials and their primary-height multisets at root order 7 differ.

**Proof.** The profiles divided by degree 3, for a=0,1,2,…, are respectively

    min(2,2^a)+min(4,2^a) = 2,4,6,6,…,
    min(3,2^a)+min(3,2^a) = 2,4,6,6,….                       (13)

Theorem A proves graph isomorphism. Formula (4) yields

    (2^6−1)/7 = 9,
    (2^12−2^6)/14 = 288,
    (2^18−2^12)/28 = 9216.

Their linear maps cannot be conjugate because their characteristic polynomials differ. ∎

**Sharpness.** Degree 18 is the smallest degree at which equal companion graphs can have different primary-exponent multisets after forgetting irreducible labels but retaining root order. This is a stronger ambiguity than choosing different primitive irreducibles, which already occurs in degree 3.

To prove sharpness, at any fixed r the profile at a=0 determines the number c of present irreducibles, and its final value determines the sum E of their positive exponents. For E≤5 the profile determines the exponent partition:

- c=1: the sum determines the sole exponent;
- c=2: the value at a=1 distinguishes the case with an exponent 1; if both exponents are ≥2, E≤5 leaves a unique partition;
- c≥3: after subtracting 1 from each exponent there are at most two units left. If there are two, the a=1 value distinguishes placing them on one exponent from placing them on two.

Hence a differing exponent partition requires E≥6. It also requires at least two distinct irreducibles of the same order. Degrees 1 and 2 each supply only one irreducible with nonzero constant term, so that order has degree d≥3. The common total register degree is at least dE≥18. Example (11) attains it. Adding other factors or a nilpotent x-factor cannot lower the bound. ∎

### Characteristic polynomial versus full module versus state spectrum

For a general F2-linear operator, the characteristic polynomial alone does not classify cycles: I_2 and the size-two unipotent Jordan block both have characteristic polynomial (x+1)^2, but their graphs have respectively four fixed points and two fixed points plus one two-cycle. Use all elementary divisors. For a direct sum of primary blocks, (5) simply sums over blocks, including repeated blocks for one p. In that wider category the stronger height-partition ambiguity already occurs in dimension 6: unipotent Jordan types (2,4) and (3,3) have identical cycle counts. The degree-18 sharpness above is specifically for cyclic companion modules, where two blocks of the same irreducible cannot occur.

For an invertible companion, the 2^n-by-2^n complex state-permutation matrix has characteristic polynomial

    ∏_(m≥1) (z^m−1)^(C_f(m)).                               (14)

Its traces recover F_f(m), hence its full complex spectrum recovers the unlabelled graph. This is a different matrix over a different field from the n-by-n companion. Formula (14) loses the F2-linear module. For arbitrary noninvertible functional graphs, the complex spectrum loses tree shapes; Theorem A's tree statement uses the special single nilpotent companion block, not a general spectral theorem.

## 4. Compatible limits and their exact supernatural invariant

A supernatural number is S=∏_ℓ ℓ^(s_ℓ), where each s_ℓ belongs to {0,1,2,…,∞}. It records divisibility: an ordinary integer d divides S when v_ℓ(d)≤s_ℓ for all primes ℓ.

Take monic polynomials with nonzero constant term and compatible divisibility

    f_1 | f_2 | f_3 | … .

Let A=lim← F2[x]/f_k, with reduction maps, and let T be multiplication by the compatible element x. This is a continuous automorphism of a compact F2-algebra. Put L_k=L(f_k). Then L_k|L_(k+1).

### Theorem C: operator and pointed-orbit clocks

The closure of the powers of T, in the topology of agreement on finite quotients, is

    G_S = lim← Z/L_k Z,
    S = lcm_k L_k = ∏_ℓ ℓ^(sup_k v_ℓ(L_k)).                 (15)

It acts on A. Its orbit through 1 is homeomorphic to G_S, with T corresponding to addition of 1. More generally, for v=(v_k) let ℓ_k(v) be the period of v_k. Then ℓ_k(v)|ℓ_(k+1)(v), and

    closure{T^n v:n∈Z} ≅ lim← Z/ℓ_k(v)Z.                  (16)

The supernatural number for this pointed orbit is lcm_k ℓ_k(v), which can be smaller than S. In particular, zero is fixed even when S has infinite prime exponents.

**Proof.** At quotient k, the powers of T form exactly a cyclic group of order L_k. Compatible residues of the time index define compatible powers on all finite quotients; the integers are dense in their inverse limit. The action is faithful because T^n1=x^n, and fixing 1 in every quotient means zero time modulo every L_k. For v, perform the identical argument on its finite cycles, pointed at v_k. Compactness makes the resulting continuous bijection onto the orbit closure a homeomorphism. ∎

This classifies the pointed orbit closure with its distinguished translation, not A as a ring, additive group, module, or whole dynamical system. An arbitrary unordered family of registers has no inverse limit until transition maps have been specified.

### Every supernatural clock occurs

Given any divisibility chain of positive integers L_k, take

    f_k=x^(L_k)−1.

Then f_k|f_(k+1), and multiplication by x has exact order L_k. This realizes every supernatural S by choosing a cofinal sequence of its finite divisors. The corresponding matrix is the cyclic shift on L_k bits. This is an existence construction, not a claim of minimal register dimension.

The universal choice L_k=k! realizes ∏_ℓ ℓ^∞. Squarefree f_k realize only odd clock orders: in characteristic two a squarefree invertible polynomial has no repeated-factor dyadic clock. Conversely, x^m−1 is squarefree for every odd m and has order m. Therefore

    lcm of all squarefree binary companion orders = ∏_(ℓ odd) ℓ^∞,
    lcm of all invertible binary companion orders = ∏_(all ℓ) ℓ^∞.  (17)

The missing 2^∞ is supplied precisely by unbounded primary thickness, for example f_k=(x+1)^(2^k).

Primitive irreducible cycles alone have lengths 2^n−1. Their supernatural lcm is the odd product in (17): every odd prime power ℓ^a divides 2^d−1 for d the multiplicative order of 2 modulo ℓ^a. These are multiplication-by-a-primitive-element clocks. Frobenius on F_(2^n), whose order is n, is a different linear operator; conflating them incorrectly deletes or inserts the prime 2.

### Direct limits are not inverse limits

The multiplicative group of F2bar is the direct limit of finite-field multiplicative groups along field inclusions. Every individual element has finite odd order; abstractly it is the direct sum of the odd-primary Prüfer groups. The inverse-limit procyclic group with supernatural order ∏_(ℓ odd)ℓ^∞ is compact and uncountable and has elements of infinite order. A shared supernatural divisibility inventory does not make these groups isomorphic. The duality between direct-limit characters and inverse-limit compact groups is a separate construction.

## 5. Exact controls and scope

Reproduce with:

    python3 verify_lfsr_tournament_clocks.py

The script has no third-party dependencies and uses exact integer/bit arithmetic. Its explicit verification gates remain active under optimization; normal and `python3 -B -O verify_lfsr_tournament_clocks.py` outputs agree byte-for-byte.

- All 2,046 monic binary polynomials of degrees 1 through 10, including zero constant terms: 1,398,100 states enumerated. Direct successor traversal agrees with the primary cycle formula and the complete preperiod distribution.
- 11,479 independent polynomial-gcd fixed-point checks. Within this universe, graph signatures and primary profiles give the same 621 classes.
- Degree-18 hostile (11): both full graphs are traversed, totaling 524,288 additional states. The counts in (12) are checked independently of the profile calculation.
- Source-adjunction towers: independent tournament canonicalization through order 5, exact F2 matrix ranks, direct state traversal through cutoff 4, and Pascal-row checks. See the companion note for the mathematical content.

The file `verification.json` records the exact universe and counts. Computation verifies stated finite universes; Theorems A–C and the degree-18 minimum have independent proofs above. No LRC(14), Collatz, Gilbreath, plane-coloring, or sphaleron-dynamics theorem is claimed.
