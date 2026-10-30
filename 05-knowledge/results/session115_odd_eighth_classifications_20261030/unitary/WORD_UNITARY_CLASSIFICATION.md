# Exactly when word-address transport preserves products

**Session 115, 30 October 2026. PROVED scoped classification.** This extends the word-address unitary in [the two-clock predecessor, §3a](../../source_address_two_clocks_20261004.md#3a-every-multiplier-has-the-same-linear-geometry-but-not-the-same-products). It is not a Collatz convergence theorem or a literature-priority claim.

## 1. The full classification

Fix odd positive multipliers p,q>=3. For a nonempty positive valuation word w=(a_1,...,a_m), let A(w)=sum a_i and s_i=a_1+...+a_i, with s_0=0. Its source address is

    r_w(q) = -sum_(i=1)^m 2^(s_(i-1)) q^(-i) modulo 2^A(w).

For each cost B the words are in bijection with the odd residues modulo 2^B. Thus all words of cost at most A, together with the empty constant label, index the character basis of

    E_A = all complex functions on Z/2^A Z.

Let U_(q->p) send the constant to itself and the character at r_w(q)/2^B to the character at r_w(p)/2^B. The predecessor proves that U is unitary, preserves the frequency filtration, and descends to a unitary U^odd on functions on the odd residues. These quadratic facts hold for every p,q; the following stronger properties do not.

**Theorem.** For A>=1, the following are equivalent:

1. p=q modulo 2^(A-1).
2. U on E_A preserves pointwise products.
3. U preserves complex conjugation, equivalently sends every real-valued function to a real-valued function.
4. U preserves nonnegative real functions.
5. U is pullback by a permutation of Z/2^A Z.

The corresponding five statements are also equivalent on the odd-residue function algebra. Whenever they hold, the operator is explicitly

    U(chi_(k/2^A)) = chi_(c k/2^A),
    U f(x) = f(c x),             c=q p^(-1) modulo 2^A.

For A=0 there are only constants and all properties hold. For A=1 the congruence modulo 1 is vacuous. For unequal p,q the last multiplicative depth is

    A_max = v_2(p-q)+1,

and the first failure is at A_0=v_2(p-q)+2. On the union of all finite character algebras, and hence for the associated global operator wherever the operations are defined, multiplicativity/real preservation/positivity/pullback holds **if and only if p=q**.

This classifies the actual finite algebras, including products that lower frequency cost. It does not treat the set of primitive frequencies alone as a multiplicatively closed algebra.

## 2. Multiplicativity: necessity and sufficiency

Every character of E_A is a power of the primitive character with frequency

    g=-q^(-1)/2^A,

whose word is (A). Its image has frequency -p^(-1)/2^A. Therefore any multiplicative U must multiply every frequency by c=q/p modulo 2^A. This necessity is valid on odd residues as well: the same primitive character generates their function algebra, now with the relation chi_g^(2^(A-1))=-1 rather than +1.

For A>=2 use the word (1,A-1). Its source numerator is -q^(-1)-2q^(-2). The difference between its actual image and the image required by multiplication by c is

    [-p^(-1)-2p^(-2)] - (q/p)[-q^(-1)-2q^(-2)]
       = 2(p-q)/(p^2 q).

It vanishes modulo 2^A exactly when p=q modulo 2^(A-1). Characters whose frequencies differ nontrivially modulo 1 are also distinct on the odd residues, since x=1 already distinguishes them. This proves necessity.

Conversely suppose 2^(A-1) divides p-q. In

    r_w(p) - (q/p) r_w(q)

its i=1 summand vanishes exactly. Every summand with i>=2 is divisible by p-q and by 2^(s_(i-1)), with s_(i-1)>=1; odd denominators are units. Hence the difference vanishes modulo 2^B for every word of cost B<=A. The map is therefore exactly frequency multiplication by c on the entire algebra. Since c is odd, x->cx permutes both all residues and odd residues. It preserves products, conjugation and positivity.

## 3. The first defect is a signed half-character

Suppose p!=q, put t=v_2(p-q)>=1 and A=t+2. Let tau_A be the full frequency permutation in numerator coordinates modulo 2^A. For every even k, earlier-depth multiplicativity gives tau_A(k)=ck. For every odd k,

    tau_A(k) = c k + 2^(A-1) 1[k=q modulo 4] modulo 2^A.       (1)

**Proof.** For a cost-A word, its i=1 address term again agrees with c times the q-address. Terms with i>=3 have prefix exponent at least 2, so their differences have valuation at least t+2=A. The i=2 term contributes only when a_1=1. In that case its difference is 2(p-q)/(p^2 q), of valuation A-1, hence exactly the nonzero half-modulus. The word starts with a_1=1 exactly when the second shortcut parity is odd: (qk+1)/2 is odd, equivalent to qk=1 modulo 4, or k=q modulo 4. This proves (1).

For odd k, exactly one of k and -k is q modulo 4. Therefore

    tau_A(-k) = -tau_A(k)+2^(A-1) modulo 2^A.                 (2)

The nonzero half-character is +1 on even seeds and -1 on odd seeds. Thus the real function chi_(k/2^A)+chi_(-k/2^A) has image, on odd x,

    chi_(tau_A(k)/2^A)(x) - chi_(-tau_A(k)/2^A)(x)
       = 2i sin(2 pi tau_A(k)x/2^A).

At x=1 this is nonzero and purely imaginary, since tau_A(k) is odd and A>=3. Thus U does not even preserve real functions, either on the full space or on the odd coset. The nonnegative function 1+cos(2 pi kx/2^A) likewise has a nonreal image at x=1, proving failure of positivity. Every permutation pullback preserves real functions, so pullback fails too. Since this first failing algebra is contained in all later E_A, none of these properties returns at a larger depth. This completes all equivalences.

A separate product witness is obtained by writing the character of (1,A-1) as a power of the character of (A). The same displayed difference shows that its actual image and its product-required image differ by the half-character. On odd seeds they are negatives.

## 4. The four odd eighths and their generalization

At p=5,q=3, t=1 and the first failure is A=3. The complete primitive-frequency permutation is

    1/8 -> 7/8,
    3/8 -> 1/8,
    5/8 -> 3/8,
    7/8 -> 5/8.

It agrees with the predecessor's quarter/eighth product hostile. The scalar required by multiplicativity is c=7 modulo8. The extra half-character occurs exactly for input numerators3,7, namely k=q modulo4. Opposite eighths receive different half-character corrections, explaining the loss of real values without a moment calculation.

For example, input cos(2 pi x/8) is real, but on odd x its image is

    i sin(2 pi 7x/8).

At x=1 this is -i/sqrt(2). Its norm is preserved by unitarity even though its pointwise algebraic meaning changes.

More generally, the first nonalgebraic shell is always the primitive odd numerators over 2^(v_2(p-q)+2), with exactly half receiving the signed half-character correction. The four eighths are the smallest member of this uniform family, rather than four universal probabilities or four universal dynamical states.

## 5. What is transported and what is lost

The source is the positive-word label set with its cost filtration; the target is the Fourier character basis. The label-preserving map transports full and odd-Haar quadratic inner products. Its failure to preserve addition of frequencies destroys products, real values and positivity at the explicitly classified depth. The extra information needed by a multiplicative consumer is the actual additive frequency law, not just the word label and cost.

The same parity labels can also define source-address permutations. Taking such a permutation as a **frequency reindexing** is different from composing functions with it on the **source variable**. The theorem does not deny the usual parity-coding conjugacies of shortcut maps; it classifies when this particular Fourier-basis unitary is itself a pointwise algebra map. No conclusion about a specified integer's orbit or Collatz convergence follows.

## Verification

`check_word_unitary.py` builds addresses in two independent ways: the closed affine sum over compositions and a recursive shortcut-parity decoder. It enumerates all odd multiplier residue pairs through depth 8, tests the congruence criterion and antipodal odd-coset relations, verifies the first-defect formula and its real-function hostile, and checks full multiplication tables through depth 5. Residue 1 represents legal multipliers congruent to 1, for example 1+2^A; the check is a finite-residue calculation, not a new global multiplier-1 claim. Explicit exception gates remain active under Python optimization. `results.json` records the actual universes and counts. The written argument covers every odd p,q and every depth; finite tests are regression controls.
