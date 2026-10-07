# Terminal proofs become guarded rays; variable-length rules cross a fixed cap

**Status:** PROVED elementary certificate transports and declared all-height
relative-coverage statements; FINITE-EXACT controls. The reset and sporadic
identities are inherited, not claimed new. Universal Collatz remains OPEN.

## Inheritance and objective

The closest mechanism is [the source-preserving selector](adaptive_boundary_selector_20261004.md)
and its learned descent cylinders, together with
[THM-4555, uniform switches at minus one](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md)
and [THM-4556, Mersenne debt states](../../01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md).
The hostile is a fixed finite word bank confronted with arbitrarily long
leading-one runs. The near miss is a new dispatcher action whose actual
edges were already present elsewhere. The least-used sidecar is the
literal partial path, now retained by the [terminal-basis compiler](collatz_terminal_basis_20261007.md).

Anchor: ground the finite missing children. Niche: lift authenticated terminal
words to larger source families. Wildcard: remove the fixed word cap from
known collision identities. Board: **original source / native residue /
strict cutoff / partial path / reusable terminal / variable run length**.
The relevant method cards are controlled forgetting with a sidecar and
re-evaluating a certificate after a fibre-changing operation.

## 1. Exact descent rays including their finite initial cut

For a nonempty positive valuation word w of length L and total A, write

    U_w(n) = (P n+B)/Q,    P=3^L, Q=2^A.

Assume Q>P. Its exact native source cylinder is

    n = r (mod 2Q),    r=(Q-B) P^(-1) (mod 2Q).

For every proper nonempty prefix i, retain its carrier(P_i,Q_i,B_i). Define

    theta=max(1, B/(Q-P), (Q_i-B_i)/P_i for proper prefixes i).

Then the exact positive-integer, nonroot, strict descent guard with no prior ROOT is

    n = r (mod 2Q),    n>theta.

**Proof.** The native congruence makes the final endpoint odd. Backwards
integrality through powers of3 forces every intermediate state integral and
odd; equivalently this follows inductively from the usual source cylinder.
All states are positive for positive n. The proper-prefix cuts exclude an
earlier1, since an odd positive state exceeds1 exactly when its numerator
exceeds its denominator. The final strict inequality is precisely
(P n+B)/Q<n. Conversely every strict first-hit descent along w satisfies
all these conditions. Thus this is an iff within the declared contracting
word, not a completeness assertion about possible words.

Let n0 be the least point of the cylinder above theta. For all t>=0,

    n=n0+2Q t,    h=(P n0+B)/Q+2P t<n.

The compiler retains the cut even when it cannot accept the least residue.
For w=(4,2), the native residue is5 mod128, but5 reaches1 before the second
letter. The first accepted source is133. Cancelling the terminal condition
would turn ROOT padding into a fictitious certificate.

## 2. A checked terminal supplies a whole ray

Suppose a supplied first-hit word w proves U_w(a)=1 for a>1. Then Q=P a+B,
so Q>P and 0<a<2Q. No earlier prefix reaches1. Therefore the full terminal
word yields the exact all-height reduction

    a+2Q t  -->  1+2P t,    t>=0.

Only t=0 has ROOT as its child. The other children are smaller but still
need their own proof. This distinction is the global induction obligation.

Keeping only the first descent prefix widens the guard. If its total is A0
and the full word's total is A, the containing cylinder is wider by2^(A-A0);
strict initial cuts remain in force. This is a ratio of cylinder measures,
not a runtime improvement or a count of newly proved integers.

The experiment consumes the independently compiled ROOT words of the37
smallest requested representatives of the missing components from the
terminal-basis package. It makes **no further orbit-discovery queries** for
these representatives. All37 give guarded descent rays; the saved output
lists them and their full/short totals. The exponent savings range from17
to129. For4591 the full cost109 becomes65; for26623,194 becomes65.
Every representative itself is already fully grounded by the finite basis.
Their arbitrary lifted children are not thereby grounded.

### Explicitly grounding an infinite subfamily inside every terminal ray

There is, however, a parameter choice with an explicit final step. For the
same supplied ROOT word w at a, put P=3^L and choose any integer s>=0. Set

    J=1+P s,
    h_s=(4^J-1)/3,
    n_s=a+4Q(4^(P s)-1)/(3P).

For every s,4^(P s)=1 mod3P. This follows by induction from4=1 mod3 and
the fact that cubing1+3^j u gives1 mod3^(j+1). Consequently
t_s=2(4^(P s)-1)/(3P) is a nonnegative integer, n_s=a+2Q t_s is in the
native cylinder, and U_w(n_s)=1+2P t_s=h_s.

For s>=1, h_s>1 and3h_s+1=4^J. Appending the **single actual valuation2J**
therefore finishes the first-hit ROOT certificate. At s=0, n_s=a and h_s=1,
so the stored word is just w; appending2 would be ROOT padding and is
explicitly forbidden. This is a fully grounded infinite generated family,
not only a conditional smaller-child rule. It is an elementary specialization
of the inherited inverse-word grammar, with no new coverage claim outside
its displayed parameterization.
Closest prior compilers are [finite seed receipts, section2](finite_seed_receipt_compiler_20261005.md)
and [frontier family completion, section6](frontier_family_compiler_20261004.md).
The present specialization starts from a checked ROOT word and gives the
exponent phase J=1+3^L s directly, without a separate phase search.

The integers can be vastly larger than a feasible explicit binary string.
The finite tuple(a,w,s) retains their exact formula and ROOT word; its final
exponent2(1+3^L s) is an ordinary stored integer. To read the source modulo
any positive m without expanding it, keep the denominator precision:

    R=4^(P s) mod(3P m),
    n_s mod m = [a+4Q(R-1)/(3P)] mod m.

R=1 mod3P makes this division exact. Reducing modulo m before dividing
would in general lose the answer. `CompletedLift`, `completed_word`, and
`completed_residue` implement this lossless symbolic certificate. Bounded
expansion independently replays16 small controls; all37 newly grounded
representatives are also lifted with exact modular checks, without expanding
their enormous generated sources. This constructs explicit infinite
completed subfamilies but does not decide whether an arbitrary input belongs
to an unrestricted ROOT basin, nor prove that their union covers all positive
integers. Membership in each displayed family **is** decidable by the supplied
`recognize_completed`: compute h=(Pn+B)/Q, require3h+1 to be a power4^J,
and require J=1 mod P. The actual lifted source then has the unique parameter
s=(J-1)/P. This check recognizes only the stated family.

## 3. A variable-length reset crosses a fixed budget obstruction

Write n=2^K t-1 with K=v2(n+1)>=2 and t odd. Put

    r=v2(3^K t-1).

If r>=2, the known reset collision supplies the smaller child h=(n-1)/2.
The actual words are

    source: 1^(K-1),r+1;
    child:  1^(K-2),2,r-1;
    join:   (3^K t-1)/2^r.

The leading-one formula U^j(n)=3^j 2^(K-j)t-1 proves the source word.
Since3^K t=1 mod4,3^(K-1)t=3 mod4, giving the child's first reset2;
one more odd step gives the same join. The only child-ROOT exception is
n=3, where the child word is empty and the source word is(1,4).
No word is padded at1. This is the D=1 identity of THM-4555, implemented
without an a priori cap on K. Literal generation and replay still cost K
steps and growing integer arithmetic; the macro is not declared free.

### An infinite family genuinely outside the frozen previous selector

Let K=1458j, j>=1, and n=2^K-1. Under the frozen B8 policy from
[the preceding routes package](collatz_uncovered_join_routes_20261007.md),
each of the eight macros consumes exactly128 leading ones. The total word
is1^1024 and the final state is3^1024*2^(K-1024)-1>n.

Here is the complete branch argument, not just a few sampled lifts. At a
checkpoint after d<=896 ones, v2(x+1)=K-d, v2(x+5)=2, and v2(11x+19)=3.
Thus the repeated12 and112 rules fail and no proper S(h)=4h+1 sibling is
present. A core-bank descent word of length<=128 cannot match because all
those actual prefixes expand. A virtual row, even if recognized, needs
K-d letters, exceeding128. The to-ROOT test fails. The fallback is therefore
exactly128 ones. At the final checkpoint K-1024>=434, there is still no
sibling. These facts prove the complete trace at every j.

Also2^1458=1 mod2187. Hence n=0 mod2187 and neither old inverse-family
guard n=91 mod162 or n=111 mod4374 fires. The entire family is PENDING in
the previous refined policy. But K is even, so lifting the exponent gives
r=2+v2(K)>=3. The uncapped reset pays h=2^(K-1)-1<n for every member.

In the small universe through32767, this rule adds22 local dispatcher
reductions beyond the earlier31. The terminal-basis package separately
shows that all of those literal edges were already retained in the merged
observation graph. This is useful global relative coverage, not additional
finite graph completion.

## 4. One reset-2 branch can also be crossed

THM-4555(vi) gives the sporadic collision

    (2,6,c) ~ (4,1,1,c+2),    c>=1.

Consequently a source with word1^(K-1),(2,6,c), K>=2, has a checked smaller
child(n-1)/2 with word1^(K-2),(4,1,1,c+2). Both routes have length K+2.
The implementation detects the actual three letters after the initial run,
retains their guard, and independently replays both complete prefixes.

Consider the odd Mersenne exponents

    K=10207+46656s=1458(7+32s)+1,    s>=0.

The preceding B8 argument still applies. Now n=1 mod2187, so both inverse
guards still fail. Odd K makes the first reset2. The next exponent is
1+v2(K+1)=6 because K=31 mod64. The sporadic rule therefore pays
2^(K-1)-1. That child's even exponent admits the ordinary reset rule,
giving a further smaller dependency2^(K-2)-1. This proves a two-rule
descent in the exponent on the entire displayed family.

The positive source2^1459-1 passes neither of these two collision rules.
This is an explicit remaining guard obligation, not evidence of divergence.
Searching further admissible collision identities or a variable-depth
cover of this residual branch is the next precise target.

## 5. Why finite terminal lifting alone cannot be the global answer

Any finite collection of fixed contracting valuation words has a maximum
length L. For M>L, the source2^M-1 starts with L ones, and every nonempty
prefix among them expands:3^j/2^j>1 with positive carry. It therefore
belongs to none of these fixed descent guards. Finite exceptions do not
repair this, since M can be chosen beyond all of them.

This obstruction applies to fixed forward words, not to all conceivable
Collatz certificates. The variable-length reset above is precisely an
operation outside that restricted class. It repairs an infinite subfamily
and leaves other branches explicit. Finite grounding and a universal
guard-cover/termination theorem remain distinct obligations.

## Reproduction

Run `python -B 04-computation/experiments/collatz_terminal_lifts_20261007.py`
and the same command with `-O`. The
[script](../../04-computation/experiments/collatz_terminal_lifts_20261007.py)
and [output](collatz_terminal_lifts_20261007.out) retain317 contracting words
over1..4 of length1..4, all odd sources3..1023 for their iff controls,
large-parameter lifts, all finite reset guards through32767, the37 supplied
representative certificates, and forged/type/ROOT-padding controls.
All-height claims above are algebraic proofs; huge finite orbit searches
are not substituted for their quantifiers.
