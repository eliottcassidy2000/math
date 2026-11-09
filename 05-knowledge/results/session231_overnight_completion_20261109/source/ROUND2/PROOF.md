# Two homogeneous clocks: exact bank rigidity and inherited fan coordinates

Status: PROVED candidate, awaiting independent audit. Chan–Pak and F2 completeness remain OPEN. No extension of interval coverage and no literature-priority claim.

## 1. Recovery before claims

The ray experiment searched 61,504 source templates against the physical order-at-most-seven bank. It found 73 parity/offset rays with offsets in [-255,255]. Almost all were elementary single-block rays. This discovery count is FINITE-EXACT, not the theorem below.

The relevant older results are stronger than a superficial search suggests:

- Session191 `session191_tournament_growth_20261104/one_apex/PROOF.md`, Sections6.1–6.2, classifies every ONE expanding homogeneous run in a strong one-apex word. A dyadic slope forces offset+1 there.
- `square_sum_spectrum_overnight_20261031e/CENTERED-FANS.md`, Sections1–6, ALREADY proves midpoint-palindrome cancellation, infinite fixed-deficit rays, strict energy ordering, an O(sqrt X) deficit bound, and the permanent deficit9 obstruction. It includes actual strong graph constructions. None of those is claimed new here.
- Session220 `selectors/PROOF.md` rules out rational finite-state selectors for the complete fixed F2 bank. Our unbounded clock is a specified arithmetic model, not a purported evasion proving universal coverage.
- The suffix law, f reversal symmetry, and square-root inverse are inherited from the prior sessions cited in the first packet.

A full tracked-Markdown search of canon and result notes for palindrome/binary/apex, expanding-run, synchronized-run and square-completion terms recovered the two principal predecessors above. The present candidate addition is the arbitrary-bank TWO-clock rigidity theorem. The half-word formulas below are an exact coordinate refinement of the inherited centered-fan family, not a new family.

## 2. The two-clock theorem, arbitrary fixed context and physical row

Let a be positive odd, b>=0, c>=1. Fix arbitrary binary words U,W,S and bits beta,gamma. Empty fixed words and leading zeros are allowed. Let v_s be the integer represented by

    U beta^s W gamma^s S,      s>=1.

Assume that for some integers j>=0 and D,

    a g(v_s)-b e(v_s)+c(v_s+1)=2^(2s+j)+D              (T)

for EVERY s>=1.

**Theorem.** Necessarily beta=gamma. If beta=gamma=1, either W consists entirely of ones (including empty), or a=1. If beta=gamma=0, either W consists entirely of zeros (including empty), or a=1.

Thus for EVERY physical rooted bank row with a>1, a successful pure-dyadic two-clock family in this model collapses to a single homogeneous expanding run. Adding any number of such rows cannot yield a genuinely separated two-clock pure-dyadic ray. This is not a theorem about arbitrary target-dependent witnesses, different clock rates, three changing blocks, or an unrestricted F2 row.

The boundary a=1 is sharp: singleton row (1,1,1) and v_s=1^s 00 1^s give 2^(2s+2)-1 for all s>=1. These are actual Hamilton-path counts, not formal states.

### Proof

Set z=2^s. Absorb the fixed suffix S into the readout. The exact suffix law shows that its effect has the form

    A g(v)+K v+D_0,

where A=a, K=c 2^|S|+a(2^|S|-1-rev(S))>0, and D_0 is constant because the parity before S is the fixed final bit gamma. The empty suffix has K=c. In particular A stays odd; no physicality claim about an artificial new row is needed.

Let P=value(U), M=2^|W|, r=value(W), R=rev_|W|(r), and C=M-1-R. Fixed empty W has M=1,r=R=C=0. Streaming the two runs gives the following exact coefficients of z^2 and z in the final readout:

    beta gamma   quadratic coefficient       linear coefficient
        1 1       K M(P+1)                   A C(P+1)-K(M-1-r)
        0 0       (A+K) M P                 (A+K)r-A P R
        1 0       (A+K) M(P+1)              -A(R+1)(P+1)-(A+K)(M-r)
        0 1       K M P                     A P(1+C)+K(r+1).

All omitted terms are independent of z. These formulas follow directly from

    append 1^s: v -> z(v+1)-1,   g -> g;
    append 0^s: v -> zv,         g -> g+(z-1)v+(v mod2).

For (T), the quadratic coefficient is 2^j and the linear coefficient is zero, because polynomial equality holds at infinitely many z. In mixed cases the displayed linear coefficient is strictly negative or strictly positive. Hence beta=gamma.

We use the elementary reversal lemma: if n>0 and A is positive odd, and A rev_k(n)/n is an integral or fractional power of two, then A=1. Strip trailing and leading zeros: n=2^t u, u odd with l bits, and rev_k(n)=2^(k-t-l) rev_l(u). Comparing odd parts forces A rev_l(u)=u. Both odd integers u and rev_l(u) have exactly l bits, so their ratio is strictly less than2. Thus A<2 and A=1.

In case11 put q=M-1-r. If q=0, W is all ones and the two runs merge. Otherwise C=rev(q)>0. The quadratic coefficient being a power of two makes K and P+1 powers of two. Linear cancellation gives

    A rev(q)(P+1)=K q,

so the reversal lemma forces A=1.

In case00, the nonzero quadratic coefficient forces P>0, and both A+K and P are powers of two. If r=0, W is all zeros and the runs merge. Otherwise cancellation gives

    A P rev(r)=(A+K)r,

and again the lemma forces A=1. This proves the theorem for arbitrarily long fixed U,W,S and every physical bank size. The endpoint coordinate b affects only the constant term and cannot repair the obstruction. QED.

The proof also covers fixed shifts in either run length by absorbing the extra copies into U,W,S. It does NOT cover distinct growing rates such as s and2s.

## 3. Exact half-word coordinates for the inherited centered fans

Let f(0)=0, f(2u)=f(u), f(2u+1)=f(u)+u+(u mod2). Thus g(v)+f(v)=v+e(v). Let w have positive length k, integer r, reversal R, and first bit t. For s>=1 set v=value(1^s w 1^s). The singleton readout is exactly

    F_(1,1,1)(v)=2^k 4^s+(r-R)2^s+d(w),
    d(w)=g(r)+R-2^k+2-t.                              (6)

For the empty middle word k=0, the separate direct formula is F=4^s+1, so d=1. Thus fixed pure-dyadic rays occur exactly when w is palindromic. This cancellation and the resulting centered fans are INHERITED; (6) is their binary-source form.

Here are closed coordinates for the constant, useful for exact target tests. Put B=2^l, 0<=q<B, x=B-1-q, p=q mod2. Interpret x and its reversal as fixed-width l-bit words (empty if l=0).

- If w=x rev_l(x), of even length2l, then

      d=1-[q^2+2f(q)+p].                              (7)

- If w=x beta rev_l(x), of odd length2l+1, then

      d=1-2[q^2+f(q)+(1-beta)(q+p)].                   (8)

Consequently every q>=0 gives actual count families

    2^m+1-q^2-2f(q)-p,
      m even, m>=2 bitlength(q)+2;

    2^m+1-2[q^2+f(q)+(1-beta)(q+p)],
      m odd, m>=2 bitlength(q)+3, beta in {0,1}.        (9)

The statements allow q=0 with bitlength(0)=0. They are sufficient ray criteria, not coverage of all offsets.

### Algebraic proof of (7)–(8)

The pair-distance formula is

    f(n)=sum_{i>j: bit_i(n)=bit_j(n)=1}
                    (2^(i-j-1)+[i=j+1]).

It immediately gives reversal invariance and the concatenation law

    f(2^k u+r)=f(u)+f(r)+u rev_k(r)+(u mod2) first_k(r).

For an l-bit complement x=B-1-q, sum the weights incident to each selected q-bit to obtain

    f(x)=B-2-q-rev_l(q)+(q mod2)+first_l(q)+f(q).       (10)

This formula is used for l>=1; l=0 is checked directly. Indeed the incident weight at position i equals 2^i+2^(l-1-i)-[i=0]-[i=l-1], and f(B-1)=B-2. Substituting (10), g=r+e-f, and

    f(x rev(x))=2f(x)+x^2+(x mod2),
    f(x beta rev(x))=2f(x)+2x^2+2beta[x+(x mod2)]

into (6) yields (7) and (8).

### Exact constant-time candidate selection, retaining the full clock

For an even exponent m and proposed offset d<=1, put T=1-d. Since 0<=f(q)<=q-1 for q>=1,

    q^2 <= q^2+2f(q)+p < (q+1)^2.

Thus the ONLY possible q is floor(sqrt(T)); test exact equality in(7), then the exponent bound in(9).

For odd m put T=(1-d)/2, an integer for odd d. Both candidate expressions in(8) lie in [q^2,(q+1)^2), so again q=floor(sqrt(T)) is the only possibility. Test the two beta choices and the exponent bound. This gives an exact test using one integer square root and O(log(2+|d|)) arithmetic operations for the recurrence for this mirrored-word family, rather than the older fixed-diameter binary search. It does not erase the exponent or diameter sidecar, and does not apply to all one-apex witnesses.

The inherited missing deficit9 is transparent: d=-9. For even exponents q=floor(sqrt10)=3 gives T=14, not10. For odd exponents q=floor(sqrt5)=2 gives T=4 or6, not5. Hence no mirrored-word member has offset-9, regardless of clock length. This is a REPROOF of CENTERED-FANS Section5, not a new global forbidden ray.

## 4. Physical graph witnesses and scope

For every accepted certificate in(9), construct the literal chain-incidence word

    BITS = 1^s w 1^s 0,

with s=(m-|w|)/2>=1. Take vertices apex x and chain vertices1,...,|BITS|. Direct each chain edge i->j for i<j; direct x->i iff bit_i=1. Its Hamilton count is g(value(BITS))=F_(1,1,1)(v), by the inherited recurrence. The first bit is1 and final bit0, so the graph is strong. This is a full vertex-level witness, with order m+2; no independently unverified rooted bank entry is needed for these positive rays.

The bank obstruction concerns exact scalar readouts of authentic one-apex sources. It introduces no fictitious physical state. Our discovery program's other positive-offset rays can use bank entries whose physical realizability is inherited from THM-4875, but those isolated elementary rays are not promoted as a new coverage advance.

## 5. What changed, and what is still missing

The experiment did not fill any new target cylinder. It led to a principled stopping theorem: a large physical bank cannot repair the two-clock pure-dyadic cancellation failure when its Hamilton coefficient exceeds1. This strictly goes beyond Session191's one-run restriction, while remaining a narrow ansatz obstruction. The palindrome construction itself is old; the new half-word coordinates expose its exact square-root candidate and preserve the clock.

An all-target approach must still prove exact F2 occupancy. For this particular ray-based plan, a genuine next step must use at least three changing blocks, unequal growing clock rates, or target-dependent middle words/contexts not fixed along the ray. None is proved sufficient here. The first packet's inverse still supplies only complete decision, not occupancy.

### Scalar-only exact membership corollary

The half-word coordinates actually remove a search over exponents for this particular family. In(9), the exponent bound gives q+1<=2^((m-2)/2) for even m, or <=2^((m-3)/2) for odd m. The square brackets above are strictly below (q+1)^2. Hence every such count satisfies

    (3/4)2^m+1 < N <= 2^m+1.

For integer N>=5 there is at most one possible m, namely ceil(log2(N-1)), computed exactly as bitlength(N-2). Determine that m, form d=N-2^m, make the sole square-root candidate test above, and check the exponent bound. This is an IFF membership test for the mirrored-word family with s>=1, with a literal strong graph constructor. It is not an IFF test for the entire one-apex spectrum. `decode.py` implements it with arbitrary-precision integers.

This does not contradict the inherited statement that diameter matters when decoding arbitrary symmetric subset energies: our output is a chosen canonical minimal-width half-word together with the uniquely compatible exponent, not the original arbitrary presentation or future-context equivalence.

### Transitive-row support control

The seven supplied bank rows with a=1 are (1,1,1) and (1,0,2^j), j=1,...,6. They add no scalar support to the singleton row. Indeed

    g(2^(j+1)v+2^(j+1)-2)=g(v)+2^j(v+1),

and singleton readout at u is g(2u). Conversely every g-value can be written with even argument by deleting trailing ones, except g(0)=1 which is already the singleton v=0 value. Thus the a=1 boundary of the bank-rigidity theorem is exactly the old one-apex support, not an unexplored enlargement hidden in the transitive bank rows.
