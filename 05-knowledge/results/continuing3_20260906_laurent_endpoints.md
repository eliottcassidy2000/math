# Actual carried-return families at endpoints 39,45,51,57,63

**Status: PROVED ANALYTICALLY + COMPLETE EXACT CERTIFICATE + INDEPENDENTLY AUDITED.**
The [independent audit](continuing3_20260906_laurent_endpoints_audit.md)
accepts the full proof with 14789 exact gates. The real-root input retains
the proved/CITED scope of THM-4436. General all-h transport remains OPEN.

Concurrent credit: incoming commit 28f41c846 already proves
[endpoint39](long_frontier_sep06_endpoint39.md). Our h=6 slice is an
independent reconstruction: all 258 coefficients and rational contents
agree exactly. The extension beyond that incoming result is endpoints
**45,51,57,63**, with 34 additional characteristic polynomials and 2912
positive shifted coefficients. The five-family statement below keeps the
shared boundary case for a common reproducible certificate.

## 1. Statement

For any h in {6,7,8,9,10}, integer g>=3h+2 with gcd(g,6h+3)=1, and arbitrary
nonzero complex coefficients alpha,beta,gamma, put N=6h+3 and

    f(u)=alpha*u^(-N)+beta*u^(2g-N)+gamma*u^(3g-N).

The first nonzero constant-term moment of f occurs at exactly g or 2g.
Both possibilities occur for every allowed pair h,g. There are exactly h
distinct coefficient phases canceling the mass-g moment; every such phase
has a nonzero mass-2g moment. In particular first detection is at most
2g<3g=N+(3g-N).

| h | Fixed negative endpoint N | Allowed g | First/doubled channel counts |
|---|---:|---|---|
| 6 | 39 | g>=20, gcd(g,39)=1 | 7 / 14 |
| 7 | 45 | g>=23, gcd(g,45)=1 | 8 / 16 |
| 8 | 51 | g>=26, gcd(g,51)=1 | 9 / 18 |
| 9 | 57 | g>=29, gcd(g,57)=1 | 10 / 20 |
| 10 | 63 | g>=32, gcd(g,63)=1 | 11 / 22 |

These are five infinite families with g unbounded, not a claim about every
trinomial whose negative endpoint is at most63. The two positive exponents
must have exactly the displayed relation to g. All signs below belong to
the specified normalized response, not a general complex raw moment.

## 2. Inheritance and the changed operation

The closest proved mechanism is the endpoint33 characteristic certificate
in `05-knowledge/results/second_20260906_laurent.md`. Its explicitly stated
next endpoint was39. The normalization and degree bound for every h were
already proved in Section5 of
`05-knowledge/results/overnight7_20260906_laurent_midpoint_transport.md`.
They are inherited inputs, not new formulas claimed here.

The root supplier is **THM-4436**,
`01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md`.
Its precise scope is a proof relative to a cited finite real-root preserver.
The canonical hostile is the free-scale carried-response model, which can
retain real-rootedness while losing the original crossing coefficient.
The corrected near miss is a positive sign scan at finitely many g. The
least-used sidecar is the inherited all-h characteristic degree bound,
which makes independent finite polynomial identity verification complete.

The board is: original charge fibres; inverse carry; the first-root quotient;
the characteristic polynomial; coefficientwise positivity after shifting
x=g-3h-1 to x-1>=0. Passing to the quotient forgets values off the first
cancellation locus. The original coefficient monomial, exact fibre maps,
and carry identity restore the interpretation of the quotient response.

The new operation is to construct the whole symbolic certificate for the
five named dimensions, not to infer signs from selected roots or g values.
All-h positivity remains open. Continuing to scan larger h alone would not
prove its induction step; this cycle stops at a coherent five-family proof.
Targeted searches of current canon/results for endpoint39, width39, the
literal supports and seven-first-channel closure found only the incoming
open question. No external literature-priority claim is made.

## 3. Literal first and doubled rows

Write x=g-3h-1>=1, tau=alpha*gamma^2/beta^3 and
X=alpha^x beta^(3h) gamma. The complete return fibres at masses g and 2g are

    (n_alpha,n_beta,n_gamma)=(x+j,3h-3j,1+2j), j=0,...,h,
    (n_alpha,n_beta,n_gamma)=(2x+e,6h-3e,2+2e), e=-1,...,2h.

Indeed the charge equation at mass g is 2n_beta+3n_gamma=6h+3;
nonnegativity and parity give exactly the first list. Doubling the equation
gives the second list. Its lower carry e=-1 is valid because 2x-1>=1.
For every positive feasible mass m, reduction modulo g gives g|Nm;
coprimality therefore gives g|m. The j=0 first channel is feasible, so the
first support return is exactly g. This step requires the gcd assumption.

Define the monic first polynomial and normalized doubled Laurent response:

    p_x(t)=sum_(j=0)^h (2h+1)! (x+h)_(h-j)
                          /[(3h-3j)!(1+2j)!] *t^j,
    q_x(t)=sum_(e=-1)^(2h) (2x+2h)_(2h-e)
                          /[(6h-3e)!(2+2e)!] *t^e.

The falling factorial (a)_r includes r factors and (a)_0=1. Literal
multinomial expansion gives exactly

    CT(f^g)=X binom(g,2h+1) p_x(tau),
    CT(f^(2g))=X^2 (2g)_(4h+2) q_x(tau).

All omitted scalars are displayed and nonzero in this range. THM-4436 with
A=2,B=3,r=0,z=1 gives h distinct strictly negative roots of p_x for every
integer x>=1. Its constant coefficient is positive; t is invertible at
every such root. Choosing beta=gamma=1 and alpha equal to one of the h
roots attains that cancellation phase. Coefficients all equal to one give
the noncanceling first-detection alternative g.

Without gcd(g,N)=1 the support-return claim is false. Uniformly for these
h, take g=3h+3, so gcd(g,N)=3. The count triple (1,h-1,1) gives first mass
h+1=g/3; reduction modulo g/3 proves no earlier mass. This is a hostile
to dropping the gcd, not a failure of the polynomial certificate.

## 4. Exact polynomial quotient and structural degree bound

Begin over Q(x)[t]/(p_x), where t is invertible. The possible denominator
from its inverse cancels in the polynomial identity

    q_(-1)(x)/p_0(x)
      = [2^(h+1)(3h)!/((6h+3)!(2h+1)!)]
           *x*product_(j=0)^(h-1)(2x+2j+1).

The even factors in (2x+2h)_(2h+1) are 2x,2(x+1),...,2(x+h);
the h factors x+1,...,x+h cancel p_0(x). The remaining factors are the
displayed odd factors. This proves the identity for every h under discussion.

Reduce the nonnegative part of q_x by monic division and substitute

    t^-1 = -sum_(j=1)^h p_j(x)t^(j-1)/p_0(x).

The resulting R_x(t)=sum_(j=0)^(h-1) R_j(x)t^j belongs to Q[x,t] and
has deg_x R_j<=2h-j. Giving both x and t weight one proves this for
monic reduction; the inverse contribution has degree h+1 in its canceled
factor and h-j-1 in p_(j+1), giving the same total degree2h-j.
The producer also checks the full polynomial identity t*q_x=t*R_x mod p_x,
including the lower carry. Invertibility is not asserted in Q[x][t]/(p_x)
before this cancellation step.

Let M_x multiply by R_x in the basis 1,t,...,t^(h-1). Entry(i,j) has
degree at most2h+j-i. A k-by-k principal determinant has total degree
at most2hk because its row and column index sums cancel. Thus, writing

    C_x(z)=det(zI-M_x)=z^h+c_1(x)z^(h-1)+...+c_h(x),

we have deg_x c_k<=2hk. The program uses rational polynomial
Faddeev-LeVerrier recurrences and verifies their terminal Cayley-Hamilton
matrix identity exactly. The finite degree bound is established before
checking values or inspecting the coefficient signs.

## 5. Complete positivity object and its consequence

For every h in {6,...,10} and k=1,...,h, the accompanying certificate gives

    c_k(x)=d_(h,k)*sum_(j=0)^(2hk) a_(h,k,j)*(x-1)^j,
    d_(h,k)>0, every a_(h,k,j) a strictly positive integer.

The full content and all integer coefficients are in
`continuing3_20260906_laurent_endpoints_certificate.json`; no coefficient
is replaced by a digest or a sign summary. Counts by h are258,399,584,819,
and1110, totaling3170. They prove c_k(x)>0 for every real x>=1.

Consequently C_x(z)>0 for every real z>=0. At any real root rho of p_x,
evaluation of the quotient multiplication operator gives
C_x(q_x(rho))=0. The real number q_x(rho) therefore cannot be nonnegative;
it is strictly negative. This implication holds for every real x>=1 and
every real first root at that parameter; no real-root theorem at nonintegral
x is assumed. For integral x, THM-4436 supplies every first root on the
negative real axis, proving the statement in Section1 through the literal
coefficient maps. The original complex scalar X^2 is nonzero, so the
normalized negative response suffices for noncancellation.

## 6. Verification and remaining boundary

The declared finite proof universe is all five symbolic quotient operators,
their full40 characteristic coefficient polynomials, and all3170 shifted
coefficients. The separate literal-control bank uses x=1,2,3,7,16 at each h:
25 indexed rows, of which18 satisfy the first-return gcd. Every nonnegative
charge triple at masses g and2g is enumerated, reconstructing675 complete
return fibres. The lower carry is retained. The uniform earlier-return
gcd hostile is checked for each h. No inherited physical filters apply to
this coefficient problem.

The producer imports SymPy and the standard library, no repository
mathematical implementation or numerical root solver. It passes4451
always-active exact gates. The symbolic algebra and complete positive
certificate, not a negative finite search, provide the universal quantifiers.
The independent literal-fibre/integer-Berkowitz reconstruction verifies
665 parameter matrices, enough by the proved degree bound to establish
every characteristic identity. Both producer and referee normal/optimized
outputs match their frozen outputs exactly.

    python -B 04-computation/continuing3_20260906_laurent_endpoints.py --certificate 04-computation/continuing3_20260906_laurent_endpoints_certificate.json
    python -B -O 04-computation/continuing3_20260906_laurent_endpoints.py

The complete JSON certificate has923826 LF bytes and SHA256
`b55bec0ba7f1063396af1a8db9be725f279fd01a403adb94b9b026b9490d8e62`.
Source SHA256: `3b5e9127cd02881afbd2a3cfe8041c990cd2a8688a4eeac86a4cb389c110eae8`.
Output SHA256: `f5ea177b5b3d3d49bd29dd677a35b8244a17ca5c38eb90e44867f5684b56c6e4`.
The manifest records all filed artifacts and the pre-promotion proof identity.

The next genuine question is a positivity-preserving recurrence or another
uniform same-root argument in h. Five endpoint certificates do
not supply such an operation. The two-anchor model's finite-phase question
is also distinct: a literal binomial family does not exhaust that larger
semialgebraic coefficient class. General two-rung noncancellation remains OPEN.
