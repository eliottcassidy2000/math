# State transport, rather than a coincidence of small numbers

## Status and purpose
PROVED algebraic identities below. ANALOGY to tournament interfaces is explicit; no modular-form theorem proves Chan–Pak here. The purpose is a reusable diagnostic: a change of presentation can preserve an observable while changing the state needed to continue it. This answers why scalar agreement, including matching a target, is insufficient for recursive construction.

## A marked substitution calculus
For any sequence a_j and scalar q define

b_n = sum_{j=0}^n binom(n+j,n-j) q^(n-j) a_j.

Writing A(x)=sum a_j x^j and B(y)=sum b_n y^n, formal summation gives

B(y)=(1-qy)^(-1) A(y/(1-qy)^2).                         (1)

Indeed set n=j+r and sum binom(2j+r,r)(qy)^r=(1-qy)^(-2j-1). This is a formal identity, requiring no convergence. On half-densities F_A(x)=sqrt(x)A(x), F_B(y)=sqrt(y)B(y), with x=y/(1-qy)^2 and coherent square roots, F_A(x(y))=F_B(y). Its marked operator transports as

Theta_x F_A = (1-qy)/(1+qy) Theta_y F_B, Theta_z=z d/dz. (2)

Thus the correct interface is (F,Theta F), not F alone. The coordinate rule is exact wherever denominators are nonzero, and formal near y=0. A higher-order marked observable needs its full jet, obtained by repeatedly applying the transported differential operator. This is analogous to, not an identification with, retaining rooted endpoint and cover coordinates under tournament substitution.

## A finite positive state implementation and its spectral escape
For a_j=b^j, (1) becomes

B(y)=(1-qy)/(1-(2q+b)y+q^2 y^2).

For nonnegative integers q,b, this has a literal nonnegative weighted-walk realization

M=[[q+b,qb],[1,q]], b_n=(M^n)_{11}.

The characteristic polynomial is X^2-(2q+b)X+q^2 and discriminant b(b+4q). For positive q,b and nonsquare discriminant both roots occur with nonzero coefficients: neither root equals q because the characteristic polynomial evaluated at q is -bq. Consequently this sequence cannot be represented as a finite sum of polynomial-in-n multiples of integer powers r^n. To see this, its reduced generating function has poles at the reciprocals of two irrational roots, whereas any such integer-mode sum has poles only at reciprocals of integers. This makes a precise test for when a clock change escapes an integer-mode construction class.

These are weighted walks, not tournaments. For q=b=1 the sequence begins1,2,5,13,34, so its even terms already forbid interpreting the entire sequence as tournament Hamiltonian-path counts. A physical tournament compiler is an additional theorem, not supplied by the matrix. Conversely, any restricted tournament construction known to have only integer exponential modes cannot realize this full sequence. This claim is conditional on that class's integer-mode theorem, not a statement about all tournaments.

## The user's modular formula: exactly what12,24,127 do
Use the standard convention U=(127+5 sqrt(645))/2. Then UU'=1 and U+U'=127. Therefore L_m=U^m+U'^m satisfies

L_0=2,L_1=127,L_{m+2}=127L_{m+1}-L_m,
L_3=127^3-3*127=2048002.

Set t=-432 U^3. Then

j=t+864+432^2/t=-432(L_3-2)=-960^3.

The cube is essential under this definition of U; the pasted formula without it uses a different convention or drops that exponent. This is an exact degree-two norm identity, not a classification of tournament cores.

Let a_k=binom(2k,k)binom(3k,k)binom(6k,3k), and define s_1B by (1) with q=-432. Then x=1/j,y=1/t satisfy x=y/(1+432y)^2. For these algebraic values

f=(1+432y)/(1-432y)=320 sqrt(645)/8127.

Transport of the linear observer5418k+263 gives

36*[5418 f (k+1/2)+(263-5418/2)]
 =24*[-3669+320 sqrt(645)(k+1/2)].

Thus the displayed24 is recovered from an exact coordinate-and-observer change of the classical level1, N43 series, not from asserting that24 or127 govern tournament realizability. The choices of i and square roots must be coherent; take sqrt(j)=i sqrt(-j),sqrt(t)=i sqrt(-t). The resulting transformed formula has denominator(-432 U^3)^(k+1/2).

Primary source for the input identity and companion framework: H.H.Chan and S.Cooper, Rational analogues of Ramanujan's series for1/pi, Mathematical Proceedings of the Cambridge Philosophical Society153(2012),361–383, DOI10.1017/S0305004112000254, Theorem2.1 and Tables1–3; accessible copy https://mrc.sdu.edu.cn/ziliao/67.pdf. The formal manipulations above are independently derived. No claim of a new Ramanujan–Sato identity is made.

## The meaningful bridge to the70 birth families
A head repair has two obligations: its scalar output and its reusable interface. Formula(2) demonstrates exactly how a scalar-preserving coordinate change can require a nontrivial interface multiplier. For a tournament operation, the analogous desired certificate must specify adjacency construction plus the entire named-port/cover transition, not just match a recurrence. A promising candidate is a core-changing compiler with a transported, rather than frozen, endpoint state; a rank must decrease on the source target while every guard survives. Neither irrational spectral modes nor a norm-one unit supplies this compiler automatically.

## Repo-history map and limiting analogies
- monster_levels_square_sum_7_21_20261022.md separates Monster/Baby Monster prime differences from unsupported exponent-cycle interpretation. Its sets{29,41,59,71} and{17,19,23,31,47} are modular-level data; the cycle{4,2,5,3,7,13} is not a tournament transition law.
- borwein_ladder_owner_numbers_vsc_twelve_20261028.md identifies{4,2,5,3,7,13}={d+1:d divides12}, with the appropriate prime subset. This suggests divisor-labelled transition types, but the labels are not yet physical operations.
- session211/square-sum-port-bridge/PROOF.md supplies a stronger structural analogy: adjoining a degree-two vertex to a pathless graph turns two-component spanning path covers with prescribed endpoints into Hamiltonian paths. The Q23/Q25 counts3/10 cannot be sent by H=2p+1 into tournaments, because they would give7/21. The useful object is the endpoint-cover interface, not the scalar count.
- session201/twelve-factorization/PROOF.md gives finite odd cores with dyadic rays and an all-height gap for the restricted twelve-form grammar. That is why a new clock/state compiler should be tested for escaping the old class rather than merely repackaging it.

The synthesis is a research direction, not a proof of unrestricted odd coverage: transport marked interfaces, then prove a physical lifting theorem, then prove a well-founded inverse has an admissible branch for every target outside7,21.
