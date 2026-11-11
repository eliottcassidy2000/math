# Two growing flanks: exact product families and a boundary

**PROVED scoped identities and obstruction; FINITE-EXACT balanced probe.**
The inherited positive core compiler is the only tournament-count input.
This does not prove a common source return, and does not rule out a changing
address or a different growing source. The related earlier obstruction is
[the finite-selector note](../../chanpak_five_core_bridges_20261110/cover/SELECTOR.md),
which treats a different quantifier: a rational selector for a positive-density
set of targets. A single dyadic ray is outside that density conclusion.

## 1. A positive family with two independent growing blocks

For a fixed core C, write d=H(C), L=d+K(1). With the uniform cut O empty,
the positive observer is(h,d,L,i,o,s,e)=(d,d,L,L-d,0,d,0).
For any integers alpha,beta>=1 and independently a,b>=1, put

    A=alpha*2^a,     B=beta*2^b-1.

The binary recurrence gives

    g(A)=A+g(alpha)-alpha+(alpha mod2),
    g(B)=g(beta-1).

Substitution into the exact chart cancels both separate growing terms:

    Phi(A,B)=L alpha beta*2^(a+b)+C,
    C=d[g(alpha)-alpha+(alpha mod2)+g(beta-1)-1].          (1)

This supplies an actual tournament for every a,b. At fixed a+b it supplies
many binary addresses with the same count; it does not claim those labelled
constructions are distinct isomorphism classes. The full-cut case is the
transposed formula. In the five-core bank, (d,L) is respectively
(1,2),(3,12),(5,32),(9,84),(9,84).

The quantity g(alpha)-alpha+(alpha mod2) is unchanged on doubling alpha.
Likewise g(2beta-1)=g(beta-1). Thus the constant in(1) depends only on
the odd parts of alpha,beta. This is a genuine arithmetic normalization,
not a permutation of the five core types.

## 2. Complete classification of independent fixed-offset cancellation

Fix a core and cut, positive alpha,beta and integer offsets t,u. Consider

    A_a=alpha*2^a+t,    B_b=beta*2^b+u

for all sufficiently large, INDEPENDENT a,b. For every fixed offset t,
the binary recursion gives eventual constants gamma,kappa with

    g(alpha*2^a+t)=gamma*2^a+kappa,    0<=gamma<=alpha.     (2)

One direct proof repeatedly halves t with floor division. A fixed t>=0
eventually becomes0, and a negative t eventually becomes-1. The recurrence
at each earlier step adds an affine function of2^a with fixed parity;
the terminal cases are the two identities in section1. The bounds on gamma
follow from1<=g(v)<=v+1 by letting a grow. This proves(2), rather than fitting
its coefficients to data.

Let lambda_A=gamma_A/alpha, lambda_B=gamma_B/beta. To have an identity

    Phi(A_a,B_b)=L alpha beta*2^(a+b)+constant             (3)

for all sufficiently large independent a,b, the coefficients of2^a and2^b
must vanish separately. Therefore

    Lu+i+d lambda_A=0,    Lt+o+d lambda_B=0.             (4)

Here0<=i,o<=K(1)=L-d, and0<=lambda_A,lambda_B<=1. Since L>d,
the only possibilities for the first equation are

    i=0,u=0,lambda_A=0; or i=L-d,u=-1,lambda_A=1.

The second equation has the corresponding two alternatives. Positivity
of the singleton-path weights means i=0 exactly when I is empty and
i=K(1) exactly when I is the whole core. Thus(3) forces a uniform cut.
For O empty it forces t=0,u=-1; for I empty it forces t=-1,u=0.
These are precisely section1 and its transpose. This completely classifies
this independent product-identity architecture for every fixed nonempty core.

The independence is essential. If a and b are tied, the separate linear
terms can have the same exponential rate and cancel each other. That larger
class is probed in section4; it is not excluded by this proof.

## 3. Why this family does not supply the +4 target ray

For the desired target102*2^k+19, the constant in(1) must equal19 under any
fixed reindexing of a+b. In the five-core bank, all cores except K1 have
d=3,5, or9, and their constants are divisible by d. Thus only K1 can work.
Matching the growing coefficient then forces the odd parts of alpha,beta
to multiply to51. The four ordered possibilities give

| alpha odd part | beta odd part | Constant C |
|---:|---:|---:|
|1|51|37|
|3|17|15|
|17|3|-5|
|51|1|-39|

None is19. Therefore no five-core identity in the complete independent
fixed-offset class(3) gives this ray. The surviving positive families are
explicit neighboring dyadic rays; the failure identifies the missing
freedom as correlated growth, changing offsets, or a new source operation.

## 4. The correlated-growth hostile probe

The script next tests A=alpha*2^m+t and B=beta*2^m+u with one common m.
For each of the90 raw charts and j=0,...,9 it takes every ordered positive
factorization L alpha beta=102*2^j and all offsets t,u in[-4096,4096].
The exact fixed-offset lemma makes g(A),g(B) affine in2^m for every m>=14.
Consequently the chart is

    102*2^j*4^m + D*2^m + C.

The search groups exact coefficients to solve D=0 without a quadratic
offset-pair scan, then tests C=19. A third evaluation independently checks
each computed affine profile. It finds no such identity in this finite box.
The respective numbers of zero-D candidates for j0,...,9 are
8,20,32,44,64,84,104,124,144,164; all fail the constant test.

This is a complete finite probe of the displayed parameter box, not an
all-offset theorem and not evidence of absent tournament witnesses. In fact
the independent inverse already supplies all25 tested target witnesses in
`ternary_plus_four.json`. What fails here is their proposed simple common
address formula. A larger unstructured scan is not justified by this negative
signal; a useful next attempt should add a changing carry/offset coordinate
or use the reset lane's literal growing source rather than repeat this box.

`python -B biflank.py` checks4500 positive identities by the inherited signed
compiler and records the exact hostile probe in `biflank.json`. The map's
source is a core, cut and two addresses; its target is the Hamilton count;
it preserves exact realization and loses the marked return state if only
the scalar is retained. The necessary sidecar is still the literal source.
