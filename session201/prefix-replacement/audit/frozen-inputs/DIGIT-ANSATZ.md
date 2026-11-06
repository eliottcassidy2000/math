# Why separate binary index changes do not give a coefficient-wise doubling rule

**PROVED, pending independent audit; separate stopping-rule sidecar.** This
records the precise limit of one attempted constructive induction. It does
not obstruct mixed-index changes, target-dependent replacements, or functional
identities holding only on special subsets of the index pairs.

Use the raw source operators O/Z/J and

    Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]].

Let g(1)=1, g(2h)=g(h), and

    g(2h-1)=g(h)+h-(h mod2), h>=2.

The unit profile is tau(u)=(g(u),u mod2,u,1), u>=2. For a chronological
middle word W, define

    Psi_W(u,v)=tau(u)^T Q M_W tau(v).

The considered middle words are empty, or begin and end with J. Write
A=g(u), B=g(v), p=u mod2 and q=v mod2. After fixing the relevant residue
classes, treat A,u,B,v as independent formal variables when testing a
coefficient-wise polynomial identity.

The binary index operations and their sibling replacement are

    O: u -> 2u,     g -> g;
    Z: u -> 2u-1,   g -> g+u-(u mod2);
    S: if u is even, u -> u-1, g -> g+u/2-(u/2 mod2);
       if u is odd,  u -> u+1, g -> g-(u+1)/2+((u+1)/2 mod2).

The sibling formulas follow directly from the recurrence. Restrict to
indices for which the intended physical unit profiles have index at least2.
Every finite separate composition, on a fixed input residue class, has form

    u'=a*u+b, A'=A+c*u+d,
    v'=a'*v+b', B'=B+c'*v+d',                    (1)

with fixed rational coefficients. The g coefficients remain exactly1, and
the new indices do not depend on g. Index exchange is also allowed.

**Proposition.** For any two middle words W,V in the stated class and k>=1,
no substitutions of form (1) make the coefficient-wise formal identity

    Psi_W(u',v') = 2^k Psi_V(u,v) + constant                  (2)

hold. This remains true with exchange of the two indices on either side.

**Proof.** For nonempty W the coefficient of AB is

    (Q M_W)_(0,0)=(M_W)_(3,0).

Apply M_W to the formal raw basis vector e_D=(1,0,0,0). Its first J gives
(1,1,1,1). Thereafter D,d remain odd and U is congruent to E modulo2 under
every O/Z/J letter. Hence the displayed AB coefficient is odd.
Separate substitutions (1) and index exchange preserve that coefficient.
If both words are nonempty, the left side of (2) therefore has odd AB
coefficient and the right side has even AB coefficient.

For the empty word the AB coefficient is0 and the expansion is

    Psi_empty=A+B+u*v-u-v+p*q-p-q+1.

A cross empty/nonempty identity is impossible by its AB coefficient. If
both words are empty, substitutions (1) preserve the individual A and B
coefficients1, while the right side of (2) has both coefficients2^k.
This covers all cases. QED.

The proposition is a formal-identity statement. It does not replace A,B by
free parameters in the target-coverage problem, where A=g(u),B=g(v) remain
mandatory. Extra relations on a selected pair of indices can collapse the
formal coefficients and are outside this argument. In particular the result
does not prove failure of the finite central-menu occupancy conjecture.
It stops an unsuitable coefficient-matching proof strategy so that the
constructive search can focus on coupled indices or genuinely different
physical replacements.
