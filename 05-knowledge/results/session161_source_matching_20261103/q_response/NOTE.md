# Equal edge jumps with opposite triangle responses

Status: proved for explicit pairs of strong Boolean two-apex tournaments. Source: `05-knowledge/results/apex_wedge_constructor_20261103/{THEORY,PROFILE}.md`, especially A2, P2, and P5, and S159's full path-cover transfer. The argument concerns an edge jump, not arbitrary prescribed absolute path counts or a global spectrum theorem.

## Family and exact response

For integers `j>=2` and `1<=q<j`, take

    P = 2^j + 2^q - 1,       R = P + 2.

Write `p` for the ordinary `j+1`-bit binary expansion of P, and write `s` for the reverse of the bitwise complement of the `j+1`-bit expansion of R. As in A2, set `w=p 01 s`, `u=p 10 s`, build the transitive base `0,...,m-1` with `m=2j+4` vertices, and add apexes `a,b` with a->i iff `w_i=1` and b->i iff `u_i=1`. Denote the two choices a->b and b->a by T+ and T-. Both are strong: w and u start in 1 and end in 0. Both have exactly `2j+6` vertices and differ only at the named edge.

Let `Delta_r=F_r(T+)-F_r(T-)`, with `F_r` the ordered r-slot path-cover count. Define

    L_h(r) = ((r+1)^h+(r-1)^h)/2,
    T_q(r) = L_q(r) - sum_{a=1}^{q-1} r^a L_{q-a}(r).

Then for every integer r>=1, indeed as a polynomial identity,

    Delta_r = 2 r^(2j+2-q) T_q(r),                         (1)
    T_q(r) = r^q + (1-r)(r+1)^q/2 + (r+1)(r-1)^q/2.       (2)

Proof: P and R have the same `j+1` binary length, and their changed tail is `0 1^q` to `1 0^(q-1) 1`. In the source weighted charge `C_r` the changed position begins at index `j-q`. P2 therefore gives

    C_r(R)-C_r(P) = r^(j-q) T_q(r),
    Delta_r = 2 r^(j+2) (C_r(R)-C_r(P)),

establishing (1). The finite geometric sum

    sum_{a=0}^q r^a L_{q-a}(r)
      = ((r+1)^(q+1)-(r-1)^(q+1))/2

proves (2). Since `T_q(1)=1`, every source pair has the *same* positive original edge jump `Delta_1=2`.

Let G be the strong directed-triangle wrapper. Its positive path-cover transfer is linear on differences, and at r=1 yields `Delta'_1=Delta_1+Delta_2`. Consequently the literally wrapped pair differs at the same named edge and has exact Hamilton-path jump

    W(j,q) = 2 + 2^(2j+2-q) (2^(q+1) - 3^q + 3).          (3)

The bracket equals 4 and 2 for q=1 and q=2. For every q>=3 it is negative: `3^3>2^4+3` and induction preserves the inequality. Thus W is positive for q=1,2 and negative for q>=3. These are unequal physical source pairs with the same source order and scalar jump, not a claim that their absolute source counts agree.

For fixed j>=4, q=1,2,3 give three strongly connected source pairs of the same order and with original jump +2, with respective full responses and wrapped jumps

    q=1:  2r^(2j+2)             -> 2+2^(2j+3)
    q=2:  2r^(2j)               -> 2+2^(2j+1)
    q=3:  2r^(2j)(2-r^2)        -> 2-2^(2j+2).

The first two reproduce the source P3 and P4 (the q=1 pair has the same profile as P3); the third realizes an equal-order, opposite-sign divergence. Formula (3) generalizes all three and demonstrates why retaining only `Delta_1` cannot predict this wrapper's jump.

## Exact recurrence along an explicit subfamily

Fix any integer c>=1 and set j=q+c, q>=1. Define

    x_q = (W(q+c,q)-2) / 2^(2c+2)
        = 2*4^q - 6^q + 3*2^q.

This is integral and obeys the exact recurrence

    x_q = 12 x_(q-1) - 44 x_(q-2) + 48 x_(q-3), q>=4,
    x_1=8, x_2=8, x_3=-64.

Indeed its characteristic roots are 4, 6, and 2; the recurrence follows from `(X-4)(X-6)(X-2)`. Each term is the normalized response of an actual strong pair on `2(q+c)+6` source vertices. This is a recurrence in the family parameter q, not a recurrence claimed for repeated G wrappers.

## Independent finite controls

Run `python3 check.py` or `python3 -O check.py`. Direct integer charge calculations check (1), (2), (3), and sign for j=2,...,14, all 1<=q<j, and r=1,...,7: 637 profile evaluations. An independent subset/endpoint Hamilton-path dynamic program builds the literal tournaments, not their profiles, and checks the q=1,2,3 instances at j=4 on 14 vertices and their wrappers on 16 vertices. The outputs are:

    q   source H(T-) H(T+)    wrapped H(G(T-)) H(G(T+))
    1       37215   37217             4783239    4785289
    2       84381   84383             8794485    8794999
    3      198021  198023            16474725   16473703

These controls are not needed for the all-parameter algebraic proof.
