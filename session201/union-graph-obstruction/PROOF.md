# The twelve-template union witness relation is not binary-automatic

Status: **PROVED / INDEPENDENTLY AUDITED**. This is a separate corollary
packet and does not alter the frozen digital-graph proof. It concerns the
three-variable witness relation, not its one-variable target support.
The twelve-template all-target occupancy assertion remains **OPEN**.

The strict-minimum argument below extends the individual-template
obstruction to the union. This packet checks all matrices and includes
the short nonregularity dependency explicitly.

## 1. Actual physical cone and strict central increase

Use source columns s=(D,E,U,d)^t, seed (1,0,2,1)^t, and

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d),
    J(s)=(D-E+U,D-E,D+U,D).

Let K(s,t)=s^t Q t, with

    Q=[0  0  0  1]
      [0  1  0 -1]
      [0  0  1 -1]
      [1 -1 -1  1].

The nonnegative physical coordinates y=(a,f,e,p)^t are defined by

    y=(D-E-d,U-D+E-d,E,d-E)^t,   s=P y,
    P=[1 0 2 1]
      [0 0 1 0]
      [1 1 2 2]
      [0 0 1 1].

The seed has y=(0,0,0,1). The generator actions in these coordinates are

    O: (a,f,e,p) -> (a+e,a+2f+e+2p,0,e+p),
    Z: (a,f,e,p) -> (2a+f+e+p,f,e+p,0),
    J: (a,f,e,p) -> (f+p,0,a+e+p,e).

Thus all actual physical sources have y>=0 and d=e+p>=1. The root size d
never decreases under these operations.

For a generator c, the exact increment matrix is

    R_c=P^t Q(M_c-I)P.

Direct integer multiplication gives

    R_O=[1 1 2 2]        R_Z=[1 1 1 1]
        [1 1 2 2]            [1 1 1 1]
        [1 1 2 2]            [2 2 2 3]
        [1 1 3 2],           [2 2 2 2],

    R_J=[1 0 2 1]
        [0 0 1 1]
        [2 1 3 3]
        [1 1 3 2].

Every entry is nonnegative, and every entry of the lower-right2-by-2 block
is at least2. Consequently, for any two actual physical sources s=P y and
t=P z,

    K(s,M_c t)-K(s,t)=y^t R_c z
       >=2(y_e+y_p)(z_e+z_p)=2d_s d_t>=2.               (1)

Iterating a chronological word W and using nondecreasing d_t yields

    K(s,M_W t)-K(s,t)>=2|W| d_s d_t.                    (2)

In particular every nonempty word increases K strictly. This is proved on
the actual nonnegative source cone, not by treating arbitrary signed
profiles as physical sources.

For actual unit sources

    tau(u)=(g(u),u mod2,u,1)^t, u>=2,

put Psi_W(u,v)=K(tau(u),M_W tau(v)). Formula (2) gives

    Psi_W(u,v)>=Psi_empty(u,v)+2|W|.                    (3)

## 2. The empty template is an automatic minimum selector, if the union were

Let C be the twelve-word menu, including the empty word, and define

    R(u,v,n) iff u,v>=2 and n=Psi_W(u,v) for some W in C.

For every fixed pair(u,v), (3) makes Psi_empty(u,v) the unique minimum
OUTPUT of the finite set of possible n. Nonempty words may share other
outputs; no uniqueness among them is asserted.

Suppose R were recognized by a finite automaton on padded, synchronized
binary representations of(u,v,n). Binary integer order is a regular relation.
Regular languages are closed under Boolean operations and existential
projection. Therefore the formula

    R(u,v,n) and not exists m [m<n and R(u,v,m)]          (4)

would define another binary-automatic relation. By the strict minimum just
proved, (4) is precisely the graph n=Psi_empty(u,v).

Now fix u=2. Since g(2)=1 and2 is even, the empty formula becomes

    Psi_empty(2,v)=g(v)+v-(v mod2), v>=2.                (5)

Section3 proves that the binary graph of (5) is not regular. This contradicts
the automaton obtained from (4). Hence the union witness relation R is NOT
binary-automatic, even though it is a union of only twelve explicit forms.

The same argument applies to any collection of O/Z/J words that includes
empty: its minimum at fixed indices is still empty by (2). This extension
concerns only the same index-preserving witness encoding.

## 3. Self-contained nonregularity dependency

The source recurrence is

    g(1)=1, g(2h)=g(h), g(2h-1)=g(h)+h-(h mod2).

It first gives g(2^r+1)=2^r+1 for r>=1 by induction. It then gives, for all
a>=2,b>=1,

    g(2^(a+b)+2^b+1)=2^(a+b)+2^b+1-2^(a-1).            (6)

For b=1, take the even parent2^a+2, strip its factor2, and apply the previous
fixed-point identity. Increasing b applies the odd-parent step to2h-1 and
preserves h-g(h), proving (6).

Use pair letters A=(1,1), Z=(0,0), B=(1,0), C=(0,1). For k,m>=1, the word

    A Z^k B C Z^(m-1) A

encodes the pair

    u=2^(k+m+2)+2^(m+1)+1,
    y=2^(k+m+2)+2^m+1.

By (6), y=g(u) if and only if k=m. Thus the infinitely many prefixes
A Z^k are pairwise distinguishable by suffixes B C Z^(m-1) A. No finite
automaton can recognize the binary graph of g, even restricted to u>=2.

Finally, if h(v)=g(v)+v-(v mod2) had a regular binary graph, intersect its
lift with the relation n=y+v-(v mod2), then project away n. Fixed-coefficient
addition and parity are finite-carry regular relations, so this would give
a regular graph y=g(v), a contradiction. Common leading zero padding does
not affect either argument. This proves the required nonregularity of (5).

The regular/synchronized vocabulary is the standard one of Carpi--Maggi,
*On synchronized sequences and their separators*, RAIRO35(6) (2001),513--524:
https://www.numdam.org/item/ITA_2001__35_6_513_0/
The local g-specific obstruction above is proved here, not attributed to
that paper.

## 4. What has and has not been ruled out

The conclusion rules out an ordinary finite automaton recognizing all exact
index-index-output witnesses in this union. It does NOT rule out a finite
automaton for the projected target set

    {n: exists u,v R(u,v,n)}.

Projection of a nonregular relation can be regular; that implication cannot
be reversed. In particular a cofinite odd target support remains possible.
No counterexample to the occupancy conjecture has been obtained. Algorithms
with unbounded integer or carry state are also not ruled out.

## 5. Portable replay

    python3 verify.py
    python3 -O verify.py

The verifier contains all scalar and matrix definitions locally. Its only
input file is the copied twelve-entry `INPUT-MENU.json`, whose hash is
reported. It verifies the exact generator matrices, cone invariance,
strict lower-right blocks, and agreement of every frozen menu matrix with
its physical word. Declared finite controls check actual sources and the
unique minimum. The all-order claims follow from the displayed identities
and the infinite separating-prefix proof, not these finite controls.
