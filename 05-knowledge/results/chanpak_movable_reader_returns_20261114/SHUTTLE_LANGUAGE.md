# The whole fixed-center binary shuttle language

Here we let a movable reader trade a binary suffix with a movable source.
For the 15-path switch center from [bilateral.md](bilateral.md), the
simple four-bit move is only the shortest member of a closed two-state
language. This classifies every shuttle of the displayed form that holds
uniformly in the outer source addresses of given parities. It does not
classify coincidences at just one pair of addresses.

Use the positive binary addresses x,y, the count N(x,y), and the positive
root response R(x,y) of the two-wing chart. Fix one of the two choices
K=41,43. Both choices have (h,f,P,Q,S,T)=(15,120,42,63,6,9).
We compare the literal wings
\[
                  (2^k X+a,Y)\quad\longleftrightarrow\quad(X,2^kY+b),
                  \qquad 0\le a,b<2^k,\quad X,Y>0.
\]
Let reverse always use exactly k bits, including leading zeros. Equality
of R and equality of the coefficient of X in the count difference are
equivalent to
\[
 f a+h\operatorname{rev}_k(b)=(2^k-1)(P+h),\qquad
 f b+h\operatorname{rev}_k(a)=(2^k-1)(Q+h).
\]
Indeed the first follows directly from the rooted response formula
as the right width increases by k; it also balances the Y coefficient
of N. For the other coefficient, appending a k-bit word to address X
adds its internal contribution plus
\(X(2^k-1-\operatorname{rev}_k(a))\) and a boundary parity term.
The same coefficient equations also preserve the other unsigned response
at the converse apex. They are necessary for all X,Y of given parities,
since either variable independently ranges through an unbounded progression.

## Paired carry, read in opposite directions

Put c=rev(b), K1=P+h, K2=Q+h. Read a and c from the low end, with
bits d,e. The second equality is then being read backwards. The state
pair (u,v) starts and ends at (0,0) and updates to
\[
         ((u+f d+h e-K1)/2,\;\;2v-f e-h d+K2).
\]
Require the first entry integral and retain only states where the first
entry lies in [−K1,f+h−K1] and the second in [−K2,f+h−K2].
These bounds follow from the ordinary forward signed-carry recurrence
being the average of its previous value and the local signed digit.
The path formulation is exact in both directions: telescoping the
scaled signed carries gives the two full integer equations.
For this profile the productive bit graph has 12 states. Its edges go
between four successive layers cyclically. Reading four bits at a time
compresses it to the following transducer (the first digit in a pair
goes into a from the low end; the second goes into b from the high end):

| old state | new state | allowed (a digit, b digit) |
|---|---|---|
| 0 | 0 | (6,9), (7,8) |
| 0 | 1 | (2,9), (3,8) |
| 1 | 0 | (6,13), (7,12) |
| 1 | 1 | (2,13), (3,12) |

The original states here are (0,0) and (−30,30). The accompanying script
constructs the bounded carry graph from the equations, prunes it
backwards from (0,0), and checks the table and four layers exactly.
Thus every nonempty admissible width is k=4r. All intermediate states
and all binary choices within an edge are free. Start and end in state 0.
This gives \(2^{2r-1}\) address pairs at width 4r.

## The remaining absolute correction is a boundary term

Let d_first,d_last be the parities of the first and last a-digits
in processing order. The exact ordinary count difference (left displayed
tournament minus right) over the entire shuttle is
\[
 \Delta N=6\bigl(p(X)+p(Y)-1+d_{\rm last}-d_{\rm first}\bigr),\qquad
 p(z)=z\bmod2 .
\]
For proof, write L_0=0, L_1(x,y)=−30x+30y. In each table edge s→t
with digits a0,b0, direct substitution of the binary word recurrence
gives, for positive x,y,
\[
 N(16x+a0,y)-N(x,16y+b0)
 = L_s(16x+a0,y)-L_t(x,16y+b0)
   +6(p(x)+p(y)-1)+12(s-t).
\]
This identity has only two linear coefficients and four residual parity
values to check per edge. Both potentials telescope. At an interior seam
the parity of x is the next a-digit parity, and the parity of y is one
minus the previous a-digit parity. Their contribution telescopes to
d_last−d_first. Thus neutrality is a single explicit boundary equality
even for arbitrarily nonperiodic sequences of intermediate states.

Every such neutral shuttle preserves the full rooted profile: R_+ agrees
by the carry equations, and R_- = −15 on both wings since the last
transitive vertex beats the apex and the old puncture has exactly 15
paths. On opposite outer parities, there are \(2^{2r-2}\) neutral pairs
at width 4r for r≥2, since only the first and last edge bits must agree.
The pair of wings produces strong tournaments and changes the size of
the transitive left section by 4r. Deleting the marked apex recovers the
unique nontrivial strong component, so these are distinct rooted classes.
The ordinary count chart and its absolute inverse still apply.
This proves a large nonperiodic supply of exact changes of source and
reader, without implying that an arbitrary prescribed N lands in the
same physical chart. The number of shuttle words is not a spectrum
density estimate.

Run the bounded certificate and arithmetic replay: `python3 replay_shuttle.py`.

## Target-directed search including the shuttle condition

The digit law also removes all branching once one physical witness is
known. From a left input (x,Y), choose r and strip its r low hexadecimal
digits, leaving X>0. Starting at state 0, each stripped digit a0
determines the next state (2 or 3 means state 1, 6 or 7 means state 0)
and its parity determines the unique b0 for that transition. Reject
other digits and a nonzero final state. The resulting b-digits appear
in big-endian order, so the new right address is 16^r Y+b. The
parity boundary test completes the exact decoder; no approximate
comparison or residue-only match is needed. For a given left input
all lengths take O(log x) arithmetic steps, by growing the digit
scan and checking the new boundary at each stage. Applying this to
the fixed-center exact absolute inverse yields a decision and witness
for this whole uniform-shuttle relation at a given absolute target.
That test is allowed to reject.

For a literal example through the nonzero internal carry, choose
(X,Y)=(1,2), r=2 and (a,b)=(98,157). The two 17-vertex graphs have
addresses (354,2) and (1,669). Both materialize the same count and
both rooted responses, verified directly by subset Hamiltonian DP.

## A canonical representative, without comparing count tables

Orient each neutral shuttle from left to right. The number of left
binary digits decreases, so the rewrite always terminates. If two
lengths r<s are available at once, the shorter valid digit scan is
literally a prefix of the longer scan, ending in the initial automaton
state 0. The remaining s−r digits therefore form another accepted scan
with exactly the updated right address. Its count difference is zero
because the two endpoint equalities subtract, and its rooted response
also agrees. Taking the shorter move then this suffix is exactly the
longer move. Hence all forks join, and every class generated by these
neutral shuttles has a unique right-normal representative. This is
a small physical rewriting system: it canonicalizes a specific family
of exact returns, not all graphs with the same (H,R+,R−).

One can find the representative by taking the longest valid initial
move. If a further move were then possible, concatenating the two
accepted scans would have been a longer valid initial move, contrary
to the choice. Thus the single growing digit scan gives the normal form
and its literal graph, rather than removing duplicates by their scalar
Hamiltonian path counts.
