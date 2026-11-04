# Independent audit: retained-order budget

**Independent analytic and primary-source audit: PASS.** The numerical and
asymptotic claims of [ORDER_BUDGET.md](ORDER_BUDGET.md) are correct. This is a
budget for permanent vertices in all-strong output chains, not a restriction
on the availability of particular Hamiltonian counts.

I opened and read the primary
[Busch paper, Theorem 1 and Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
directly. For n>=3 it gives the congruence-class minimum and in particular the
weaker uniform inequality H>=5^((n−1)/3) used here. The n=23 minimum is exactly
9*5^6. This is a CITED input, not a newly proved extremal theorem.

From H_(i+1)=3H_i+d_i with d_i<=4, induction gives

    H_i<=3^i(H_0+2)−2 <3^i(H_0+2).

If every old physical name remains and at least two C2 names plus e_i total
extra names have been retained, n_i>=n_0+2i+e_i. Combining the two inequalities
with the cited lower bound and taking log base5 gives exactly

    e_i < 1−n_0+3log_5(H_0+2)+(3log_5(3)−2)i.

The strict sign is justified by the retained −2. Dividing by i proves the
limsup bound with no assumption that e_i/i converges. Independently evaluating
the coefficient gives `0.0478185834579558854...`; the exact logarithmic
expression is the theorem, and the displayed decimal is only explanatory.

For the actual (n_0,H_0)=(11,1429), four steps adding three permanent names
each would give n_4>=23 and

    H_4<=3^4*1431−2=115909 <140625=9*5^6.

These are exact integer comparisons. A larger output order only strengthens
the cited lower bound. A finite earlier prefix with this expenditure is not
ruled out by the asymptotic statement; the explicit four-step threshold pays
the constant term.

The argument allows arbitrary old-arc changes while charging them no permanent
vertex cost. It charges only vertices still present at the completed outputs;
temporary helpers later deleted do not contribute to e_i. Thus it does not
obstruct a zero-net-birth rearming step, and it does not prove that any such step
exists. Weak outputs would lose the cited strong-tournament bound and are
explicitly outside the theorem.

The pre-promotion primary read for this audit was 2,125 bytes, SHA256
`013fad111b4ec20b65c179c57d564f5394d460faff65d1277fc42280d3aae492`.
The harmless wording clarification about uncharged old-arc changes was sent
to the producer; no mathematical correction is required. There is no separate
source program for this short deduction and none is needed to replace its
exact inequalities.
