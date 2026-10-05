# Later checkpoints remove a ternary guard and bypass an infinite power-of-three family

2026-10-04. **PROVED:** complete finite candidate bounds at each supplied
pre-descent checkpoint; the exact two-branch source guard; a new all-height
dyadic reduction family; its wider common-future cuts; and its power-of-three
exponent subfamily. **FINITE-EXACT:** the explicitly bounded target scan and
independent controls. **OPEN:** universal coverage and the free-parameter
children's root certificates. No novelty claim beyond the named compared banks.

## 1. Recovery and the important boundary

The [early-reroute construction](collatz_early_reroute_20261004.md) improves
the sibling clock and covers new critical states, but proves that its entire
initial-ones grammar misses every positive power of3. Other inherited rules
leave the narrow exponent class \(a\equiv3\pmod8\) as a useful next target
among \(n=3^a\). This is a schema-specific survivor, not a claim that every
other repository rule fails there.

Closest mechanisms are the [composed source-paid reroute](collatz_complement_routing_20261004.md),
the [branch-toll rank](collatz_branch_toll_rank_20261004.md), the
[exact family lift and child cuts](collatz_join_shields_and_lifts_20261004.md),
and [THM-4512, coefficient-descent cylinders](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
Canonical hostiles remain7,27,703. The corrected near miss is to impose the
full prefix's ternary branch phase when the inverse depth is shorter than
that prefix. The missing coordinate is how much ternary precision this
particular inverse word actually consumes.

The live board is **immutable source / checked prefix / inverse depth /
truncated ternary guard / carry / smaller child / proof provenance**.
Anchor: escape the proved early-reroute obstruction. Niche: compile a supplied
prefix into exact candidate guards. Wildcard: intersect a successful dyadic
family with powers of3 by an exact exponent lift.

The initial target scan below gave no earlier reroutes. An independent
boundary audit supplied the \(n=155,\ell=4,r=10\) control. Lifting that
actual diagram, then cutting to an earlier common future and widening its
terminal valuation, produced the positive result. These simplifications are
inherited proof operations, not newly invented Collatz identities.

## 2. A complete finite search at every supplied pre-descent checkpoint

Let a checked actual word \(w\) have length \(r\), total valuation \(A\), and
carry \(B\), so
\[
x=U_w(n)=\frac{3^rn+B}{2^A}.
\]
Keep the original positive odd source \(n>1\) fixed. At a checkpoint before
ordinary first descent, \(x\ge n\). For integers \(k,\ell\ge1\), test
\[
3^\ell\mid S^k(x)+1,\qquad
h=2^{\ell-1}\frac{S^k(x)+1}{3^\ell}-1,\qquad 0<h<n.
\tag{1}
\]
Here \(S(x)=4x+1\). The actual common-future diagram is
\[
n\xrightarrow{(w,a)}J
\xleftarrow{(1^{\ell-1},2,a+2k-2)}h,\qquad a=v_2(3x+1).
\tag{2}
\]
The child's penultimate node is \(S^{k-1}(x)\). In a pre-descent scan this
is at least \(n>1\), so neither word hits1 prematurely.

**PROVED finite completeness.** Every successful (1) satisfies
\[
2\le\ell\le\lfloor\log_2n\rfloor,\qquad
1\le k\le\lfloor\ell/2\rfloor,\qquad \ell\ge D(k),
\tag{3}
\]
where \(D(k)=\min\{d:3^d>2^{d+2k-1}\}\).
Indeed \(S^k(x)+1\) has 2-adic valuation exactly one, so \(v_2(h+1)=\ell\).
Positivity and \(h<n\) give the depth bound. Write
\[
h+1=\rho(x+c_k),\qquad
\rho=\frac{2^{\ell+2k-1}}{3^\ell},\quad c_k=(1+2/4^k)/3>0.
\]
If \(\rho\ge1\), then \(h-n>-1\), since \(x\ge n\); its integral difference
cannot be negative. Thus \(\ell\ge D(k)>2k-1\), proving all bounds.

The scan retains every successful candidate strictly before ordinary first
descent. When the next actual value is below the original source, it reports
that edge separately as ORDINARY_DESCENT and does not enumerate reroutes from
that final checkpoint. The negative scan result therefore concerns joins
strictly before the prioritized ordinary-descent exit; it makes no failure
claim about alternative reroutes at that exit. A forward horizon remains necessary:
the theorem makes each checkpoint finite, not the entire forward process.
A budget exit reports PENDING with the actual checked word and frontier.

## 3. The exact guard has two branches

Put \(H_w(k)=4^k(3B+2^A)+2^{A+1}\). Then
\[
S^k(x)+1=\frac{4^k3^{r+1}n+H_w(k)}{3\cdot2^A}.
\tag{4}
\]

**PROVED guard decomposition.**

* If \(\ell\le r\), the extra condition is only
  \(H_w(k)\equiv0\pmod{3^{\ell+1}}\). There is no additional ternary condition
  on the source beyond its already checked binary word.
* If \(\ell>r\), first require \(H_w(k)\equiv0\pmod{3^{r+1}}\), then require
  \[
  4^kn+H_w(k)/3^{r+1}\equiv0\pmod{3^{\ell-r}}.
  \]

In the second case the preliminary condition is equivalent to
\[
k\equiv\kappa_w\pmod{3^r},\qquad
4^{\kappa_w}\equiv-\frac{2^{A+1}}{3B+2^A}\pmod{3^{r+1}}.
\tag{5}
\]
The right side is a principal unit, and4 generates those units: the elementary
identity \(v_3(4^{3^j}-1)=j+1\) proves its order. There is one next ternary
digit at every precision. For \(w=(1,2)\), this gives \(\kappa_w=1\pmod9\),
recovering the earlier after-reset construction. For an all-ones word it
gives zero. Formula (5) must not be imposed when \(\ell<r\).

The source-relative affine child is
\[
h+1=\lambda n+\beta,\qquad
\lambda=\frac{2^{\ell+2k-A-1}}{3^{\ell-r}},\qquad
\beta=\frac{2^{\ell-1}H_w(k)}{2^A3^{\ell+1}}>0.
\tag{6}
\]
An integer child cannot be smaller if \(\lambda\ge1\). When \(\lambda<1\),
payment is exactly
\[
n>\frac{\beta-1}{1-\lambda}.
\tag{7}
\]
The compiler intersects the exact dyadic cylinder for \((w,a)\) with the
appropriate ternary guard, then starts its arithmetic family above (7).
It also removes a possible padded-root least member as described in section7.

## 4. Discovery, earlier common future, and the wider final family

The discovery word is
\[
w=(1,2,1,1,1,2,3,1,1,2),\quad
(r,A,B)=(10,15,120749),\quad \ell=4,\ k=1.
\]
Here \(H_w(1)=1645596\) has ternary valuation five. The truncated guard passes;
the full-depth phase would require \(k=10288\pmod{59049}\) and wrongly reject1.

Fixing the next source valuation to one gives
\[
n=155+131072t,\qquad h=111+93312t=\frac{729n+669}{1024}.
\]
The words are \((w,1)\) and \((1,1,1,2,1)\). Every source prefix grows,
but the final child is smaller. Their paths already meet earlier, after
seven source edges and one child edge. Cutting there widens the family to
\(n=155+4096t,\ h=111+2916t\).

The final useful union keeps only the first six source valuations
\[
v=(1,2,1,1,1,2),\qquad
x=U_v(n)=\frac{729n+925}{256}.
\]
On the exact child-odd guard
\[
\boxed{\,n_t=155+2048t,\qquad h_t=111+1458t
       =\frac{729n_t+669}{1024},\qquad t\ge0\,}
\tag{8}
\]
one has \(x=4h_t+1\) and \(x\equiv5\pmod8\). Therefore, with the actual
\(a=v_2(3x+1)\ge3\),
\[
n_t\xrightarrow{(v,a)}J_t
\xleftarrow{(a-2)}h_t,\qquad 0<h_t<n_t.
\tag{9}
\]
The inequality is \(295n_t>669\), true at the least source and thereafter.
The six source prefix coefficients are
\[
3/2,\ 9/8,\ 27/16,\ 81/32,\ 243/64,\ 729/256,
\]
all greater than one with positive carries. This uses the inherited sibling
identity \(U(4h+1)=U(h)\) at a later, source-paid checkpoint.

The union guard2048 is minimal for this fixed six-letter prefix together with
an odd quarter-child: the exact prefix alone has modulus512; requiring
\(x=5\pmod8\) fixes two more binary digits. It is64 times broader than the
initial discovery cylinder.

The efficacy split matters:

* For even \(t=2u\), \(a=3\). The family
  \(n=155+4096u,\ h=111+2916u\) has all seven source steps above \(n\);
  the last coefficient is \(2187/2048>1\). This is a genuine bypass before
  the source's first ordinary descent.
* For odd \(t\), \(a\ge4\), and
  \(J_t\le(2187n_t+3031)/4096<n_t\).
  This half already first-descends at its seventh edge.

Both halves have the smaller child in (8). All sources have root-distance
valuation \(v_2(n-1)=1\); any positive odd \(h<n\) has lower incoming graph
rank, since \(E(h)\le3(h-1)^2/4<3(n-1)^2/4=E(n)\).
The whole progression is not critical (155 itself is2 mod3); its members
outside the strengthened origin bank are critical.

### Comparison with the named earlier banks

Either dyadic half fixes first valuations \(1,2,1,1\), excluding immediate
descent, the first-reset-at-least-three switch and binary16. The only compatible
binary-debt boundary would require the fourth exponent six. The maximal initial
one-run is exactly one.

The strengthened original-source bank, including the inherited small inverse
rules, has density at most
\(\delta_0=0.4458180914749392\), proved with directed rounding in the
early-reroute note. CRT gives the same relative occupancy inside either dyadic
half of (8). Other initial-one positions have only position one here, with
indices \(k=3,6,\ldots\). Their depths start at \(D(3)=9\) and increase by at
least ten, because \((3/2)^{10}<4^3\). Their occupancy is at most
\[
\epsilon_1=\frac{3^{-9}}{1-3^{-10}}.
\]
The previous after-reset run-one row \(k=1\) is already the origin bank's
\(n=2\pmod3\) class. Its remaining rows \(k=10,19,\ldots\) start at precision31
and increase by at least30, with total occupancy at most
\[
\epsilon_2=\frac{3^{-31}}{1-3^{-30}}.
\]
Larger inverse depths only refine these guards. Each bank has its proved
height cutoff: below \(X\) only \(O(\log X)\) indices can apply. Counting residue
intersections contributes \(O(\log X)\) endpoint errors; finite truncations and
the geometric tails therefore justify natural-density limits.

Subtracting these earlier domains leaves a disjoint new dispatch class occupying
at least
\[
1-\delta_0-\epsilon_1-\epsilon_2>0.55413
\tag{10}
\]
of the **growing seven-step half**. Its density among positive odd integers is
at least the left side of (10) divided by2048. This is32 times the same lower
bound on the original discovery cylinder. The wider union has an additional
half already handled by ordinary seventh-step descent. This comparison does
not exhaust all repository constructions, and the new children are not assumed
to be grounded.

## 5. An infinite repair inside the power-of-three obstruction

For \(b\ge3\), powers of3 modulo \(2^b\) are exactly the odd residues1 or3
modulo8; the order is \(2^{b-2}\). This follows from
\(v_2(3^{2^j}-1)=j+2\), \(j\ge1\). An exact binary digit lift gives
\[
3^a\equiv155\pmod{2048}
\quad\Longleftrightarrow\quad a\equiv483\pmod{512}.
\tag{11}
\]
Thus every \(n=3^{483+512t}\), \(t\ge0\), has the smaller child in (8).
The nontrivial growing-join half is
\[
a\equiv483\pmod{1024},
\]
corresponding to the source cylinder155 modulo4096. All these exponents are
\(3\pmod8\), inside the previous narrow survivor, and every source was proved
outside the entire initial-ones reroute grammar.

The least source \(3^{483}\) has766 bits and is independently replayed. For
comparison, the original discovery cylinder selected
\(a=27107\pmod{32768}\), whose least source has42,964 bits; that earlier
literal control is also retained. Cutting the shared future and retaining the
variable terminal valuation materially widens the exponent family.

These are proved smaller dependencies for every indicated exponent. They do not
prove that all their children reach1 or that all exponents \(a=3\pmod8\) are covered.
The ternary depth changes in a controlled way: every child in (8) is
\(3(37+486t)\), so \(v_3(h)=1\). In particular, a covered tower source
\(3^a\) moves from arbitrarily large ternary valuation to exactly one,
while also decreasing the original integer. This does not assert that
the next selector rule applies to every resulting cofactor.

## 6. Bounded target probe and costs

The initial fixed targets were \(3^a\) for \(a=3,11,\ldots,67\), plus7 and703.
At each checkpoint the scanner tests every candidate permitted by (3), comparing
the direct integer guard with the independently compiled prefix-carry guard.

All eleven sources reach ordinary first descent within the declared horizon128.
Their depths are
\[
37,7,5,4,6,5,7,4,13;\qquad4,51.
\]
No earlier reroute succeeds in this fixed target list. The work is143 literal
odd edges and36,803 candidate guard tests. These negative results are preserved;
the new infinite family comes from the short-inverse-depth boundary witness.

Other fixed controls: five discovery-family parameters, including one million;
128 members of the widened2048 family, checking its efficacy split; three
symbolic exponent lifts and the two literal large sources;96 full-phase
comparisons for four words and indices1..24; recovery of the earlier
\((1,2),\ell=33,k=10\) family; and a703 scan capped at eight edges that retains
its PENDING prefix and frontier. No oracle child or known root route is a
search input.

## 7. First-hit correction and reproduction

The generic arithmetic compiler initially admitted a valid affine diagram
with root padding:
\[
\operatorname{compile\_family}((4),1,1,2):
\quad 5\xrightarrow{(4,2)}1
\xleftarrow{(2,2)}1.
\]
It is a common-future identity but not a strict first-hit word. The corrected
API checks both least-member paths before export. If either contains1 before
its last vertex, it advances one common period and records that removal.
Every fixed-word intermediate increases strictly with the parameter, so one
increment raises all possible root prefixes above1; they cannot reappear
later. The regression now returns source133, child33 and records the removed
root member. The displayed new families require no such removal.

The [script](../../04-computation/experiments/collatz_checkpoint_reroute_20261004.py)
exports exact guard compilation, first-hit-safe fixed-word family compilation,
the widened quarter-child rule and bounded scanning. Outputs are checked
smaller-child obligations or ordinary descent segments. A rooted certificate
still needs grounded child evidence.

Run:

    python -B 04-computation/experiments/collatz_checkpoint_reroute_20261004.py
    python -O -B 04-computation/experiments/collatz_checkpoint_reroute_20261004.py

The [saved output](collatz_checkpoint_reroute_20261004.out) records exact constants,
work, controls and limits. All checks remain active under optimized Python.
The next obligation is to organize further later-checkpoint families while
retaining original-source payment; universal coverage remains open.

