# Weighted prime debt and the two-event boundary form

**PROVED elementary extensions of the inherited block action.** The scalar
entry theorem and the unrestricted Hamilton-path spectrum remain OPEN. The
prime-count asymmetry theorem, the singleton/D minimum-cost inverse, and the
unweighted p-singular bound are inherited, not new here.

## 1. Source, operation and retained clock

Use the actual marked parent (A,x,S) and state (F,q,g,h) of overnight packet h:
F=H(A+x), h=H(A), g counts parent paths ending in S, and q counts ordered
spanning covers (P,Q) with Q nonempty and first(Q) in S. Append an actual
nonempty tournament B after A, with specified cut S_B. Its profile is
(H,E,P,Q,s,t,f), where H=H(B), E counts ordered spanning two-path covers
allowing empty slots, P counts covers whose first slot is nonempty and ends
outside S_B, Q counts those whose second slot starts in S_B, s counts full
paths starting in S_B, t counts those ending outside S_B, and f=H(B+x).

The inherited exact action is

    F'=HF+Pq+(H-s)g+(f-H)h,
    q'=Eq+Qh, g'=(H-t)h, h'=Hh.                         (1)

Every letter retains B's adjacency, cut, position and cost |B|. The source
labels and the integer clock are not replaced by these scalar coordinates.
In particular an inverse arithmetic solution still needs an actual word.

The current baseline is `spectrum_descent_jacobian_20261102/DESCENT.md`:
it already gives minimum-cost full-state membership in the singleton/D
language. Its admissible state (9,7,3,5) is globally unrealizable despite
attainability of scalar9. We do not claim another solution to that inverse.
The present results prune and explain candidate branches.

## 2. Weighted valuation gate

For a prime p and each block occurrence i, set

    w_i=min(v_p(H_i),v_p(E_i)).

Let W=sum w_i and let w_(1)>=w_(2) be the two largest weights, padded with
zeros if necessary. For every actual integral initial state and every such
word, without restrictions on the intervening blocks,

    v_p(F_final) >= W-w_(1)-w_(2).                       (2)

In particular a target with v_p(N)=e rejects every word for which the right
side exceeds e. For a prime target N different from p this right side must
vanish. When all positive weights are1, (2) is the **already proved**
r-2 theorem in Session124 `PROOF.md`, last section. The weighted statement
uses the sizes of individual prime powers as well.

**Proof.** Order the coordinates as h,(q,g),F. In the matrix of (1), every
nonzero entry is either a diagonal loop or raises this three-level order.
There are no transitions between q and g; the g diagonal is zero. A nonzero
term in a matrix-product entry therefore uses at most two off-diagonal
transitions. Every other block contributes a diagonal H or E and hence at
least w_i powers of p. A term can omit the valuation payment at no more than
two positions, whose combined weights are at most w_(1)+w_(2). Each term,
and hence their integer sum, is divisible by the stated power. Multiplying
by an integral input cannot lower it. QED.

This is a divisibility lower bound, not an equality: p-adic cancellation can
increase it. Exact min-plus propagation of the valuations of **all entries**
of (1) gives a possibly stronger lower bound. An inverse routine may retain
the maximum of these authenticated bounds. Replacing a block by a coarser
or finer factorization can change the bound; retain earlier certificates
instead of assuming every compression preserves sharpness.

For a nontrivial weighted control, the strong nine-vertex circulant with
steps{8,2,6,5} has H3159=3^5*13 and E48924=3^4*604, hence weight4 at3.
Three such blocks give the bound v_3(F)>=4, whereas counting merely the
three singular events gives1. Cuts1,3,7 (bitmasks) give an actual28-vertex
word with output9418116918213 and v_3=10. The output is checked by exact
profile transport; no full order28 subset-DP run is claimed.

**The two unpaid transitions are necessary.** Let C_j denote a directed
triangle block whose cut has j vertices. Starting at (1,0,0,1), C_1 C_1 gives

    (F,q,g,h)=(41,45,3,9).

For every r>=1, append r empty-cut triangles C_0. The exact actual graph has

    H=45*12^r-3^r=3^r(45*4^r-1),     |V|=7+3r.           (3)

Thus v_3(H)=r, exactly the bound (r+2)-2. The graph is strong: each parent
triangle is strong, the marked vertex reaches the first component, every
later component reaches the mark, and the parent components point forward.
The endpoint r=0 has H41, not44; the first appended block consumes g=3.
The general empty-cut diagonal law is inherited from AMPLIFIERS.md. Formula
(3) is its sharp valuation application, with the initial endpoint retained.

## 3. Two prime-symmetric events: adjacency versus separation

Here is finer information than the number of singular blocks. Suppose B_1
and B_2 each have p vertices and are invariant under a p-cycle, p an odd
prime. For i=1,2 define

    H_i=p*u_i,  E_i=p*v_i,  d_i=v_i-u_i.

Free action on ordered path covers gives integrality; cyclic symmetry makes
endpoint and cover-start counts uniform. For a cut of size a_i,

    Q_i=a_i*d_i, P_i=(p-a_i)*d_i,
    s_i=a_i*u_i, H_i-t_i=a_i*u_i.                        (4)

Let h_0 be the parent count just before B_1. Work modulo p.

* If the two blocks are consecutive, their output satisfies

      F_2 = -a_1 a_2 (u_1 u_2+d_1 d_2) h_0 (mod p).     (5)

* If a nonempty actual block word W separates them, write E(W) for the
  product of its E values. Then

      F_2 = -a_1 a_2 d_1 d_2 E(W) h_0 (mod p).           (6)

After B_2, q,g,h are all0 modulo p. Any subsequent suffix multiplies this
F residue by its product of H values. Thus (5)-(6), including prefix and
suffix factors, are complete necessary p-residue tests for the two-event
word, independent of the earlier F and of all cuts inside W.

**Proof.** After B_1, modulo p, h_1=0, q_1=a_1 d_1 h_0 and
g_1=a_1 u_1 h_0. If the next letter is B_2, substitute these three values
in (1) to obtain (5). If at least one intervening letter occurs, its update
sets g=0, while q is multiplied by E. All remaining intervening letters
preserve h=g=0. The second singular letter then gives (6) and annihilates
q,g,h. The suffix statement follows immediately from (1). QED.

This is an exact map from actual two-event words to a bilinear form over
F_p. Consecutive events retain the vector (u,d); a nonempty intervening word
projects away the endpoint term u_1 u_2. It can therefore destroy a residue
**or remove a cancellation**. The retained sidecar is adjacency/separation,
the cuts, h_0, E(W), and the suffix H product—not the scalar output alone.

For repeated identical blocks, the adjacent coefficient is u^2+d^2. If
p=3 mod4 and (u,d) is nonzero modulo p, it cannot vanish: a nonzero ratio
squaring to -1 would contradict Fermat's little theorem. For p=1 mod4,
isotropic directions can occur. This is the specific finite-field connection;
it is not a claim that all prime behavior reduces to a tournament diagram.

| block | p | (u,d) modulo p | adjacent u^2+d^2 | separated d^2 |
|---|---:|---:|---:|---:|
| directed triangle |3|(1,0)|1|0|
| cyclic regular five |5|(3,1)|0|1|
| circular seven, H175 |7|(4,0)|2|0|
| Paley seven, H189 |7|(6,2)|5|4|

Consequently, two triangle blocks separated by **any** nonempty block word
force divisibility by3; this extends the inherited singleton-separation rule.
Two circular-seven blocks have the analogous7 obstruction. For regular-five
blocks, adjacency forces divisibility by5 while separation need not do so.
That last statement concerns the5-residue only; other prime obstructions may
still prohibit a prime output. Neither a nonzero residue nor avoidance of
(2) proves an actual desired integer is reached.

The dyadic-cover blocks D=(H5,E32) and B19=(H19,E256) have w_p=0 for every
odd p. This explains their escape from the inherited singular-block debt.
It does not overcome the proved sparsity of a fixed empty-cut amplifier bank.

## 4. A stronger actual-source terminal budget

If S is nonempty, append one zero singleton (A->z->x). Its output N obeys

    N >= (|S|+2) H(A),  in particular H(A)<=floor(N/3).   (7)

For each Hamiltonian word of A, at least |S|+2 distinct extensions preserve
that word: insert the adjacent pair z,x immediately before each S vertex;
append z,x after the whole word; and insert x at the start of any run of S
and put z last. There is at least one such run. All displayed words are
legal, the families have different positions of z, and removing x,z
recovers the old word. Orphan paths only increase N. This proves (7).

This is a budget for the stated positive one-zero terminal source, not for
arbitrary scalar witnesses. It can be combined with the inverse cost bound,
the positive coordinate guards and the weighted valuation test. Universal
selection of an actual state above each target N remains OPEN.

## 5. Reproduction and boundary

`python3 valuation.py` reconstructs314 actual cut profiles by literal
permutation/cut enumeration, checks12500 word/prime valuation inequalities,
5376 two-event formulas,122 separate full-adjacency counts,54 terminal
bounds and the exact sharp family through r50. All assertions remain active
under Python -O. The all-word conclusions follow from the proofs, not the
finite tests. Same-agent review with independent counting paths is used;
no independent subagent or formal verification is claimed.
