# Open-middle endpoint-two repair

**PROVED.** Let n>=9. Every closed blue marked word P=(a,M,z)
with endpoint scores (n-3,2), an open middle M, and exceptional port
u other than the last two middle vertices has a literal three-black
detour. The complementary penultimate-port branch is treated in
[probe.md](probe.md), and the last-port branches in
[last_port.md](last_port.md) and [lift.md](lift.md).

Here open means first(M)->last(M); it is a chord condition, not a
statement that the middle tournament is nonstrong. As in
[theory.md](theory.md), u->a is the unique extra incoming middle arc,
and z->v with v=r_P(u) is its reflected partner.

Write P=(a,b,W,z), x=first(W), c=last(W). Since u is not last or
penultimate, v is neither b nor x. In particular b->z and x->z.
Because M is open, b->c. Use the inherited guarded cyclic ledger

    Q1=(b,W,z,a),
    Q2=(z,W,b,a),
    Q3=(a,z,W,b).

The guards are z->a, b->c, x->z. The directedness and exact path-edge
XOR are the guarded three-move identity in
[kernels_games_detours_20261020f/detours.md, Section 2](../kernels_games_detours_20261020f/detours.md).
These are genuine K moves reaching B=K_P(T) on the same labelled set.

The endpoint score excesses (sum minus n-1) at their successive sources
are d_T(b)-2, -1, and d_T(b)-3. Vertex b beats the three distinct vertices
x,c,z, so Q1 and Q2 are black. If d_T(b)>3, Q3 is black as well.

Suppose d_T(b)=3. At B, which is also an endpoint of the Q3 move,

    N+_B(a)={b,z,u},     N-_B(b)={a,c,z}.

Reversal of Q3 pairs a with b and z with c. If Q3 were blue, it would
send the first displayed neighborhood to the second, forcing
r_Q3(u)=z. Thus u=c, contrary to the hypothesis. This covers the score
equality boundary without discarding the actual port.

The proof uses neither strongness of B nor a module assumption on M.
The order restriction here is the common synthesis bound; the local
argument also applies at smaller orders whenever its distinct-vertex
and port hypotheses hold.
