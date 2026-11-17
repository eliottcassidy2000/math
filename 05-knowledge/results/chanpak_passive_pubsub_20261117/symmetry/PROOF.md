# One physical two-arc switch with endlessly nested visibility and exact +2

A single hidden port change can remain undetectable to one observer, remain visible to a richer observer, and be exposed as an exact unit repair by growing the same subscriber. This is a literal nested tournament tower, not a claim of Chan–Pak completeness.

## Rooted observation quotients

Write (h,R,S) for a tournament with ordinary Hamilton count h, weighted full-or-optional path covers ending in its designated port (R), and full paths starting at that port (S). The [ordinary migration subscriber](https://github.com/eliottcassidy2000/math/blob/f9b24f62064508d7ceca42a337ebaadb91f0d22b/05-knowledge/results/chanpak_bank_fusion_20261117/migration/PROOF.md#the-three-coordinate-root-transfer) with binary word u reads

    h' = F(u) h + (u+1)R − [u even]S

(F(0)=1; F(2v+1)=F(v), F(2v)=F(v)+v+[v odd]). Its rooted successor coordinates R',S' are integer linear functions of h,R,S; specifically R'=2^k R + (1+RC_k(u))h, S'=h for positive standard binary words and S'=h−S for u=0.

Consequently the exact quotient for observers all of whose first migration words are odd is (h,R), even if the later words are unrestricted. The first odd read and its whole successor state depend only on h,R; conversely the reads at u=1 and u=3 are h+2R and h+4R, which determine both. The hidden S direction is therefore **erased** by fresh odd migration of the port. By contrast the full subscribers distinguish all three coordinates; h from the empty read, then R from u=1 and S from u=0 suffice. Even the same u=0 subscriber applied repeatedly is observable by counts at times 0,1,2: these are h, h+R−S, h+2R. These claims concern the rooted boundary states, not tournament isomorphism.

## The literal primitive pair and growing the same subscriber

Take one five-vertex tournament with outgoing masks (20,1,2,7,14). The port r=0 has (h,R,S,e,L,P₂)=(9,11,1,1,11,84), while r=3 has (9,11,3,1,23,84), where e counts full ends and L weights optional covers starting at the port. Roots are structurally distinct (their outdegrees are 2 and 3).

For any positive u, use a new subscriber s with transitive spine B of standard binary word u (bit 1 means s beats that vertex), all B→A, and s beating all A except the chosen r, which beats s. Both graphs are strong. Their only arc differences are the two s contacts incident to old 0 and 3. Their counts are:

    odd u:   9F(u)+11(u+1) on both tracks;
    even u: 9F(u)+11(u+1)−3 and that value +2.

Now append more vertices at the end of **the same B**, between old B and A, retaining s and all its old arcs. Appending one bit gives induced extensions of both graphs, so we can alternate hiding (append 1) with revealing (append 0) indefinitely. This is distinct from migrating the port at each step: the latter would consume the hidden state, while here the same old graph behind the subscriber is preserved.

For example for u the word (10)^k the even shell gives the adjacent odd values

    (34·4^k+26)/3+6k−3 and +2.

Between consecutive even shells, appending 1 gives equal counts. Alternatively u=2^k−1 gives common 9+11·2^k, and appending a 0 gives 31·2^k−5 and +2. At u=12 the exact two values 239,241 are both primes; at u=424 they are 8191,8193. No prime infinitude is inferred.

## Rich observers distinguish every odd hidden shell

Let P_q be the optional q-slot polynomial. Let E_q,O_q count selected nonempty paths ending at the old root with even/odd vertex count, weighting the complement by P_q. Take differences port 3 minus port 0. Then

    Delta E_q = −q³,    Delta O_q=q².

For E_q in order five only degrees 3 and 1 occur: the first is the indegree difference −1, and at q=1 the signed endpoint identity and equal R give difference −1, so the q coefficient is 0. For O_q the possible degrees are 4,2,0. Singleton leading terms and full-end constants cancel; at q=1 the difference is +1. This proves both polynomial formulas.

For the subscriber stage S=B+s, let U_e,U_o be the corresponding weighted counts of even and odd selected paths STARTING at s. The only crossing cycles in S→A except r→s use this sole return arc, so each uses one run in S and one in A. Their total length is odd iff the two parities differ; all leftover ordinary packs lie separately inside the two blocks. The pure terms cancel since the old graph is the same. Thus

    Delta P_q = 2q(U_e Delta O_q + U_o Delta E_q)
              = 2q³(U_e−q U_o).

For i counted from zero in a k-bit word, d=k−1−i. A first spine selection at a bit-1 position i gives factor q^i; subsequent positions are freely selected or left in the complement. Hence

    U_e+U_o = q^k + sum_(bit i is 1) q^i(q+1)^d,
    U_e−U_o = −q^k + sum_(bit i is 1) q^i(q−1)^d.

When u is odd, at q=1 the indicated difference Delta P_1 vanishes. At every q>1, the singleton-stage contribution to U_e−qU_o is −q^(k+1), while the last bit contributes +q^(k−1). Every other bit's contribution is nonpositive, since (1−q)(q+1)^d+(1+q)(q−1)^d≤0 for d≥1. Thus **every** richer slot read q>1 detects the switch on every odd shell, even where the whole ordinary Hamilton count hides it. For u=2^k−1, the explicit detection P₂ difference is −8(3^k+3).

## Any prescribed odd phase, not every absolute value

At an even read write u=2a, F(u)=1+2b. Its lower member is 17+22a+18b. For an odd target N, put q₀=(N−17)/2. The congruences a=5q₀ and b=−6q₀ give exactly the required linear phase because 22·5−18·6=2. Apply the pinned [protected-address compiler](https://github.com/eliottcassidy2000/math/blob/1266194f74c870905649eab3ccb75974bda8d783/05-knowledge/results/chanpak_wheel_square_selector_20261117/CONGRUENCE.md) with modulus `2M`, desired `u≡10q₀` and `F(u)≡1−12q₀`, any fixed prefix starting in 1, and suffix 0. The requested address is even, so this suffix is compatible with the dyadic prescription; the requested F is odd, as required. The lifted integer read is `N mod 2M`, hence `N mod M`. It therefore gives a physical exact +2 pair in every odd lower residue modulo any M. This relies on actual joint address closure, not arbitrary independent numerical substitution into F. Exact equality to a prescribed N is not guaranteed; for example direct bounded evaluation u≤N/11 rejects 65537 in either track of this fixed family.

## An intentionally different reflection reader

Both old ports still have e=1 and can be used as the two alternative marked inside flanks in the Farey construction, against any even ordinary right flank. Each window keeps the reflected +2. But this reader does **not** respect the (h,R) quotient of odd migrations: it also sees L=11 versus 23. Its windows differ by `12(147+256w)−10[u odd]`, where u is the left binary prefix and w is the even right address. The difference follows from the same P₂,R numerics, local apex increment L−1, and signed end-mask difference 2. This illustrates why the quotient is strategy-dependent.

Replay: `python verify.py` (also with `-O`). No random search premise; the graph and root coordinates are checked by literal path DP.
