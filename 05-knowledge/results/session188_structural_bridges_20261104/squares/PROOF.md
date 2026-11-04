# Parity inventory in growing square-sum graphs

Let G_n have vertices 1,...,n and an edge between distinct vertices whose sum is a square. Write O_n and E_n for the numbers of odd-odd and even-even edges, respectively. Then, for all n>=1,

    O_n-E_n = floor(sqrt((n-1)/2))    if n is odd,
              floor(sqrt((n+1)/4))    if n is even.

The integer formula can be evaluated with isqrt((n-1)//2) or isqrt((n+1)//4). The newly added vertex n creates exactly floor(sqrt(2n-1))-floor(sqrt(n)) edges. The square-root formula for births counts k with n<k^2<=2n-1; its upper bound excludes the forbidden loop n+n.

**Proof of parity inventory.** Same-parity edges have even-square sum (2k)^2. Each such edge has the unique form {2k^2-t, 2k^2+t} with integer 1<=t<=min(2k^2-1,n-2k^2). The edge is odd-odd for odd t and even-even for even t. Thus each k contributes one to O_n-E_n exactly when its available radius is positive and odd, otherwise zero. If n is odd, every positive available radius is odd: both candidate bounds are odd. There are floor(sqrt((n-1)/2)) such k. If n is even, the unsaturated radius n-2k^2 is even, whereas the saturated radius 2k^2-1 is odd. Saturation occurs exactly when 4k^2<=n+1, yielding floor(sqrt((n+1)/4)). QED.

**Path consequence.** Let n>=2 and let P be any Hamiltonian path, with p_oo and p_ee the numbers of odd-odd and even-even edges used, and let e_o,e_e count its endpoint incidences of odd/even parity. Degree counting gives

    p_oo-p_ee = (number odd vertices)-(number even vertices)-(e_o-e_e)/2.

Consequently p_oo-p_ee lies in {0,1,2} for odd n and {-1,0,1} for even n. Every Hamiltonian path therefore omits at least max(0,O_n-E_n-2) odd-odd edges for odd n, or max(0,O_n-E_n-1) for even n. For a Hamiltonian cycle, p_oo-p_ee equals the odd-even vertex count, giving omission bounds max(0,O_n-E_n-1) for odd n and O_n-E_n for even n. These bounds concern the unused inventory, not the existence of a path or cycle.

The sharp fluctuation between adjacent n is structural: all active even-square families contribute one at odd n, whereas only saturated ones do at even n. For example the inventory difference at n=25 is 3 and at n=26 is 2, despite the latter graph containing all the former edges.

`verify.py` independently enumerates all edges for every n<=500, comparing the displayed formula. Both normal and optimized Python runs give the attached receipt. This is an elementary identity; no literature novelty is claimed, and it does not prove the square-sum Hamiltonicity problem.

**Order 25 illustration.** Literal depth-first enumeration in `path25.py` recovers all 20 directed Hamiltonian paths. Vertex 18 is the unique degree-one vertex, hence an endpoint of each path. Their (p_oo,p_ee) distribution is (5,4):2, (6,4):4, (6,5):6, (7,5):4, (7,6):2, (8,6):2. The gap is always 1 or 2, as forced by odd n and the even endpoint 18. Since G_25 has O_25=9 and E_25=6, the all-order inventory lower bound says at least one odd-odd edge is unused; the distribution shows the 20 paths omit between 1 and 4 odd-odd edges. This finite enumeration is a checked illustration, not a proof about other n.
