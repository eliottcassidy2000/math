# Exact joint wheels and six primes at a single flip Jacobian

**FINITE-EXACT** for the group generators and the six primes; **PROVED** for lifting the wheels to arbitrary powers. See the [main circle transfer](FILLER-CIRCLES.md). Use the eight-site third-slot frame with \(h=11D+42Q+20R\) and \(d=-3D-6Q\). All bank indices below are zero-based, as in the [literal bank](third-slot-bank.py).

**An 82-vertex wheel.** Take quotas \(3,1,1,3\) of types \(2,6,11,13\). The ratios reduce modulo 77 to \(25,59,16,39\). Their bounded exponent vectors (64 choices) generate all 60 units jointly, not just in separate projections. Therefore \(h\) runs through every phase with \(h-35\) a unit modulo 77, at one fixed physical order and fixed \(D,Q,d\). The raw starts are automatically \(1\bmod3\) and \(3\bmod5\), as well as odd. This does not imply that every target in a legal phase occurs.

**Two larger joint wheels.** The generator lists are:
- primes \(7,13,17,19,29,37,41,47\), bank types \(3,4,6,8,9,10,12,13\);
- primes \(13,23,31,37,41,47\), bank types \(1,2,3,6,7,13\).

The [checker](filler-crt-check.py) verifies generation on the simultaneous product of the listed prime squares. It computes discrete-log rows in the cyclic unit groups and checks spanning modulo every prime \(q\) dividing their order, in the coordinates where \(q\) divides that cyclic order. A finite abelian subgroup is the whole group when it maps onto \(G/qG\) for every such prime \(q\). For every higher prime power the latter quotient stabilizes: the cyclic order at \(p^e\) is \(p^{e-1}(p-1)\), and reduction to \(p^2\) induces the same relevant mod-\(q\) log maps for \(e\ge2\). Projection handles any primes used only to the first power. Therefore these exact same types generate at arbitrary prescribed positive powers **simultaneously**. At any chosen modulus take a common quota at least the unit group's exponent minus one, allowing all generator exponents to be chosen within the physical quota. Only unit \(R\) phases are claimed.

The six-prime set also admits ordinary small-prime avoidance: take the quota a multiple of four so \(D,Q\equiv1\bmod5\); \(D\) is a unit mod 3 while \(R\) is zero there. Hence the cofactor is prime to 2,3,5. Each of the three affine reads \(h,H_{10},H_{01}\) excludes at most one unit \(R\) phase at each listed prime, because the slopes \(20,116,44\) are units there. There is a joint phase that keeps all three reads units, as in the earlier same-clock wheel. The eight-prime list alone does not protect against small primes outside that wheel.

**A small literal prime grid.** Use one factor each of types \(0,8,11\) in that order. Every choice lives on 37 vertices and has \(D=7683910325,\ Q=36380429405056,\ d=-218305628161311\). Of the eight starts, six are prime, as certified by exact recursive Pratt certificates:

| Choice (0 = first member) | Start | Result |
|---|---:|---|
| 000 | 460648142166497107 | prime |
| 001 | 461036408494562707 | prime |
| 010 | 461046213093554707 | prime |
| 011 | 461434816060212307 | divisible by 109 |
| 100 | 461025215748156307 | prime |
| 101 | 461413800957885907 | prime |
| 110 | 461423613609357907 | prime |
| 111 | 461812535734159507 | divisible by 17 |

For each prime the same marked pair gives its exact sevenfold endpoint and both edit orders increase. The counts use the proved arbitrary-module identity rather than a subset census on 37 vertices. This is an existence example with no claim of minimum order or infinitely many prime starts.

Reproduce the compact wheel and the circle checks with `python3 filler-circle-check.py`; reproduce the group lifts and six primality certificates with `python3 filler-crt-check.py`.
