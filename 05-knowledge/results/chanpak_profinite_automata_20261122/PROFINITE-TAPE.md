# A literal unit-seam family is dense in every odd profinite residue

**Proved.** Use the singleton-receiver chart of the existing Fibonacci perfect-target construction. Each positive address \(y\) supplies two actual strong tournaments differing in one marked arc, with ordinary Hamilton counts
\[
 L(y),L(y)+2,\qquad L(y)=43+81y+25f(y)+25[y\text{ odd}],
 \quad f(0)=0,
\]
where \(f(2y+1)=f(y)\) and \(f(2y)=f(y)+y+[y\text{ odd}]\). The formula also holds at zero, with a forced singleton that plays no role in the positive-address argument. The formula and the literal common arc are inherited from the repository's [Fibonacci receiver packet](../chanpak_e61_fibonacci_target_20261122/README.md). This note changes neither endpoint nor the marked arc.

For every modulus \(q\), every odd integer \(N\), and any prescribed finite binary prefix and suffix (the prefix starting in 1), there are arbitrarily large positive addresses in that cylinder for which
\[
 (L(y),L(y)+2)\equiv(N,N+2)\pmod q .
\]
Nevertheless this family does not contain the pair \(45,47\) (its minimum lower output is 43 at zero, and all positive addresses give at least 149). Stronger, the globally excluded count 7 shows that the pair \(7,9\) cannot be the two actual counts of any tournaments, although it passes every congruence in this one literal strong seam family. There is no contradiction: the witnesses to the congruences have unbounded size.

## 1. The odd-modulus tape group

An append digit \(b'\) acts on \((y,f,b)\), where \(b\) is the previous last bit, by
\[
 (y,f,b)\mapsto(2y+b',\ f+(1-b')(y+b),\ b').
\]
Over \(\mathbb Z/m\) for odd \(m\), any fixed tape from a specified incoming last bit is an invertible affine map on \((y,f)\). A loop tape that ends in its original last bit acts as a finite permutation, so every inverse of a loop is implementable by positive repetition. Reading words is always chronological; products of maps below are composition, rightmost map first.

Work first with incoming last bit 1, use \(u=y+1\), and take loops “1” and “01”. Their maps are
\[
 A:(u,f)\mapsto(2u,f),\qquad B:(u,f)\mapsto(4u-2,f+u).
\]
Set \(D=BA^{-2}\), \(E=ADA^{-1}\), \(C=ED^{-2}\). Direct composition gives
\[
 D:(u-2,f+u/4),\quad E:(u-4,f+u/8),\quad
 C:(u,f-3u/8-1).
\]
The map \(C(ACA^{-1})^{-2}\) is exactly \((u,f+1)\); this holds also modulo powers of 3. The map \(D\) translates \(u\) by a unit. Thus one first reaches any desired \(u\) by iterations of \(D\), then any desired \(f\) by iterations of the pure translation. The action of last-bit-1 loops on the odd-modulus state is transitive.

For last bit 0, let \(T\) append 1 from bit 0 and \(U\) append 0 from bit 1. Given any loop \(M\) at bit 1, \(UMT\) and \(UT\) are actual loops at bit 0. Their quotient \((UMT)(UT)^{-1}=UMU^{-1}\) proves last-bit-0 transitivity as well.

## 2. Lift the whole odd tape action through a dyadic reset

Factor a general modulus as \(q=2^k m\), \(m\) odd, and assume \(k\ge1\) (for odd \(q\) the first section alone suffices). The quantity \(f(y)\) is always even: an appended zero adds \(y+[y\text{ odd}]\), an even integer, and an appended one adds 0. After a long zero run, the dyadic state has \(y=0\bmod2^k\), last bit 0, and \(f\) in the even subgroup.

Choose \(t\ge k\) with \(2^t\equiv1\bmod m\). Appending \(0^t\) from last bit 0 is the identity on the odd state, since its exact formula is \((y,f)\mapsto(2^t y,f+(2^t-1)y)\). Hence any last-bit-0 loop may be padded by \(0^t\), preserving its odd action and returning the dyadic state to the reset. On the resulting finite state space
\[
 (y\bmod m,\ f\bmod m,\ f\bmod2^k\text{ even}),
\]
each padded loop acts by a permutation: its odd action is invertible, and its dyadic action adds a fixed even constant to \(f\) because the incoming dyadic \(y\) and last bit are fixed. Its inverse is again realized by positive repetition, now on the *joint* state space.

Take the two equal-length padded loops “\(110\,0^t\)” and “\(000\,0^t\)”. From the reset, their dyadic \(f\) increments are respectively \(-2\) and 0: the first has \(y=6,f=4\) after “110” and gains \(-6\bmod2^k\) in the padding. Over the odd modulus they have identical \(y\)-slope \(2^{t+3}\). Consequently their quotient has the affine unipotent form \((y,f)\mapsto(y+\alpha,f+\beta y+\gamma)\), which has order dividing \(m\): its \(m\)th power adds \(m\alpha\) to \(y\), and \(m\beta y+m\gamma+\binom m2\alpha\beta\) to \(f\), all zero modulo odd \(m\). Raising the joint quotient to \(m\) therefore fixes every odd state and translates the dyadic \(f\) by \(-2m\). This generates every even residue modulo \(2^k\). Combined with the transitive odd projection, **the joint state space is transitive**; both components can be prescribed independently by physical words.

## 3. Profinitely universal, even inside an arbitrary cylinder

Begin with any positive leading prefix and append \(0^t\) to reach the dyadic reset (for arbitrary incoming last bit this initial padding need not preserve the odd state; transitivity starts from whichever state is reached). Use the joint state action, then append the prescribed suffix \(r\) of length \(\ell\). On the odd state this last append is a bijection. Since \(\gcd(81,25)=1\), the final affine output \(43+81y+25f+25[y\text{ odd}]\) hits every residue modulo \(m\), with its last bit fixed. On the dyadic factor the last suffix fixes the final \(y\bmod2^k\) starting from the reset \(y=0\); it translates the freely chosen \(f\bmod2^k\) by a fixed even amount. Because 25 is odd, this gives every odd output residue modulo \(2^k\). CRT gives the desired \(N\bmod q\). One can insert arbitrarily many copies of a suitable positive power of any padded loop acting identically on the joint finite state; the addresses then grow unboundedly without changing either modulus or cylinder.

Thus the closure of the literal pair set in \(\widehat{\mathbb Z}^2\) is exactly
\[
 \{(u,u+2):u\equiv1\bmod2\}.
\]
The inclusion from left to right follows because \(f\) is even and \(81y+25[y\text{ odd}]\) is even; the other inclusion is the simultaneous arbitrary-modulus construction. Already one fixed first/last-bit cylinder has the same closure. Since \(L(y)\ge43+81y\), any finite target still has only finitely many possible addresses and many have none. Modulo \(q\) the append rule is a finite automaton, but no intersection or selection relying solely on congruence attainability (even with finitely many forced initial and terminal contact bits) can distinguish a missing integer pair from this family. Exact absolute inverse methods remain indispensable.

## 4. All congruence fibers have the same limiting frequency

For each fixed \(q\), choose the free digits between any fixed prefix and suffix independently and fairly. As their number tends to infinity, the output residues tend to the uniform law on the allowed (odd) residues modulo \(q\), exponentially for each fixed \(q\).

Indeed, use the full modular state \((y\bmod q,f\bmod q,b)\); when \(q\) is even require \(f\) even and \(b=y\bmod2\), and when \(q\) is odd permit both \(b\). The two digit actions make a finite strongly connected multigraph. For even \(q\), every state can reach the reset by zeros, the joint action is transitive there, and any output dyadic \(y\) value can be set by appending its \(k\) padded bits. That last append is bijective on the odd coordinates and translates the even \(f\) subgroup, so every full state is reachable. For odd \(q\), use the transitive loop and a final bit to reach the desired last-bit chamber.

Every state has two incoming edges counted with multiplicity. At even \(q\) the output last bit determines the appended digit and \(2y+b'=y'\bmod q\) has exactly two solutions for the incoming \(y\); each determines the preceding bit and the unique even \(f\) residue. At odd \(q\) the incoming \(y\) is unique and each of the two preceding last bits gives one incoming state. Since outdegree is also two, the fair transition matrix has uniform stationary measure. The all-zero state has a zero self-loop, so convergence is exponential at each fixed modulus. Applying the fixed suffix to uniform states leaves \(f\) uniform in its even dyadic subgroup independently of other relevant inputs, and gives uniform odd coordinates conditional on the incoming last bit and the dyadic residues, because the fixed append is invertible modulo \(m\). The output linear form then has the stated uniform law. This strengthening still gives no exact hit for any finite target absent from the absolute chart.

## 5. Complete local test for every Fibonacci receiver chart

The tape proof applies unchanged to any exact chart \(K_j+cy+hf(y)\) with \(h\) odd and outputs odd, as for the normalized receivers with \(h=25H(V)\). Set \(g=\gcd(c,h)\), necessarily odd. In the chamber of last bit \(j\), its profinite image is exactly
\[
 \{u\in\widehat{\mathbb Z}: u\equiv1\bmod2,\ u\equiv K_j\bmod g\}.
\]
Necessity is the integer linear form. Sufficiency follows from the transitive odd state, since the image of \((y,f)\mapsto cy+hf\) in \(\mathbb Z/m\) is precisely the subgroup generated by \(g\); the dyadic reset independently supplies every odd residue as before. Any prescribed suffix of that parity and any valid leading prefix are allowed. With no prescribed terminal bit, take the union of the two chambers. In particular, no refinement involving higher prime powers or additional primes can reject a target that survives this gcd test in one chamber.

For the bank of 610 normalized receivers through order nine at \(E_{61}-1\), 372 charts survive this complete local test in both parities, another 24 only in parity zero, and 214 in neither. In particular the singleton and the four-vertex chart among the three order-123 height survivors have \(g=1\), so no modulus rejects their unrestricted parity chambers. The remaining cyclic-triangle chart has \((c,h,K_0,K_1)=(411,75,367,392)\); its two output residues modulo 3 are 1 and 2, while \(E_{61}-1\equiv0\pmod3\). Thus the height table actually reduces to **two** order-123 candidates by the complete local test. This does not settle *fixed order*: the local witnesses constructed here can have arbitrarily many free bits and the height filter constrains their absolute lengths. The positive result is an exact limit on what unbounded congruence tests can certify, paired with the already-existing finite-width inverse question.
