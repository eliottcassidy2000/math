# The seam reader has two incompatible automata clocks

Use the fixed strong [RESET-PORT](../chanpak_mirimanoff_unit_lifts_20261122/RESET-PORT.md) chart, with lower count \(L(y)=149+162y+25f(y)\), companion \(L(y)+2\), and
\(f(2y+1)=f(y),\ f(2y)=f(y)+y+[y\text{ odd}],\ f(0)=0\).
The subscript \(y\) is the integer address, not the address length.

## Exactly which residue columns are eventually periodic

**Theorem.** \(L(y)\bmod M\) is 2-automatic for every \(M\ge1\), and eventually periodic exactly when \(M\mid50\).

The reader states \((y\bmod 2M,f(y)/2\bmod M)\) are finite. The append rule respects them, proving 2-automaticity. Since \(f\) is even, the claim when \(M\mid50\) is immediate.

Conversely assume an eventual period \(D\ge1\). Take \(B=2^b>2D\), \(u=B/2-D\), \(v=B/2\), and let \(\operatorname{rev}\) reverse \(b\) padded bits. Use the exact inherited suffix formula
\[
f(Bz+w)=f(z)+(B-1-\operatorname{rev}(w))z+f(w)+[w<B/2][z\text{ odd}].
\]
The difference in \(L\) at \(Bz+v\) and \(Bz+u\) is \(K+25dz-25[z\text{ odd}]\), where \(d=\operatorname{rev}(u)-\operatorname{rev}(v)\). It vanishes modulo \(M\) for all large \(z\), so its alternating first differences \(25(d-1)\) and \(25(d+1)\) both vanish. Hence \(M\mid 50\).

For an odd prime \(\ell\), apply [Cobham's theorem](https://research.ibm.com/publications/on-the-base-dependence-of-sets-of-numbers-recognizable-by-finite-automata) to the finitely many symbol fibers: \(\ell\)-automaticity of this same sequence holds exactly for \(\ell=5\). This excludes 13 of the 14 odd Ogg primes. In the [Rowland–Yassawi model](https://arxiv.org/abs/1209.6008), a column of a linear cellular automaton with memory over a finite field of characteristic \(\ell\) and eventually periodic initial conditions is \(\ell\)-automatic. Thus no such CA in an odd characteristic other than 5 can produce the full residue sequence in ordinary address time, even over an extension field. This does not forbid nonlinear CA, other initial data, or nontrivial sampling schedules. Nor does it obstruct a **binary** controller from compiling addresses for any prime.

## Exact state cost for odd moduli

**Theorem.** For odd \(M\), put \(d=M/\gcd(M,50)\). The minimum size of a most-significant-bit-first binary DFA with output for \(L(y)\bmod M\), allowing leading zeroes, is \(M\) if \(d=1\) and \(2Md\) otherwise.

Use \((x,e,z)=(y\bmod M,y\bmod2,f(y)/2\bmod d)\), forgetting \(e,z\) if \(d=1\). On appending bit \(b\), \(z\) increases by \((1-b)(y+e)/2\bmod d\); two is a unit there. All states are reachable by projecting the fully reachable \(2M^2\)-state arithmetic machine of RESET-PORT. Let \(O_w\) be the output after appending \(w\). Independently of the original \(e,z\),
\[
O_1-O_{\emptyset}=162(x+1),\qquad O_{00}-O_0=374x .
\]
The second identity uses the zero step at the even address \(2y\). Since \(\gcd(162,374)=2\), these recover \(x\bmod M\). Then \(O_{\emptyset}\) recovers \(z\bmod d\) by canceling \(\gcd(50,M)\). For \(d>1\) the two possible values of \(e\) give readings \(O_0\) separated by \(25\bmod M\), which is nonzero. Thus continuations of length at most two distinguish all reachable states.

For the squarefree product \(m\) of the odd Ogg primes the exact state count is \(2m^2/5=262104763919816946870742191364383210\). This is not a lower bound of that size on arithmetic memory: a modular tape machine stores \(O(\log m)\) bits. The two clocks are different: the DFA advances in binary address digits, while the excluded CA advances in the integer address itself.
