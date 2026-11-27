# Cullen and Woodall: reflected odd runs and a two-sign digit machine

This is the intended Cullen–Woodall continuation of [the guarded machines note](../guarded_moore_pollock_collatz_cullen_machines_20261127/README.md). Wilson and Woodall are different; the Wilson checksum in that note is separate arithmetic. Here use the classical notation of [Keller's study](https://cir.nii.ac.jp/crid/1362262944614601984):
\[
C_m=m2^m+1,\qquad W_m=m2^m-1,\qquad Y_m=m2^m,\qquad m\ge1 .
\]

## 1. Exact pairing of the signed Collatz branches

Run the halved \(3x+1\) map on \(W_m\) and the halved \(3x-1\) map on \(C_m\), only as long as both inputs are odd. For any even center \(Y\),
\[
\frac{3(Y-1)+1}{2}=\frac{3Y}{2}-1,\qquad
\frac{3(Y+1)-1}{2}=\frac{3Y}{2}+1.
\]
Thus both maps share one evolving center and keep the endpoints exactly two apart. After \(j\) legal odd steps the center is \(3^jm2^{m-j}\), interpreted as an integer since \(j\le m+\nu_2(m)\). Precisely
\[
m+\nu_2(m)
\]
odd steps are legal before both endpoints first become even. Every two steps multiply the common center by \(9/4\), giving \(\lfloor(m+\nu_2(m))/2\rfloor\) full pairs and possibly one final odd step. At that first even state exactly one endpoint is divisible by four, so the common odd-step law does not extend to the subsequent division runs.

For \(m=3\) the Woodall number is \(23\) and the Cullen number is \(25\). Their synchronized states are
\[
(23,25)\to(35,37)\to(53,55)\to(80,82).
\]
The fact that 23 is also the prime at which the modular \(9/4\) clock stalls for one level is a numerical coincidence here, not a consequence of the signed-run proof. The statement is independent of whether either input is prime.

## 2. One digit controller, two disjoint targets

For odd \(p\), let \(d=\operatorname{ord}_p(2)\). At a fixed phase \(m\bmod d\), divisibility of Cullen by \(p^r\) asks for \(Y_m=-1\bmod p^r\), whereas Woodall asks for \(Y_m=+1\bmod p^r\). The same machine solves either sign; at each next precision the two sets are disjoint. In particular the two numbers are coprime at every fixed positive index (their difference is 2 and both are odd).

**Common-clock theorem.** Pick a finite set of distinct odd primes, one sign \(\epsilon_p\in\{-1,1\}\) and precision \(r_p\ge1\) for each. Choose a common clock \(L\) divisible by each \(\operatorname{ord}_p(2)\), and a phase \(j\in[0,L)\), and restrict to \(m=j+Lt\). Write \(s_p=\nu_p(L)\). There is a simultaneous solution if and only if the frozen digits pass:
\[
j2^j\equiv\epsilon_p\pmod {p^{\min(s_p,r_p)}}\quad\text{for each }p.
\]
When this holds, the permitted \(t\)'s form exactly one residue class modulo \(\prod_p p^{\max(0,r_p-s_p)}\).

Indeed \(2^L=1\bmod p^{s_p+1}\) by binomial lifting. Varying \(t\) by \(h p^b\) at any stage changes \(Y_m\) by
\[
h\,p^{s_p+b}(L/p^{s_p})2^j\pmod {p^{s_p+b+1}}.
\]
The coefficient stripped of its \(p\)-power is a unit, independent of the digit already written. The first \(s_p\) digits are frozen, each later digit has exactly one choice, and the remaining classes in \(t\) combine by CRT. Taking one prime and \(L=d\) recovers the ordinary single-prime controller for either sign. A phase obstruction is real when another prime's clock imports a factor \(p\) into \(L\); treating all local choices as independent without first testing the phase is incorrect.

Example: with \(L=6\), require \(27\mid C_m\) and \(49\mid W_m\). In phase \(m=0\bmod6\) the Cullen condition already fails mod 3. In phase \(m=2\bmod6\) both conditions reduce to the single class \(m=1262\bmod2646\). Exchanging the signs does not work in that phase, since the Woodall target at 3 is frozen out. The modular existence statement produces progressions of factor certificates, not prime certificates or bounded absolute witnesses for unrelated equations.

Run `node replay.js` from this directory for the signed Collatz checks and a brute-force comparison of the common-clock inverse over 612 prime/phase/precision/sign choices.
