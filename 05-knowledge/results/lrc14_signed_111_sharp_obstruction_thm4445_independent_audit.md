# Independent audit of THM-4445

**Verdict: PASS.** No constant, equality locus, carrier statement, parity
leader, endpoint policy, or finite cutoff needs mathematical correction.

The referee is a third implementation. It imports neither the candidate nor
a repository computation engine. A modular one-coordinate kernel solver,
separate ray compiler, and generic rational capped-line integrator agree on
all 1,901 typed rows through height 223.

It independently reproduces:

- the unique sorted sign cone \(a+b=c\);
- the complete ray \(k(1,1,-1)\) with strict roof and deleted-third address;
- continuum physical bulk \(9/98\) and network range
  \([9/98,39/392]\);
- the two-sided error \(4/(7c)\);
- sole row at or below \(6/77\): \((1,4,5)\), value \(1/28\);
- nonexception lower endpoint \(31/392\), uniquely \((1,7,8)\);
- upper endpoint \(6/55\), uniquely \((1,10,11)\);
- every parity-sheet maximum and first strict upper cutoff.

The raw arithmetic first makes the lower bound exceed \(6/77\) at \(c=42\).
Since \(3\mid42\), the theorem's \(c=43\) is exactly the first admissible
height. The analogous lower cutoff for \(31/392\) is 45; the head through
44 and the tail therefore meet without a gap.

The ordinary and optimized runs have the same normalized transcript and
execute 417,669 checks. This proves a local obstruction classification, not
body entry, synchronization, or LRC(14).
