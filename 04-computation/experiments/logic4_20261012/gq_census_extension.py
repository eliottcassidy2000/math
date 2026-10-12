"""Extend the G_q automorphism census to prime powers 83 <= q <= 128 (HYP-9256 evidence)."""
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
import sys
qs=[q for q in range(82,129) if factor_pp(q)]
print("q  p^e   |V|   deg  |Aut|          |Aff|=q^2(q+1)2e  ratio  orbits(Aut)  |Stab_0|  orbits(Stab_0)")
for q in qs:
    adj,p,e,mu,elems,idx,add=gq(q)
    g,orb,out=nauty_grpsize(adj); aff=q*q*(q+1)*2*e
    gs,orbs=nauty_stab_orbits(adj,0)
    G=float(g) if g else float('nan')
    print("%-3d %d^%d  %5d  %3d  %-14s %-15d %-6.4g %-12s %-9s %s"%(q,p,e,q*q,q+1,g,aff,G/aff,orb,gs,orbs)); sys.stdout.flush()
