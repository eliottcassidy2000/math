"""
KEY: tight => COVERING (multiple of every q in {2..n-1}), the q-witness at t=1/q (THM-523).
Averaged lenses miss this. Verify tight sets cover; AP=minimal cover (q itself); doubling replaces q by mq.
"""
def frac(x): x=x%1.0; return min(x,1-x)
def covers_q(S,q): return any(v%q==0 for v in S)
def covering_report(S,n):
    miss=[q for q in range(2,n) if not covers_q(S,q)]
    return miss
n=14
print("tight => contains a multiple of every q in {2..n-1} (else q-witness t=1/q gives M>=1/q>1/n):")
for name,S,tight in [("AP {1..13}",list(range(1,n)),True),("GW {..11,13,24}",list(range(1,12))+[13,24],True),
        ("swap12->26",[v for v in range(1,n) if v!=12]+[26],False),
        ("drop 7 (no mult of 7)",[v for v in range(1,n) if v!=7],False)]:
    miss=covering_report(S,n)
    print(f"  {name:>22} (tight={tight}): missing multiples of q = {miss if miss else 'NONE (covering)'}")
print()
print("the DOUBLING as covering-preserving: GW replaces runner 12 (covers q=12) by 24=2*12 (STILL covers q=12):")
for q in [12]:
    print(f"  q={q}: AP uses {q} itself; GW uses 24={24//q}*{q} (a MULTIPLE of {q}, still covers q). Both cover q.")
print("  => the doubling operad k->2k preserves the covering (2k is a multiple of k) -- that's WHY it's a valid")
print("     retiling. The finer rational conditions (Jacobsthal) gate WHICH q can use a multiple instead of q.")
print()
print("SYNTHESIS of the lens-hunt:")
print("  * AVERAGED lenses (int m^2 moments; additive energy) = DEAD ENDS (AP+1 max energy but non-tight;")
print("    non-tight 12->26 has lower int m^2). They average away the pointwise covering.")
print("  * FOURIER-POSITIVITY (sum f(vt)>=1) is the RIGHT framing but LOCALIZES to the combinatorial covering")
print("    at rationals t=1/q -- it IS the q-witness (THM-523), not a genuine analytic min.")
print("  * The tight-skeleton = COMBINATORIAL COVERING (mult of every q) refined by finer Farey conditions.")
print("    Lowness: covering forces a mult of every q; minimal = AP (q itself); sporadics = bounded gated")
print("    multiple-swaps (doubling). The right lens is COVERING/COMBINATORIAL, not averaged/analytic.")
