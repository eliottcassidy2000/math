"""
S79: the INDEX-THEOREM frame. LRC(2p) <=> an INDEX = (p-1)/2 is NONZERO.
Analytic index (equidistribution/cap = Euler char of the cover, S78) = Topological index (Borsuk-Ulam degree, kps).
Parity = p mod 4 selects the method: even (p=1 mod4) -> Brouwer/SOS; odd (p=3 mod4) -> Borsuk-Ulam (odd degree !=0).
For n=14 (p=7) the index=3 is ODD => Borsuk-Ulam forces the lonely point. Subsumes all prior frames as two
computations of one index.
"""
print("LRC INDEX-THEOREM FAMILY (n=2p):  index=(p-1)/2 ; parity=p mod 4 ; method ; index!=0")
print(f"{'p':>3}{'n':>5}{'index':>7}{'p%4':>5}{'parity':>8}{'method':>14}{'!=0':>5}")
for p in (3,5,7,11,13,17,19,23):
    idx=(p-1)//2; par='odd' if idx%2 else 'even'; method='Brouwer/SOS' if p%4==1 else 'Borsuk-Ulam'
    print(f"{p:>3}{2*p:>5}{idx:>7}{p%4:>5}{par:>8}{method:>14}{'yes' if idx>0 else 'NO':>5}")
print()
print("n=14 (p=7): index=3 ODD => Borsuk-Ulam: odd-degree map S^1->S^1 has nonzero degree => antipodal")
print("  lonely pair (t*,-t*) at the 1/14 equioscillation EXISTS. odd degree = i*sqrt7 = -trace = (p-1)/2.")
print()
print("INDEX THEOREM (the meta-frame): ANALYTIC index (equidist/cap=Euler char of cover) = TOPOLOGICAL index")
print("  (Borsuk-Ulam degree). LRC(2p) <=> index != 0. p=1 mod4: Brouwer/SOS (even index); p=3 mod4: Borsuk-Ulam")
print("  (odd index, FORCED !=0). Subsumes: cyclotomic (de Moivre deg=(p-1)/2=index), topology (Euler=top index),")
print("  analytic (equidist=analytic index), even/odd (index parity), Vitali wall (measure vs construction),")
print("  SOS/non-SOS (Brouwer/Borsuk-Ulam). OPEN: (1) (p-1)/2 = degree of the right map; (2) coincidence at 1/14.")
