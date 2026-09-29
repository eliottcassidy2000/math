"""gamma_J(h) = c_J/(h^-3/2 e^-hI) across h from the jseries outputs: convergence of the fixed-J limits."""
import re, sys
hs=[int(x) for x in sys.argv[1:]] or [200,300,400,600,800,1200,1600]
tab={}
for h in hs:
    rows={}
    for line in open(f'jseries_h{h}.out'):
        m=re.match(r'\s+(\d+): ([\d.e+-]+) ([\d.e+-]+) ([\d.e+-]+) (\S+)\s+([\d.]+) ([+-][\d.]+)', line)
        if m: rows[int(m.group(1))]=(float(m.group(6)), float(m.group(7)), float(m.group(4)))
    tab[h]=rows
print("gamma_J(h) = c_J/(h^-3/2 e^-hI): modulus (argument) by h; last column |mu_hat_J(1)|")
print("  J   " + "   ".join(f"h={h:5d}" for h in hs))
for J in list(range(0,13))+[15,17,20,25,30,40,50]:
    cells=[]
    for h in hs:
        r=tab[h].get(J); cells.append(f"{r[0]:6.3f}({r[1]:+.2f})" if r else "   -   ")
    print(f" {J:3d}  "+"  ".join(cells)+f"   {tab[hs[-1]].get(J,(0,0,0))[2]:.2e}")
