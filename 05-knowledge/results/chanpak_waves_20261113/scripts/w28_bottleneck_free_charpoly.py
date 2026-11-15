import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ribbon_intrep_span2.py')).read().split("# test")[0])
from sympy import Matrix, symbols, factor, Poly, nroots, factor_list
x=symbols('x')
idx={'0':0,'1':1,'2':2,'12':3}
def Mw(b):
    M=Matrix.eye(dim)
    for t in b.split('.'): M=Matrix(MATS[idx[t]].tolist())*M
    return M
for b in ['0.2.1.2','0.12.2.2','1.2','0.2']:
    cp=Mw(b).charpoly(x).as_expr()
    fl=factor_list(cp)
    print(b, fl)
    for fac,mult in fl[1]:
        if Poly(fac,x).degree()>1:
            rts=[r for r in nroots(fac) if abs(r.as_real_imag()[1])<1e-12]
            print('   real roots', rts, ' 4th roots of max', max(abs(complex(r)) for r in nroots(fac))**0.25)
