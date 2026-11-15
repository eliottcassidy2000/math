import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ribbon_intrep_span2.py')).read().split("# test")[0])
from sympy import Matrix, symbols, factor
x=symbols('x')
def nm(a): return ''.join(map(str,sorted(a))) or '0'
def Mword(u):
    M=np.identity(dim,dtype=object)
    for li in u: M=MATS[li].dot(M)
    return M
idx={'0':0,'1':1,'2':2,'12':3}
for bs in ['1.2','2.1','12.0.2','1.2.2','0.12.0','2.0.12','12.0','0.2.1','1.12.0','2.12.0']:
    u=[idx[t] for t in bs.split('.')]
    M=Mword(u); MM=Matrix(M.tolist())
    cp=factor(MM.charpoly(x).as_expr())
    ranks=[]
    P=np.identity(dim,dtype=object)
    for j in range(1,9):
        P=M.dot(P); ranks.append(Matrix(P.tolist()).rank())
    print(bs, 'charpoly', cp, 'ranks of powers', ranks)
