import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ribbon_intrep_span2.py')).read().split("# test")[0])
from sympy import Matrix
def nm(a): return ''.join(map(str,sorted(a))) or '0'
res={}
for n in range(1,7):
    cnt=0; r1=[]
    for u in itertools.product(range(4),repeat=n):
        M=np.identity(dim,dtype=object)
        for li in u: M=MATS[li].dot(M)
        r=Matrix(M.tolist()).rank()
        res[u]=r
        if r==1: r1.append('.'.join(nm(L[i]) for i in u))
    ranks=sorted(set(res[u] for u in itertools.product(range(4),repeat=n)))
    print('len',n,'rank values',ranks,'rank-1 words',len(r1), r1[:30])
