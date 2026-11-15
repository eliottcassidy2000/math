import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, sys
sys.path.insert(0, RIB)
from ribbon import letters
from sympy import Matrix
w=int(sys.argv[1])
R=os.path.join(RIB,f'span{w}_intrep.txt')
tok=open(R).read().split(); pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=[nxt() for _ in range(dim)]; fin=[nxt() for _ in range(dim)]
MATS={}
for _ in range(nl):
    li=nxt(); MATS[li]=Matrix([[nxt() for _ in range(dim)] for _ in range(dim)])
def mask(tokn): return sum(1<<(int(c)-1) for c in tokn if c in '123')
def Mw(block):
    M=Matrix.eye(dim)
    for t in block.split('.'): M=MATS[mask(t)]*M
    return M
for b in sys.argv[2:]:
    M=Mw(b); P=Matrix.eye(dim); ranks=[]
    for j in range(1,9):
        P=M*P; ranks.append(P.rank())
    tr=M.trace()
    print(b,'trace',tr,'ranks of powers',ranks)
