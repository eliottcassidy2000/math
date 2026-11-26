"""Exact uncentered Gram and frozen three-star spectral counts."""
from collections import Counter
from nested_core_check import orbit_chart
from slot_interface_check import read_board, BRANCHES
from antipodal_options import compressed
import sys

P,_=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
F=(56,0,28)
O=[h for h in range(57) if h not in F]
d=56
C=[set() for _ in range(168)]
B=[set() for _ in C]
for i,h in enumerate(F):
    for x in range(d):
        for j,k in enumerate(F):
            if i!=j:C[i*d+x].add(j*d+P[h,k][x])
        for j,k in enumerate(O):B[i*d+x].add(j*d+P[h,k][x])
for i,X in enumerate(B):
    for j,Y in enumerate(B):
        assert len(X & Y)==d*(i==j)-len(C[i]&C[j])-(j in C[i])+(i//d!=j//d)
assert P[56,0]==tuple(range(56))==P[56,28]
assert P[0,28]==tuple((x+28)%56 for x in range(56))
Ctypes={2:27,-1:54,1:56,-2:28}
Gtypes=Counter()
for e,m in Ctypes.items():Gtypes[d-e*e-e]+=m
assert sum(Gtypes.values())==165 and Gtypes=={50:27,56:54,54:84}
n=57*55
assert n==3135 and n*8//15==1672
print("Frozen Gram eigenvalue multiplicities:",sorted(Gtypes.items()))
print("Paired packets 165+165; free spectrum",1672-165,(n-1672)-165)

def action_halves(board, ys):
    _,D=read_board(board)
    _,cap,_=compressed(board)
    assert ys[0]==7 and len(ys)==53
    slots=[D[h][y] for h,y in enumerate(ys) if h]
    E=coarse=imbalance=0
    for g in range(3):
        C=Counter(p[g] for p in slots)
        for r,mask in enumerate(cap[g]):
            a,b=C[r],C[r+28]
            assert mask&1 or not a
            assert mask&2 or not b
            if mask==3:
                E+=(a-1)**2+(b-1)**2
                coarse+=(a+b-2)**2
                imbalance+=(a-b)**2
            elif mask:
                E+=(a+b-1)**2
                coarse+=2*(a+b-1)**2
    assert 2*E==coarse+imbalance
    return E,coarse,imbalance

if len(sys.argv)>1:
    board=sys.argv[1]
    r,_=read_board(board)
    w=[{13:16,33:19}.get(h,z) for h,z in zip(BRANCHES,r)]
    a=action_halves(board,w)
    assert a==(6,6,6)
    print("Marked baseline: action, twice coarse, twice phase",a)
    for path in sys.argv[2:]:
        w=list(map(int,open(path).read().split()))
        assert w
        print(path,action_halves(board,w))
