"""Greatest-fixed-set finite quotients, rather than seeded address expansion."""
from fractions import Fraction as F
import json

def need(b,m):
    if not b:raise ArithmeticError(m)
CHARTS=[(1,1,1),(4,6,179),(7,11,16177),(7,11,12289),(7,11,9697)]

def main():
    counts=[]
    for L in range(18):
        M=1<<(L+1);S=set(range(1,M,2));maps=[(1<<D,C,pow(3**J,-1,M)) for J,D,C in CHARTS]
        while True:
            T={(a*y-b)*u%M for y in S for a,b,u in maps}
            need(T<=S,'descending full-set attractor')
            if T==S:break
            S=T
        counts.append(len(S))
    initial=[1,2,3,4,6,9,13,18,24,31,40]
    need(counts[:11]==initial,'initial occupied cells')
    for i in range(11,18):need(counts[i]==counts[i-1]+counts[i-6]+3*counts[i-11],'occupied recurrence')
    a=F(1357705,1000000);b=F(1357706,1000000)
    f=lambda x:x**11-x**10-x**5-3
    need(f(a)<0<f(b),'growth bracket')
    need(a**1000>2**441 and b**1000<2**442,'dimension bracket')
    need(all(F(n)<=F(17,8)*a**i for i,n in enumerate(initial)),'upper constant')
    # Integer counting uses L=floor(log2 X), not an infinite-address count.
    for X in [1,2,3,4,7,8,15,16,100,1000]:
        L=X.bit_length()-1;M=1<<(L+1)
        need(M>X and 1<<L<=X,'at most one bounded integer per occupied cell')
    print(json.dumps({'fixed_set_cell_counts':counts,'verified_precision_levels':18,
        'lambda_bracket':[str(a),str(b)],'dimension_bracket':['441/1000','442/1000'],
        'upper_constant':'17/8'},sort_keys=True))
if __name__=='__main__':main()
