#!/usr/bin/env python3
"""Independent tiny exhaustive checks of the pair-state dictionary."""
from itertools import product, permutations, combinations
from collections import Counter
from fractions import Fraction
import json

def ham(a):
    n=len(a)
    return sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
def decode(a,mode):
    n=len(a)
    return [[0 if i==j else (a[i][j] if mode=='ordinary' else (1-a[j][i] if mode=='swap' else a[i][j]*(1-a[j][i]))) for j in range(n)] for i in range(n)]
def matrices(n,states):
    pairs=list(combinations(range(n),2))
    for word in product(states,repeat=len(pairs)):
        a=[[0]*n for _ in range(n)]
        for (i,j),(x,y) in zip(pairs,word): a[i][j]=x;a[j][i]=y
        yield a

def main():
    out={'status':'PASS','all_four_state_checks':[],'ternary_checks':[]}
    for n in range(1,5):
        universe=0; ordinary=Counter(); swapped=Counter(); fibers=Counter(); total_cancel=0
        for a in matrices(n,[(0,0),(1,0),(0,1),(1,1)]):
            s=decode(a,'swap');c=decode(a,'cancel');h=ham(a);hs=ham(s)
            assert decode(s,'swap')==a
            assert h%2==hs%2
            assert all(a[i][j]-a[j][i]==s[i][j]-s[j][i] for i in range(n) for j in range(n))
            ordinary[h]+=1;swapped[hs]+=1;fibers[tuple(map(tuple,c))]+=1;total_cancel+=ham(c);universe+=1
        assert ordinary==swapped
        for c,count in fibers.items():
            missing=sum(c[i][j]==c[j][i]==0 for i,j in combinations(range(n),2))
            assert count==2**missing
        import math
        assert Fraction(total_cancel,universe)==Fraction(math.factorial(n),4**(n-1))
        out['all_four_state_checks'].append({'n':n,'universe':universe,'cancellation_image_size':len(fibers)})
        total_h=total_c=universe=0
        for a in matrices(n,[(1,0),(0,1),(1,1)]):
            c=decode(a,'cancel');assert decode(a,'swap')==c
            ha,hc=ham(a),ham(c);assert ha>=hc and ha%2==hc%2
            total_h+=ha;total_c+=hc;universe+=1
        assert Fraction(total_h,universe)==Fraction(math.factorial(n)*2**(n-1),3**(n-1))
        assert Fraction(total_c,universe)==Fraction(math.factorial(n),3**(n-1))
        out['ternary_checks'].append({'n':n,'universe':universe,'ordinary_sum':total_h,'decoded_sum':total_c})
    a=[[0,1,0,0],[1,0,0,1],[1,0,0,0],[0,0,0,0]]
    assert ham(a)==1 and ham(decode(a,'cancel'))==0
    out['cancellation_parity_hostile']={'ordinary_rows':[''.join(map(str,r)) for r in a],'ordinary_H':ham(a),'cancellation_H':ham(decode(a,'cancel')),'full_swap_H':ham(decode(a,'swap'))}
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
