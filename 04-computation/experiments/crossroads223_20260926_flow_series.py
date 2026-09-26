"""Exact fringe series for the P2 prefix optimum. Standard library only."""
from array import array
from fractions import Fraction
from functools import lru_cache
import json


def require(test,message):
    if not test:raise AssertionError(message)


@lru_cache(None)
def direct_tree(h,i):
    if h==0:return 0,1,0
    if i%2==0:
        x,y,_=direct_tree(h-1,3*i//2)
        zero,one=y,1+x
        return zero,one,min(zero,one)-min(x,y)
    a,b,_=direct_tree(h-1,(3*i-1)//2)
    x,y,_=direct_tree(h-1,(3*i+1)//2)
    zero,one=a+min(x,y),1+min(a,b)+y
    return zero,one,min(zero,one)-min(a,b)-min(x,y)


def main():
    differences=array('i',[1]);partial=Fraction(0);direct_checks=0
    for h in range(1,21):
        modulus=len(differences);new=array('i');total=0
        for i in range(2*modulus):
            if i%2==0:
                x=differences[(3*i//2)%modulus]
                delta,toll=1-x,int(x>=1)
            else:
                x=differences[((3*i-1)//2)%modulus]
                y=differences[((3*i+1)//2)%modulus]
                delta=1+min(x,0)+max(y,0)
                toll=min(max(-x,0),1+max(y,0))
            require(abs(delta)<=h+1 and 0<=toll<=h,'universal depth bound')
            if h<=7:
                a,b,t=direct_tree(h,i)
                require((b-a,t)==(delta,toll),'direct full-tree audit')
                direct_checks+=1
            total+=toll;new.append(delta)
        differences=new
        require(sum(differences)==len(differences),'mean difference identity')
        partial+=Fraction(total,3**(h+1))
        print(json.dumps({'depth':h,'residues':2**h,'exact_toll_sum':total,'partial_decimal':float(partial)}))
    tail=Fraction(23)*Fraction(2,3)**21
    print(json.dumps({'independent_full_tree_cases':direct_checks,'partial_exact':str(partial),
                      'tail_bound_exact':str(tail),'limit_interval_exact':[str(partial),str(partial+tail)],
                      'limit_interval_decimal':[float(partial),float(partial+tail)]}))
    print('ALL CHECKS PASSED')


if __name__=='__main__':main()
