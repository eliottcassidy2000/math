from collections import Counter
from itertools import combinations
from math import isqrt
from fractions import Fraction as F

def sq(n):
    A=[set() for _ in range(n+1)]
    labels=Counter()
    for i,j in combinations(range(1,n+1),2):
        x=i+j
        if isqrt(x)**2==x:
            labels[x]+=1; A[i].add(j); A[j].add(i)
    # Count four cycles from opposite pairs.
    c4=sum(len(A[i]&A[j])*(len(A[i]&A[j])-1)//2
           for i,j in combinations(range(1,n+1),2))//2
    c3=sum(len(A[i]&A[j]) for i,j in combinations(range(1,n+1),2) if j in A[i])//3
    # Counts of all traversals, with one cycle anchor and reversal quotient.
    def dfs(root,last,rem,want_cycle):
        if not rem:return int(not want_cycle or root in A[last])
        if any(not (A[x]&(rem|{last,root} if want_cycle else rem|{last})) for x in rem):
            return 0
        return sum(dfs(root,x,rem-{x},want_cycle) for x in A[last]&rem)
    V=set(range(1,n+1))
    h=sum(dfs(i,i,V-{i},False) for i in V)//2
    c=dfs(1,1,V-{1},True)//2
    return [sum(labels.values()),sum(v for k,v in labels.items() if k%2),
            sum(v for k,v in labels.items() if not k%2),c3,c4,h,c]

def algebra():
    assert 125*25-25**2==2500 and 25*(F(225,4))<2500
    assert 25+2*32*25==1625 and 1+32+32*31+25*24==1625
    assert 25*(5*7-25)>4*F(225,4)
    for c,rank_max in [(30,117),(120,105),(210,93)]:
        assert 121-(2*c)//15==rank_max
    assert 2*120+100>2*25*6
    assert [(i,160-i) for i in range(20,151,30) if 20<=160-i<=150]==[
        (20,140),(50,110),(80,80),(110,50),(140,20)]
    Q={i*i%25 for i in range(25)}
    assert Q=={0,1,4,6,9,11,14,16,19,21,24}

if __name__=='__main__':
    algebra()
    result={n:sq(n) for n in (25,32)}
    print('n: [edges, cross parity, same parity, triangles, squares, paths, cycles]')
    assert result=={25:[32,17,15,0,0,10,0],32:[46,24,22,1,0,392,1]}
    print(result)
