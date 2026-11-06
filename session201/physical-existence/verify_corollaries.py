"""Exact follow-up corollary checks; standard library and explicit exceptions."""
import hashlib,json,math
from fractions import Fraction
from itertools import product
from pathlib import Path
from verify import T,B,K,step,state,dagger,F,graph,path_cycle_counts,need

def main():
    rays=[(1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1)]
    C=[[K(step(v,'J'),w)for w in rays]for v in rays]
    A=[step(v,'J')[2]for v in rays]
    need(C==[[2,1,4,3],[1,1,2,2],[4,2,7,6],[3,2,6,5]],'J pairing matrix')
    R=[[16*C[i][j]-7*A[i]*A[j]for j in range(4)]for i in range(4)]
    need(R==[[4,2,8,6],[2,9,4,11],[8,4,0,12],[6,11,12,17]],'J residual matrix')
    need(all(x>=0 for row in R for x in row),'J residual positivity')
    words=rectangles=nonintegral=Jcross=0
    for seed in [T,B]:
        for le in range(7):
            for letters in product('OZJ',repeat=le):
                word=''.join(letters);left=[T];right=[seed]
                for c in word:left.append(step(left[-1],c))
                for c in word[::-1]:right.append(step(right[-1],dagger(c)))
                right=right[::-1];n=K(T,right[0]);words+=1
                for denom in [1,2,3]:
                    for numerator in range(2*denom,12*denom+1):
                        L=Fraction(numerator,denom);R=max(Fraction(seed[2]),3*(n-1)/L)
                        need(L>=2 and R>=seed[2]and L*R>=3*(n-1),'valid rational rectangle budgets')
                        i=max(i for i,s in enumerate(left)if s[2]<=L)
                        need(left[i][2]<=L and right[i][2]<=R,'last admissible-left cut meets rectangle')
                        rectangles+=1;nonintegral+=L.denominator>1 or R.denominator>1
                cross=next((i for i,(v,w)in enumerate(zip(left,right))if v[2]>=w[2]),None)
                if cross not in [None,0]and word[cross-1]=='J':
                    H=min(max(left[cross-1][2],right[cross-1][2]),max(left[cross][2],right[cross][2]))
                    need(H<=math.isqrt((16*n)//7),'J-specific crossing cutoff')
                    Jcross+=1
    sharp=[]
    for k in range(1,21):
        d=2**k+1;v=state('Z'*k+'J');n=K(step(v,'J'),v)
        need(v==(2*d-1,d-1,2*d,d),'J sharp physical source')
        need(n==7*d*d-2*d and step(v,'J')[2]==4*d-1,'J sharp formulas')
        w='Z'*k+'JJJ'+'O'*k
        need(F(state(w))==n,'J sharp literal source word')
        heights=[max(state(w[:i])[2],state(dagger(w[i:]))[2])for i in range(len(w)+1)]
        need(min(heights)==4*d-1,'J sharp best cut over whole word')
        sharp.append(dict(k=k,target=n,best_cut_height=4*d-1))
    literal_word='ZJJJOJ';rows=graph(literal_word)
    need(path_cycle_counts(rows)==(57,1),'J sharp literal graph at k1')
    return dict(status='PASS',complete_word_length_cap=6,words=words,rectangle_checks=rectangles,
                noninteger_rectangle_checks=nonintegral,J_crossings=Jcross,J_matrix=C,J_sharp_family=sharp,
                literal_graph=dict(word=literal_word,vertices=len(rows),H=57,C=1),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Exact finite checks support the all-order proved corollaries; no universal coverage or global optimality claim.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
