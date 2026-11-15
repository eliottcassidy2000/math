"""Independent finite controls for the exact continuation-horizon theorem."""
from itertools import combinations, permutations, product
from math import factorial
from collections import Counter
import json, random

def ordered_adjacency(n, code):
    A = [[False]*n for _ in range(n)]
    for bit,(a,b) in enumerate(combinations(range(n),2)):
        A[a][b] = bool(code >> bit & 1)
        A[b][a] = not A[a][b]
    return A

def query(k, matching):
    # vertex 0 is the reader hole
    n=2*k+1
    A=[[False]*n for _ in range(n)]
    for a,b in combinations(range(n),2):
        yes= (b<=k) if a==0 else (not matching or b==a+k) if a<=k<b else True
        A[a][b]=yes
        A[b][a]=not yes
    return A

def coeffs(A):
    q=len(A)-1
    vals=[0]*(q//2+2)
    for word in permutations(range(1,q+1)):
        for k in range(1,len(vals)):
            for slots in combinations(range(q+1),k):
                P=[]; S=set(slots)
                for i in range(q+1):
                    if i in S: P.append(0)
                    if i<q: P.append(word[i])
                if all(A[a][b] for a,b in zip(P,P[1:])): vals[k]+=1
    return vals

def permanent_bound(A):
    n=len(A); m=(n-1)//2
    out=[i for i in range(1,n) if A[0][i]]
    inc=[i for i in range(1,n) if A[i][0]]
    if len(out)!=m: return 0
    return factorial(m)*sum(all(A[a][b] for a,b in zip(out,p)) for p in permutations(inc))

def paths(A):
    return sum(all(A[a][b] for a,b in zip(w,w[1:])) for w in permutations(range(len(A))))

def insert(A, D):
    # first |D| vertices replace hole; every pair compiled literally
    d=len(D); n=len(A)+d-1
    B = [[False]*n for _ in range(n)]
    for a,b in combinations(range(n),2):
        z = D[a][b] if b<d else A[max(0,a-d+1)][b-d+1]
        B[a][b]=z; B[b][a]=not z
    return B

def covers(A, k):
    n=len(A)
    return sum(all(i in c or A[w[i]][w[i+1]] for i in range(n-1))
               for w in permutations(range(n)) for c in combinations(range(n-1),k-1))

def ensure(guard, *details):
    if not guard: raise AssertionError(details)

def run():
    checked=Counter()
    seed=random.Random(4856)
    for q in range(5):
        n=q+1
        for code in range(1<<(n*(n-1)//2)):
            A=ordered_adjacency(n,code); vals=coeffs(A)
            ensure(all(v%factorial(i-1)==0 for i,v in enumerate(vals) if i),q,code,vals)
            checked['factorial_contexts']+=1
            if q%2==0:
                ensure(vals[-1]==permanent_bound(A), q,code, vals)
                checked['permanent_contexts']+=1
            if seed.randrange(40)==0:
                d=3
                D=ordered_adjacency(d,seed.randrange(8))
                ensure(paths(insert(A,D))==sum(v*covers(D,i) for i,v in enumerate(vals) if i),q,code)
                checked['literal_substitutions']+=1
    for m in range(4):
        for match in (True,False):
            A=query(m,match); vals=coeffs(A)
            ensure(vals[-1]==factorial(m)**(1 if match else 2),m, match, vals)
            checked['separating_queries']+=1
    for modulus in (2,3,4,6,8,9,12,24,36,72):
        for codes in product(range(8),repeat=2):
            X,Y=(ordered_adjacency(3,k) for k in codes)
            actual=all((paths(insert(query(m,True),X))-paths(insert(query(m,True),Y)))%modulus==0 for m in range(3))
            expected=all((factorial(k-1)*(covers(X,k)-covers(Y,k)))%modulus==0 for k in range(1,4))
            ensure(actual==expected,modulus,codes)
            checked['modular_comparisons']+=1
    for d in range(1,5):
        for code in range(1<<(d*(d-1)//2)):
            D=ordered_adjacency(d,code)
            E=insert(query(1,True),D)
            W={k:factorial(k-1)*covers(D,k) for k in range(1,d+1)}
            for j in range(1,d+3):
                expected = (W.get(j+1,0) + j*(j+2)*W.get(j,0)
                      + j*(j-1)*(2*j-1)*W.get(j-1,0)
                      + j*(j-1)**2*(j-2)*W.get(j-2,0))
                ensure(factorial(j-1)*covers(E,j)==expected,d,code,j)
                checked['band_formula_grades']+=1
    return checked

if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True))
