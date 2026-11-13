#!/usr/bin/env python3
"""Independent algebra/template audit; does not import producer matrix code."""
from hashlib import sha256
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def need(q, message):
    if not q:
        raise RuntimeError(message)

def skew(n, entries):
    out = [[0]*n for _ in range(n)]
    for (i,j), x in zip(combinations(range(n),2), entries):
        out[i][j], out[j][i] = x, -x
    return out

def determinant(A):
    n = len(A)
    total = 0
    for p in permutations(range(n)):
        term = (-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        for i,j in enumerate(p):
            term *= A[i][j]
        total += term
    return total

def identity(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]

def multiply(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]

def transpose(A):
    return [list(row) for row in zip(*A)]

def mutate(A, k):
    n = len(A)
    return [[-A[i][j] if i==k or j==k else
             A[i][j]+max(A[i][k],0)*max(A[k][j],0)
             -max(-A[i][k],0)*max(-A[k][j],0)
             for j in range(n)] for i in range(n)]

def ordinal(*blocks):
    sizes = [len(b) for b in blocks]
    n = sum(sizes)
    result = [[0]*n for _ in range(n)]
    offset = 0
    for block in blocks:
        size = len(block)
        for i in range(size):
            for j in range(size):
                result[offset+i][offset+j] = block[i][j]
            for j in range(offset+size, n):
                result[offset+i][j], result[j][offset+i] = 1,-1
        offset += size
    return result

def from_rows(rows):
    n=len(rows)
    need(all(not (rows[i] >> i & 1) for i in range(n)), 'no loops')
    need(all((rows[i] >> j & 1)+(rows[j] >> i & 1)==1
             for i,j in combinations(range(n),2)), 'complete orientations')
    return [[0 if i==j else 1 if rows[i] >> j & 1 else -1
             for j in range(n)] for i in range(n)]

def birth(rows, root):
    n=len(rows)
    return tuple(row | ((1 << n) if i==root else 0)
                 for i,row in enumerate(rows)) + (((1 << n)-1) ^ (1 << root),)

def rooted_literal(rows):
    n=len(rows);full=(1 << n)-1
    @lru_cache(None)
    def paths(mask):
        labels=tuple(i for i in range(n) if mask >> i & 1)
        return tuple(p for p in permutations(labels)
                     if all(rows[p[j]] >> p[j+1] & 1 for j in range(len(p)-1)))
    H=len(paths(full))
    roots=[r for r in range(n) if sum(rows[i] >> r & 1 for i in range(n))==1]
    result=[]
    for root in roots:
        q=len(paths(full ^ (1 << root)))
        s=sum(p[0]==root for p in paths(full))
        F=sum(sum(p[-1]==root for p in paths(mask))*len(paths(full ^ mask))
              for mask in range(1,full+1) if mask >> root & 1)
        need(F-H==q-s,'independent endpoint-cover/deletion defect agreement')
        result.append({'root':root,'H':H,'deletion':q,'start':s,'end_two_cover':F,'defect':F-H})
    return result

def main():
    mutation_checks = 0
    for entries in product(range(-2,3), repeat=3):
        A=skew(3,entries)
        for k in range(3):
            E=identity(3)
            E[k][k]=-1
            for i in range(3):
                if i != k:
                    E[i][k]=max(A[i][k],0)
            need(determinant(E)==-1 and multiply(E,E)==identity(3), 'integer involution')
            need(mutate(A,k)==multiply(multiply(E,A),transpose(E)), 'FZ/congruence equality')
            mutation_checks += 1
    C=skew(3,(1,-1,1))
    basis_checks = 0
    for n in range(4):
        for entries in product(range(-2,3), repeat=n*(n-1)//2):
            R=skew(n,entries); A=ordinal(C,R); size=n+3
            vectors=[]
            vectors.append([int(i==0) for i in range(size)])
            vectors.append([int(i==1) for i in range(size)])
            vectors.append([int(i<3) for i in range(size)])
            for r in range(n):
                v=[0]*size;v[0]=1;v[1]=-1;v[r+3]=1
                vectors.append(v)
            P=transpose(vectors)
            need(determinant(P)==1,'triangle basis unimodular')
            paired=[[sum(v[i]*A[i][j]*w[j] for i in range(size) for j in range(size))
                     for w in vectors] for v in vectors]
            expected=[[0]*size for _ in range(size)]
            expected[0][1]=1;expected[1][0]=-1
            for i in range(n):
                expected[2][i+3]=3;expected[i+3][2]=-3
                for j in range(n):
                    expected[i+3][j+3]=R[i][j]
            need(paired==expected,'complete triangle pairing')
            need(determinant(A)==9*determinant(ordinal([[0]],R)),'factor-nine determinant')
            basis_checks += 1
    templates=[]
    for a in range(4):
        b=3-a
        for left in product((-1,1),repeat=a*(a-1)//2):
            for right in product((-1,1),repeat=b*(b-1)//2):
                for orientation in (-1,1):
                    triangle=[[orientation*x for x in row] for row in C]
                    A=ordinal(skew(a,left),triangle,skew(b,right))
                    templates.append({'left_size':a,'det':determinant(A)})
    need(len(templates)==40 and {t['det'] for t in templates}=={9,81}, 'all triangle templates')
    controls=[]
    for rows in ((52,33,34,39,46,0),(58,60,57,16,32,8)):
        A=from_rows(rows)
        H=sum(all(A[p[i]][p[i+1]]==1 for i in range(5))
              for p in permutations(range(6)))
        roots=[j for j in range(6) if sum(A[i][j]==1 for i in range(6))==1]
        controls.append({'rows':rows,'H':H,'det':determinant(A),'indegree_one':roots})
    need([(r['H'],r['det'],r['indegree_one']) for r in controls]
         ==[(9,25,[3,4]),(9,81,[0,1,2])], 'actual H9 separation')
    G=(10,12,1,4)
    first=birth(G,1)
    prime_A=birth(first,4)
    prime_B=(58,60,1,52,36,4)
    need(first==(10,28,1,4,13) and prime_A==(10,28,1,4,45,15),'literal two births')
    G_profiles=rooted_literal(G)
    need(G_profiles[1]=={'root':1,'H':5,'deletion':3,'start':2,
                        'end_two_cover':6,'defect':1},'H5 ready predecessor')
    prime_controls=[]
    for rows in (prime_A,prime_B):
        prime_controls.append({'rows':rows,'det':determinant(from_rows(rows)),
                               'root_profiles':rooted_literal(rows)})
    need([r['det'] for r in prime_controls]==[49,1],'prime determinant separation')
    need([[(p['root'],p['H'],p['deletion'],p['start'],p['end_two_cover'],p['defect'])
           for p in r['root_profiles']] for r in prime_controls]
         ==[[(4,17,5,5,17,0),(5,17,9,6,20,3)],
            [(0,17,1,1,17,0),(1,17,9,8,18,1)]],'full prime rooted counts')
    A=skew(3,(1,1,1));B=skew(3,(-2,1,1))
    need(all((A[i][j]-B[i][j])%3==0 for i in range(3) for j in range(3)), 'same residues')
    need((mutate(A,1)[0][2],mutate(B,1)[0][2])==(2,1),'residue mutation hostile')
    sources=['05-knowledge/results/mutation_collatz_chanpak_20261112/matrix/PROOF.md',
             '05-knowledge/results/mutation_collatz_chanpak_20261112/matrix/check.py']
    print(json.dumps({
        'status':'PASS; independent algebra and complete small ordinal templates',
        'mutation_checks':mutation_checks,'basis_checks':basis_checks,
        'ordinal_templates':len(templates),'ordinal_determinants':[9,81],
        'literal_H9_controls':controls,
        'prime_H17_controls':prime_controls,'H5_predecessor_profiles':G_profiles,
        'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
    },indent=2,sort_keys=True))

if __name__=='__main__':
    main()
