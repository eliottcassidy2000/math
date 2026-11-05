#!/usr/bin/env python3
"""Independent literal path/cut audit of the typed gluing output.
No subset-cover DP is used for the independent expected profile.
"""
from collections import Counter
from itertools import permutations,product
import hashlib,json,pathlib
import verify as v
HERE=pathlib.Path(__file__).resolve().parent

def literal(a,colors,k):
    n=len(a);out=Counter()
    if not n:return {(0,)*(k*k):1}
    for word in permutations(range(n)):
        for cuts in range(1<<(n-1)):
            paths=[];start=0;legal=True
            for i in range(n-1):
                if cuts>>i&1:
                    paths.append(word[start:i+1]);start=i+1
                elif not(a[word[i]]>>word[i+1]&1):legal=False;break
            if not legal:continue
            paths.append(word[start:])
            # Lexicographic list order picks precisely one ordering of path components.
            if paths!=sorted(paths):continue
            sig=[0]*(k*k)
            for path in paths:sig[colors[path[0]]*k+colors[path[-1]]]+=1
            out[tuple(sig)]+=1
    return dict(out)

def main():
    counts=[]
    for n in range(1,5):
        count=0
        for mask in range(1<<(n*(n-1)//2)):
            a=v.tournament(n,mask)
            for colors in product(range(2),repeat=n):
                v.require(literal(a,colors,2)==v.path_profile(a,colors,2),'literal profile mismatch')
                count+=1
        counts.append({'n':n,'profiles':count})
    selected=[]
    for word in ('120','1201','12','2'):
        a=v.two_apex(word);colors=(0,1)+(2,)*len(word)
        lp=literal(a,colors,3)
        v.require(lp==v.path_profile(a,colors,3),'selected literal mismatch')
        selected.append({'word':word,'H':v.profile_h(lp),'profile_terms':len(lp)})
    # Literal verification of the kernel on indistinguishable-type fragments:
    # collect canonical path covers on three named fragments with endpoint tags.
    # Type construction: A fragments (0,1),(1,0), B fragment (0,0).
    # C all A->B, hence no cycle; choose root reversal separately as hostile.
    for C in product((0,1),repeat=4):
        ma=(0,1,1,0);mb=(1,0,0,0)
        types=[(0,0,1),(0,1,0),(1,0,0)];a=[0]*3
        for i,(side,s,t) in enumerate(types):
            for j,(other,u,w) in enumerate(types):
                if side!=other and (C[t*2+u] if side==0 else not C[u*2+t]):a[i]|=1<<j
        out=Counter()
        for word in permutations(range(3)):
            for cuts in range(4):
                paths=[];start=0;legal=True
                for i in range(2):
                    if cuts>>i&1:paths.append(word[start:i+1]);start=i+1
                    elif not(a[word[i]]>>word[i+1]&1):legal=False;break
                if not legal:continue
                paths.append(word[start:])
                if paths!=sorted(paths):continue
                sig=[0]*4
                for p in paths:sig[types[p[0]][1]*2+types[p[-1]][2]]+=1
                out[tuple(sig)]+=1
        v.require(dict(out)==v.kernel(ma,mb,2,C),'kernel literal mismatch')
    # Degree-only hostile: the triangle has one directed cycle cover but zero empty-path term.
    tri=v.alternating(1);p=literal(tri,(0,0,0),1)
    v.require(p=={(1,):3,(2,):3,(3,):1},'triangle cycle-safe profile hostile')
    receipt={'status':'PASS','all_binary_colored_tournaments':counts,'selected_two_apex':selected,'literal_fragment_kernels':16,
      'cycle_hostile':{'rows':list(tri),'pathcover_coefficients':{'1':3,'2':3,'3':1},'discarded_directed_cycle_covers':1},
      'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
      'implementation_sha256':hashlib.sha256((HERE/'verify.py').read_bytes()).hexdigest()}
    (HERE/'literal_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
