#!/usr/bin/env python3
"""Independent permutation/cut and subset-convolution audit of check.py."""
import itertools as it
import json
import math
from fractions import Fraction
from pathlib import Path
from check import need, marked_covers, coupled, compressed, slot, ordinal, profile, add_cut, cut_polynomial


def permutation_covers(rows):
    n=len(rows)
    C=[0]*(n+1);S=[0]*4;E=[0]*4
    for perm in it.permutations(range(n)):
        bad=sum(1<<i for i in range(n-1) if not (rows[perm[i]]>>perm[i+1]&1))
        for cuts in range(1<<(n-1)):
            if cuts & bad != bad:
                continue
            j=cuts.bit_count()+1
            C[j]+=1
            if j<=3:
                S[j]+=perm[0]==0
                firstend=(cuts & -cuts).bit_length()-1 if cuts else n-1
                E[j]+=perm[firstend]==0
    return C,S,E


def raw_h(rows):
    # Suffix-start recursion, independent of endpoint dynamic programming.
    from functools import lru_cache
    n=len(rows)
    @lru_cache(None)
    def walk(mask,u):
        if mask==0:
            return 1
        options=rows[u]&mask
        total=0
        while options:
            bit=options & -options
            total+=walk(mask^bit,bit.bit_length()-1)
            options^=bit
        return total
    return sum(walk(((1<<n)-1)^(1<<u),u) for u in range(n))


def polynomial(C):
    ans=[Fraction(0)]*len(C)
    basis=[Fraction(1)]
    for j,c in enumerate(C):
        for k,b in enumerate(basis):
            ans[k]+=c*b
        nextbasis=[Fraction(0)]*(len(basis)+1)
        for k,b in enumerate(basis):
            nextbasis[k]-=j*b/(j+1)
            nextbasis[k+1]+=b/(j+1)
        basis=nextbasis
    return ans


def main():
    count=0; Hamilton={}; profile_sets={}; samples={}
    for n in range(1,6):
        pairs=list(it.combinations(range(n),2));Hamilton[n]=set();profile_sets[n]=set()
        for bits in range(1<<len(pairs)):
            rows=[0]*n
            for k,(i,j) in enumerate(pairs):
                u,v=(i,j) if bits>>k&1 else (j,i)
                rows[u]|=1<<v
            C,S,E=permutation_covers(rows)
            observed=marked_covers(tuple(rows))
            padded=(C+[0]*4)[:4]
            need(observed==(padded,S,E),('marked permutation check',rows,observed,(padded,S,E)))
            P=profile(rows)
            need(C==[math.factorial(j)*x for j,x in enumerate(P)],'all covers independent')
            power=polynomial(C)
            need(all(x.denominator==1 and x>=0 and (not x or k%2==n%2) for k,x in enumerate(power)), 'parity and integrality')
            c=math.gcd(*C)
            vals=[sum(C[j]*math.comb(r,j) for j in range(min(r,n)+1)) for r in range(n+1)]
            need(c==math.gcd(*vals) and math.factorial(n)%c==0,'fixed content')
            Hamilton[n].add(C[1]);profile_sets[n].add((C[1],vals[2] if n>=2 else 2,c))
            count+=1
    A=(8,1,3,6);B=(24,1,3,6,14)
    need(polynomial(permutation_covers(A)[0])==[0,0,4,0,1], 'A polynomial')
    need(polynomial(permutation_covers(B)[0])==[0,2,0,8,0,1], 'B polynomial')
    for a,b in [(A,B),(B,A),(A,A),(B,B)]:
        rows=coupled([a,b]);actual=raw_h(rows)
        need(actual==compressed([a,b])[0], 'suffix recursion reader audit')
        samples[str((len(a),len(b)))]=actual
    R=json.loads((Path(__file__).parent/'receipt.json').read_text())
    need(R['census_count']==count==1099,'receipt census')
    need({str(n):sorted(profile_sets[n]) for n in profile_sets}==
         {n:sorted(tuple(v) for v in values) for n,values in R['small_profiles_H_F2_content'].items()},'receipt profiles')
    parent=ordinal(A,B)
    *_,ev=cut_polynomial(parent)
    cut_values=[]
    for mask in range(512):
        observed=raw_h(add_cut(parent,mask))
        need(observed==ev(mask),('all H55 cuts',mask))
        cut_values.append(observed)
    need(math.gcd(*cut_values)==1 and len(set(cut_values))==354 and all(v not in cut_values for v in (111,167,169)), 'primitive but noncovering')
    # Modulus is a literal value range, not a probabilistic distribution.
    height_controls=0
    for bound in range(0,13):
        M=math.factorial(7)
        for n in range(-bound,bound+1):
            candidates=[a for a in range(-bound,bound+1) if a%M==n%M]
            need(candidates==[n], 'height-selected source')
            height_controls+=1
    print(json.dumps({'status':'independent exact algorithms; same agent review',
      'permutation_cut_tournaments':count,'H_through5':{n:sorted(v) for n,v in Hamilton.items()},
      'suffix_recursion_reader_counts':samples,'height_unique_controls':height_controls,
      'H55_all_cuts_independent':512,'H55_response_gcd':math.gcd(*cut_values),
      'all_checks':'PASS'},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
