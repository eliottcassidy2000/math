#!/usr/bin/env python3
"""Separate implementations by the same agent; no imports of other checkers."""
from fractions import Fraction
from itertools import product, permutations
from math import comb, gcd, lcm, factorial
from pathlib import Path
import json


def check(ok, msg):
    if not ok: raise RuntimeError(msg)


def brute_profile(rows):
    n=len(rows); bank=[set() for _ in range(n+1)]
    for perm in permutations(range(n)):
        for cut in range(1<<(n-1)):
            paths=[]; part=[perm[0]]; legal=True
            for j in range(n-1):
                if cut>>j&1:
                    paths.append(tuple(part));part=[perm[j+1]]
                elif rows[perm[j]]>>perm[j+1]&1:
                    part.append(perm[j+1])
                else:
                    legal=False;break
            if legal:
                paths.append(tuple(part));bank[len(paths)].add(tuple(sorted(paths)))
    return [len(b) for b in bank]


def cyclic_words(profiles):
    ordered=[[factorial(i)*c for i,c in enumerate(P)] for P in profiles]
    caps=[len(p)-1 for p in profiles]
    result=0
    for start in range(3):
        visits=[0]*3
        for length in range(1,sum(caps)+1):
            visits[(start+length-1)%3]+=1
            if any(v>c for v,c in zip(visits,caps)):break
            if min(visits)>0:
                term=1
                for p,v in zip(ordered,visits):term*=p[v]
                result+=term
    return result


def main():
    counts={}
    orbit_cases=0
    for m in range(1,51):
        factors=[];remaining=m;p=2
        while p<=remaining:
            e=0
            while remaining%p==0:remaining//=p;e+=1
            if e:factors.append((p,e))
            p+=1
        for multiplier in range(1,13):
            for source in range(m):
                seen={};y=source
                while y not in seen:
                    seen[y]=len(seen);y=multiplier*y%m
                actual=(seen[y],len(seen)-seen[y])
                tail=0;period=1
                for p,e in factors:
                    v=0;z=source
                    while v<e and z%p==0:v+=1;z//=p
                    if multiplier%p==0:
                        c=0;z=multiplier
                        while z%p==0:c+=1;z//=p
                        tail=max(tail,(e-v+c-1)//c)
                    else:
                        modulus=p**(e-v);order=1;z=multiplier%modulus
                        while z!=1%modulus:z=z*multiplier%modulus;order+=1
                        period=lcm(period,order)
                check(actual==(tail,period),('CRT orbit',m,multiplier,source,actual,tail,period))
                orbit_cases+=1
    counts['CRT_multiplication_orbits']=orbit_cases
    automata=0
    for B in range(2,9):
        for m in range(1,25):
            u=m
            while gcd(u,B)>1:u//=gcd(u,B)
            s=m//u;K=0
            while B**K<m or B**K%s:K+=1
            actual=[];theory=[]
            for r in range(m):
                bits=[]
                for k in range(K+1):
                    bits.extend((B**k*r+t)%m==0 for t in range(B**k))
                actual.append(tuple(bits))
                theory.append((r%u,tuple(((-B**k*r)%m if (-B**k*r)%m<B**k else None)
                                         for k in range(K))))
            for r in range(m):
                for t in range(m):
                    check((actual[r]==actual[t])==(theory[r]==theory[t]),('suffix',B,m,r,t))
            automata+=1
    counts['literal_suffix_automata']=automata

    cases=0
    for k in range(8):
        totals={M:0 for M in range(1,31)}
        for word in product((0,1),repeat=k):
            slope=Fraction(1);intercept=Fraction(0)
            for b in word:
                slope*=Fraction(3 if b else 1,2)
                intercept=((3 if b else 1)*intercept+b)/2
            residues=[]
            for n in range(2**k):
                y=slope*n+intercept
                if y.denominator==1:residues.append(n)
            check(len(residues)==1,('native uniqueness',word))
            r=residues[0]
            for M in totals:
                outputs={int(slope*(r+2**k*t)+intercept)%M for t in range(M)}
                totals[M]+=len(outputs)
                check(len(outputs)==M//gcd(M,3**sum(word)),('phase',k,word,M))
                cases+=1
        for M,value in totals.items():
            check(value==sum(comb(k,o)*M//gcd(M,3**o) for o in range(k+1)),('sum',k,M))
    counts['rational_word_phase_cases']=cases

    minimum_cases=0
    for M0 in (1,2,3,4,6,9):
        for o in range(5):
            expected=M0*3**o if M0%3==0 else M0
            for M in range(1,expected+1):
                seen={};valid=True
                for t in range(lcm(M,M0)):
                    y=(3**o*t)%M
                    if y in seen and seen[y]!=t%M0:
                        valid=False;break
                    seen[y]=t%M0
                if valid:
                    check(M==expected,('adaptive minimum',M0,o,M,expected))
                    break
            minimum_cases+=1
    counts['exhaustive_minimal_modulus_cases']=minimum_cases

    seeds=[(16,1,3,7,14),(20,1,2,7,14)]
    seed_profiles=[brute_profile(s) for s in seeds]
    check(seed_profiles[0]==seed_profiles[1]==[0,9,33,31,10,1],'seed profile')
    for perm in permutations(range(5)):
        check(not all(((seeds[0][i]>>j)&1)==((seeds[1][perm[i]]>>perm[j])&1)
                      for i in range(5) for j in range(5)),'seed isomorphism')
    queries=[]
    for s in seeds:
        rows=[v+(32 if i else 0) for i,v in enumerate(s)]+[1]
        queries.append(brute_profile(rows))
    cert=json.loads(Path(__file__).with_name('tournament_certificate.json').read_text())
    check(queries==cert['query_profiles'],'query profiles')
    for rec in cert['small_wrappers']:
        r=rec['r'];tt=tuple(sum(1<<j for j in range(i+1,r)) for i in range(r))
        pp=brute_profile(tt)
        check(cyclic_words([seed_profiles[0],seed_profiles[1],pp])==rec['H'],('wrapper',r))
        for bit in (0,1):
            check(cyclic_words([queries[bit],seed_profiles[0],pp])==rec['reads'][bit],('read',r,bit))
    def literal_H(rows):
        N=len(rows)
        def walk(mask,last):
            if mask==(1<<N)-1:return 1
            choices=rows[last]&~mask;result=0
            while choices:
                bit=choices&-choices;choices-=bit
                result+=walk(mask|bit,bit.bit_length()-1)
            return result
        return sum(walk(1<<v,v) for v in range(N))
    for first in cert['coupled_one_vertex_reads'][:4]:
        check(literal_H(first['rows'])==first['H'],'nonmemoized coupled read audit')
    seven=[row|32 for row in seeds[0]]+[64,31]
    check(literal_H(seven)==93,'sharpness component')
    check(all(30431%p for p in range(2,175)),'trial division prime')
    counts['brute_permutation_profiles']=8
    counts['cyclic_run_word_H_reproductions']=12
    counts['nonmemoized_coupled_graphs']=4
    counts['sharpness_component_paths']=93
    print(json.dumps({'status':'PASS: separate same-agent implementation',
                      'counts':counts,
                      'seed_profiles':seed_profiles,'query_profiles':queries},indent=2,sort_keys=True))


if __name__=='__main__':main()
