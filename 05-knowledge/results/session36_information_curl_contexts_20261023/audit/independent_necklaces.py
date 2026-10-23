#!/usr/bin/env python3
"""Labelled isomorphism orbits then cyclic SCC orbits, without repo imports."""
from independent_audit import mat,codeof,curl,components,rotate,H,need
import math,itertools,json
from fractions import Fraction
from pathlib import Path
OUT=Path(__file__).resolve().parent

def relabel(z,p):return [[z[i][j] for j in p] for i in p]

def phi(n):return sum(math.gcd(i,n)==1 for i in range(1,n+1))

def logordinary(a):
    # (log A)' A = A', exact rational formal coefficients.
    b=[Fraction(0)]*len(a)
    for n in range(1,len(a)):
        b[n]=a[n]-sum(Fraction(j,n)*b[j]*a[n-j] for j in range(1,n))
    return b

def main():
    rows=[];a=[1];strong=[0];c=[0];weighted=[1];stabilizers=0;Htypes=[[]];Hcurl=[[]]
    for n in range(1,7):
        perms=list(itertools.permutations(range(n)));remaining=set(range(1<<(n*(n-1)//2)));representative={};reps=[]
        while remaining:
            b=min(remaining);z=mat(n,b);orb={codeof(relabel(z,p)) for p in perms}
            for x in orb:representative[x]=b
            remaining-=orb;reps.append(b)
        states={};numstrong=0
        for b in reps:
            z=mat(n,b);blocks=components(z);numstrong+=len(blocks)==1
            state=min(representative[codeof(rotate(z,blocks,k))] for k in range(len(blocks)))
            if state in states:need(states[state]==H(z),('unlabelled H',n,b))
            states[state]=H(z)
            word=[];autprod=1
            for block in blocks:
                sub=[[z[i][j] for j in block] for i in block]
                orbit=[codeof(relabel(sub,p)) for p in itertools.permutations(range(len(block)))]
                word.append((len(block),min(orbit)));autprod*=sum(x==codeof(sub) for x in orbit)
            rotations=sum(word==word[j:]+word[:j] for j in range(len(word)))
            actual=sum(curl(relabel(z,p))==curl(z) for p in perms)
            need(actual==rotations*autprod,('curl stabilizer',n,b,actual,rotations,autprod));stabilizers+=1
        Htypes.append([H(mat(n,b)) for b in reps]);Hcurl.append(list(states.values()))
        a.append(len(reps));strong.append(numstrong);c.append(len(states));weighted.append(sum(Htypes[-1]))
        row=dict(n=n,tournament_isotypes=len(reps),strong_isotypes=numstrong,curl_metric_isotypes=len(states),sum_H_over_curl_isotypes=sum(states.values()))
        rows.append(row);print(row,flush=True)
    loga=logordinary(a)
    expected=[0]+[sum(Fraction(phi(d),d)*loga[n//d] for d in range(1,n+1) if n%d==0) for n in range(1,7)]
    need(expected==c,('unlabelled necklace series',expected,c))
    # Confirm T=SEQ(S) at isomorphism-type level as well.
    for n in range(1,7):need(a[n]==sum(strong[k]*a[n-k] for k in range(1,n+1)),('ordinary SCC sequence',n))
    weighted_checks=[];naive_failures=[]
    for power in (0,1,2):
        logs={d:logordinary([1]+[sum(h**(d*power) for h in Htypes[n]) for n in range(1,7)]) for d in range(1,7)}
        for n in range(1,7):
            predicted=sum(Fraction(phi(d),d)*logs[d][n//d] for d in range(1,n+1) if n%d==0)
            actual=sum(h**power for h in Hcurl[n])
            need(predicted==actual,('Adams weighted transform',n,power))
            weighted_checks.append(dict(n=n,power=power,value=actual))
            naive=sum(Fraction(phi(d),d)*logs[1][n//d] for d in range(1,n+1) if n%d==0)
            if naive!=actual:naive_failures.append(dict(n=n,power=power,actual=actual,naive=str(naive)))
    report=dict(status='PASS',rows=rows,unweighted_necklace_coefficients=c[1:],stabilizers_checked=stabilizers,Adams_checks=weighted_checks,naive_failures=naive_failures,formula='sum_{d>=1} phi(d)/d log A(z^d), A ordinary unlabelled tournament series')
    (OUT/'necklace_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
