#!/usr/bin/env python3
"""Fixed-core Gold strengthening of the independently supplied four-collar kernel."""
from itertools import combinations, product
from fractions import Fraction
from pathlib import Path
import json
from gold_terminal_chain import close,count,extensions,table,gold,deficits,require
CORE=(0,0,1,0,8,24,3,5,189)
PORTS=(0,6)
PAIR=(2,4)

def insert(A,u,B,a,b,C,v):
    na,nb=len(A),len(B);p=list(A)
    p += [(x<<na)|((1<<na)-1) for x in B]
    p += [(x<<(na+nb))|((1<<(na+nb))-1) for x in C]
    if A:p[na+a]&=~(1<<u)
    if C:p[na+nb+v]&=~(1<<(na+b))
    return tuple(p)

def flags(e,f):
    v=[[[0,0],[0,0]] for k in range(4)]
    for w in extensions(CORE):
        ff=(True,w[0]==PORTS[0],w[-1]==PORTS[1],w[0]==PORTS[0] and w[-1]==PORTS[1])
        i=0 if w.index(e[0])<w.index(e[1]) else 1
        j=0 if w.index(f[0])<w.index(f[1]) else 1
        for k,yes in enumerate(ff):v[k][i][j]+=int(yes)
    return v

def corner(v,r,s):return [[v[0][i][j]+r*v[1][i][j]+s*v[2][i][j]+r*s*v[3][i][j] for j in range(2)] for i in range(2)]

def hand_flags():
    v=[[[0,0],[0,0]] for _ in range(4)]
    for q in range(2,6):
        for t in range(1,q):
            multiplicity=6-q
            for j in range(8):
                i=0 if q<=3 else 1
                k=0 if j<q else 1
                w=8-t if j<t else 8-j
                v[0][i][k]+=multiplicity*w
                v[2][i][k]+=multiplicity
                if t==1 and j>=1:
                    v[1][i][k]+=multiplicity*w
                    v[3][i][k]+=multiplicity
    return v

def table_pred(p,e,f):
    result=[]
    for a in (True,False):
        row=[]
        for b in (True,False):
            pp=list(p)
            for x,y in (e if a else e[::-1], f if b else f[::-1]):pp[y]|=1<<x
            row.append(count(tuple(pp)))
        result.append(row)
    return result

def main():
    found=[]
    for f in combinations(range(9),2):
        v=flags(PAIR,f);cc=[corner(v,r,s) for r,s in product((0,1),repeat=2)]
        if all(gold(t) for t in cc):
            slack=min(Fraction(x,sum(sum(t,[]))) for t in cc for x in sum(deficits(t),[]))
            found.append((slack,f,v,cc))
    found.sort(reverse=True)
    require(found,'no common Gold certificate')
    best=found[0];slack,f,v,cc=best
    require(f==(1,2) and v==hand_flags(),'hand shuffle count')
    # Recount the exact four collar sources independently with ideal DP.
    verified=[]
    for r,s in product((0,1),repeat=2):
        p=insert((0,) if r else (),0,CORE,*PORTS,(0,) if s else (),0)
        rel=[(i,j) for j,mask in enumerate(p) for i in range(len(p)) if mask>>i&1]
        real=table(len(p),rel,tuple(x+r for x in PAIR),tuple(x+r for x in f))
        require(real==corner(v,r,s),'corner DP mismatch')
        verified.append(dict(r=r,s=s,table=real,deficits=deficits(real)))
    contexts=[(),(0,),(0,1),(0,0),(0,0,1),(0,0,1,3)]
    context_count=0
    for A,C in product(contexts,repeat=2):
        umax=[i for i in range(len(A)) if not any(ds>>i&1 for ds in A)] or [0]
        vmin=[i for i,ds in enumerate(C) if ds==0] or [0]
        for u,z in product(umax,vmin):
            ea,ec=count(A),count(C)
            r=Fraction(sum(w[-1]==u for w in extensions(A)),ea) if A else Fraction(0)
            s=Fraction(sum(w[0]==z for w in extensions(C)),ec) if C else Fraction(0)
            p=insert(A,u,CORE,*PORTS,C,z)
            got=table_pred(p,tuple(x+len(A) for x in PAIR),tuple(x+len(A) for x in f))
            expected=[[ea*ec*x for x in row] for row in corner(v,r,s)]
            require(got==expected and gold(got),'context table')
            require(min(sum(deficits(got),[]))*slack.denominator>=slack.numerator*sum(sum(got,[])),'uniform context slack')
            context_count+=1
    bridge=[];p=();low_slave_controls=[]
    for m in range(1,13):
        p=CORE if not p else insert(p,len(p)-9+PORTS[1],CORE,*PORTS,(),0)
        H=count(p)
        if m<=3:
            low=[]
            for x,y in combinations(range(len(p)),2):
                if p[x]>>y&1 or p[y]>>x&1:continue
                for a,b in ((x,y),(y,x)):
                    slaves=sum(((p[z]>>a&1) and not((p[z]>>b&1) or (p[b]>>z&1))) or ((p[b]>>z&1) and not((p[z]>>a&1) or (p[a]>>z&1))) for z in range(len(p)) if z not in (a,b))
                    if slaves<=1:
                        pp=list(p);pp[b]|=1<<a
                        if 2*count(tuple(pp))>=H:low.append((a,b,slaves))
            require(not low,'unexpected low-slave majority certificate')
            low_slave_controls.append(dict(blocks=m,majority_low_slave_pairs=low))
        margins=[]
        for j in range(m):
            got=table_pred(p,tuple(9*j+x for x in PAIR),tuple(9*j+x for x in f))
            require(sum(sum(got,[]))==H and gold(got),'bridge Gold')
            margin=Fraction(min(sum(deficits(got),[])),H)
            require(margin>=slack,'bridge uniform slack')
            margins.append(str(margin))
        bridge.append(dict(blocks=m,vertices=len(p),extensions=H,block_normalized_min_slacks=margins))
    def relation(i,j):return -1 if CORE[j]>>i&1 else (1 if CORE[i]>>j&1 else 0)
    modules=[]
    for mask in range(1,(1<<9)-1):
        S=[i for i in range(9) if mask>>i&1]
        if len(S)==1:continue
        if all(len({relation(x,y) for x in S})==1 for y in range(9) if y not in S):
            nonchain=any(relation(x,y)==0 for x,y in combinations(S,2))
            require(not nonchain,'nonchain proper module')
            modules.append(S)
    core_low=[]
    for a,b in product(range(9),repeat=2):
        if a==b or CORE[a]>>b&1 or CORE[b]>>a&1:continue
        k=sum(((CORE[z]>>a&1) and not((CORE[z]>>b&1) or (CORE[b]>>z&1))) or ((CORE[b]>>z&1) and not((CORE[z]>>a&1) or (CORE[a]>>z&1))) for z in range(9) if z not in (a,b))
        if k>1:continue
        V=flags((a,b),(a,b));nums=[]
        for r,s in product((0,1),repeat=2):
            t=corner(V,r,s);H=sum(sum(t,[]));num=sum(t[0]);nums.append(num)
            require(83*num<=35*H,'core low-slave corner bound')
            pp=insert((0,) if r else (),0,CORE,*PORTS,(0,) if s else (),0)
            real=table_pred(pp,(a+r,b+r),(a+r,b+r))
            require(real==t,'independent low-slave corner count')
        core_low.append(dict(pair=[a,b],slaves=k,corner_numerators=nums))
    require(len(core_low)==16,'low-slave list completeness')
    out={'core_low_slave_corner_certificate':core_low,'low_slave_controls':low_slave_controls,'proper_modules':modules,'hand_count_verified':True,'arbitrary_context_checks':context_count,'bridge_checks':bridge,'core':CORE,'ports':PORTS,'first':PAIR,'second':f,'flag_tables':v,'corners':verified,'uniform_normalized_slack':str(slack),'common_second_pairs':[list(t[1]) for t in found]}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
