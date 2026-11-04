#!/usr/bin/env python3
"""Independent exact formal-profile blindness audit; no tournament realization."""
from pathlib import Path
from fractions import Fraction
from itertools import product
from math import comb,lcm
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'coarse_positive_obstruction'

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def kernel(a,c,reverse):
    rows=[0]*5
    for i in range(3):
        rows[i]|=1<<((i+(2 if reverse else 1))%3)
        for mask,p in [(a,3),(c,4)]:
            if mask>>i&1:rows[i]|=1<<p
            else:rows[p]|=1<<i
    rows[3]|=1<<4
    # State is (old multiplicities, used-port mask, last label).
    levels=[{} for _ in range(12)]
    for v in range(5):
        e=[0,0,0];ports=0
        if v<3:e[v]=1
        else:ports=1<<(v-3)
        levels[1][tuple(e)+(ports,v)]=1
    series={}
    for length in range(1,12):
        for state,z in levels[length].items():
            e=state[:3];ports,last=state[3:]
            if ports==3:series[e]=series.get(e,0)+z
            if length==11:continue
            for v in range(5):
                if not rows[last]>>v&1:continue
                nxt=list(e);pp=ports
                if v<3:
                    if nxt[v]==3:continue
                    nxt[v]+=1
                else:
                    bit=1<<(v-3)
                    if pp&bit:continue
                    pp|=bit
                key=tuple(nxt)+(pp,v);levels[length+1][key]=levels[length+1].get(key,0)+z
    return series

def numerator(a,c,reverse):
    series=kernel(a,c,reverse);out={}
    for e in product(range(4),repeat=3):
        value=0
        for k,m in enumerate((1,-3,3,-1)):
            if min(e)>=k:value+=m*series.get(tuple(i-k for i in e),0)
        if value:out[e]=value
    return out

def third_coefficients(poly,A,B,cap):
    out=[0]*(cap+1)
    for (a,b,c),coefficient in poly.items():
        for k in range(min(len(A)-a,len(B)-b,cap+1-c)):
            out[c+k]+=coefficient*comb(k+2,2)*A[a+k]*B[b+k]
    return out

def main():
    files=['PROOF.md','verify.py','profile_data.json','normal.json','optimized.json']
    pins={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in files}
    expected={'PROOF.md':'d4528e9b0067e453b62e5933d48838932d66af16063b8e45930a76c947397b23','verify.py':'ab180af0ac4806204c80fd05927880354e492e1f21560fbc132ddf0cbd5596df','profile_data.json':'b6fc86843d75f2870748625af3da7fe7de65b7144a70bace5c6eb648f2229fd9','normal.json':'d368f10ca611cbf699c5bb652d244304214d3f22ed838d367126a37290db126a','optimized.json':'d368f10ca611cbf699c5bb652d244304214d3f22ed838d367126a37290db126a'}
    need(pins==expected,'reviewed source pins');data=json.loads((SRC/'profile_data.json').read_text());receipt=json.loads((SRC/'normal.json').read_text())
    A=data['A'];B=data['B'];C=[Fraction(x) for x in data['C']];Ci=data['C_integer'];D=data['C_scaling_denominator']
    need(A==[0,1,8,2,5,2,9,5,4,7,2] and B==[0,1,9,1,8,9,3,5,6,1,8],'explicit integer sequences')
    need(len(C)==14 and all(q[0]==0 and q[1]==1 and all(x>0 for x in q[1:]) for q in [A,B,C]),'all normalized rational entries positive')
    need(D>0 and len(Ci)==14 and all(Fraction(n,D)==r for n,r in zip(Ci,C)),'every denominator-cleared entry')
    common=1
    for x in C:common=lcm(common,x.denominator)
    need(D%common==0,'positive common denominator')
    need(B[2]%2 and A[3]%6,'factorial divisibility obstruction')
    values=[];matrix=[]
    for index,(a,c) in enumerate(product(range(8),repeat=2)):
        vectors=[third_coefficients(numerator(a,c,r),A,B,13) for r in (False,True)]
        rational=[sum(Fraction(n)*x for n,x in zip(row,C)) for row in vectors]
        integer=[sum(n*x for n,x in zip(row,Ci)) for row in vectors]
        need(rational[0]==rational[1] and rational[0]>0,'independent all64 rational equality')
        need(integer[0]==integer[1] and integer[0]>0 and all(x==D*y for x,y in zip(integer,rational)),'all integer equalities and scaling')
        old=receipt['records'][index];need(old['block_query']==[a,c] and integer==old['integer_profile_counts'] and [str(x) for x in rational]==old['normalized_counts'],'all retained exact values')
        row=[f-g for f,g in zip(*vectors)];need(sum(Fraction(x)*y for x,y in zip(row,C))==0,'full coefficient row annihilates exact profile')
        matrix.append(row);values.append({'query':[a,c],'integer_values':integer,'rational_values':[str(x) for x in rational]})
    for flags in ([],['-O']):
        replay=subprocess.check_output([sys.executable,*flags,str(SRC/'verify.py')]);need(replay==(SRC/'normal.json').read_bytes(),'paired producer replay')
    need(pins=={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in files},'frozen source')
    return dict(status='PASS',method='Independent bounded-port automaton followed by full numerator recovery and exact third-profile coefficient extraction',complete_rational_entries_checked=len(C),common_denominator=D,minimal_common_denominator=common,queries=64,exact_formal_values=128,all64_contrast_rows=matrix,all64_recomputed_values=values,source_sha256=pins,scope='Strict positivity is insufficient for formal cover-profile separation; the explicitly unrealizable data are not a tournament or named-response collision.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
