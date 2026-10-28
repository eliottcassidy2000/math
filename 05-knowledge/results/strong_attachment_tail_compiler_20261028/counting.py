"""Exact contraction of the terminal 1-cycle stage; preserves inherited profiles."""
from collections import defaultdict
from fractions import Fraction
from math import factorial, gcd, comb
from pathlib import Path
import argparse
import importlib.util
import json
import sys
from time import perf_counter

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py'
spec=importlib.util.spec_from_file_location('inherited_local',BASE)
inherited=importlib.util.module_from_spec(spec); spec.loader.exec_module(inherited)

def need(ok,message):
    if not ok: raise RuntimeError(message)

def scales(bound):
    facts=[factorial(i) for i in range(bound+1)]
    twos=[(f&-f).bit_length()-1 for f in facts]
    return facts,twos,[f>>v for f,v in zip(facts,twos)]

def prefix(bound,*,core=False,stop=1,progress=False):
    width=max(1,bound.bit_length()); digit=(1<<width)-1
    masks=[(1<<(width*c))-1 for c in range(bound//2+2)]
    facts,twos,odds=scales(bound)
    states={(0,0):1}; rows=[]; precomputed=0
    for ell in range(bound if bound%2 else bound-1,stop,-2):
        idx=(ell-1)//2
        row=sum(gcd(ell,j)<<(width*i) for i,j in enumerate(range(1,ell,2)))
        tables={}
        for degree,_ in states:
            if degree in tables: continue
            entries=[]
            for k in range((bound-degree)//ell+1):
                h=ell**k*odds[k]
                base=k*(ell-1)//2+ell*k*(k-1)//2-twos[k]
                for excess in (range((bound-degree)//ell-k+1) if core and k else (0,)):
                    nd=degree+ell*(k+excess)
                    factor,remainder=divmod(odds[nd],odds[degree]*h)
                    need(not remainder,'prefix transfer integer')
                    if core and excess: factor*=(-1)**excess*comb(k+excess-1,excess)
                    count=max(0,(min(ell-2,bound-nd)-1)//2+1)
                    entries.append((k,nd,base,factor,masks[count],k*row));precomputed+=1
            tables[degree]=entries
        new=defaultdict(int); transitions=0
        for (degree,response),weight in states.items():
            charge=(response>>(width*idx))&digit
            for k,nd,base,factor,mask,offset in tables[degree]:
                new[nd,(response+offset)&mask]+=(weight*factor)<<(base+k*charge)
                transitions+=1
        states={key:value for key,value in new.items() if value}
        rows.append(dict(length=ell,states=len(states),transitions=transitions))
        if progress and (ell<=9 or ell%20==1):
            print(json.dumps(dict(stage=ell,states=len(states))),file=sys.stderr,flush=True)
    return states,dict(rows=rows,transitions=sum(r['transitions'] for r in rows),
                       precomputed_transfers=precomputed,peak_states=max([1]+[r['states'] for r in rows]))

def evaluate(coefficients,power):
    """Sparse Horner at 2**power; coefficients in descending exponent order."""
    answer=0; previous=coefficients[0][0]
    for exponent,weight in coefficients:
        answer=(answer<<(power*(previous-exponent)))+weight; previous=exponent
    return answer<<(power*previous)

def contracted(bound,*,core=False,progress=False):
    need(type(bound) is int and bound>=0,'nonnegative integer bound')
    states,stats=prefix(bound,core=core,progress=progress)
    facts,twos,odds=scales(bound)
    bydegree=defaultdict(list)
    for (degree,response),weight in states.items(): bydegree[degree].append((response,weight))
    totals=[0]*(bound+1); evaluations=0; products=0
    for degree,coefficients in sorted(bydegree.items()):
        coefficients.sort(reverse=True)
        totals[degree]+=evaluate(coefficients,0); evaluations+=1
        for k in range(1,bound-degree+1):
            value=evaluate(coefficients,k); evaluations+=1
            base=k*(k-1)//2-twos[k]
            value <<= base
            for remaining in (range(k,bound-degree+1) if core else (k,)):
                n=degree+remaining
                factor,remainder=divmod(odds[n],odds[degree]*odds[k])
                need(not remainder,'tail transfer integer')
                if core and remaining>k: factor*=(-1)**(remaining-k)*comb(remaining-1,k-1)
                totals[n]+=factor*value; products+=1
    answers=[]
    for n,total in enumerate(totals):
        value,remainder=divmod(total,odds[n]); need(not remainder,'final integer')
        answers.append(value)
    stats.update(tail_input_states=len(states),tail_polynomials=len(bydegree),tail_evaluations=evaluations,
                 tail_products=products,mode='pair-free' if core else 'ordinary',bound=bound)
    return answers,stats

def contracted13(bound,*,progress=False):
    """Sum the complete {1,3} tail with one moment per total future degree."""
    states,stats=prefix(bound,stop=3,progress=progress)
    facts,twos,odds=scales(bound); width=max(1,bound.bit_length());digit=(1<<width)-1
    polynomials=defaultdict(list);short=defaultdict(list)
    for (degree,response),weight in states.items():
        c=response&digit
        if bound-degree<3:
            short[degree].append((c,weight));continue
        r3=(response>>width)&digit
        need((3*c-r3)%2==0,'three-tail parity')
        q=(3*c-r3)//2; b=c-q
        need(q>=0 and b>=0,'three-tail divisibility counts')
        polynomials[degree,q].append((b,weight<<twos[degree]))
    coeffs=[]
    for S in range(bound+1):
        row=[]
        for k in range(S//3+1):
            r=S-3*k
            orbit,remainder=divmod(facts[S],3**k*facts[k]*facts[r])
            need(not remainder,'tail permutation count')
            e=k+3*k*(k-1)//2+r*(r-1)//2+k*r
            row.append((S-2*k,orbit<<e))
        coeffs.append(row)
    kernels={};total=[0]*(bound+1);products=0;evaluations=0
    for (degree,q),poly in sorted(polynomials.items()):
        poly.sort(reverse=True)
        for S in range(bound-degree+1):
            key=S,q
            if key not in kernels:kernels[key]=evaluate(coeffs[S],q)
            value=evaluate(poly,S);evaluations+=1
            total[degree+S]+=comb(degree+S,degree)*kernels[key]*value;products+=1
    for degree,poly in short.items():
        poly.sort(reverse=True)
        for r in range(bound-degree+1):
            value=evaluate(poly,r)<<twos[degree];evaluations+=1
            total[degree+r]+=comb(degree+r,degree)*value<<(r*(r-1)//2);products+=1
    values=[]
    for n,numerator in enumerate(total):
        value,rem=divmod(numerator,facts[n]);need(not rem,'three-tail final integer');values.append(value)
    stats.update(tail_input_states=len(states),tail_polynomials=len(polynomials)+len(short),
                 tail_products=products,tail_evaluations=evaluations,kernel_evaluations=len(kernels),
                 mode='ordinary-13',bound=bound)
    return values,stats

def main():
    if hasattr(sys,'set_int_max_str_digits'): sys.set_int_max_str_digits(0)
    parser=argparse.ArgumentParser();parser.add_argument('--bound',type=int,default=100)
    parser.add_argument('--core',action='store_true');parser.add_argument('--baseline',action='store_true')
    parser.add_argument('--tail13',action='store_true')
    parser.add_argument('--output');parser.add_argument('--progress',action='store_true');args=parser.parse_args()
    start=perf_counter()
    if args.baseline:
        values,stats=inherited.local_dp(args.bound,core=args.core)
        need(all(v.denominator==1 for v in values),'baseline integer');values=list(map(int,values))
    elif args.tail13:
        need(not args.core,'three-tail implementation is ordinary only')
        values,stats=contracted13(args.bound,progress=args.progress)
    else:values,stats=contracted(args.bound,core=args.core,progress=args.progress)
    result=dict(values=values,statistics=stats,seconds=perf_counter()-start)
    text=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else:print(text,end='')

if __name__=='__main__':main()
