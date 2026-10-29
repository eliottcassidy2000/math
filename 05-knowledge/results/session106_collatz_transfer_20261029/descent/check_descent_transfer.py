#!/usr/bin/env python3
"""Standalone exact atlas and scoped hostile checks; no assert statements."""
import argparse, hashlib, json
from pathlib import Path

def require(p, msg):
    if not p: raise ValueError(msg)

def v2(n):
    require(type(n) is int and n>0, 'valuation domain')
    return (n & -n).bit_length()-1

def step(n,b=3):
    require(type(n) is int and n>0 and n%2==1, 'positive odd source required')
    a=v2(b*n+1)
    return (b*n+1)>>a,a

def carrier(word,b=3):
    A,B,C=1,1,0
    for a in word:
        require(type(a) is int and a>0, 'invalid letter')
        A,B,C=b*A,B<<a,b*C+B
    return A,B,C

def certificate(word,b=3):
    require(bool(word), 'empty certificate')
    A,B,C=carrier(word,b)
    require(B>A, 'noncontracting slope')
    mod=2*B
    residue=((B-C)*pow(A,-1,mod))%mod
    smallest=residue if residue>1 else residue+mod
    require((B-A)*smallest>C, 'unpaid cylinder')
    return {'word':list(word),'A':A,'B':B,'C':C,'modulus':mod,'residue':residue,'smallest_nonterminal':smallest}

def verify_cert(c,b=3):
    word=c['word']; require(c==certificate(word,b),'forged carrier or guard')
    # Independent symbolic affine propagation on x = smallest + modulus*t.
    u=c['smallest_nonterminal']; v=c['modulus']
    for a in word:
        numerator=b*u+1; coefficient=b*v
        require(numerator%(1<<a)==0 and coefficient%(1<<a)==0,'division')
        u=numerator>>a; v=coefficient>>a
        require(u%2==1 and v%2==0,'exact valuation for every t')
    require(u<c['smallest_nonterminal'] and v<c['modulus'],'uniform strict descent')
    # Check exact physical clocks independently for three literal members.
    for t in (0,1,17):
        n=c['smallest_nonterminal']+t*c['modulus']; x=n
        for a in word:
            x,actual=step(x,b); require(a==actual,'replay guard')
        require(x<n,'payment replay')

def build_atlas(bits):
    modulus=1<<bits; atlas={}; unresolved=[]
    for r in range(1,modulus,2):
        n=r if r>1 else r+modulus
        x=n; word=[]; A,B,C=1,1,0; D=0; accepted=False
        while True:
            y,a=step(x)
            if D+a+1>bits: break
            word.append(a); A,B,C=3*A,B<<a,3*C+B; D+=a; x=y
            if B>A and (B-A)*n>C:
                try: c=certificate(word)
                except ValueError: continue
                verify_cert(c)
                atlas[tuple(word)]=c; accepted=True
                break
        else: raise RuntimeError('unreachable')
        if not accepted: unresolved.append(r)
    certs=sorted(atlas.values(),key=lambda c:(c['modulus'],c['residue']))
    # Independently enumerate coverage; cylinders must be disjoint (first payment).
    covered=set()
    for c in certs:
        for r in range(c['residue'],modulus,c['modulus']):
            require(r not in covered,'overlapping certificates'); covered.add(r)
    require(covered.isdisjoint(unresolved),'mixed status')
    require(covered|set(unresolved)==set(range(1,modulus,2)),'incomplete partition')
    require(modulus-1 in unresolved,'lost all-one hostile')
    return {'bits':bits,'odd_residue_count':modulus//2,'covered_residue_count':len(covered),'unresolved_residue_count':len(unresolved),'certificate_count':len(certs),'max_word_length':max(map(lambda c:len(c['word']),certs)),'certificates':certs,'unresolved_residues':unresolved}

def hostiles():
    for k in range(1,257):
        n=(1<<(k+1))-1; x=n
        for j in range(1,k+1):
            x,a=step(x)
            require(a==1,'ones guard')
            require(x==pow(3,j)*(1<<(k+1-j))-1,'ones endpoint')
            require(x>n,'ones source expansion')
            # Every one is a leaf, so Q=1 returns to zero each step.
        require(x==2*pow(3,k)-1,'ones last endpoint')
    for b,word,n in [(5,[1,1,5],13),(5,[5,1,1],83)]:
        A,B,C=carrier(word,b); x=n
        for a in word:
            x,actual=step(x,b); require(actual==a,'5 guard')
        require(x==n and (B-A)*n-C==0,'5 cycle payment')
    rejects=0
    for w in ([],[0],[1],[1,1,5]):
        try: certificate(w,5)
        except ValueError: rejects+=1
        else: raise ValueError('accepted malformed or unpaid hostile')
    c=certificate([2]); bad=dict(c,C=c['C']+2)
    try: verify_cert(bad)
    except ValueError: rejects+=1
    else: raise ValueError('accepted corrupt carry')
    return {'ones_lengths_checked':256,'five_cycles_checked':2,'rejections_checked':rejects}

def finite_descent(limit):
    max_steps=0; witness=None
    # Exactly bounded experiment, not an all-integer certificate.
    for n in range(3,limit+1,2):
        x=n; count=0
        while x>=n and count<10000:
            x,a=step(x); count+=1
        require(x<n,'finite range descent failed or work cap exceeded')
        if count>max_steps: max_steps,witness=count,n
    return {'maximum_source':limit,'odd_sources_checked':(limit-1)//2,'work_cap_per_source':10000,'max_first_descent_steps':max_steps,'first_maximizer':witness}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',default=None); p.add_argument('--atlas',default=None); args=p.parse_args()
    all_atlases=[build_atlas(b) for b in (8,12,16)]
    summary={'status':'PASS','scope':'Infinite arithmetic cylinders proved; finite exhaustive atlas and bounded orbit checks only. No universal descent claim.','atlas_summaries':[{k:v for k,v in a.items() if k not in ('certificates','unresolved_residues')} for a in all_atlases],'hostiles':hostiles(),'finite_descent':finite_descent(100001)}
    full=json.dumps(all_atlases[-1],sort_keys=True,indent=2)+'\n'
    summary['atlas_sha256']=hashlib.sha256(full.encode()).hexdigest()
    if args.atlas: Path(args.atlas).write_text(full)
    out=json.dumps(summary,sort_keys=True,indent=2)+'\n'
    if args.output: Path(args.output).write_text(out)
    print(out,end='')
if __name__=='__main__': main()
