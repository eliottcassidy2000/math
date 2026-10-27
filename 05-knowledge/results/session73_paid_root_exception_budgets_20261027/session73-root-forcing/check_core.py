#!/usr/bin/env python3
"""Exact native phase/payment controls. Python and -O must agree."""
from fractions import Fraction as F
from math import factorial
from itertools import product
import json, hashlib

checks = 0
def need(cond, message):
    global checks
    checks += 1
    if not cond:
        raise RuntimeError(message)
def ratio(x):
    x=F(x);return [x.numerator,x.denominator]
def step(n):
    z=3*n+1;a=(z&-z).bit_length()-1
    return z>>a,a

def replay(n, limit=10000):
    need(type(n) is int and n>0 and n%2==1,'source type')
    word=[]; states=[n]
    while states[-1]!=1 and len(word)<limit:
        m,a=step(states[-1]);word.append(a);states.append(m)
    need(states[-1]==1,'unproved cutoff source: '+str(n))
    return states,word

def weight_from_word(word):
    if not word:return F(1)
    K=sum((a-1)//2 for a in word);tau=len(word)
    return F(2*factorial(K)*factorial(tau),factorial(tau+K+1))
def array_weight(L,K):
    return F(2*factorial(K)*factorial(L+1),factorial(L+K+2))
def siblings(n,j):return ((3*n+1)*4**j-1)//3

def main():
    records={};occupation={};source_bill=F(0);source_clock=F(0)
    phase_cases=0;word_digest=hashlib.sha256()
    for n in range(1,1002,2):
        states,word=replay(n);w=weight_from_word(word)
        sign=(-1)**sum(word)
        word_digest.update((str(n)+':'+','.join(map(str,word))+';').encode())
        records[n]=(states,word,w,sign)
        for t,(x,a) in enumerate(zip(states[:-1],word)):
            e=(-1)**sum(word[t:]);ep=(-1)**sum(word[t+1:]);sigma=(-1)**a
            b=F(3*x,3*x+1);d=1-b*b;c=sigma*b
            need(ep==sigma*e,'shortcut parity')
            force=ep-c*e;residual=d*e-c*force
            need(ep*ep+residual*residual/d==e*e+force*force/d,'signed passive identity')
            need(force*force/d==F(1,6*x+1),'phase-aware cost')
            raw=ep-b*e
            need(raw*raw/d==(F(6*x+1) if a%2 else F(1,6*x+1)),'unsigned phase price')
            phase_cases+=1
        if n%3==0:
            path_bill=sum((F(1,6*x+1) for x in states[:-1]),F(0))
            source_bill+=w*path_bill
            source_clock+=w*(len(word)+sum(word))
            for x in states[:-1]:occupation[x]=occupation.get(x,F(0))+w
    vertex_bill=sum((mass/F(6*x+1) for x,mass in occupation.items()),F(0))
    need(source_bill==vertex_bill,'source-preserving occupation regrouping')
    for x,mass in occupation.items():
        states,word=replay(x)
        need(mass<=weight_from_word(word),'finite occupation below native W')
    need(source_bill<F(13,57),'native uniform paid bound')

    # Analytic incoming sibling tail, including the exceptional killed root row.
    row_cases=0
    for n in range(1,102,2):
        states,word,w,sign=records[n]
        for J in range(0,12):
            if n==1:
                row=sum((F(2,(j+1)*(j+2)) for j in range(1,J+1)),F(0))+F(2,J+2)
                need(row==1,'root inverse row')
            elif n%3==0:
                need(w-w==0,'leaf row')
            else:
                L=len(word)-1;K=sum((a-1)//2 for a in word)
                row=sum((array_weight(L+1,K+j) for j in range(J)),F(0))+array_weight(L,K+J)
                need(row==w,'exact infinite inverse row with tail')
            row_cases+=1
    # Complete finite-value grid for the native contraction, infinite rows paid above.
    targets=[1,3,5,7,13];contraction_cases=0
    for vals in product(range(-2,3),repeat=len(targets)):
        norm=sum((records[n][2]*v*v for n,v in zip(targets,vals)),F(0))
        image=sum((records[n][2]*v*v for n,v in zip(targets,vals) if n%3),F(0))
        defect=sum((records[n][2]*v*v for n,v in zip(targets,vals) if n%3==0),F(0))
        need(norm-image==defect,'native contraction defect')
        contraction_cases+=1

    divergence=[];clock_rows=[];harmonic=F(0);clock_partial=F(0)
    for k in range(20):
        n=(10*64**k-1)//3
        states,word=replay(n);w=weight_from_word(word)
        need(word==[6*k+1,4],'odd-phase native family word')
        need(n%3==0 and w==F(4,(3*k+2)*(3*k+3)*(3*k+4)),'odd-phase family mass')
        b=F(3*n,3*n+1);d=1-b*b
        need(w*(1+b)**2/d==w*(20*64**k-1),'actual signed ROOT divergent invoice')
        need(w*(1-b)**2/d==w/F(6*n+1),'phase repair family')
        divergence.append({'k':k,'n':n,'lambda':ratio(w),'unsigned_paid':ratio(w*(6*n+1)),'phase_paid':ratio(w/F(6*n+1))})
        n=(64**(k+1)-1)//3
        states,word=replay(n);w=weight_from_word(word);clock=len(word)+sum(word)
        need(word==[6*k+6] and clock==6*k+7,'root ray literal clock')
        need(w==F(2,(3*k+3)*(3*k+4)),'root ray mass')
        need(w*clock>=F(1,k+1),'harmonic clock lower bound')
        harmonic+=F(1,k+1);clock_partial+=w*clock
        clock_rows.append({'k':k,'n':n,'lambda':ratio(w),'standard_clock':clock,'shortcut_clock':sum(word),'snapshot_frontier_paid':ratio(w/(1-F(3*n,3*n+1)**2))})
    # All depth horizons in this declared finite test: exact infinite family formula.
    target=1;snapshot=[]
    for horizon in range(13):
        st,wo=replay(target)
        need(len(wo)==horizon and target%3!=0,'horizon target')
        amin=2 if target%3==1 else 1;base=(2**amin*target-1)//3;j0=(-base)%3
        family=[]
        for k in range(8):
            j=j0+3*k;n=siblings(base,j);states,word=replay(n)
            need(n%3==0 and states[1]==target,'named inverse frontier port')
            need(len(word)==horizon+1,'frontier outside exact horizon')
            w=weight_from_word(word);b=F(3*n,3*n+1);d=1-b*b
            # Target is present but this new predecessor is zero in the horizon profile.
            ep=(-1)**sum(wo);e=0;c=(-1)**word[0]*b
            need(w*(ep-c*e)**2/d==w/d,'snapshot phase-aware frontier invoice')
            family.append({'k':k,'source':n,'lambda':ratio(w),'paid':ratio(w/d)})
        snapshot.append({'horizon':horizon,'target':target,'families':family})
        candidates=[]
        for a in range(1,8):
            z=2**a*target-1
            if z%3==0:
                p=z//3
                if p>target and p%3 and p%2:candidates.append(p)
        need(bool(candidates),'unit inverse extension')
        target=min(candidates)
    # Clock convention hostile and λ-only quotient hostile.
    need(records[21][3]==1 and (-1)**(len(records[21][1])+sum(records[21][1]))==-1,'two clocks really differ')
    need(F(1,6)>0,'root delta is null under lambda but pullback at21 is not')
    result={'status':'FINITE-EXACT controls; all-order claims require report proofs','checks':checks,'sources_odd_1_to_1001':501,'phase_transition_cases':phase_cases,'authenticated_leaf_sources':len([n for n in records if n%3==0]),'occupation_vertices':len(occupation),'finite_leaf_paid_bill':ratio(source_bill),'finite_leaf_standard_clock_bill':ratio(source_clock),'incoming_row_tail_cases':row_cases,'contraction_vector_cases':contraction_cases,'source_word_sha256':word_digest.hexdigest(),'unsigned_gauge_divergence':divergence,'literal_clock_harmonic_control':{'partial_clock_bill':ratio(clock_partial),'harmonic_lower_bound':ratio(harmonic),'rows':clock_rows},'snapshot_frontier_families':snapshot,'scope':{'eventual_phase_profile_uniform_bound':'13/57','native_spatial_current_uniform_bound':'OPEN','all_positive_source_ROOT':'OPEN','finite_zero_extension_total_forcing':'diverges, even with native phase','guarded_certified_path_payment':'bounded by13/57 with original lambda masses'}}
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
