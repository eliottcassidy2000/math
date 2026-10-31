#!/usr/bin/env python3
"""Exact positive renewal price of every mixed periodic singleton word."""
import json,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from audit_port_words import adjacency,count,word_ports

def act(state,bits):
    t,q,g,h=state
    for b in bits:
        if b=='0':t,q,g=t+q+g,2*q,0
        else:q,g=2*q+h,h
    return t,q,g,h

def profile(bits):
    n=len(bits);B=1<<n;last=int(bits[-1]);z=act((0,0,last,1),bits);o=act((0,1,last,1),bits)
    gamma,beta=z[:2];alpha=o[0]-gamma
    return B,alpha,beta,gamma,F(gamma)-F(alpha*beta,B-1)

def renewal(bits):
    n=len(bits);den=(1<<n)-1;price=F(0);runs=[]
    for j in range(n):
        if bits[j]=='0' and bits[(j-1)%n]=='1':
            r=1
            while bits[(j+r)%n]=='0':r+=1
            end=(j+r)%n;rot=bits[end:]+bits[:end];u=F(int(rot,2),den)
            local=1-u*(1-F(1,1<<r));price+=local;runs.append([r,local])
    return price,runs

def main():
    words=0;minimums=[];formula_controls=0;adjacency_controls=0
    for L in range(1,15):
        best=None;minwords=[]
        for mask in range(1<<L):
            bits=format(mask,f'0{L}b');B,alpha,beta,gamma,delta=profile(bits)
            mixed='0' in bits and '1' in bits
            if (delta>0)!=mixed:raise RuntimeError('positive iff mixed')
            if mixed:
                energy,runs=renewal(bits)
                if energy!=delta or any(v<=0 for r,v in runs):raise RuntimeError('renewal equality')
                sharp=F((1<<(L//2))+(1<<((L+1)//2))-2,(1<<L)-1)
                R=len(runs)
                if delta<R*sharp:raise RuntimeError('run-count refined lower bound')
                if (delta==sharp)!=(R==1 and abs(bits.count('0')-bits.count('1'))<=1):
                    raise RuntimeError('sharp equality classification')
                corr=[sum(bits[j]=='0' and bits[(j+t)%L]=='1' for j in range(L)) for t in range(L)]
                if any(corr[t]>R*min(t,L-t) for t in range(1,L)):
                    raise RuntimeError('cyclic transition union bound')
                if delta!=R-sum((F(corr[t]*(1<<(L-1-t)),B-1) for t in range(1,L)),F(0)):
                    raise RuntimeError('correlation energy')
                if best is None or delta<best:best=delta;minwords=[bits]
                elif delta==best:minwords.append(bits)
            if L<=7:
                for h in [1,5]:
                    state=(7,3,int(bits[-1])*h,h)
                    for r in range(5):
                        drift=F(alpha)*(3+F(h*beta,B-1))*F(B**r-1,B-1)+h*delta*r
                        expected=7+drift
                        if act(state,bits*r)[0]!=expected:raise RuntimeError('period formula')
                        formula_controls+=1
            if L<=4:
                for r in [1,2,3]:
                    prefix=[('S',int(bits[-1]))];word=prefix+[('S',int(b)) for b in bits*r]
                    state=act((1,0,0,1),bits[-1]);actual=act(state,bits*r)[0]
                    if count(adjacency(word),word_ports(word))!=actual:raise RuntimeError('physical tournament count')
                    adjacency_controls+=1
            words+=1
        minimums.append({'length':L,'minimum':str(best),'examples':minwords[:4],'minimizer_words':len(minwords)})
    paid=[]
    # D(mask11) gives(F,q,g,h)=(9,18,4,5); S0 gives(31,36,0,5).
    for r in range(9):
        word=[('D',11),('S',0)]+[('S',int(b)) for b in '1100'*r]
        expected=(1<<(4*r+5))+2*r-1
        if act((31,36,0,5),'1100'*r)[0]!=expected or count(adjacency(word),word_ports(word))!=expected:
            raise RuntimeError('actual tuned positive-price family')
        paid.append({'repetitions':r,'vertices':len(adjacency(word)),'raw_H':expected})
    print(json.dumps({'status':'passed','complete_binary_words':words,'through_length':14,
                      'state_formula_controls':formula_controls,'actual_adjacency_controls':adjacency_controls,
                      'finite_minimum_drift':minimums,
                      'tuned_native_family':paid,
                      'scope':'All-size positivity and the sharp minimum/equality cases follow the cyclic-correlation proof. Exact controls cover all binary words through14.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
