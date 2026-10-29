"""Exact ordered valuation cylinders and signed accelerated Collatz controls."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import json

CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def valuation(n):
    if n==0: raise ValueError('pole')
    return (abs(n)&-abs(n)).bit_length()-1

def carry(word):
    a=c=0
    for e in word: c=3*c+2**a; a+=e
    return a,c

def follows(n,word,sign):
    for e in word:
        x=3*n+sign
        if not x or valuation(x)!=e: return False
        n=x//2**e
    return True

def cycle(start,sign):
    seen={}; nodes=[]; word=[]; n=start
    while n not in seen:
        seen[n]=len(nodes); nodes.append(n)
        x=3*n+sign; e=valuation(x); word.append(e); n=x//2**e
    need(seen[n]==0,'control is a cycle')
    a,c=carry(word); d=2**a-3**len(word)
    need(Fraction(sign*c,d)==start,'signed cycle equation')
    need(follows(-start,word,-sign),'negation conjugacy')
    return {'sign':sign,'nodes':nodes,'ordered_valuations':word,'total_binary_clock':a,
            'ternary_clock':len(word),'carry':c,'denominator':d,'affine_slope':str(Fraction(3**len(word),2**a))}

def main():
    # All ordered positive words with length<=5 and total clock<=12.
    records=0; residue_checks=0
    for length in range(1,6):
        for word in product(range(1,13),repeat=length):
            if sum(word)>12: continue
            a,c=carry(word); modulus=2**(a+1)
            residues=[]
            for sign in [1,-1]:
                x=((2**a-sign*c)*pow(3**length,-1,modulus))%modulus
                need(x%2==1 and follows(x,word,sign),'exact valuation cylinder')
                # Independent exhaustive residue check for total clock <=8.
                if a<=8:
                    actual=[y for y in range(1,modulus,2) if follows(y,word,sign)]
                    need(actual==[x],'independent cylinder uniqueness'); residue_checks+=modulus//2
                residues.append(x)
            need((sum(residues))%modulus==0,'source negation retains ordered clocks')
            records+=1
    examples=[cycle(1,1),cycle(1,-1),cycle(5,-1),cycle(17,-1)]
    a,c=carry((1,2)); aa,cc=carry((2,1))
    need(a==aa and c!=cc,'same scalar clocks lose carry')
    # Same slope, distinct translations: port differences alone lose the source.
    need((a,c,cc)==(3,5,7),'ordered carry witness')
    result={'status':'FINITE-EXACT controls for proved identities','valuation_words':records,
            'literal_odd_residue_tests':residue_checks,'cycles':examples,
            'same_clock_different_carry':{'words':[[1,2],[2,1]],'carries':[c,cc]},'checks':CHECKS}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
