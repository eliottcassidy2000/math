"""Exact full-observer phase repair, with primitive mod-11 source guards."""
from itertools import product
import json
import sympy as s
from effects import HERE,need,data
from algebra import power

STATES=(0,1,3,4,5,6,7,None)
def step(r,letter):
    a,c=(61,9) if letter=='3' else (421,33)
    if r is None:return 0
    den=(a-2*r)%11
    return None if den==0 else c*pow(den,-1,11)%11
def advance(states,word):
    states=tuple(states)
    for letter,count in word:
        for _ in range(count):
            need(letter!='5' or 7 not in states,'every primitive C5 guard')
            states=tuple(step(r,letter) for r in states)
    return states
def neutral(j,e):
    a,b=12*7**(e-1),48*7**(e-1)
    return [('3',j),('5',b),('3',a-j)]
def program(states,e):
    j=next(j for j in range(4) if 7 not in advance(states,[('3',j)]))
    word=neutral(j,e)
    if j==3:word+=neutral(1,e)
    word+=[('5',1)]
    return j,word
def matrix_word(word,ops,mod):
    M=s.eye(6)
    for letter,count in word:M=(M*power(ops[letter],count,mod)).applyfunc(lambda x:x%mod)
    return M
def main():
    ops={k:s.Matrix(v) for k,v in json.loads((HERE/'effects.json').read_text())['natural_matrices'].items()}
    tests=[];examples=[]
    for e in range(1,4):
        mod=7**e;cache={}
        for r in range(1,4):
            histogram={};maximum=0
            for states in product(STATES,repeat=r):
                j,word=program(states,e);key=tuple(word)
                if key not in cache:
                    cache[key]=matrix_word(word,ops,mod)
                need(cache[key]==ops['5'].applyfunc(lambda x:x%mod),'full six-coordinate phase action')
                end=advance(states,word);need(len(set(end))==1,'shared controller synchronization')
                cost=sum((4 if k=='3' else 8)*n for k,n in word)
                maximum=max(maximum,cost);histogram[j]=histogram.get(j,0)+1
            bound=(432 if r<=2 else 864)*7**(e-1)+8
            need(maximum<=bound,'constructive native correction bound')
            tests.append({'precision':e,'number_of_sources':r,'tuples':len(STATES)**r,'j_histogram':histogram,'maximum_native_moves':maximum,'bound':bound})
        if e==1:
            for states in ((7,),(0,7),(0,7,3)):
                j,word=program(states,e);examples.append({'states':states,'j':j,'word':word,'end':advance(states,word)})
    # Each member of the four-cycle prevents the first C5 after any C3 prefix.
    hostile=(0,7,3,None)
    need(all(7 in tuple(step(r,'3') for r in hostile) for _ in [0]),'hostile cycle invariant')
    for j in range(12):need(7 in advance(hostile,[('3',j)]),'four-source first-C5 obstruction')
    need(ops['5'][:2,:2].det()%7==3 and all(pow(int(ops['3'][:2,:2].det()),j,7)!=3 for j in range(6)),'C3 powers cannot replace C5 phase action')
    # Actual positive integer guards, starting at the prior H497 primitive entry.
    # The pinned k2 source enters via C5: (6729,2167) -> (2828575,222057).
    b,t=2828575,222057
    import math
    need(math.gcd(b,t)==1 and b%3 and t%3==0,'actual primitive source entry')
    initial=[b,t];j,word=program((t*pow(b,-1,11)%11,),1)
    for letter,count in word:
        for _ in range(count):
            b,t=((61*b-2*t,9*b) if letter=='3' else (421*b-2*t,33*b))
            need(math.gcd(b,t)==1 and b>t>0 and b%2 and t%2 and b%3 and t%3==0,'exact integer guard along legal program')
    out={'status':'FINITE_EXACT','exhaustive_controller_tests':tests,'examples':examples,'hostile_four_states':hostile,
         'integer_replay':{'initial':initial,'word':word,'final_decimal_digits':len(str(b)),'final_mod7':[b%7,t%7]}}
    (HERE/'legalizer.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: all 584 source tuples at each of three precisions, full phase identities and primitive integer replay')
if __name__=='__main__':main()
