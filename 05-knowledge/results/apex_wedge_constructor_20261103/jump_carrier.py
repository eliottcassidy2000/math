"""Six integral coordinates for the edge-jump series, with a unit Hankel minor."""
from itertools import product
from pathlib import Path
import json
import sympy as s
from mechanisms import *
HERE=Path(__file__).resolve().parent
INITIAL=s.Matrix([1,0,0,0,0,0])
OUTPUT=s.Matrix([[0,-1,1,-1,1,1]])
PREFIXES=('', '2', '1', '20', '10', '201')
SUFFIXES=('', '0', '2', '00', '02', '200')

def transition(letter):
    a=letter&1;b=(letter>>1)&1
    return s.Matrix([[1,0,0,0,0,0],[0,2,0,1,0,0],[0,0,2,0,1,0],
                     [a,0,0,0,0,0],[b,0,0,0,0,0],[0,b,-a,0,0,1]])

def state(word):
    # z=(constant, floor(Qw/2), floor(Qu/2), last_w, last_u, area).
    M=N=a=b=A=0
    for letter in word:
        w=letter&1;u=(letter>>1)&1
        M,N,A,a,b=2*M+a,2*N+b,A+u*M-w*N,w,u
    return s.Matrix([1,M,N,a,b,A])

def evaluate(word):
    _,M,N,a,b,A=state(word)
    return int(A+N+b-M-a)

def decode(values,observation):
    z=observation.inv()*s.Matrix(values)
    need(all(v.q==1 for v in z),'integer jump observation lattice')
    constant,M,N,a,b,A=map(int,z)
    need(constant==1 and M>=0 and N>=0 and a in (0,1) and b in (0,1),'actual word guard')
    qw,qu=2*M+a,2*N+b
    n=max(qw.bit_length(),qu.bit_length())
    word=tuple(((qw>>j)&1)+2*((qu>>j)&1) for j in reversed(range(n)))
    need(state(word)==z,'actual area replay')
    return word

def main():
    transitions=[transition(d) for d in range(4)]
    def mat(word):
        result=s.eye(6)
        for d in word:result=transitions[int(d)]*result
        return result
    R=s.Matrix.hstack(*(mat(p)*INITIAL for p in PREFIXES))
    O=s.Matrix.vstack(*(OUTPUT*mat(q) for q in SUFFIXES))
    H=O*R
    need(R.det()==O.det()==-1 and H.det()==1,'unit reachable/observable/Hankel minor')
    need(all(v.q==1 for v in O.inv()),'integral observer inverse')
    checked=0
    for n in range(7):
        for word in product(range(4),repeat=n):
            z=state(word)
            need(z==mat(word)*INITIAL,'exact update identity')
            w=tuple(d&1 for d in word);u=tuple((d>>1)&1 for d in word)
            need(evaluate(word)==jump(w,u)//2,'area versus contraction')
            values=[evaluate(word+tuple(map(int,q))) for q in SUFFIXES]
            decoded=decode(values,O)
            k=0
            while k<len(word) and word[k]==0:k+=1
            need(decoded==word[k:],'all-size normalized source decoder finite controls')
            checked+=1
    for n in (25,100,1000):
        word=tuple((j*j+3*j+1)%4 for j in range(n))
        values=[evaluate(word+tuple(map(int,q))) for q in SUFFIXES]
        need(decode(values,O)==word,'long source decode')
        w=tuple(d&1 for d in word);u=tuple((d>>1)&1 for d in word)
        need(2*evaluate(word)==gap_h(w,u)-gap_h(w,u,False),'long literal gap evaluation')
    # A legal integral linear state is not automatically an actual source word.
    false=INITIAL.copy();false[-1]=1
    bad=O*false
    try:decode(list(bad),O)
    except RuntimeError:pass
    else:raise RuntimeError('unrealized area must be rejected')
    out={'status':'PROVED_FINITE_LINEAR_CERTIFICATE','linear_dimension':6,
         'coordinate_order':['1','floor(Qw/2)','floor(Qu/2)','last_w','last_u','area'],
         'transition_matrices':[m.tolist() for m in transitions],'output':list(OUTPUT),
         'prefixes':PREFIXES,'suffixes':SUFFIXES,'hankel':H.tolist(),'hankel_determinant':1,
         'observation_inverse':O.inv().tolist(),'exhaustive_words_through_6':checked,
         'long_lengths':[25,100,1000],'scope':'all common chain append words; edge-jump observations; normalized words modulo leading 00'}
    (HERE/'jump_carrier.json').write_text(json.dumps(out,indent=2,default=int)+'\n',encoding='utf8')
    print('PASS: sharp integral dimension six; actual-word recovery;',checked,'words')
if __name__=='__main__':main()
