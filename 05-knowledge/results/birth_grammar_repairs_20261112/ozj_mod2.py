"""Complete finite matrix-semigroup automaton modulo2, with operation flags."""
from common import *
from ozj_periodic import M
from collections import Counter
import sympy as sp

def flatten(A):return tuple(v%2 for row in A for v in row)
def step(A,label):
    return tuple(sum(M[label][i][k]*A[4*k+j] for k in range(4))%2 for i in range(4) for j in range(4))

def run():
    root=(flatten([[int(i==j) for j in range(4)] for i in range(4)]),0);states=[root];ids={root:0};words=[''];edges=[]
    for idx,(A,flags) in enumerate(states):
        row=[]
        for label in ['O','Z','J']:
            B=step(A,label);f=flags|(1 if label=='J' else 2);key=(B,f)
            if key not in ids:ids[key]=len(states);states.append(key);words.append(words[idx]+label)
            row.append(ids[key])
        edges.append(row)
    x=sp.Symbol('x');hist=Counter();records=[]
    for idx,(A,flags) in enumerate(states):
        coeff=[int(c)%2 for c in sp.Matrix(4,4,A).charpoly(x).all_coeffs()];hist[flags,tuple(coeff)]+=1
        records.append({'id':idx,'matrix':list(A),'flags':flags,'first_word':words[idx],'charpoly_mod2':coeff,'transitions_O_Z_J':edges[idx]})
    out={'scope':'Exhaustive closure of the exact4x4 O/Z/J matrices modulo2, with flags recording presence of J and presence of O or Z. This finite closure covers all words; its characteristic-polynomial consequences are proved separately.','states':len(states),'histogram':[{'flags':f,'charpoly':list(c),'states':n} for (f,c),n in hist.items()],'records':records}
    save('ozj_mod2.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
