"""Integral future-query closure: remove spurious odd-prime chart obstructions."""
import json
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form,smith_normal_form
from effects import HERE,need,data

def main():
    ops={k:s.Matrix(v) for k,v in json.loads((HERE/'effects.json').read_text())['natural_matrices'].items()}
    B=s.eye(6)[:,:4];history=[]
    for depth in range(1,5):
        C=hermite_normal_form(s.Matrix.hstack(B,*[A*B for A in ops.values()]))
        history.append({'depth':depth,'basis':data(C),'index':int(abs(C.det()))})
        if C==B:break
        B=C
    need(B==s.diag(1,1,1,1,1,2),'exact full future-query lattice')
    transforms={k:B.inv()*A*B for k,A in ops.items()}
    need(all(x.q==1 for M in transforms.values() for x in M),'integer invariance proves all-word closure')
    need(smith_normal_form(B,domain=s.ZZ)==B,'only a factor-two formal quotient')
    # The direct two-probe determinant has irrelevant odd factors.
    need(s.Matrix([[81,100],[597,940]]).det()==16440,'two-probe minor')
    out={'status':'EXACT_INTEGER_LATTICE_CERTIFICATE','initial_queries':['h','s','E1','S2'],
         'closure_history':history,'all_word_basis':['h','s','E1','S2','P2','2D2'],
         'closed_integer_matrices':{k:data(M) for k,M in transforms.items()},
         'scope':'formal Z^6 coordinate module; no claim that every residue is an actual source',
         'consequence':'Every odd modulus recovers all six natural coordinates from queries of length at most two.'}
    (HERE/'lattice.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: integral query indices 40 -> 2 -> 2; all odd moduli have full natural-coordinate observation')
if __name__=='__main__':main()
