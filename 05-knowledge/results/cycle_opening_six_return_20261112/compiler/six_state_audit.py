"""Exact linear-state minimality and periodic-program obstruction controls.
Run from any directory; requires SymPy. No floating-point calculations.
"""
from pathlib import Path
import itertools, json, sys
import sympy as sp
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f

def require(ok,message):
    if not ok: raise RuntimeError(message)

ALPHABET={'E':(), 'B':(1,), 'C':(0,1), 'D':(1,1)}
BASIS=[tuple(int(i==j) for i in range(6)) for j in range(6)]
M={k:sp.Matrix.hstack(*[sp.Matrix(f.advance(v,w)) for v in BASIS]) for k,w in ALPHABET.items()}
SEED=sp.Matrix(f.advance(f.SINGLE,(1,)))
ELL=sp.Matrix([[1,0,0,0,1,0]])

def program(word):
    matrix=sp.eye(6)
    for letter in word: matrix=M[letter]*matrix
    return matrix

def rows_literal(word):
    T,root=f.a.funnel((0,),[0],(1,))
    for letter in word: T,root=f.a.funnel(T,[root],ALPHABET[letter])
    return T,root

def direct_state(T,root):
    p,e,s=f.a.bin_values(T,root,3)
    return sp.Matrix((p[1],e[1],s[2],p[2],e[2],s[3]))

def main():
    # These short legal words were found by a breadth-first rational-rank probe.
    columns=['','E','B','C','D','EE']
    observers=['','E','B','D','EE','EB']
    reach=sp.Matrix.hstack(*[program(w)*SEED for w in columns])
    observe=sp.Matrix.vstack(*[ELL*program(w) for w in observers])
    hankel=observe*reach
    difference_words=['E','B','C','D','EE','EB']
    differences=sp.Matrix.hstack(*[program(w)*SEED-SEED for w in difference_words])
    difference_hankel=observe*differences
    require(difference_hankel.det()==811008,'affine six-state lower bound')
    require(reach.det()!=0 and observe.det()!=0,'six reachable and observable coordinates')
    require(hankel.det()==552960,'exact Hankel determinant')
    # Independent subset dynamic programming on literal tournaments, including
    # the actual ready child rather than merely the compiled linear readout.
    controls=[]
    for word in sorted(set(columns+observers)):
        T,r=rows_literal(word)
        direct=direct_state(T,r)
        require(direct==program(word)*SEED,'literal six-coordinate state '+word)
        child,_=f.a.birth(T,r)
        count=f.a.tables(child)[0][-1]
        require(count+1==(ELL*direct)[0],'literal full ready-child count '+word)
        require(f.a.strong(T),'strong legal parent '+word)
        controls.append({'word':word,'order':len(T),'parent':list(T),'root':r,
                         'state':list(map(int,direct)),'ready_count':count})
        f.a.tables.cache_clear()
    # Minimality uses all 36 concatenations as genuine future observations.
    for future in observers:
        for past in columns:
            require((ELL*program(past+future)*SEED)[0]==
                    hankel[observers.index(future),columns.index(past)],'Hankel order')
    # Similarity makes the upper blocks nonnegative: x=(e,h-e,s2-h).
    Q=sp.Matrix([[0,1,0],[1,-1,0],[-1,0,1]])
    upper={k:Q*M[k][:3,:3]*Q.inv() for k in 'EB'}
    require(upper['E']==sp.Matrix([[0,1,0],[1,1,1],[1,0,0]]),'E upper cone')
    require(upper['B']==sp.Matrix([[1,1,0],[2,2,2],[0,0,0]]),'B upper cone')
    require(M['E'][3:,3:].det()==2 and M['B'][3:,3:].det()==12,'lower determinant budget')
    lower_B=M['B'][3:,3:]
    require(all(int(x)%3==0 for x in lower_B**3),'B lower nilpotent modulo3')
    periodic=[]
    for length in range(1,9):
        mixed=0; pure=0
        for letters in itertools.product('EB',repeat=length):
            word=''.join(letters); mat=program(word); eigen=3**length
            lower_det=(mat[3:,3:]-eigen*sp.eye(3)).det()
            require(lower_det!=0,'lower target eigenvalue absent '+word)
            if 'E' in word:
                require((mat-eigen*sp.eye(6)).det()!=0,'full target eigenvalue absent '+word)
                mixed+=1
            else:
                require((mat[:3,:3]-eigen*sp.eye(3)).det()==0,'inherited parent3mode')
                pure+=1
        periodic.append({'length':length,'mixed_programs':mixed,'pure_B_programs':pure})
    # A single fixed tail inventory is a finite control, not the theorem.
    fixed=[]
    for length in range(7):
        words=[()] if length==0 else [bits+(1,) for bits in itertools.product((0,1),repeat=length-1)]
        hits=[]
        for word in words:
            mat=sp.Matrix.hstack(*[sp.Matrix(f.advance(v,word)) for v in BASIS])
            if (mat-3*sp.eye(6)).det()==0: hits.append(list(word))
        fixed.append({'tail_length':length,'tails':len(words),'eigen3_tails':hits})
    data={'status':'PASS','scope':'exact symbolic certificates plus finite controls',
          'seed':list(map(int,SEED)),'alphabet':{k:list(v) for k,v in ALPHABET.items()},
          'matrices':{k:[list(map(int,row)) for row in v.tolist()] for k,v in M.items()},
          'reach_words':columns,'observer_words':observers,
          'reach_determinant':int(reach.det()),'observer_determinant':int(observe.det()),
          'difference_words':difference_words,'difference_hankel_determinant':int(difference_hankel.det()),
          'difference_hankel':[[int(x) for x in row] for row in difference_hankel.tolist()],
          'hankel_determinant':int(hankel.det()),'hankel':[[int(x) for x in row] for row in hankel.tolist()],
          'literal_controls':controls,'periodic_controls':periodic,'fixed_tail_controls':fixed}
    (HERE/'six_state_audit.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','hankel_determinant':int(hankel.det()),'literal_controls':len(controls),'periodic_programs':sum(x['mixed_programs']+x['pure_B_programs'] for x in periodic)}))

if __name__=='__main__': main()
