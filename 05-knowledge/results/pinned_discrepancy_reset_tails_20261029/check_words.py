"""Independent integer-cylinder and coherent-attachment controls."""
from fractions import Fraction as F
from itertools import product
import json
import numpy as np
from check_resets import reset, require


def compile_word(word):
    A, B, D = 1, 0, 1
    for a in word:
        A, B, D = 3*A, 3*B+((1 << (a-1))-1)*D, (1 << a)*D
    return A, B, D, (-B*pow(A, -1, D)) % D


def literal_word(j, length):
    result = []
    for _ in range(length):
        n = 2*j-1
        t = 3*n+1
        a = (t & -t).bit_length()-1
        result.append(a)
        j = ((t >> a)+1)//2
    return tuple(result), j


def word_controls():
    alphabet = [(k,a) for k in range(3) for a in range(2,5)]
    rows = []
    tested = 0
    for b in range(1,4):
        coeffs = {F(1): F(0), F(3,4): F(0)}
        for blocks in product(alphabet, repeat=b):
            word = tuple(a for k,last in blocks for a in [1]*k+[last])
            A,B,D,r = compile_word(word)
            for t in (1,2,17):
                j = r+D*t
                found,end = literal_word(j,len(word))
                require(found==word and end*D==A*j+B, 'literal source/cylinder/endpoint')
                changed,_ = literal_word(j+D//2,len(word))
                require(changed!=word, 'one erased terminal bit must not preserve the word')
            J = 10**18
            threshold = (D*J-B)//A
            first = r+D*max(0,(threshold-r)//D+1)
            require(first>0 and A*first+B>D*J, 'first AP tail source')
            # Decreasing-function integral bounds, with w(j) between (j+1)^-2 and j^-2.
            lower = F(J,D*(first+1))
            upper = F(J,first*first)+F(J,D*first)
            leading = F(A,D*D)
            require(lower<=upper and lower-F(1,10**10)<=leading<=upper+F(1,10**10),
                    'fixed-word leading coefficient survives exact tail bounds')
            for z in coeffs:
                coeffs[z] += z**len(word)*leading
            tested += 1
        for z,actual in coeffs.items():
            one = sum(z**(k+1)*F(3**(k+1),4**(k+a)) for k,a in alphabet)
            require(actual==one**b, 'finite word-bank coefficient tensor product')
            rows.append({'blocks':b,'z':str(z),'finite_bank_coefficient':str(actual),
                         'all_word_lower_coefficient':str((z/(4-3*z))**b)})
    native = 0
    for source in range(2,10002):
        j=source; word=[]
        for b in range(1,5):
            if j==1: break
            k,_,a,j=reset(j)
            word.extend([1]*k+[a])
            A,B,D,r=compile_word(word)
            require(source % D==r and A*source+B==D*j, 'repeated native reset source preservation')
            native += 1
    # A single positively weighted source already forces a finite-scale tail cost.
    j=14; receipt=[]
    for b in range(10):
        k,_,a,j=reset(j)
        receipt.append({'block':b+1,'k':k,'a':a,'endpoint_j':j})
    require(j==456 and F(j-1,14*15)==F(13,6),'original source27 tail obstruction')
    return {'bounded_bank_words':tested,'native_block_prefixes':native,'coefficients':rows,
            'source27_receipt':receipt,'all_cutoff_envelope_constant_lower_bound':'13/6'}


def ramified_controls():
    p=5
    points=list(product(range(p),repeat=2))
    def g(v):
        x,y=v
        return ((2*x+3*y)%p,(3*x+5*y)%p)
    left=set(points)-{(0,0)}
    symbol={(0,0):0}
    orbit_pairs=0
    while left:
        v=min(left); orbit=set(); t=v
        while t not in orbit:
            orbit.add(t); t=g(t)
        negative={((-x)%p,(-y)%p) for x,y in orbit}
        require(not orbit & negative,'odd ramified antipodal packets')
        for t in orbit: symbol[t]=1
        for t in negative: symbol[t]=-1
        left-=orbit|negative
        orbit_pairs += 1
    A=np.array([[symbol[((v[0]-u[0])%p,(v[1]-u[1])%p)] for v in points] for u in points],dtype=np.int64)
    require(np.all(A+A.T==0) and np.all(A.sum(axis=1)==0),'regular kernel tournament')
    # Order each Z/25 pair by its Z/5 quotient, then its actual kernel coordinate.
    lifts=[(x+p*s,y+p*t) for x,y in points for s,t in points]
    def lifted_symbol(dx,dy):
        if dx%p or dy%p: return symbol[(dx%p,dy%p)]
        return symbol[(dx//p,dy//p)]
    K=np.array([[lifted_symbol((v[0]-u[0])%25,(v[1]-u[1])%25) for v in lifts]
                for u in lifts],dtype=np.int64)
    tensor=np.kron(A,np.ones((25,25),dtype=np.int64))+np.kron(np.eye(25,dtype=np.int64),A)
    require(np.array_equal(K,tensor),'actual ramified difference matrix equals lexicographic product')
    P=np.kron(np.eye(25,dtype=np.int64),np.ones((25,1),dtype=np.int64))
    require(np.array_equal(K@P,25*P@A),'fixed quotient mode amplifies by 25')
    # An integer basis for the perpendicular fibre space has the kernel action.
    V=np.vstack([np.eye(24,dtype=np.int64),-np.ones((1,24),dtype=np.int64)])
    E=np.kron(np.eye(25,dtype=np.int64),V)
    require(np.array_equal(K@E,np.kron(np.eye(25,dtype=np.int64),A@V)),
            'zero-sum fibre space retains only the kernel action')
    return {'p':p,'base_vertices':25,'lift_vertices':625,'orbit_pairs':orbit_pairs,
            'quotient_mode_multiplier':25,'literal_tensor_identity':True,
            'perpendicular_fibre_identity':True}


if __name__=='__main__':
    print(json.dumps({'words':word_controls(),'ramified':ramified_controls()},indent=2,sort_keys=True))
