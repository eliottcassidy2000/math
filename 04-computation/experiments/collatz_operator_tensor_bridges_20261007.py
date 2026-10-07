"""Exact local interfaces for the accepted operator/tensor paper premises.

This checker does not re-prove the papers, infer any universal orbit result,
submit a prize solution, or import actual moment measurements.
"""
from fractions import Fraction as F
from itertools import product
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def natural(n):
    need(type(n) is int and n >= 0, "exact natural required")


def affine(word):
    need(type(word) is tuple, "valuation tuple required")
    p,q,b=1,1,0
    for a in word:
        need(type(a) is int and a>=1, "positive integer valuation required")
        p,q,b=3*p,q*(1<<a),3*b+q
    return p,q,b


def odd_step(n, sign=1):
    need(type(n) is int and n%2==1, "odd integer required")
    need(type(sign) is int and sign in (-1,1), "sign required")
    y=3*n+sign
    a=0
    while y%2==0:
        y//=2; a+=1
    return y,a


def replay(n, word, sign=1):
    affine(word)
    states=[n]
    for a in word:
        n,actual=odd_step(n,sign)
        need(actual==a, "literal valuation guard failed")
        states.append(n)
    return tuple(states)


def mm(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def kron(a,b):
    return tuple(tuple(a[i][j]*b[k][l] for j in range(len(a[0])) for l in range(len(b[0])))
                 for i in range(len(a)) for k in range(len(b)))


def matrix(word,sign=1):
    p,q,b=affine(word)
    return ((F(p,q),F(sign*b,q)),(F(0),F(1)))


def inv(a):
    need(len(a)==2 and len(a[0])==2,"two by two matrix required")
    d=a[0][0]*a[1][1]-a[0][1]*a[1][0]
    need(d!=0,"nonsingular matrix required")
    return ((a[1][1]/d,-a[0][1]/d),(-a[1][0]/d,a[0][0]/d))


def comm(a,b):
    return mm(mm(mm(a,b),inv(a)),inv(b))


def free_reduce(w):
    out=[]
    for c in w:
        if out and out[-1]==c.swapcase():out.pop()
        else:out.append(c)
    return ''.join(out)


def inverse_word(w):
    return ''.join(c.swapcase() for c in reversed(w))


def c5_uncut(sizes, counts):
    return sum(counts[i]*counts[(i+1)%5]
               +(sizes[i]-counts[i])*(sizes[(i+1)%5]-counts[(i+1)%5])
               for i in range(5))


def main():
    checks=0
    def check(ok,label):
        nonlocal checks
        need(ok,label); checks+=1
    words=[w for r in range(5) for w in product(range(1,5),repeat=r)]
    sign_matrix=((F(-1),F(0)),(F(0),F(1)))
    for w in words:
        p,q,b=affine(w)
        residue=((q-b)*pow(p,-1,2*q))%(2*q)
        x=residue+2*q
        first=replay(x,w);second=replay(x+2*q,w)
        check(first[-1]==F(p*x+b,q),"affine carry matches literal orbit")
        check(second[-1]-first[-1]==2*p,"paired same-word separation")
        check(replay(-x,w,-1)==tuple(-v for v in first),"signed conjugacy including guards")
        m=matrix(w)
        check(mm(mm(sign_matrix,m),sign_matrix)==matrix(w,-1),"marked sign-conjugacy matrix")
        tensor=kron(m,m)
        check(tensor[0][0]==m[0][0]**2 and tensor[0][1]==m[0][0]*m[0][1],
              "full tensor retains marked slope and carry")
        check(mm(inv(m),m)==((1,0),(0,1)),"all unmarked nonsingular matrices restrict to identity")
    sample=words[:21]
    for u,v in product(sample,repeat=2):
        a,b=matrix(u),matrix(v)
        check(kron(mm(a,b),mm(a,b))==mm(kron(a,a),kron(b,b)),"ordered tensor functor")
    for n in range(1,1000,2):
        _,a=odd_step(n,1);_,b=odd_step(n,-1)
        check(min(a,b)==1 and max(a,b)>=2,"positive plus/minus initial guards differ")
    one=matrix((1,));square=mm(one,one)
    lam=one[0][0];j=lam+1/lam
    check(square[0][0]==F(9,4),"valuation-one slope square")
    check(lam**2+lam**-2==j*j-2==F(97,36),"reciprocal trace doubling")
    check(replay(7,(1,1))==(7,11,17) and replay(15,(1,1))==(15,23,35),
          "two supplied positive orbits share exactly the displayed segment")
    check(odd_step(17)[1]!=odd_step(35)[1],"shared-word exponent stops when guards split")
    check(replay(7,(1,1,2,3))[-1]==replay(3,(1,))[-1]==5,"different-word common future")
    check(odd_step(5,-1)==(7,1) and odd_step(7,-1)==(5,2),"positive minus-cycle hostile")
    check(matrix((1,))[0][0]!=matrix((2,))[0][0],"mutual restriction loses orbit slope")

    # Finite phase and weight identities from MM Proposition 3.1.  These
    # test the exact character-orthogonality index, not floating roots of unity.
    separation=0
    for size in range(1,7):
        modulus=5*size
        survivors=[]
        for g,h,u,v in product(range(1,size+1),repeat=4):
            phase=u-v+2*(g-h)
            check((phase%modulus==0)==(phase==0),"no wrapped Fourier alias")
            if phase==0:
                weight=g*g+h*u-h*h-h*v
                check(weight==(g-h)**2>=0,"legwise weights total squared label mismatch")
                if weight==0: survivors.append((g,h,u,v))
            separation+=1
        check(survivors==[(h,h,u,u) for h in range(1,size+1) for u in range(1,size+1)],
              "M blocks with a retained M-dimensional dot-product index")

    # Formal affine group, with all residue and positivity guards forgotten.
    d=((F(2),F(0)),(F(0),F(1)))
    e=((F(2,3),F(-1,3)),(F(0),F(1)))
    c=comm(d,e);conjugate=mm(mm(d,c),inv(d))
    check(c==((1,F(-1,3)),(0,1)),"affine commutator is translation")
    check(comm(c,conjugate)==((1,0),(0,1)),"commuting translation commutators")
    cw='DEde';vw='D'+cw+'d'
    relation=free_reduce(cw+vw+inverse_word(cw)+inverse_word(vw))
    check(bool(relation),"affine relation is nontrivial in the free group")
    for k in range(1,25):
        check(F(3,2)**k>1,"unbounded powers forbid unitary similarity")

    # Complete small blowup universe: coordinatewise affine rounding reduces
    # arbitrary split parts to monochromatic parts, so the symbolic minimum
    # is valid for all positive sizes; the census is only a control.
    blowups=0
    for sizes in product((1,2,3),repeat=5):
        exact=min(c5_uncut(sizes,c) for c in product(*(range(a+1) for a in sizes)))
        claimed=min(sizes[i]*sizes[(i+1)%5] for i in range(5))
        check(exact==claimed,"complete small C5 blowup bipartization")
        blowups+=1
    for bad in (lambda:affine([1]),lambda:affine((True,)),lambda:affine((0,)),
                lambda:odd_step(True),lambda:odd_step(3,True)):
        try:bad()
        except ValueError:check(True,"typed input hostile")
        else:raise ValueError("hostile accepted")
    print(json.dumps({"checks":checks,"valuation_words":len(words),"tensor_compositions":len(sample)**2,
      "shared_label_tuples":separation,"c5_blowup_size_vectors":blowups,
      "affine_kernel_free_word":relation,"slope_square":"9/4","reciprocal_trace_square":"97/36",
      "status":"Exact interface controls; cited paper theorems accepted as premises",
      "scope":"No universal orbit theorem, source floor, or prize submission."},indent=2,sort_keys=True))


if __name__=='__main__':main()
