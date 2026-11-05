"""Injective alternative three-observable encoding of physical source responses."""
from dependencies.verify_weak_start_one import oracle

def signature(s):
 D,E,U,d=s;return D+2*d,E+d,U-d

def recover(A,B,C):
 """Only a candidate until membership is checked. Includes the singleton."""
 if min(A,B)<1 or C<0:return None
 if B%2:d=B if A>=3*B else A-B-C
 else:d=B//2 if 2*A<=3*B+2*C else A-B-C
 D,E,U=A-2*d,B-d,C+d
 if D<1 or d<1 or E<0 or E>D or U<D or d>D or D%2==0 or d%2==0:return None
 return D,E,U,d

def member(A,B,C):
 """Exact membership, returning the full quadruple and a literal source recipe."""
 s=recover(A,B,C)
 if s is None:return None
 D,E,U,d=s;witnesses=oracle(D,E,U)
 if d not in witnesses:return None
 return s,witnesses[d]

def n_value(s,t):
 D,E,U,d=s;F,G,V,f=t
 return (D-E)*(F-G)+D*V+U*F+3*U*V-d*(G+V)-f*(E+U)

def n_signature(x,y):
 A,B,C=x;P,Q,R=y
 return (A-B)*(P-Q)+A*R+C*P+3*C*R

def m_value(s,t):
 D,E,U,d=s;F,G,V,f=t
 return D*F+2*(D*V+U*F)+6*U*V-E*V-U*G-d*(G+V)-f*(E+U)
