"""Physical C2 families from any rooted source with exactly one starting path."""
from core import birth,glue,gadget,constants,require

def seeds(V,b):
 # U has a sole out-neighbor b. In F, a is the first and y the last clone.
 U,a=glue(V,{b},[0]);F,y=glue(U,{b},[0]);M,_=glue(F,{a},[0])
 # The transitive Dcell has entry n (last) and exit n+1 (first).
 W,_=glue(V,{b},[0,5,1]);P,r=glue(W,{len(V)+1},[0])
 return {'N':(F,a),'M':(M,y),'P':(P,r),'Q':(F,a)}

def extend_retained(T,r,v):
 require(v>=1,'positive suffix parameter')
 T=list(T)
 if v==1:return T,r
 for c in bin(v-1)[2:][::-1]:T,_=birth(T,[r],int(c),1)
 return T,r

def construct(kind,V,r,v):
 require(kind in ('N','M','P','Q'),'family name')
 require(v>=(2 if kind=='N'else 1),'suffix domain')
 require(kind=='P'or len(V)>=2,'adjacent source must be nontrivial')
 T,a=seeds(V,r)[kind]
 return glue(T,{a},gadget(v))if kind=='Q'else extend_retained(T,a,v)

def scalar(kind,s,v):
 D,E,U,d=s;g,_,_=constants(v);t=v%2
 if kind=='N':return (D-E+U)*g+(D+3*U-d)*v-U-E-t*(D-E+d)
 if kind=='M':return (D+2*U)*g+(2*D-E+6*U-d)*v-E-U-t*(U+d)
 if kind=='P':return (D+4*U)*g+(4*D-E+6*U)*v-t*U
 require(kind=='Q','family name')
 return (D-E+3*U-2*d)*g+(2*D-E+4*U-d)*v-t*(U-d)

def seed_profiles(s):
 D,E,U,d=s
 return {'N':(D-E+3*U-2*d,U-d,2*D-E+4*U-d,D-E+U),
         'M':(3*D-2*E+6*U-2*d,D+U-d,3*D-E+8*U-d,D+2*U),
         'P':(5*D-E+9*U,D+3*U,5*D-E+10*U,D+4*U)}

def parameter_limit(kind,s,target):
 D,E,U,d=s
 if kind=='N':return (target+E+d)//(D+3*U-d)
 if kind=='M':return (target-D+E+d)//(2*D-E+6*U-d)
 if kind=='P':return (target-D-3*U)//(4*D-E+6*U)
 return target//(2*D-E+4*U-d)
