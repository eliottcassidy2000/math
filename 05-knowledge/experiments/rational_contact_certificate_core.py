"""Standard-library-only exact polynomial models for quadratic contact proofs."""
from fractions import Fraction as F
from math import lcm

class Poly:
    def __init__(self,n,terms=None):
        self.n=n;self.terms={tuple(k):F(v) for k,v in (terms or {}).items() if v}
    @classmethod
    def scalar(cls,n,x):return cls(n,{(0,)*n:F(x)})
    def _coerce(self,x):
        if isinstance(x,Poly):
            if x.n!=self.n:raise ValueError('Polynomial variable-count mismatch')
            return x
        return Poly.scalar(self.n,x)
    def __add__(self,x):
        x=self._coerce(x);out=self.terms.copy()
        for m,c in x.terms.items():out[m]=out.get(m,F(0))+c
        return Poly(self.n,out)
    __radd__=__add__
    def __neg__(self):return Poly(self.n,{m:-c for m,c in self.terms.items()})
    def __sub__(self,x):return self+-self._coerce(x)
    def __rsub__(self,x):return self._coerce(x)+-self
    def __mul__(self,x):
        x=self._coerce(x);out={}
        for m,c in self.terms.items():
            for k,d in x.terms.items():
                mon=tuple(a+b for a,b in zip(m,k));out[mon]=out.get(mon,F(0))+c*d
        return Poly(self.n,out)
    __rmul__=__mul__
    def __truediv__(self,x):return self*F(1,x)
    def __pow__(self,k):
        if not isinstance(k,int) or k<0:raise ValueError('Nonnegative integer exponent required')
        p=Poly.scalar(self.n,1);b=self
        while k:
            if k%2:p=p*b
            b=b*b;k//=2
        return p
    def __eq__(self,x):return self.terms==self._coerce(x).terms
    def __bool__(self):return bool(self.terms)
    def coeff(self,index,k):
        return Poly(self.n-1,{m[:index]+m[index+1:]:c for m,c in self.terms.items() if m[index]==k})
    def degree(self,index):return max((m[index] for m in self.terms),default=-1)
    def clear_denoms(self):
        d=1
        for c in self.terms.values():d=lcm(d,c.denominator)
        return self*d
    def substitute(self,values):
        out={}
        for mon,c in self.terms.items():
            new=list(mon)
            for i,x in values.items():c*=F(x)**mon[i];new[i]=0
            new=tuple(new);out[new]=out.get(new,F(0))+c
        return Poly(self.n,out)

def variables(n):
    ans=[]
    for i in range(n):
        m=[0]*n;m[i]=1;ans.append(Poly(n,{tuple(m):F(1)}))
    return tuple(ans)

def decode(data,n=5):return Poly(n,{tuple(mon):F(c) for mon,c in data})

def lift_t(p,k):return Poly(p.n+1,{(k,)+m:c for m,c in p.terms.items()})

def quotient_t2plus1(p):
    d=p.degree(0);rows=[p.coeff(0,j) for j in range(d+1)]
    quotient=Poly(p.n)
    for j in range(d,1,-1):
        q=rows[j];quotient+=lift_t(q,j-2);rows[j]=Poly(p.n-1);rows[j-2]-=q
    if rows[0] or rows[1]:raise RuntimeError('Expected exact divisibility by t²+1')
    return quotient

def system(pattern):
    """Return displayed necessary equations in order (a,b,d,v,m)."""
    t,a,b,d,v,m=variables(6);q=a*a+b*b
    if pattern in ('adjacent_plus','adjacent_minus'):
        if pattern=='adjacent_plus':
            c=-m*(q-1)/4;rho=1+m*a;RS=t*(t-1);RT=t+1
        else:
            c=m*a/2;rho=1+m*(q-1)/2;RS=t*(t+1);RT=t-1
        QA=t*t-2*c*t+rho;QB=t*t-2*a*t+q
        X=t*t-(a+c)*t+a*c-b*d;Y=-(b+d)*t+a*d+b*c
        E=4*(t*t+1)*(X*X-Y*Y+v*X*Y)+m*m*RT*QB*QB-4*RS*QA*QA
        if E.degree(0)>5:raise RuntimeError('Leading-degree cancellation failed')
        out=[E.coeff(0,j).clear_denoms() for j in range(6)]
        out += [(16*(d*d-rho+c*c)).coeff(0,0).clear_denoms(),(v*v-2*m*m+4).coeff(0,0)]
    elif pattern=='opposite':
        c=m*(q-1+2*a)/8;rho=1+m*(q-1-2*a)/4
        QA=t*t-2*a*t+q;QB=t*t-2*c*t+rho
        X=t*t-(a+c)*t+a*c-b*d;Y=-(b+d)*t+a*d+b*c
        N=4*(t*t-1)*QB*QB-m*m*t*QA*QA
        J=quotient_t2plus1(N);E=J-4*(X*X-Y*Y+v*X*Y)
        if E.degree(0)!=3:raise RuntimeError('Opposite polynomial has wrong degree')
        out=[16*E.coeff(0,j) for j in [3,2,1,0]]
        out += [(64*(d*d-rho+c*c)).coeff(0,0),(v*v-2*m*m+4).coeff(0,0)]
    else:raise ValueError(pattern)
    return out
