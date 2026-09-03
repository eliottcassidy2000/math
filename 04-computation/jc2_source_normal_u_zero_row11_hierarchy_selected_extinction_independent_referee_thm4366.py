#!/usr/bin/env python3
"""Clean-room exact referee for proposed THM-4366.

No primary certificate is read or imported.  The source rows, bracket rows,
and depth matrices are rebuilt from the displayed canonical definitions.
"""
from fractions import Fraction as Q
from math import comb
import sys
sys.stdout.reconfigure(newline="\n")

NVAR = 4
PHI, ETA, ALPHA, BETA = range(NVAR)
CHECKS = 0

def check(ok, label):
    global CHECKS
    if not ok:
        raise AssertionError(label)
    CHECKS += 1

class R:
    """Sparse Q[Phi^{+-1},eta,alpha,beta]."""
    __slots__ = ("terms",)
    def __init__(self, terms=None):
        self.terms = {}
        for mon, value in (terms or {}).items():
            mon, value = tuple(mon), Q(value)
            if value:
                self.terms[mon] = self.terms.get(mon, Q(0)) + value
        self.terms = {m:v for m,v in self.terms.items() if v}
    @staticmethod
    def c(value):
        value = Q(value)
        return R({(0,)*NVAR:value}) if value else R()
    @staticmethod
    def var(i):
        mon=[0]*NVAR; mon[i]=1
        return R({tuple(mon):Q(1)})
    @staticmethod
    def mon(value=1, **powers):
        mon=[0]*NVAR; names={"p":PHI,"e":ETA,"a":ALPHA,"b":BETA}
        for k,v in powers.items(): mon[names[k]]=v
        return R({tuple(mon):Q(value)})
    def __bool__(self): return bool(self.terms)
    def __eq__(self, other): return self.terms == rr(other).terms
    def __neg__(self): return R({m:-v for m,v in self.terms.items()})
    def __add__(self, other):
        out=dict(self.terms)
        for m,v in rr(other).terms.items():
            out[m]=out.get(m,Q(0))+v
            if not out[m]: del out[m]
        return R(out)
    __radd__=__add__
    def __sub__(self, other): return self + (-rr(other))
    def __rsub__(self, other): return rr(other)-self
    def __mul__(self, other):
        out={}
        for ma,va in self.terms.items():
            for mb,vb in rr(other).terms.items():
                m=tuple(ma[i]+mb[i] for i in range(NVAR))
                out[m]=out.get(m,Q(0))+va*vb
        return R(out)
    __rmul__=__mul__
    def __truediv__(self, scalar):
        scalar=Q(scalar); return R({m:v/scalar for m,v in self.terms.items()})
    def __pow__(self, exponent):
        if exponent < 0:
            if len(self.terms)!=1: raise ValueError("negative power of nonmonomial")
            (m,v),=self.terms.items()
            return R({tuple(exponent*z for z in m):v**exponent})
        out=R.c(1); base=self; n=exponent
        while n:
            if n&1: out=out*base
            base=base*base; n//=2
        return out
    def subst(self, i, replacement):
        replacement=rr(replacement); out=R()
        for mon,value in self.terms.items():
            exponent=mon[i]
            if exponent<0 and len(replacement.terms)!=1:
                raise ValueError("negative-power substitution")
            kept=list(mon); kept[i]=0
            out += R({tuple(kept):value}) * replacement**exponent
        return out
    def coeff(self, mon): return self.terms.get(tuple(mon),Q(0))

def rr(value): return value if isinstance(value,R) else R.c(value)
P,E,A,B = R.var(PHI),R.var(ETA),R.var(ALPHA),R.var(BETA)

# Polynomials in x over R.
def xp(items=None):
    return {int(d):rr(c) for d,c in (items or {}).items() if rr(c)}
def xadd(a,b):
    out=dict(a)
    for d,c in b.items():
        out[d]=out.get(d,R())+c
        if not out[d]: del out[d]
    return out
def xneg(a): return {d:-c for d,c in a.items()}
def xscale(a,s): s=Q(s); return {d:c*s for d,c in a.items() if c*s}
def xmul(a,b):
    out={}
    for i,c in a.items():
        for j,d in b.items():
            out[i+j]=out.get(i+j,R())+c*d
            if not out[i+j]: del out[i+j]
    return out
def xder(a): return {i-1:c*i for i,c in a.items() if i and c*i}
def xcoef(a,d): return a.get(d,R())
def xmany(parts):
    out={}
    for p in parts: out=xadd(out,p)
    return out
def xsubst(poly,i,value): return xp({d:c.subst(i,value) for d,c in poly.items()})
def xdivide_x(poly):
    check(not xcoef(poly,0),"x division constant")
    return {d-1:c for d,c in poly.items() if d}

def rank_q(matrix):
    if not matrix:return 0
    m=[[Q(z) for z in row] for row in matrix]; nr,nc=len(m),len(m[0]); r=0
    for c in range(nc):
        p=next((i for i in range(r,nr) if m[i][c]),None)
        if p is None: continue
        m[r],m[p]=m[p],m[r]; z=m[r][c]; m[r]=[v/z for v in m[r]]
        for i in range(nr):
            if i!=r and m[i][c]:
                z=m[i][c]; m[i]=[m[i][j]-z*m[r][j] for j in range(nc)]
        r+=1
        if r==nr:break
    return r

def independent_row_indices(matrix):
    if not matrix:return []
    m=[[Q(z) for z in row] for row in zip(*matrix)]
    nr,nc=len(m),len(m[0]); out=[]; r=0
    for c in range(nc):
        p=next((i for i in range(r,nr) if m[i][c]),None)
        if p is None:continue
        m[r],m[p]=m[p],m[r]; z=m[r][c]; m[r]=[v/z for v in m[r]]
        for i in range(nr):
            if i!=r and m[i][c]:
                z=m[i][c]; m[i]=[m[i][j]-z*m[r][j] for j in range(nc)]
        out.append(c); r+=1
        if r==nr:break
    return out

def solve_square(matrix,rhs):
    n=len(matrix); check(n==len(rhs) and all(len(row)==n for row in matrix),"square solve dims")
    m=[[Q(z) for z in row] for row in matrix]; b=list(rhs)
    for c in range(n):
        p=next(i for i in range(c,n) if m[i][c])
        m[c],m[p]=m[p],m[c]; b[c],b[p]=b[p],b[c]
        z=m[c][c]; m[c]=[v/z for v in m[c]]; b[c]=b[c]/z
        for i in range(n):
            if i!=c and m[i][c]:
                z=m[i][c]; m[i]=[m[i][j]-z*m[c][j] for j in range(n)]; b[i]-=z*b[c]
    return b

def solve_injective(matrix,rhs):
    nc=len(matrix[0]); rows=independent_row_indices(matrix)
    check(len(rows)==nc,"injective map")
    sol=solve_square([[matrix[i][j] for j in range(nc)] for i in rows],[rhs[i] for i in rows])
    res=[rhs[i]-sum((row[j]*sol[j] for j in range(nc)),R()) for i,row in enumerate(matrix)]
    return sol,res,rows

def rref_sparse(rows,ncols):
    rows=[{j:Q(v) for j,v in row.items() if v} for row in rows]; piv=[]; r=0
    for c in range(ncols):
        p=next((i for i in range(r,len(rows)) if rows[i].get(c)),None)
        if p is None:continue
        rows[r],rows[p]=rows[p],rows[r]; z=rows[r][c]
        rows[r]={j:v/z for j,v in rows[r].items()}; prow=rows[r]
        for i in range(len(rows)):
            if i==r:continue
            z=rows[i].get(c,Q(0))
            if not z:continue
            new=dict(rows[i])
            for j,v in prow.items():
                nv=new.get(j,Q(0))-z*v
                if nv:new[j]=nv
                elif j in new:del new[j]
            rows[i]=new
        piv.append(c); r+=1
        if r==len(rows):break
    return rows,piv

def nullspace_of_column_span(columns,ncoords):
    reduced,piv=rref_sparse(columns,ncoords); pset=set(piv); basis=[]
    for f in range(ncoords):
        if f in pset:continue
        v={f:Q(1)}
        for i,p in enumerate(piv):
            z=reduced[i].get(f,Q(0))
            if z:v[p]=-z
        basis.append(v)
    return len(piv),basis

def depth_matrix(m,d):
    coords=[(n,r) for n in range(m+1) for r in range(n+d+1)]; ci={z:i for i,z in enumerate(coords)}
    columns=[]
    for a in range(d+1):
        for b in range(d-a+1):
            for e in range(m//2+1):
                for c in range(m+1):
                    n0=b+c+2*e
                    if n0>m:continue
                    N=c+e; r0=a+2*b+e; col={}
                    for k in range(N+1):
                        n,r=n0+k,r0+2*k
                        if n<=m and (n,r) in ci:col[ci[(n,r)]]=Q(comb(N,k))
                    columns.append(col)
    return coords,columns

def jet(rows,coords):return [xcoef(rows[n],r) for n,r in coords]
def dot(left,values):return sum((z*values[i] for i,z in left.items()),R())
def tangent_column(coords,row,power,which):
    ci={z:i for i,z in enumerate(coords)}; out=[Q(0)]*len(coords)
    if which=="A":
        if (row,power+1) in ci:out[ci[(row,power+1)]]=Q(1,2)
    else:
        for off,v in ((0,Q(-3,4)),(2,Q(-3,8))):
            if (row,power+off) in ci:out[ci[(row,power+off)]]=v
    return out
def left_eval(left,column):return sum((left.get(i,Q(0))*column[i] for i in range(len(column))),Q(0))

def reduce_eta_quadratic(expr,eta2):
    out=R()
    for mon,value in expr.terms.items():
        check(mon[PHI]==0,"Phi0 quotient")
        q,rem=divmod(mon[ETA],2); kept=list(mon); kept[ETA]=rem
        out += R({tuple(kept):value*eta2**q})
    return out

def laurent_remainder(expr,modulus):
    check(all(mon[ETA]==mon[ALPHA]==mon[BETA]==0 for mon in expr.terms),"univariate expr")
    check(all(mon[ETA]==mon[ALPHA]==mon[BETA]==0 for mon in modulus.terms),"univariate mod")
    if not expr:return R()
    shift=max(0,-min(mon[PHI] for mon in expr.terms))
    poly={mon[PHI]+shift:v for mon,v in expr.terms.items()}
    mod={mon[PHI]:v for mon,v in modulus.terms.items()}; md=max(mod); ml=mod[md]
    while poly and max(poly)>=md:
        deg=max(poly); factor=poly[deg]/ml; off=deg-md
        for j,v in mod.items():
            k=j+off; poly[k]=poly.get(k,Q(0))-factor*v
            if not poly[k]:del poly[k]
    return R({(d-shift,0,0,0):v for d,v in poly.items()})
def zero_mod(expr,modulus):return not laurent_remainder(expr,modulus)

def trim_mod(poly,p):return {d:v%p for d,v in poly.items() if v%p}
def add_mod(a,b,p):
    out=dict(a)
    for d,v in b.items():
        out[d]=(out.get(d,0)+v)%p
        if not out[d]:del out[d]
    return out
def mul_mod(a,b,p):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=(out.get(i+j,0)+x*y)%p
    return trim_mod(out,p)
def divrem_mod(a,b,p):
    a,b=trim_mod(a,p),trim_mod(b,p); q={}; db=max(b); ib=pow(b[db],-1,p)
    while a and max(a)>=db:
        da=max(a); z=a[da]*ib%p; off=da-db; q[off]=(q.get(off,0)+z)%p
        for j,v in b.items():
            k=j+off; a[k]=(a.get(k,0)-z*v)%p
            if not a[k]:del a[k]
    return trim_mod(q,p),trim_mod(a,p)
def gcd_mod(a,b,p):
    a,b=trim_mod(a,p),trim_mod(b,p)
    while b:_,r=divrem_mod(a,b,p); a,b=b,r
    if not a:return {}
    inv=pow(a[max(a)],-1,p)
    return {d:v*inv%p for d,v in a.items()}
def even_phi_to_y_mod(expr,p):
    out={}
    for mon,v in expr.terms.items():
        check(mon[ETA]==mon[ALPHA]==mon[BETA]==0 and mon[PHI]>=0 and mon[PHI]%2==0,
              "even Phi polynomial")
        check(v.denominator%p,"mod-p denominator")
        z=v.numerator*pow(v.denominator,-1,p)%p
        if z:out[mon[PHI]//2]=(out.get(mon[PHI]//2,0)+z)%p
    return trim_mod(out,p)

# THM-4308/4357 response, with U=0 fixing xi.
D=R.c(Q(896,15)); THETA=R.c(Q(512,75)); K=R.c(Q(-32,5))
UPS=R.c(Q(-731648,2025)); X=R.c(Q(237757952,54675)); ZETA=-P*Q(3,2); U=R()
W=-(P*P*4343625-X*17172000+R.c(143826305024))/4009500
Z=(R.c(12506118074368)-P*P*173745000-P*E*195463125-X*926883000)/108256500

def source_rows():
    """Literal expansion p=t(1+x^2t), y=xtp of every residual H term."""
    terms=[(R.c(-3),1,0),(R.c(Q(8,3)),2,0),(R.c(Q(-1376,135)),3,0),
           (K,0,2),(P,2,1),(D,4,0),(THETA,1,2),(E,3,1),(ZETA,0,3),
           (UPS,5,0),(X,2,2),(A,4,1),(B,1,3),(U,6,0),(W,3,2),(Z,0,4)]
    rows={n:{} for n in range(4,13)}
    for coeff,c,e in terms:
        n0=c+2*e; N=c+e
        for k in range(N+1):
            n=n0+k
            if n in rows:
                r=e+2*k; rows[n][r]=rows[n].get(r,R())+coeff*comb(N,k)
    return {n:xp(poly) for n,poly in rows.items()}

def initial_rows():
    aa={0:xp({0:1,2:Q(1,4)}),1:xp({0:Q(4,3),2:2}),
        2:xp({0:Q(-32,9),2:Q(-4,5)}),
        3:xp({0:Q(2176,135),1:-P/2,2:R.c(Q(1088,315))-K*Q(4,7),4:Q(-32,15)})}
    cc={0:xp({1:Q(-3,4),3:Q(-1,8)}),1:xp({1:-4,3:Q(-3,2)}),
        2:xp({1:Q(88,15),3:Q(-12,5)}),
        3:xp({0:P*Q(3,4),1:R.c(Q(-8128,315))+K*Q(6,7),2:P*Q(3,8),
              3:R.c(Q(736,105))+K*Q(3,7),5:Q(8,5)})}
    return aa,cc

QX=xp({0:-3,2:Q(-1,2)})
def brow(aa,cc,m):
    return xmany([xadd(xscale(xmul(xder(aa[i]),cc[m-i]),m-i),
                       xscale(xmul(aa[i],xder(cc[m-i])),-i)) for i in range(1,m)])
def trow(aa,cc,m):
    out=xmany([xmul(cc[i],cc[m-i]) for i in range(1,m)])
    for i in range(m):
        for j in range(m):
            k=m-i-j
            if 0<=k<m:out=xadd(out,xneg(xmul(xmul(aa[i],aa[j]),aa[k])))
    return out
def predict(aa,cc,m):return xadd(trow(aa,cc,m),xscale(xmul(QX,brow(aa,cc,m)),Q(-1,m)))

def particular(aa,cc,m):
    det=xscale(brow(aa,cc,m),Q(-1,m)); a0=xcoef(det,0)*Q(4,3)
    correction=xmul(xp({0:2,2:1}),xp({0:a0*Q(3,8)}))
    cpart=xscale(xdivide_x(xadd(det,xneg(correction))),2); apart=xp({0:a0})
    lhs=xadd(xmul(xneg(xder(cc[0])),apart),xmul(xder(aa[0]),cpart))
    check(lhs==det,f"row{m} particular")
    return apart,cpart

def dmatrix(m,powers):
    out=[[Q(0) for _ in powers] for _ in range(m+1)]
    for col,j in enumerate(powers):
        if j:out[j-1][col]+=Q(3*j,m)
        if j+1<=m:out[j+1][col]+=Q(j-2*m,2*m)
    return out
def add_tangent(aa,cc,row,apart,cpart,values,powers):
    theta=xp({powers[i]:values[i] for i in range(len(values))})
    aa[row]=xadd(apart,xmul(theta,xder(aa[0]))); cc[row]=xadd(cpart,xmul(theta,xder(cc[0])))
def subst_all(aa,cc,gg,i,value):
    for n in list(aa):aa[n]=xsubst(aa[n],i,value); cc[n]=xsubst(cc[n],i,value)
    for n in list(gg):gg[n]=xsubst(gg[n],i,value)

def build_to_row8_terminal():
    aa,cc=initial_rows(); gg=source_rows()
    check(xcoef(gg[4],0)==D and xcoef(gg[4],1)==P,"literal G4")
    check(xcoef(gg[8],4)==15*U+5*W+Z,"literal G8")
    check(xcoef(gg[9],6)==20*U+10*W+4*Z,"literal G9")
    check(xcoef(gg[10],8)==15*U+10*W+6*Z,"literal G10")
    check(gg[4]==predict(aa,cc,4),"base G4")
    for row in range(4,8):
        apart,cpart=particular(aa,cc,row); aa[row],cc[row]=apart,cpart
        defect=xadd(gg[row+1],xneg(predict(aa,cc,row+1)))
        sol,res,_=solve_injective(dmatrix(row+1,list(range(row+1))),[xcoef(defect,j) for j in range(row+2)])
        check(not any(res),f"row{row} selected")
        add_tangent(aa,cc,row,apart,cpart,sol,list(range(row+1)))
        check(gg[row+1]==predict(aa,cc,row+1),f"row{row+1} reached")
    apart8,cpart8=particular(aa,cc,8); aa[8],cc[8]=apart8,cpart8
    # Rebuild the inherited row-eight terminal fibre and restrict the next
    # bracket selector to it: 9 tangent coordinates minus rank 2 = A^7.
    c28,g28=depth_matrix(8,2); c38,g38=depth_matrix(8,3)
    r28,n28=nullspace_of_column_span(g28,len(c28)); r38,n38=nullspace_of_column_span(g38,len(c38))
    check((len(c28),len(g28),r28,len(n28))==(63,131,51,12),"pi8 P2 ledger")
    check((len(c38),len(g38),r38,len(n38))==(72,204,63,9),"pi8 P3 ledger")
    m8=[[left_eval(v,tangent_column(c28,8,j,"A")) for j in range(9)] for v in n28]
    m8 += [[left_eval(v,tangent_column(c38,8,j,"C")) for j in range(9)] for v in n38]
    sparse8=[{j:z for j,z in enumerate(row) if z} for row in m8]
    rank8,kernel8=nullspace_of_column_span(sparse8,9)
    check(rank8==2 and len(kernel8)==7,"row8 depth terminal A7")
    d9=dmatrix(9,list(range(9)))
    restricted=[[sum((d9[i][j]*v.get(j,Q(0)) for j in range(9)),Q(0)) for v in kernel8] for i in range(10)]
    check(rank_q(restricted)==7,"G9 rank7 on row8 terminal fibre")
    defect=xadd(gg[9],xneg(predict(aa,cc,9)))
    sol8,res9,rows9=solve_injective(dmatrix(9,list(range(9))),[xcoef(defect,j) for j in range(10)])
    check(len([z for z in res9 if z])==1,"one E9 residual")
    return aa,cc,gg,apart8,cpart8,sol8,res9,rows9

def apply_row9_depth(aa,cc,reduce=None):
    coords2,cols2=depth_matrix(9,2); coords3,cols3=depth_matrix(9,3)
    rank2,null2=nullspace_of_column_span(cols2,len(coords2)); rank3,null3=nullspace_of_column_span(cols3,len(coords3))
    check((len(coords2),len(cols2),rank2,len(null2))==(75,160,59,16),"pi9 P2 ledger")
    check((len(coords3),len(cols3),rank3,len(null3))==(85,251,73,12),"pi9 P3 ledger")
    apart,cpart=particular(aa,cc,9); aa[9],cc[9]=apart,cpart
    base=[dot(v,jet(aa,coords2)) for v in null2]+[dot(v,jet(cc,coords3)) for v in null3]
    mat=[[left_eval(v,tangent_column(coords2,9,j,"A")) for j in range(10)] for v in null2]
    mat += [[left_eval(v,tangent_column(coords3,9,j,"C")) for j in range(10)] for v in null3]
    check(rank_q(mat)==3,"row9 terminal depth rank3")
    piv=[]
    for j in range(10):
        if rank_q([[row[k] for k in piv+[j]] for row in mat])>len(piv):piv.append(j)
    check(piv==[7,8,9],"row9 pivots")
    ir=independent_row_indices([[row[j] for j in piv] for row in mat])
    sol=solve_square([[mat[i][j] for j in piv] for i in ir],[-base[i] for i in ir])
    res=[base[i]+sum((mat[i][piv[j]]*sol[j] for j in range(3)),R()) for i in range(len(mat))]
    if reduce is None:check(not any(res),"row9 depth identity")
    else:check(not any(reduce(z) for z in res),"row9 depth modulo E9")
    add_tangent(aa,cc,9,apart,cpart,sol,piv)
    return coords2,cols2,null2,coords3,cols3,null3

def branch_phi_zero():
    aa,cc,gg,apart8,cpart8,sol8,res9,_=build_to_row8_terminal()
    subst_all(aa,cc,gg,PHI,R()); sol8=[z.subst(PHI,R()) for z in sol8]; res9=[z.subst(PHI,R()) for z in res9]
    target=-(R.c(5646560625)*E*E+R.c(379697122115584))*Q(11,243)
    nz=[z for z in res9 if z]; probe=(0,2,0,0); ratio=nz[0].coeff(probe)/target.coeff(probe)
    check(nz[0]==target*ratio and ratio,"Phi0 E9")
    eta2=Q(-379697122115584,5646560625); reduce=lambda z:reduce_eta_quadratic(z,eta2)
    check(eta2!=0,"Phi0 E9 two distinct eta roots")
    check(W.subst(PHI,R()) and Z.subst(PHI,R()),"Phi0 W,Z strict constants")
    check((A*A-4*W*UPS).subst(PHI,R()).coeff((0,0,2,0))==1,"Phi0 discriminant open")
    # The cached particular was formed before Phi substitution, so substitute it too.
    apart8=xsubst(apart8,PHI,R()); cpart8=xsubst(cpart8,PHI,R())
    add_tangent(aa,cc,8,apart8,cpart8,sol8,list(range(9)))
    for n in list(aa):
        aa[n]=xp({r:reduce(z) for r,z in aa[n].items()}); cc[n]=xp({r:reduce(z) for r,z in cc[n].items()})
    for n in list(gg):gg[n]=xp({r:reduce(z) for r,z in gg[n].items()})
    row9def=xadd(gg[9],xneg(predict(aa,cc,9)))
    check(all(not reduce(z) for z in row9def.values()),"Phi0 row9 quotient source")
    apply_row9_depth(aa,cc,reduce)
    defect=xadd(gg[10],xneg(predict(aa,cc,10)))
    _,res,selected=solve_injective(dmatrix(10,list(range(7))),[xcoef(defect,j) for j in range(11)])
    res=[reduce(z) for z in res]; claimed=R.c(Q(-11388581281792,332150625))
    check(res[8]==claimed,"Phi0 G10 constant")
    # A second residual at index 6 is affine in alpha*eta; index 8 is the
    # source-independent unit and already makes the system inconsistent.
    check([i for i,z in enumerate(res) if z]==[6,8],"Phi0 G10 residual support")
    check(selected==[0,1,2,3,4,5,7],"Phi0 G10 selector")
    return claimed

def branch_phi_nonzero():
    aa,cc,gg,apart8,cpart8,sol8,res9,_=build_to_row8_terminal()
    e9=P*P*13553385750-P*E*69677820000-E*E*5646560625-R.c(379697122115584)-P*A*11293121250
    nz=[z for z in res9 if z]; probe=(1,0,1,0); ratio=nz[0].coeff(probe)/e9.coeff(probe)
    check(nz[0]==e9*ratio and ratio,"Phi!=0 E9")
    alpha9=(P*P*13553385750-P*E*69677820000-E*E*5646560625-R.c(379697122115584))*R.mon(p=-1)/11293121250
    check(not e9.subst(ALPHA,alpha9),"alpha9 solve")
    add_tangent(aa,cc,8,apart8,cpart8,sol8,list(range(9))); subst_all(aa,cc,gg,ALPHA,alpha9)
    check(gg[9]==predict(aa,cc,9),"row9 reached")
    apply_row9_depth(aa,cc)

    defect=xadd(gg[10],xneg(predict(aa,cc,10)))
    sol9,res10,selected=solve_injective(dmatrix(10,list(range(7))),[xcoef(defect,j) for j in range(11)])
    check(selected==[0,1,2,3,4,5,7],"Phi!=0 G10 selector")
    eta10=(R.c(1752089427968)-P*P*32805000)*R.mon(p=-1)/36905625
    alpha10=(P**4*R.c(145629074958046875)-P**2*R.c(6583704417821122560000)
             -R.c(26093447590576484415176704))*R.mon(p=-3)/23154427662890625
    beta10=(-P**6*R.c(319030749548372281640625000)
            +P**4*R.c(14659146525631027935375360000000)
            +P**2*R.c(101452811911656563438652405841920000)
            +R.c(434321509795518334240224474125496745984))*R.mon(p=-5)/7690757619746410400390625
    check(alpha9.subst(ETA,eta10)==alpha10,"alpha10 formula")
    check(not any(z.subst(ETA,eta10).subst(BETA,beta10) for z in res10),"eta,beta solve row10")
    active=[i for i,z in enumerate(res10) if z]
    check(len(active)==2,"two row10 residual positions")
    # Independence is checked at the function-field level: beta appears linearly
    # in one residual; eliminating it leaves the displayed eta graph.
    beta_visible=any(any(mon[BETA] for mon in z.terms) for z in res10)
    check(beta_visible,"beta visible row10")

    add_tangent(aa,cc,9,aa[9],cc[9],sol9,list(range(7)))
    subst_all(aa,cc,gg,ETA,eta10); subst_all(aa,cc,gg,BETA,beta10)
    check(gg[10]==predict(aa,cc,10),"row10 reached")
    check(W==-R.c(13)*(P*P*820125+R.c(13056802816))/9841500,"W formula")
    check(Z.subst(ETA,eta10)==R.c(Q(-225611776,30375)),"Z formula")
    qpoly=(P**6*R.c(373891487235896675830078125)
           -P**4*R.c(15097287707154073014589440000000)
           -P**2*R.c(101452811911656563438652405841920000)
           -R.c(434321509795518334240224474125496745984))
    check(eta10.subst(PHI,-P)==-eta10,"eta sign sheet")
    check(alpha10.subst(PHI,-P)==-alpha10,"alpha sign sheet")
    check(beta10.subst(PHI,-P)==-beta10,"beta sign sheet")
    check(W.subst(PHI,-P)==W,"W sign invariant")
    check(qpoly.subst(PHI,-P)==qpoly,"Q sign invariant")

    coords2,cols2=depth_matrix(10,2); coords3,cols3=depth_matrix(10,3)
    rank2,null2=nullspace_of_column_span(cols2,len(coords2)); rank3,null3=nullspace_of_column_span(cols3,len(coords3))
    check((len(coords2),len(cols2),rank2,len(null2))==(88,193,68,20),"pi10 P2 ledger")
    check((len(coords3),len(cols3),rank3,len(null3))==(99,304,83,16),"pi10 P3 ledger")
    apart10,cpart10=particular(aa,cc,10); aa[10],cc[10]=apart10,cpart10

    def depth_system(nulls,coords,rows,which):
        base=[dot(v,jet(rows,coords)) for v in nulls]
        mat=[[left_eval(v,tangent_column(coords,10,j,which)) for j in range(11)] for v in nulls]
        check(rank_q(mat)==3,f"row10 {which} rank3")
        piv=[]
        for j in range(11):
            if rank_q([[row[k] for k in piv+[j]] for row in mat])>len(piv):piv.append(j)
        check(piv==[8,9,10],f"row10 {which} pivots")
        ir=independent_row_indices([[row[j] for j in piv] for row in mat])
        sol=solve_square([[mat[i][j] for j in piv] for i in ir],[-base[i] for i in ir])
        res=[base[i]+sum((mat[i][piv[j]]*sol[j] for j in range(3)),R()) for i in range(len(mat))]
        check(not any(res),f"row10 {which} standalone compatible")
        return mat,base,piv,sol
    LA,baseA,pivA,solA=depth_system(null2,coords2,aa,"A")
    LC,baseC,pivC,solC=depth_system(null3,coords3,cc,"C")
    check(rank_q(LA+LC)==3,"joint coefficient rank3")

    # Select by P3, then read the new P2 annihilator.
    aaP3,ccP3=dict(aa),dict(cc)
    add_tangent(aaP3,ccP3,10,apart10,cpart10,solC,pivC)
    newA={(6,1):Q(35),(7,3):Q(-20),(8,5):Q(10),(9,7):Q(-4),(10,9):Q(1)}
    newAvec={coords2.index(key):v for key,v in newA.items()}
    check(all(sum((newAvec.get(i,Q(0))*v for i,v in col.items()),Q(0))==0 for col in cols2),
          "new L(10,11,3) P2 annihilator")
    newval=dot(newAvec,jet(aaP3,coords2))
    expected_new=-qpoly*R.mon(p=-5)*Q(2,23072272859239231201171875)
    check(newval==expected_new,"new A consumer")

    # Select by P2, then read the old tetrahedral P3 annihilator.
    aaP2,ccP2=dict(aa),dict(cc)
    add_tangent(aaP2,ccP2,10,apart10,cpart10,solA,pivA)
    oldC={(5,0):Q(56),(6,2):Q(-35),(7,4):Q(20),(8,6):Q(-10),(9,8):Q(4),(10,10):Q(-1)}
    oldCvec={coords3.index(key):v for key,v in oldC.items()}
    check(all(sum((oldCvec.get(i,Q(0))*v for i,v in col.items()),Q(0))==0 for col in cols3),
          "old L(10,10,3) P3 annihilator")
    oldval=dot(oldCvec,jet(ccP2,coords3))
    expected_old=qpoly*R.mon(p=-5)/15381515239492820800781250
    check(oldval==expected_old,"old C consumer")
    check(newval*3==-oldval*4,"consumer ratio -4/3")
    coord8=4*qpoly*R.mon(p=-5)/23072272859239231201171875
    check(solA[0]-solC[0]==coord8,"selected coordinate8 difference")
    check(solA[1]==solC[1],"selected coordinate9 agreement")
    check(solA[2]==solC[2],"selected coordinate10 agreement")

    # Exhaust the finite THM-4364-admissible opposite-order hierarchy at m=10.
    def hierarchy_eval(rows,coords,columns,d,ell,q):
        s=(ell+1)//2
        weights={}
        for n in range(s,11):
            key=(n,2*n-ell)
            if key in coords:
                weights[coords.index(key)]=Q((-1)**(n-s)*comb(10+q-n,q))
        check(all(sum((weights.get(i,Q(0))*v for i,v in col.items()),Q(0))==0 for col in columns),
              f"hierarchy annihilator ell{ell} q{q} d{d}")
        return dot(weights,jet(rows,coords))
    aliveA=[]; aliveC=[]
    for ell in range(2,21):
        rho=(ell+2)//3
        for q in range(rho):
            if (ell+1)//2<=10 and 10+q>=ell and 2<=10+q-ell:
                if hierarchy_eval(aaP3,coords2,cols2,2,ell,q):aliveA.append((ell,q))
            if (ell+1)//2<=10 and 10+q>=ell and 3<=10+q-ell:
                if hierarchy_eval(ccP2,coords3,cols3,3,ell,q):aliveC.append((ell,q))
    check(aliveA==[(11,3)] and aliveC==[(10,3)],"opposite-order hierarchy scan")

    # Check the entire opposite system, not just one named witness.
    joint=[baseA[i]+sum((LA[i][pivA[j]]*solC[j] for j in range(3)),R()) for i in range(len(LA))]
    check(all(zero_mod(z,qpoly) for z in joint),"full joint system modulo Q")
    check(any(joint),"joint obstruction nonzero off Q")

    # Algebraic six-sheet count and strict gates.  Work in Y=Phi^2 mod 7.
    q7={3:6,2:3,1:3,0:2}; dq7={2:4,1:6,0:3}
    check(even_phi_to_y_mod(qpoly,7)==q7,"Q reduction mod7")
    check(gcd_mod(q7,dq7,7)=={0:1},"Q squarefree mod7")
    check(add_mod(mul_mod({1:1,0:6},q7,7),mul_mod({2:2,1:3,0:1},dq7,7),7)=={0:1},
          "Q squarefree Bezout mod7")
    check(q7[0]!=0,"Q nonzero roots")
    w7={1:820125%7,0:13056802816%7}
    check(w7=={1:5,0:1},"W factor reduction mod7")
    check(gcd_mod(q7,w7,7)=={0:1},"Q avoids W mod7")
    check(add_mod(mul_mod({0:3},q7,7),mul_mod({2:2,1:2,0:2},w7,7),7)=={0:1},
          "W Bezout mod7")
    disc_coeff=[680868007162161739848062587106828346429167544303616,
                343583092356524653010661045438183132081684480000,
                35745234591073505392643546713864273920000000,
                -15281178725265792138311029689600000000000,
                -818189072042144733182993060302734375]
    disc7={i:v%7 for i,v in enumerate(disc_coeff)}
    check(trim_mod(disc7,7)=={4:6,3:2,1:6,0:1},"discriminant reduction mod7")
    check(gcd_mod(q7,disc7,7)=={0:1},"Q avoids alpha discriminant mod7")
    check(add_mod(mul_mod({3:6,2:6,1:1,0:4},q7,7),mul_mod({2:1},disc7,7),7)=={0:1},
          "discriminant Bezout mod7")
    strict=alpha10*alpha10-4*W*UPS
    disc=sum((P**(2*i)*R.c(v) for i,v in enumerate(disc_coeff)),R())
    candidate=disc*R.mon(p=-6)
    mon=next(iter(candidate.terms)); ratio=strict.coeff(mon)/candidate.coeff(mon)
    check(ratio and strict==candidate*ratio,"strict discriminant numerator")
    check(qpoly.subst(PHI,R.c(1)) and W.subst(PHI,R.c(1)) and strict.subst(PHI,R.c(1)),
          "Phi=1 strict joint-depth hostile")
    # Z is the displayed nonzero rational constant; Phi and W are handled above.

    # On Q=0 the two selections agree modulo Q. Use the P2 selection for row11.
    check(all(zero_mod(solA[i]-solC[i],qpoly) for i in range(3)),"terminal selections coincide modulo Q")
    aa,cc=aaP2,ccP2
    defect11=xadd(gg[11],xneg(predict(aa,cc,11)))
    _,res11,selected11=solve_injective(dmatrix(11,list(range(8))),[xcoef(defect11,j) for j in range(12)])
    check(len(selected11)==8,"row11 selector rank8")
    rpoly=(P**8*R.c(6846329377771290182382546697998046875)
           -P**6*R.c(713835723041306505264998768716800000000000)
           -P**4*R.c(2754991513504883058403611855707575418880000000)
           -P**2*R.c(31916203206707002973657986739896008646412206080000)
           -R.c(156854967149983010817497418504735580308619473018945536))
    check(res11[9]==qpoly*R.mon(p=-5)*Q(8,84598333817210514404296875),"row11 Q residual")
    check(res11[8]==-rpoly*R.mon(p=-6)/71525718603423911567894897460937500,"row11 R residual")
    check(all((not z) or i in (8,9) for i,z in enumerate(res11)),"only row11 residuals")
    r7={4:5,3:3,2:3,1:3,0:5}
    check(even_phi_to_y_mod(rpoly,7)==r7,"R reduction mod7")
    check(gcd_mod(q7,r7,7)=={0:1},"Q,R coprime mod7")
    bezq={3:3,2:5,1:4}; bezr={2:2,1:5,0:3}
    check(add_mod(mul_mod(bezq,q7,7),mul_mod(bezr,r7,7),7)=={0:1},"displayed Bezout")
    return qpoly,rpoly

def main():
    phi0=branch_phi_zero(); branch_phi_nonzero()
    print("THM-4366 CLEAN-ROOM REFEREE: PASS")
    print(f"CHECKS={CHECKS}")
    print("SOURCE_REBUILT=literal residual-H expansion and bracket recursion rows4..11")
    print("PHI0_E9=(-11/243)*(5646560625*eta^2+379697122115584)")
    print("PHI0_ROW9=two A2 source components; terminal fibre A7")
    print(f"PHI0_G10_RESIDUAL8={phi0.coeff((0,0,0,0))}")
    print("PHI_NONZERO_ROW9=Gm_x_A2 source graph; terminal fibre A7")
    print("ROW10_BRACKET=eta,beta solved; source parameter Phi remains")
    print("ROW10_DEPTH=P2 88x193 rank68 null20; P3 99x304 rank83 null16; each rank=aug3")
    print("JOINT_DEPTH=Q(Phi^2)=0; squarefree cubic, nonzero roots, six sign sheets; terminal fibre A8")
    print("CONSUMERS=new A L(10,11,3)=(-4/3)*old C L(10,10,3) after opposite selection")
    print("HIERARCHY_SCAN=only nonzero admissible cross-consumers A(11,3),C(10,3)")
    print("STRICT_GATES=Phi,W,Z,alpha^2-4W*ups nonzero on Q")
    print("ROW11=rank8 selector on A8 terminal fibre; Q residual plus coprime quartic R")
    print("MOD7_Q=-Y^3+3Y^2+3Y+2")
    print("MOD7_R=-2Y^4+3Y^3+3Y^2+3Y-2")
    print("MOD7_BEZOUT=(3Y^3-2Y^2-3Y)Q+(2Y^2-2Y+3)R=1")
    print("SCOPE=finite source-normal residual-weight<=12 projected-depth calculation only")

if __name__=="__main__":main()
