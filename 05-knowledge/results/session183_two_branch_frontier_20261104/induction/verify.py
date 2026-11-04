"""No dependencies; independent physical joints and rational escape controls."""
import json
from pathlib import Path
from fractions import Fraction
OUT=Path(__file__).resolve().parent

def require(x,msg):
    if not x: raise RuntimeError(msg)

def joint_pair(f):
    ports=[i for i in range(f.bit_length()) if f>>i&1]
    return sum((1<<(k-i-1))+(k==i+1) for i in ports for k in ports if i<k)

def joint_scan(f):
    j=g=e=0
    for i in range(f.bit_length()):
        b=f>>i&1
        j,g,e=j+b*(g+e),2*g+b,b
    return j

def literal_joint_matrix(m):
    # TT paths on a given subset have exactly one order. Neither slot empty.
    out=[[0]*m for _ in range(m)]
    full=(1<<m)-1
    for left in range(1,full):
        right=full^left
        endpoint=left.bit_length()-1
        start=(right&-right).bit_length()-1
        out[endpoint][start]+=1
    return out

def reverse_fraction(u):
    return sum((Fraction(1,1<<i) for i in range(u.bit_length()) if u>>i&1),Fraction())

def test_product(u):
    m=u.bit_length();bits=[i for i in range(m) if u>>i&1]
    w=len(bits);a=sum(bool(u>>(i+1)&1) for i in bits)
    B=sum((Fraction(1,1<<(k-i)) for i in bits for k in bits if i<k),Fraction())
    R=reverse_fraction(u);J=joint_scan(u)
    require(J-u*R/2==a-(w+B)/2,'exact reversed product')
    require(abs(Fraction(J,u)-R/2)<=Fraction(m,u),'rational error')

def main():
    physical=0;directed=0
    for m in range(1,9):
        matrix=literal_joint_matrix(m)
        for i in range(m):
            for k in range(m):
                expected=(1<<(i-k-1)) if i>k else int(k==i+1)
                if i==k:expected=0
                require(matrix[i][k]==expected,'physical directed joint');directed+=1
        for f in range(1<<(m-1),1<<m):
            actual=sum(matrix[i][k] for i in range(m) for k in range(m) if f>>i&1 and f>>k&1)
            require(actual==joint_pair(f)==joint_scan(f),'literal aggregate joint')
            test_product(f);physical+=1
    escape=0
    for u in range(1,1<<15,2):
        J=joint_scan(u);K=joint_scan(3*u);m=u.bit_length();mp=(3*u).bit_length()
        require(J==joint_pair(u) and K==joint_pair(3*u),'independent joint evaluators')
        require(8*abs(K-3*J)>=3*u-8*(3*m+mp),'sharp escape estimate')
        require(8*abs(K-3*J)>=3*u-8*(4*m+2),'loose escape estimate')
        if u>=128:require(abs(K-3*J)>6,'bounded alphabet cutoff')
        escape+=1
    hits=[]
    for u in range(1,128,2):
        d=3*joint_scan(u)-joint_scan(3*u)-2
        if d in (2,4):hits.append({'odd_mask':u,'digit':d})
    require(hits==[{'odd_mask':7,'digit':4}],'finite below-cutoff census')
    shifts=0
    for f in range(1,256):
        for k in range(-4,5):
            ff=Fraction(3*f,1)*(Fraction(2)**(-k))
            if ff.denominator!=1:continue
            ff=int(ff);d=3*joint_scan(f)-joint_scan(3*f)-2
            for t in (max(1,1-k),max(1,1-k)+2):
                H=1-joint_scan(f)+(1<<t)*f
                Hp=1-joint_scan(ff)+(1<<(t+k))*ff
                require(Hp==3*H+d,'uniform shifted chart')
            require((ff//(ff&-ff))==3*(f//(f&-f)),'odd part growth');shifts+=1
    orbit=[];u=7
    for n in range(64):
        m=u.bit_length();test_product(u)
        R=reverse_fraction(u)/2;Rp=reverse_fraction(3*u)/2
        require(abs(R-Rp)>=Fraction(1,8),'mod8 separated cylinders')
        require(8*abs(joint_scan(3*u)-3*joint_scan(u))>=3*u-8*(4*m+2),'large exact escape')
        if n in (0,1,2,8,32,63):
            orbit.append({'n':n,'odd_mask':u,'joint':joint_scan(u),'digit':3*joint_scan(u)-joint_scan(3*u)-2})
        u*=3
    result={'physical_masks':physical,'individual_joint_checks':directed,'escape_odd_masks':escape,'shifted_chart_checks':shifts,'bounded_cutoff_odd_masks':64,'hits':hits,'large_orbit_samples':orbit,'claim':'finite controls of the all-order proof, not evidence of closure'}
    s=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (OUT/'verification.json').write_text(s);print(s,end='')
if __name__=='__main__':main()
