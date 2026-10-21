"""Consume frozen exclusions and check their exact generating-function consequences."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json

HERE=Path(__file__).resolve().parent
PINS={
 'deficit_six_centers_20261021b/combined_receipt.json':'9a1bdb64fdc717cb3af1aa9c2c65eef5b59978786f161969fbe1147411d8e048',
 'deficit_twelve_20261021c/receipt.json':'fe54716a7b6e6fac30a41e118e9a59e2e80cdea0da6e7602932cfa0697779da2',
 'deficit_fourteen_20261021c/receipt.json':'2b150fd9d93ea708cd201aef0eaa41fa5d5205daed5b1dba41757cf0da3c329d',
 'center_transitivity_20261021c/receipt.json':'1cd59dbfdc5ff5c296fc69aa7f5327a41991e237ba2c291d991697e2792ad0c8',
 'center_transitivity_20261021c/receipt35.json':'360ac9b31f78935cda6699c5dc08693f75d0a02079e777d1d31be720f5e47be1',
 'center_transitivity_20261021c/receipt15.json':'4697b30544a6e55a8bc4d81e2d4b1bffd61fe00cf6eb5e37be80c33dae720280',
 'tail_geometry_20261021c/receipt.json':'fe7d8190a61e3e883f1b82ad38d6d11a9abeebc9da0833aafa2eb5fa9c16aae9',
}
CHECKS=0


def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise ValueError(msg)


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c


def add(*polys):
    c=[F(0)]*max(map(len,polys))
    for p in polys:
        for i,x in enumerate(p): c[i]+=x
    while len(c)>1 and not c[-1]: c.pop()
    return c


def shift(p,n): return [F(0)]*n+list(p)


def run():
    for path,h in PINS.items():
        raw=(HERE.parent/path).read_bytes().replace(b'\r\n',b'\n')
        need(sha256(raw).hexdigest()==h,'frozen receipt '+path)
    initial=[12,14,15,18,20,21,24,25,28,35]
    removed=[12,14,25,35]
    remaining=[x for x in initial if x not in removed]
    need(remaining==[15,18,20,21,24,28],'six remaining coefficients')
    p=[1,-1,0,-1]; u=[1,0,-1,0,0,0,-1]
    p2=mul(p,p); u2=mul(u,u); D=mul(p2,u2)
    need(len(D)==19 and D[0]==1,'degree18 denominator')
    # Formal top coefficient is a hostile algebra control, not a center claim.
    R=[F(int(i in (5,6,7,28))) for i in range(29)]
    numer=add(shift([x/2 for x in add(u2,[-x for x in mul(p2,u)])],5),
              mul(mul(p2,u),R),shift(mul(p2,[1,1,0,1]),10))
    need(len(numer)==41 and numer[-1]==-1,'degree40 numerator can occur for formal top coefficient')
    coeff=[]
    for n in range(81):
        coeff.append((numer[n] if n<len(numer) else 0)-sum(D[j]*coeff[n-j] for j in range(1,min(n,len(D)-1)+1)))
    need(all(x.denominator==1 for x in coeff),'integral formal coefficient control')
    residual=[sum(D[j]*coeff[n-j] for j in range(min(n,len(D)-1)+1)) for n in range(len(coeff))]
    need(residual[40]==-1 and all(x==0 for x in residual[41:]),'onset41 is sufficient and degree-only bound can be sharp')
    # Independent product of the three inherited rational-series summands.
    A=[1]
    for n in range(1,81): A.append(A[n-1]+(A[n-3] if n>=3 else 0))
    A2=[F(A[n//2] if n%2==0 else 0) for n in range(81)]
    direct=add(shift([x/2 for x in add(mul(A,A),[-x for x in A2])],5),mul(A2,R),
               shift(mul(mul(A2,A2),[1,1,0,1]),10))
    need(direct[:81]==coeff,'independent formal summand expansion')
    return {'status':'PROVED consequences consuming separately audited finite exclusions; six center coefficients OPEN',
            'pins':PINS,'removed_orders':removed,'remaining_orders':remaining,'center_degree_at_most':28,
            'known_trivalent_coefficients_orders_0_through_14':[int(x) for x in coeff[:15]],
            'denominator_coefficients':[int(x) for x in D],'minimal_eventual_order_inherited':18,
            'remaining_order15_must_be_non_vertex_transitive':True,
            'safe_onset':41,'formal_z28_control_is_not_an_existence_claim':True,'checks':CHECKS}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run();data=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
    out=HERE/'combined_receipt.json'
    if a.write: out.write_bytes(data)
    else: need(out.read_bytes().replace(b'\r\n',b'\n')==data,'saved combined receipt')
    print(data.decode())
