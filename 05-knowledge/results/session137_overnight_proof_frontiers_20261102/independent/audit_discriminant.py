"""Independent exact SymPy determinants and rational solve; not Bareiss probe."""
import json
import sympy as sp
from tournament_kernel import strong,require

def invariants(rows):
    n=len(rows)
    a=sp.Matrix([[int(rows[i]>>j&1) for j in range(n)] for i in range(n)])
    m=sp.eye(n)+a-a.T
    x=m.inv()*sp.ones(n,1)
    s=(sp.ones(1,n)*x)[0]
    D=m.det(method='domain-ge')/sp.Integer(2)**(n-1)
    u=(sp.eye(n)+a).det(method='domain-ge');v=(-1)**n*a.det(method='domain-ge')
    require(u==D*(1+s)/2 and v==D*(1-s)/2,'coordinate identity')
    require(s==(x.T*x)[0],'inverse response norm')
    return {'n':n,'u':int(u),'v':int(v),'disc':str(D),'s':str(s),'response':[str(y) for y in x],
            'strong':strong(tuple(tuple(int(a[i,j]) for j in range(n)) for i in range(n)))}

def run():
    rows=[30,8,122,80,2,91,19];n=len(rows)
    joined=[r|(((1<<n)-1)<<n) for r in rows]+[r<<n for r in rows]
    factor,product=invariants(rows),invariants(joined)
    require(factor['strong'] and factor['disc']=='6' and factor['s']=='2/3','factor witness')
    require(product['disc']=='26' and product['s']=='12/13','join witness')
    c3=[2,4,1];rows=[];controls=[]
    for r in range(1,5):
        old=len(rows);rows=[x|(((1<<3)-1)<<old) for x in rows]+[x<<old for x in c3]
        data=invariants(rows);controls.append(data)
        require(sp.Rational(data['disc'])==2**r+(-1)**r,'cycle join disc')
        require(sp.Rational(data['s'])==sp.Rational(2**r-(-1)**r,2**r+(-1)**r),'cycle join response')
    print(json.dumps({'factor_rows':[30,8,122,80,2,91,19],'factor':factor,'ordinal_square':product,
                      'disc_product':36,'supermultiplicativity':'REFUTED','C3_ordinal_powers':controls},sort_keys=True,indent=2))
if __name__=='__main__':run()
