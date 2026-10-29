"""Direct point counts, not Hasse polynomial arithmetic."""
import json

def require(b,m):
    if not b:raise ArithmeticError(m)

def trace(p,a,b):
    sq={i*i%p for i in range(1,p)}
    return -sum(0 if (v:=(x*x*x+a*x+b)%p)==0 else 1 if v in sq else -1 for x in range(p))

def rational_j(p):
    js=[]
    for j in range(p):
        if j==0:a,b=0,1
        elif j==1728%p:a,b=1,0
        else:
            k=j*pow(1728-j,-1,p)%p;a,b=3*k,2*k
        require((4*a**3+27*b*b)%p!=0,'smooth curve')
        if trace(p,a,b)==0:js.append(j)
    return js

def main():
    ans=[]
    for p in [37,71]:
        js=rational_j(p);fibres={}
        sq={i*i%p for i in range(1,p)}
        for l in range(2,p):
            t=-sum(0 if (v:=x*(x-1)*(x-l)%p)==0 else 1 if v in sq else -1 for x in range(p))
            if t==0:
                j=256*pow(1-l+l*l,3,p)*pow(l*l*(1-l)**2,-1,p)%p
                fibres.setdefault(j,[]).append(l)
        ans.append({'p':p,'supersingular_rational_j':js,'rational_legendre_fibres':fibres})
    require(ans[0]['supersingular_rational_j']==[8],'37 rational j')
    require(len(ans[1]['supersingular_rational_j'])==7,'71 seven j')
    require({k:len(v) for k,v in ans[1]['rational_legendre_fibres'].items()}=={17:6,24:3,40:6,48:6},'71 marked fibres')
    print(json.dumps(ans,sort_keys=True))
if __name__=='__main__':main()
