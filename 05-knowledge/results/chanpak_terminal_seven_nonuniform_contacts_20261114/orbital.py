from functools import cache
from math import comb

def rows7(s):
    return tuple(sum(1<<((i+d)%7) for d in s) for i in range(7))

def induced_table(rows):
    n=len(rows)
    @cache
    def end(s,k):
        if s==1<<k:return 1
        t=s^(1<<k)
        return sum(end(t,j) for j,arcs in enumerate(rows) if t>>j&1 and arcs>>k&1)
    return [sum(end(s,j) for j in range(n) if s>>j&1) if s else 1 for s in range(1<<n)]

def color_poly_value(path,m):
    n=(len(path)-1).bit_length()
    F=[1]+[0]*((1<<n)-1)
    vals=[]
    for _ in range(m):
        res=[0]*(1<<n)
        for S in range(1<<n):
            T=S
            while True:
                res[S]+=path[T]*F[S^T]
                if T==0:break
                T=(T-1)&S
        F=res
        vals.append(F[-1])
    return vals

if __name__=='__main__':
    for ds in ((1,2,3),(1,2,4)):
        tab=induced_table(rows7(ds))
        V=[0]+color_poly_value(tab,7)
        Q=[sum((-1)**(r-j)*comb(r,j)*V[j] for j in range(r+1))//7 for r in range(1,8)]
        print(ds, V, Q)

def glued_paley(differences):
    rows=[]
    for u in range(14):
        i=u%7
        first=u<7
        same=[(0 if first else 7)+(i+d)%7 for d in (1,2,4)]
        cross=[j for j in range(7,14) if (j-7-i)%7 in differences] if first else [
               j for j in range(7) if (i-j)%7 not in differences]
        rows.append(sum(1<<j for j in same+cross))
    return rows

def prime(n):
    from math import isqrt
    return n>=2 and all(n%i for i in range(2,isqrt(n)+1))

def triangle_modules(rows,r):
    n=len(rows)
    b=[a+sum(1<<j for j in range(n,n+r)) for a in rows]
    for i in range(n,n+r):
        b.append(sum(1<<j for j in range(i+1,n+r+1)))
    b.append((1<<n)-1)
    return b

if __name__=='__main__':
    expected=((25,234,954,2160,2880,2160,720),(27,234,954,2160,2880,2160,720))
    for k,U in enumerate(((1,2,3),(1,2,4))):
        tab=induced_table(rows7(U))
        V=[0]+color_poly_value(tab,7)
        assert V==[m**7+28*m**5+112*m**3+(34+14*k)*m for m in range(8)]
        for r in range(1,8):
            assert sum((-1)**(r-j)*comb(r,j)*V[j] for j in range(r+1))==7*expected[k][r-1]
    assert prime(827) and induced_table(triangle_modules(list(rows7((1,2,3))),2))[-1]==7*827
    U=glued_paley({0,1})
    tab=induced_table(U)
    h=tab[-1]
    f2=sum(x*y for x,y in zip(tab,tab[::-1]))
    assert (h,f2)==(15144269,400589504) and prime(h//7)
    assert sorted(x.bit_count() for x in U)==[5]*7+[8]*7
    rot=[(i//7)*7+(i+1)%7 for i in range(14)]
    assert all((U[i]>>j&1)==(U[rot[i]]>>rot[j]&1) for i in range(14) for j in range(14))
    r=2
    H=induced_table(triangle_modules(U,r))[-1]
    assert H==7*(55063605*(1<<r)-50736671)==7*169517749 and prime(H//7)
    print('two coupled Paley orbits, order-7 automorphism and both exact prime cofactors checked')
