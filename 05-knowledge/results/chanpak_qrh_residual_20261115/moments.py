from fractions import Fraction as Q
from complement import p

def control():
    count=0
    for s in range(1,9):
        g=1<<s
        for u in (1,2,19):
            for v in (1,2,37):
                mean=p(u)+p(v)+Q((g-1)*(u+v+1)+(u&1)+(v&1),2)
                c=[Q((u-v)*(1<<j)+(j==0)*((u&1)-(v&1)),2)+Q(j==s-1,6) for j in range(s)]
                var=sum(z*z for z in c)+Q((1<<(2*s-2))+6*s-7,9)
                data=[]
                for b in range(g):
                    x=g*u+b
                    y=g*v+g-1-b
                    T=p(x)+p(y)+Q((x&1)+2*(y&1),3)
                    X=[1-2*((b>>(s-1-j))&1) for j in range(s)]
                    v2=sum((1<<max(0,j-i-2))*X[i]*X[j] for j in range(s) for i in range(j))
                    assert T==mean+sum(a*z for a,z in zip(c,X))-v2
                    data.append(T)
                assert sum(data)/g==mean
                assert sum((z-mean)**2 for z in data)/g==var
                count+=1
    print('exact low-phase Walsh expansions and two moments verified on',count,'families')
if __name__=='__main__':
    control()
