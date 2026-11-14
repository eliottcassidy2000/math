if __name__=='__main__':
    for n in range(25281,2000000,3292):
        u,v=47,56
        k=0
        while True:
            U,V=u+4*v,2*u+4*v
            if (U+2*V)**2>n: break
            u,v=U,V
            k+=1
        a=n*pow(u,-1,2*v)%(2*v)
        b=(n-u*a)//v
        assert a&1 and b%2==0 and b>=2*a and a*a<=n and b*b<324*n
        x,y=a,b
        for _ in range(k):
            x,y=x+2*y,4*x+4*y
        assert 47*x+56*y==n
    print('target-selected logarithmic transfer depth and square-root source bound check')
