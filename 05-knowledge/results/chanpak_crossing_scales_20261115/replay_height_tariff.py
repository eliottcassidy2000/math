from replay_chromatic_interface import signature

def main():
    for k in range(1,13):
        got=set()
        q=1<<k
        g=k+2
        for b in range(q):
            d,c=signature((1,2),(1,2),f'{b:0{k}b}')
            n=(d<<g)+c
            if not 0<c<(1<<g) or not n<=15*q*q or n in got:
                raise AssertionError((k,b,n))
            got.add(n)
        if len(got)!=q:raise AssertionError(k)
    print('exact lower-bound witnesses for each suffix width 1..12')
    loops=0
    for u in ((1,2),(5,32),(19,256),(345,5640)):
        for w in ((1,2),(5,32),(73,1088)):
            d=(u[1]-u[0])*w[1]
            for k in range(1,8):
                for b in range(1<<k):
                    D,c=signature(u,w,f'{b:0{k}b}')
                    if D<d*(1<<k) or 2*c<d*(b+1):
                        raise AssertionError((u,w,k,b))
                    loops+=1
    print('literal profile-sheet floor controls:',loops)
if __name__=='__main__':main()
