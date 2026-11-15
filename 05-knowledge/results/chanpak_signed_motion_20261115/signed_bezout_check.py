"""Independent permutations for local profiles, subset DP for global H."""
from oracle import block, tournament, lift, path_dp

def tt(n):
    return tournament(n,(1<<(n*(n-1)//2))-1)

def verify():
    rows=[(9,{0,1},39,15,6),(9,{1,2},41,17,8),(6,{0,1,2},42,18,10),(6,{1,2,3},44,20,12)]
    for a in range(1,5):
        actual=[]
        for code,V,E,L,J in rows:
            d=tournament(5,code)
            assert block(d,V)[:5]==(5,64,E,J,L)
            H=path_dp(lift([tt(a),d,tt(a)],[set(),V,set()]))
            assert H==5+64*4**a+24*2**a-J
            actual.append(H)
        assert [x-y for x,y in zip(actual,actual[1:])]==[2]*3
        print(a,sorted(actual))
    for n in range(4,9,2):
        for a in range(1,4):
            p=n//2
            for V,J in [(set(),0),({p-1,p},2),({p-2,p+1},4)]:
                d=tt(n)
                h,q,E,j,L,*_=block(d,V)
                assert E==L and j==J
                assert path_dp(lift([tt(a),d,tt(a)],[set(),V,set()]))==1+2**(n+2*a)-J
    print('verified 43 literal strong tournaments')
if __name__=='__main__': verify()


def check_moving_q():
    rows=[(7,229304,{1,3,4,6},336,178,124,60,1),
          (6,4329,{0,1,2,5},168,91,125,62,2),
          (7,2065790,{0,2,5},336,198,129,64,1),
          (6,4329,{0,2,3,5},168,95,127,66,2),
          (6,4329,{0,3,4,5},168,99,129,68,2)]
    for z in (0,1):
        counts=[]
        for n,code,V,q,E,L,J,r in rows:
            D=tournament(n,code)
            assert block(D,V)[:5]==(9,q,E,J,L)
            a=z+r; b=z+3
            assert q*2**(r+3)==5376 and E*2**r-8*L==-636
            H=path_dp(lift([tt(a),D,tt(b)],[set(),V,set()]))
            assert H==5376*4**z-636*2**z+9-J
            counts.append(H)
        assert [x-y for x,y in zip(counts,counts[1:])]==[2]*4
        print(z,sorted(counts))
if __name__=='__main__': check_moving_q()

def check_all_moving():
    for u in (0,1,2):
        for n,c,V,h,q,E,L,J,r in [(4,63,{0,1},1,16,3,12,2,2),
                                  (5,4,{0,4},5,64,19,37,8,1)]:
            D=tournament(n,c)
            assert block(D,V)[:5]==(h,q,E,J,L)
            a=u+r
            assert q*4**r==256 and (E-L)*2**r==-36
            count=path_dp(lift([tt(a),D,tt(a)],[set(),V,set()]))
            assert count==256*4**u-36*2**u+h-J
    print('six independent checks with all raw coordinates moving')
if __name__=='__main__': check_all_moving()
