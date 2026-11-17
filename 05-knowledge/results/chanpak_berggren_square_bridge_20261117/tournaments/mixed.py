from replay import M,F,rc,act,grow,read_singleton,tab,ck
from math import isqrt
def main():
    v=act(M('0101'),1,2);ck(v,(15,38))
    for j in range(1,41):
        A=1<<j;z='0'+'1'*j+'001'
        m=M(z);s=act(m,1,2)
        ck(m,(3*A-1,8*A-6,6*A+2,16*A))
        ck(s,(19*A-13,38*A+2))
        h,R=act(M(rc(z)),*v)
        p=247*A-17
        ck(h+2*R,1729*A-119);ck(h+2*R,7*p)
        ck(p%7 in (1,5,6),True)
        ck((h+2*R+119)//1729,A)
    a,r=grow((2,0),1,'0101')
    b,r=grow(a,r,rc('011001'))
    q=read_singleton(b,r)
    ck(len(q),16);ck(tab(q)[0][-1],6797)
    ck(all(971%d for d in range(2,isqrt(971)+1)),True)
    print('40 mixed norm formulas, exact v7=1, and literal H=6797 with prime quotient 971.')
if __name__=='__main__':main()
