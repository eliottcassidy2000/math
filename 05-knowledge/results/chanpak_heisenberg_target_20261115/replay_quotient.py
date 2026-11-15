from heisenberg import dp,literal
def count_h(w):
    a,h,prev=0,1,0
    for bit in w:
        v=int(bit=='A')
        if not v:h+=a+prev
        a=2*a+v
        prev=v
    return h
for w,p in [('ABB',(5,32)),('ABABB',(19,256))]:
    if dp(literal(w))!=p:raise AssertionError(w)
    print(w,p)
for e in range(5,13):
    q,delta=1<<(18*e),475**e
    k=delta.bit_length()
    g=k+1+18*e
    hi=q*(1<<k)+delta
    seen=set()
    for i in range(2*e+1):
        for j in range(e+1):
            t=5**i*19**j
            r=(1<<k)-delta//t
            v=f'{r:0{k}b}'[::-1].translate(str.maketrans('01','BA'))
            b=int(v.translate(str.maketrans('BA','01')),2)
            a=count_h(v)
            lo=q*(b+1)+t*a
            af=count_h('A'+'B'*g+v)
            N=q*((1<<(g+k))+b+1)+t*af
            if N!=(hi<<g)+lo or not 0<t*a<q or not 0<lo<1<<g:raise AssertionError((e,i,j))
            bb,ll=divmod(lo,q); bb-=1
            rr=int(f'{bb:0{k}b}'[::-1],2)
            d=(1<<k)-rr
            if bb!=b or delta%d or delta//d!=t or ll!=t*a or lo in seen:raise AssertionError(N)
            seen.add(lo)
    print('common-quotient checksum count',e,len(seen))
