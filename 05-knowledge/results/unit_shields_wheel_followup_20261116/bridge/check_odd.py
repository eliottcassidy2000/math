from itertools import product
from check import prof, rc, spine_word, actual_h, require, check_strong

def guard(head):
    h,P=prof(head)
    t=(2*P-3*h-1)//2
    require(t>=1,'tail')
    w=head+'10'*t
    return w,t

def invert(w):
    t=0
    while w.endswith('10'):
        w=w[:-2]
        t+=1
    if len(w)<2 or not (w.startswith('1') and w.endswith('1')):
        return None
    h,P=prof(w)
    return w if t==(2*P-3*h-1)//2 else None

def check():
    u0,u1='111010','111110'
    require(prof(u0)==(39,59) and prof(u1)==(33,63),'profiles')
    require(sum(a!=b for a,b in zip(u0,u1))==1,'one edit')
    wanted={'101':(1,'10110',23,1941),'11':(2,'111010',59,5001),
            '1011':(7,'101110101010101010',191147,16247481)}
    for head,(t,w,Q,H) in wanted.items():
        require(guard(head)==(w,t),'example guard')
        require(prof(w)[1]==Q and 85*Q-14==H,'example score')
    checked=0
    for n in range(2,9):
        for m in product('01',repeat=n-2):
            head='1'+''.join(m)+'1'
            w,t=guard(head)
            h,Q=prof(w)
            require(2*Q-3*h==1 and invert(w)==head,'guard inverse')
            for i,u in enumerate([u0,u1]):
                bits=spine_word(u,rc(w))
                require(actual_h(bits)==85*Q-14+2*i,'literal count')
                if n<=3: check_strong(bits)
                checked+=1
    finite=0
    for n in range(2,17):
        for m in product('01',repeat=n-2):
            w='1'+''.join(m)+'0'
            h,P=prof(w)
            require((2*P-3*h==1)==(invert(w) is not None),'complete finite replay')
            finite+=1
    print(f'{checked} literal three-run counts; exact inverse on {finite} strong words.')

if __name__=='__main__':
    check()
