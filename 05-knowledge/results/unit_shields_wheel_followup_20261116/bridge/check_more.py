from check import prof, spine_word, actual_h, check_strong, rc, require

R = [
    ('1010010',91,'1101001111100','1101100110010',
     5,4004667,6674445,865551),
    ('11010010011100',12725,'110101001001001000','111011010110110010',
     12,18415108821,30691848035,5046246339),
    ('11101010011000',14761,'101010011110000100','110001101011010110',
     10,15233834805,25389724675,4649307361),
]

def check_more():
    checks=0
    for y,p,u0,u1,edits,r3,r5,b in R:
        k=len(y)
        a=2**k
        hy,ay=prof(y)
        require(y.startswith('1') and y.endswith('0'),'suffix strong')
        require(p==a-int(y[::-1],2) and p*ay-a*(hy-2)==1,'suffix unit')
        h0,P0=prof(u0);h1,P1=prof(u1)
        require((h1-h0,P1-P0)==(-2*a,2*p),'left change')
        require(edits==sum(x!=z for x,z in zip(u0,u1)),'edit number')
        for t in (1,2,3,4,5,9,12,27,67,100):
            for pref,m,coef in [('11',3,r3),('101',5,r5)]:
                v=rc(pref+'0'*t+y)
                D=2**t
                for i,u in enumerate([u0,u1]):
                    s=spine_word(u,v)
                    require(actual_h(s)==coef*D+b+2*i,'literal count')
                    if t==1:
                        check_strong(s)
                    checks+=1
    print(checks,'additional direct physical three-run counts passed.')

if __name__=='__main__':
    check_more()
