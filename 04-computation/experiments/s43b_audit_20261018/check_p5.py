from math import comb
FAM = {
 'A1': (lambda k: comb(2*k,k)**3, [-512,-64,-8,1,16,64,256,4096]),
 'A2': (lambda k: comb(2*k,k)**2*comb(3*k,k), [-8640,-1728,-192,-27,8,64,108,216,1458,3375]),
 'A3': (lambda k: comb(2*k,k)**2*comb(4*k,2*k), [-82944,-12288,-3969,-1024,-144,81,256,648,2304,20736]),
 'A4': (lambda k: comb(2*k,k)*comb(3*k,k)*comb(6*k,3*k), [-32768,-3375,1728,8000,54000]),
}
def ok(Af,m,p):
    pp=p*p
    if m%p==0: return None
    w=pow(m%pp,-1,pp); S0=S1=0; pw=1
    for k in range(2*p):
        t=Af(k)%pp*pw%pp
        if k<p: S0+=t
        else: S1+=t
        pw=pw*w%pp
    return (S1-Af(1)*w*S0)%pp==0
for name,(Af,ms) in FAM.items():
    print(name, 'p=5:', {m: ok(Af,m,5) for m in ms}, ' p=3:', {m: ok(Af,m,3) for m in ms})
