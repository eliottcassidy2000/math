from profile_tools import canon, extend, strong, profile

def need(ok):
    if not ok:
        raise RuntimeError('grade replay failed')

def main():
    expected=[(3,12),(5,32),(9,84),(15,216),(25,548),(45,1416)]
    T={(0,)}; sizes=[1]
    actual=[]
    for n in range(2,9):
        T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
        sizes.append(len(T))
        if n>=3:
            V=[profile(A) for A in T if strong(A)]
            w=(min(h for h,q in V),min(q for h,q in V))
            need(w in V);actual.append(w)
    need(sizes==[1,1,2,4,12,56,456,6880] and actual==expected)
    need([h+99*q for h,q in actual]==[1191,3173,8325,21399,54277,140229])
    print('all strong profile minima through order eight and six exact grade activations')
if __name__=='__main__': main()
