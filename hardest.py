from lrc import *

def search(m, cap, topN=10):
    results=[]
    for vs in combinations(range(1,cap+1), m):
        if reduce(gcd, vs)!=1: continue
        M,t=M_exact(list(vs))
        val=M*(m+1)
        results.append((val, vs, M, t))
    results.sort(key=lambda r:(r[0].numerator*1.0/r[0].denominator, r[1]))
    return results

if __name__=='__main__':
    import sys
    m=int(sys.argv[1]); cap=int(sys.argv[2])
    res=search(m,cap)
    print(f'=== m={m}, speeds<= {cap}, total primitive sets={len(res)} ===')
    print(f'min M*(m+1)={float(res[0][0]):.6f}')
    print(f'{"M*(m+1)":>12} {"speeds":<22} {"M":>10} {"argmax t":>10} {"denom(t)":>8} {"gap":>14}')
    for val,vs,M,t in res[:10]:
        gap=M-F(1,m+1)
        print(f'{str(val):>12} {str(list(vs)):<22} {str(M):>10} {str(t):>10} {t.denominator:>8} {str(gap):>14}')
