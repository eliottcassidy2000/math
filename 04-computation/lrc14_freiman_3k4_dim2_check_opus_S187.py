# Confirm the |A+^A|=13 unbounded-spread outliers are Freiman dimension-2 GAPs (unions of equal-spaced runs)
def rss(A): return len(set(A[i]+A[j] for i in range(len(A)) for j in range(i+1,len(A))))
outliers=[(0,1,7,8,9,15,16),(0,1,8,9,10,17,18),(0,2,7,9,11,16,18),(0,3,5,8,11,13,16)]
for A in outliers:
    diffs=[A[i+1]-A[i] for i in range(6)]
    # detect 2-D GAP: A = {p + q*s : ...}; check if A ⊆ {i + j*D} small grid
    found=None
    for D in range(4,10):
        # residues mod D
        res=sorted(set(a%D for a in A))
        # is A a subset of a small #residues each forming a short run?
        cols={r:sorted(a for a in A if a%D==r) for r in res}
        widths=[ (max(c)-min(c))//D+1 for c in cols.values()]
        if len(res)<=3 and max(widths)<=3:
            found=(D,res,widths); break
    print(f"{A} |A+^A|={rss(A)} spread={A[-1]} diffs={diffs} -> 2D-GAP mod {found[0] if found else '?'}: {len(found[1]) if found else '?'} cols" if found else f"{A} |A+^A|={rss(A)} spread={A[-1]} diffs={diffs} -> (structured runs)")
