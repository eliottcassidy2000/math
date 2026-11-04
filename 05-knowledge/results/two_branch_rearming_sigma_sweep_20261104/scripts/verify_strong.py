def strong(rows):
    n=len(rows)
    def reach(adj):
        seen=1; st=[0]
        while st:
            v=st.pop(); m=adj[v]&~seen
            while m:
                b=m&-m; m^=b; seen|=b; st.append(b.bit_length()-1)
        return seen==(1<<n)-1
    rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
    return reach(rows) and reach(rev)
