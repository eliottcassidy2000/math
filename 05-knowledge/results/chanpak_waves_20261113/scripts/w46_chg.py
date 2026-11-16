import subprocess, sys
def Hs(parent_lists):
    lines = [f"{len(p)+1} " + " ".join(map(str, p)) for p in parent_lists]
    out = subprocess.run(['./w46_treeH'], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout.split()
    return [int(x) for x in out]
def charges(par):          # par = parents of vertices 1..n-1
    n = len(par) + 1
    res = Hs([par] + [par + [v] for v in range(n)])
    H = res[0]
    return H, [res[1 + v] - 2 * H + 1 for v in range(n)]
if __name__ == '__main__':
    for m in range(2, 11):
        par = list(range(m - 1))           # path 0-1-2-...-(m-1)
        H, a = charges(par)
        print(f"path m={m}: H={H} charges by vertex 0..m-1: {a}")
