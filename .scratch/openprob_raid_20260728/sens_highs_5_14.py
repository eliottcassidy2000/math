import highspy, numpy as np, time
n, d = 14, 5
N = 1 << n
pc = [bin(x).count("1") for x in range(N)]
h = highspy.Highs()
h.setOptionValue("output_flag", True)
h.setOptionValue("time_limit", 14000.0)
lower = np.zeros(N); upper = np.ones(N)
lower[0] = upper[0] = 0.0
for i in range(n):
    lower[1 << i] = upper[1 << i] = 1.0
h.addVars(N, lower, upper)
h.changeColsIntegrality(N, np.arange(N, dtype=np.int32),
                        np.array([highspy.HighsVarType.kInteger]*N))
starts, idxs, vals = [0], [], []
rows = 0
t0 = time.time()
for S in range(N):
    if pc[S] <= d: continue
    T = S
    while True:
        idxs.append(T)
        vals.append(-1.0 if ((pc[S]-pc[T]) & 1) else 1.0)
        if T == 0: break
        T = (T-1) & S
    starts.append(len(idxs))
    rows += 1
print("built rows", rows, "nnz", len(idxs), f"{time.time()-t0:.0f}s", flush=True)
h.addRows(rows, np.zeros(rows), np.zeros(rows), len(idxs),
          np.array(starts[:-1], dtype=np.int32), np.array(idxs, dtype=np.int32),
          np.array(vals))
h.run()
st = h.getModelStatus()
print("FINAL (5,14) status:", st, flush=True)
if str(st).endswith("kOptimal"):
    sol = h.getSolution()
    x = np.rint(np.array(sol.col_value)).astype(int)
    print("truth table ones:", [i for i in range(N) if x[i]])
