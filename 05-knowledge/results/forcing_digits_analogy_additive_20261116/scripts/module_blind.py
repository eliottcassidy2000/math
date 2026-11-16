#!/usr/bin/env python3
"""module_blind.py: history-blind extensions.  A gadget Q (a tournament on {tau} u W) is applied to T by substituting T
for the module vertex tau (every w in W beats all of T or loses to all of T according to the arc between w and tau).
Claim: H(Q[T]) is a function of H(T) alone (over all T) iff tau lies on no directed cycle of Q, and then
H(Q[T]) = H(T) * H(Q - tau) ... more precisely H(Q[T]) = H(Q with tau as a count-H(T) module) = c * H(T) with c = H(Q)
(ordinal composition).  Test: for every gadget Q of order <= 4 with tau = vertex 0, evaluate H(Q[T]) on the two
order-4 tournaments with H = 3 (3-cycle + source, 3-cycle + sink) and the two order-5 tournaments with H = 5, and
on TT_1 (H = 1) vs the 3-cycle: if tau is on a cycle of Q the values differ for some equal-H pair."""
import itertools, subprocess
def ham(adj):
    n = len(adj); N = 1 << n; dp = [[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v] = 1
    for m in range(1, N):
        for l in range(n):
            a = dp[m][l]
            if not a: continue
            for v in range(n):
                if not (m >> v) & 1 and adj[l][v]: dp[m | (1<<v)][v] += a
    return sum(dp[N-1])
def substitute(Q, T):
    # Q adjacency on {0..q-1}, tau = 0; T adjacency on t vertices; result on T + (q-1)
    q = len(Q); t = len(T); n = t + q - 1
    adj = [[0]*n for _ in range(n)]
    for i in range(t):
        for j in range(t): adj[i][j] = T[i][j]
    for w in range(1, q):
        W = t + w - 1
        for i in range(t):
            if Q[0][w]: adj[i][W] = 1
            else: adj[W][i] = 1
        for w2 in range(1, q):
            if w2 != w and Q[w][w2]: adj[W][t + w2 - 1] = 1
    return adj
def parse(s, n):
    adj = [[0]*n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i+1, n):
            if s[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    return adj
def on_cycle(Q, v=0):
    q = len(Q); reach = {v}; frontier = [v]
    while frontier:
        x = frontier.pop(); 
        for y in range(q):
            if Q[x][y] and y not in reach: reach.add(y); frontier.append(y)
    back = {v}; frontier = [v]
    while frontier:
        x = frontier.pop()
        for y in range(q):
            if Q[y][x] and y not in back: back.add(y); frontier.append(y)
    return len((reach & back) - {v}) > 0
# test tournaments grouped by H
tests = {}
for n in (1, 3, 4, 5):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split() if n > 1 else ['']
    for s in out:
        T = parse(s, n) if n > 1 else [[0]]
        tests.setdefault(ham(T), []).append(T)
pairs = [(H, Ts) for H, Ts in tests.items() if len(Ts) >= 2]
print("equal-H test groups:", [(H, len(Ts)) for H, Ts in pairs])
viol = 0; ok = 0; total = 0
for q in (2, 3, 4):
    out = subprocess.run(['gentourng', '-q', str(q)], capture_output=True, text=True).stdout.split()
    for s in out:
        Qc = parse(s, q)
        for perm in itertools.permutations(range(q)):   # every choice of tau (vertex 0 after relabelling)
            Q = [[Qc[perm[i]][perm[j]] for j in range(q)] for i in range(q)]
            total += 1
            cyc = on_cycle(Q)
            blind = True
            for H, Ts in pairs:
                vals = {ham(substitute(Q, T)) for T in Ts}
                if len(vals) > 1: blind = False; break
            mult = all(ham(substitute(Q, T)) == H * ham([[Q[i][j] for j in range(1, q)] for i in range(1, q)]) if q > 1 else True for H, Ts in tests.items() for T in Ts[:1]) if not cyc else None
            if blind == cyc: viol += 1
            else: ok += 1
print(f"gadgets (order 2-4, every tau): {total}; 'count depends only on H' <=> 'tau on no cycle of Q': holds for {ok}, fails for {viol}")
