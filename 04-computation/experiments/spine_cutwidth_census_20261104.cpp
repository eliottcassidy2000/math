// spine_cutwidth_census_20261104.cpp -- opus-2026-11-04-S79
// All-length census of Hamiltonian-path counts of CONNECTED spine tournaments of cutwidth <= w (every internal cut
// crossed by 1..w backward arcs), values <= X, by a depth-first search over spine vertices that extends the fragment
// transfer of spine_cutwidth_hamilton_20261104.cpp one vertex at a time (the transfer map of the prefix is shared by
// all continuations).  Pruning: H(prefix) <= H(T) (insertion lemma), so a prefix whose complete-path entry exceeds X is
// abandoned.  A word ends exactly when no backward arc is pending (then the tournament is complete and is recorded);
// continuing past an empty cut would be an ordered join, and products are added in post-processing (closure).
// Usage: census w X [maxn]   -> prints attained odd values (connected), closure under products, misses, witnesses.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
using namespace std;
typedef unsigned long long u64;
struct Frag { int s; u64 e; bool operator<(const Frag& o) const { return s != o.s ? s < o.s : e < o.e; } };
typedef vector<Frag> State;
static const int GSTART = -1; static const u64 GEND = ~0ULL;      // end sentinel
struct Arc { int head, tail; };                                     // tail -> head, head < tail; tail = -1 while pending
static int W; static u64 X; static int MAXN;
static vector<Arc> D; static vector<unsigned char> att; static vector<string> witness; static u64 nodes = 0, recorded = 0;
static map<State, u64> step(const map<State, u64>& cur, int v, u64 headMaskV, u64 tailMaskV) {
    map<State, u64> nxt;
    vector<int> headArcsV, tailArcsV; for (int k = 0; k < (int)D.size(); k++) { if (headMaskV >> k & 1) headArcsV.push_back(k); if (tailMaskV >> k & 1) tailArcsV.push_back(k); }
    for (auto& kv : cur) { const State& st = kv.first; u64 c = kv.second; int m = (int)st.size();
        bool hasGS = false, hasGE = false; for (auto& f : st) { if (f.s == GSTART) hasGS = true; if (f.e == GEND) hasGE = true; }
        vector<int> pred; if (!hasGS) pred.push_back(-1);
        for (int p = 0; p < m; p++) if (st[p].e != GEND && (st[p].e & tailMaskV) == 0) pred.push_back(p);
        for (int k : headArcsV) pred.push_back(1000 + k);
        vector<int> succ; if (!hasGE) succ.push_back(-2); succ.push_back(-3);
        for (int q = 0; q < m; q++) for (int k : tailArcsV) if (st[q].s == k) succ.push_back(q);
        for (int pc : pred) for (int sc : succ) {
            if (pc >= 0 && pc < 1000 && sc == pc) continue;
            Frag nf; nf.s = (pc == -1) ? GSTART : (pc >= 1000 ? pc - 1000 : st[pc].s);
            nf.e = (sc == -2) ? GEND : (sc == -3 ? headMaskV : st[sc].e);
            State ns; for (int i = 0; i < m; i++) if (i != pc && i != sc) ns.push_back(st[i]); ns.push_back(nf);
            bool ok = true; for (auto& f : ns) if (f.s >= 0 && D[f.s].tail == v) { ok = false; break; }   // an arc ending at v entered nobody's fragment
            if (!ok) continue;
            for (auto& f : ns) if (f.e != GEND) f.e &= ~tailMaskV;
            sort(ns.begin(), ns.end()); nxt[ns] += c; } }
    return nxt; }
static u64 completeEntry(const map<State, u64>& cur) { for (auto& kv : cur) { const State& st = kv.first; if (st.size() == 1 && st[0].s == GSTART && st[0].e == GEND) return kv.second; } return 0; }
static string describe() { string s; for (auto& a : D) s += "(" + to_string(a.tail) + ">" + to_string(a.head) + ")"; return s; }
// DFS: place vertex v given the transfer map of the prefix; pending = arcs with tail == -1.
static void dfs(int v, const map<State, u64>& cur) {
    nodes++;
    vector<int> pending; for (int k = 0; k < (int)D.size(); k++) if (D[k].tail < 0) pending.push_back(k);
    int c = (int)pending.size();
    // choose subset T of pending arcs to terminate at v (tails = v), and b new arcs with head v
    for (int sub = 0; sub < (1 << c); sub++) { int t = __builtin_popcount(sub);
        for (int b = 0; c - t + b <= W; b++) {
            if (v == 0 && b == 0) continue;                      // first vertex must start an arc (n >= 2 connected); n = 1 handled by caller
            int after = c - t + b;
            u64 tailMaskV = 0; for (int i = 0; i < c; i++) if (sub >> i & 1) { D[pending[i]].tail = v; tailMaskV |= 1ULL << pending[i]; }
            u64 headMaskV = 0; size_t before = D.size(); for (int j = 0; j < b; j++) { headMaskV |= 1ULL << D.size(); D.push_back({v, -1}); }
            if (D.size() <= 62) {
                map<State, u64> nxt = step(cur, v, headMaskV, tailMaskV);
                u64 h = completeEntry(nxt);
                if (after == 0) { recorded++; if (h <= X && !att[h]) { att[h] = 1; witness[h] = to_string(v + 1) + ":" + describe(); } }
                else if (h <= X && v + 1 < MAXN) dfs(v + 1, nxt);
            }
            D.resize(before); for (int i = 0; i < c; i++) if (sub >> i & 1) D[pending[i]].tail = -1; } } }
int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: census w X [maxn]\n"); return 1; }
    W = atoi(argv[1]); X = strtoull(argv[2], 0, 10); MAXN = argc > 3 ? atoi(argv[3]) : 60;
    att.assign(X + 1, 0); witness.assign(X + 1, ""); att[1] = 1; witness[1] = "1:transitive";
    map<State, u64> start; start[State()] = 1; dfs(0, start);
    u64 cntc = 0; for (u64 v = 1; v <= X; v += 2) if (att[v]) cntc++;
    printf("cutwidth <= %d, X = %llu, maxn = %d: %llu DFS nodes, %llu complete connected words; connected odd values attained: %llu of %llu\n", W, X, MAXN, nodes, recorded, cntc, (X + 1) / 2);
    // multiplicative closure (ordered joins of connected pieces)
    vector<unsigned char> cl(att.begin(), att.end()); for (u64 a = 3; a <= X; a += 2) if (cl[a]) for (u64 b = 3; a * b <= X; b += 2) if (att[b]) cl[a * b] = 1;
    u64 cntcl = 0; vector<u64> miss; for (u64 v = 1; v <= X; v += 2) { if (cl[v]) cntcl++; else miss.push_back(v); }
    printf("with ordered joins: %llu of %llu odd values <= %llu attained; %zu misses; first misses:", cntcl, (X + 1) / 2, X, miss.size());
    for (size_t i = 0; i < miss.size() && i < 80; i++) printf(" %llu", miss[i]); printf("%s\n", miss.size() > 80 ? " ..." : "");
    printf("misses per dyadic scale:"); for (int k = 2; (1ULL << k) <= X; k++) { u64 lo = 1ULL << k, hi = min(X, (1ULL << (k + 1)) - 1), cc = 0, tot = 0; for (u64 v = lo | 1; v <= hi; v += 2) { tot++; if (!cl[v]) cc++; } printf(" [2^%d:%llu/%llu]", k, cc, tot); } printf("\n");
    printf("connected witnesses (odd values <= 99, first found; format n:(tail>head)...):\n"); for (u64 v = 1; v <= 99 && v <= X; v += 2) if (att[v]) printf("  %llu = %s\n", v, witness[v].c_str());
    return 0; }
