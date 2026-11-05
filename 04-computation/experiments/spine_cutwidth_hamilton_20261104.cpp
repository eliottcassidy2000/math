// spine_cutwidth_hamilton_20261104.cpp -- opus-2026-11-04-S79
// Hamiltonian-path counts of spine tournaments by a bounded-cutwidth fragment transfer.
// A spine tournament on vertices 0 < 1 < ... < n-1 has every arc forward (i -> j for i < j) except a set D of backward
// arcs (j -> i, j > i).  A Hamiltonian path is a permutation whose descents are backward arcs of D and whose ascents are
// not reversed pairs.  The transfer processes the spine left to right and keeps the restriction of the (unknown) path to
// the prefix: a set of vertex-disjoint fragments, each with a start type (GLOBAL START, or "entered later through the
// backward arc k", a crossing arc whose head is the fragment's first vertex) and an end type (GLOBAL END, or OPEN at a
// prefix vertex u, which can be followed later only by vertices v with (v, u) not in D; since all forward arcs exist, u
// matters only through the set of crossing backward arcs with head u).  With at most w backward arcs crossing any cut,
// the number of fragments is at most w + 1 and the state space is finite: this is the standard path-cover dynamic
// programme for digraphs of bounded cutwidth (classical technique), specialised to spine tournaments.
// Modes: "check n trials" (random D at order n, compare with the subset DP; exhaustive over all D for n <= 7),
//        "chain s1 s2 ..." (cutwidth-1 chain of consecutive backward arcs with the given spans, H printed),
//        "census w nmax X" (every D of cutwidth <= w on at most nmax spine vertices with H <= X; prunes by the
//        insertion lemma H(prefix) <= H(T); prints the attained odd values <= X and the first misses).
// Build: g++ -O2 -std=c++17 -o sch spine_cutwidth_hamilton_20261104.cpp
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
using namespace std;
typedef unsigned long long u64;

// ---------- subset DP (brute force reference) ----------
static u64 hamiltonSubsetDP(int n, const vector<vector<int>>& out) {
    vector<u64> dp((size_t)(1 << n) * n, 0); int full = (1 << n) - 1;
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { u64 c = dp[(size_t)S * n + v]; if (!c) continue;
        for (int w : out[v]) if (!(S >> w & 1)) dp[(size_t)(S | (1 << w)) * n + w] += c; }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * n + v]; return h; }
static vector<vector<int>> spineOut(int n, const vector<pair<int,int>>& D) {   // D: (j, i) with j > i meaning j -> i
    vector<vector<char>> back(n, vector<char>(n, 0)); for (auto& a : D) back[a.first][a.second] = 1;
    vector<vector<int>> out(n); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (back[j][i]) out[j].push_back(i); else out[i].push_back(j); } return out; }

// ---------- fragment transfer ----------
// A fragment: start type s (-1 = global start, else index k of the crossing backward arc (tail_k -> head) that will enter it),
// end type e (-2 = global end; else a bitmask over the current crossing arcs: the arcs whose head is the end vertex u,
// i.e. the future tails forbidden as successors of u).  Fragments are kept as a sorted vector for canonical states.
struct Frag { int s; int e; bool operator<(const Frag& o) const { return s != o.s ? s < o.s : e < o.e; } bool operator==(const Frag& o) const { return s == o.s && e == o.e; } };
typedef vector<Frag> State;
static u64 hamiltonCutwidth(int n, const vector<pair<int,int>>& D) {
    // crossing arcs at the cut after vertex t: arcs (j -> i) with i <= t < j.  Arc ids are global indices into D.
    // Represent end masks over global arc ids (D small in practice: use a 64-bit mask; require |D| <= 60).
    if (D.size() > 60) { fprintf(stderr, "too many arcs\n"); exit(1); }
    vector<vector<int>> headArcs(n), tailArcs(n);      // arcs by head i and by tail j
    for (int k = 0; k < (int)D.size(); k++) { headArcs[D[k].second].push_back(k); tailArcs[D[k].first].push_back(k); }
    map<State, u64> cur; cur[State()] = 1;             // empty prefix
    for (int v = 0; v < n; v++) {
        u64 headMaskV = 0; for (int k : headArcs[v]) headMaskV |= 1ULL << k;          // arcs whose head is v (their tails are > v: future)
        u64 tailMaskV = 0; for (int k : tailArcs[v]) tailMaskV |= 1ULL << k;          // arcs whose tail is v (their heads are < v: in the prefix)
        map<State, u64> nxt;
        for (auto& kv : cur) { const State& st = kv.first; u64 c = kv.second;
            // Options for placing v.  v's predecessor in the final path: none (global start), a prefix vertex u = end of an
            // OPEN fragment with v not forbidden ((v, u) not in D: i.e. no arc with tail v and head u, i.e. (fragment end mask & tailMaskV) == 0),
            // or a future vertex via a backward arc into v (pending arc k with head v: k in headMaskV).
            // v's successor: none (global end), a prefix vertex u' = start of a fragment whose pending arc has tail v (k in tailMaskV, the fragment's s == k),
            // or a future vertex (open end at v with mask headMaskV).
            int m = (int)st.size();
            // enumerate predecessor choice: -1 global start (if no fragment has s == -1), fragment index p (open end compatible), or pending arc k in headMaskV
            vector<int> predChoices; // encode: -1 = global start; 0..m-1 = extend fragment p; 1000+k = new fragment entered later by arc k
            bool hasGS = false; for (auto& f : st) if (f.s == -1) hasGS = true;
            if (!hasGS) predChoices.push_back(-1);
            for (int p = 0; p < m; p++) if (st[p].e != -2 && (((u64)st[p].e) & tailMaskV) == 0) predChoices.push_back(p);
            for (int k : headArcs[v]) predChoices.push_back(1000 + k);
            // successor choice: -2 global end (if no fragment has e == -2), fragment index q with s == k for some k in tailArcs[v] (v then precedes that fragment), or open (-3)
            bool hasGE = false; for (auto& f : st) if (f.e == -2) hasGE = true;
            vector<int> succChoices; if (!hasGE) succChoices.push_back(-2); succChoices.push_back(-3);
            for (int q = 0; q < m; q++) for (int k : tailArcs[v]) if (st[q].s == k) succChoices.push_back(q);
            for (int pc : predChoices) for (int sc : succChoices) {
                if (pc >= 0 && pc < 1000 && sc >= 0 && sc == pc) continue;       // cannot close a fragment onto itself
                State ns; // build: remove merged fragments, add the combined one
                Frag nf; // new fragment containing v
                // start part
                if (pc == -1) nf.s = -1; else if (pc >= 1000) nf.s = pc - 1000; else nf.s = st[pc].s;
                // end part
                if (sc == -2) nf.e = -2; else if (sc == -3) nf.e = (int)headMaskV; else nf.e = st[sc].e;
                for (int i = 0; i < m; i++) { if (i == pc || i == sc) continue; ns.push_back(st[i]); }
                ns.push_back(nf); sort(ns.begin(), ns.end());
                // validity: fragments whose pending arc tail is v but were not used must now be impossible (the arc's tail is placed): drop states containing such fragments
                bool ok = true; for (auto& f : ns) { if (f.s >= 0 && (tailMaskV >> f.s & 1) && !(sc >= 0 && st[sc].s == f.s && false)) { /* pending arc with tail v unused */ } }
                // check explicitly: any remaining fragment with s == k where tail(k) == v is dead
                for (auto& f : ns) if (f.s >= 0 && D[f.s].first == v) ok = false;
                if (!ok) continue;
                // open-end masks must only refer to arcs whose tails are still in the future: arcs with tail v are now placed, clear them
                for (auto& f : ns) if (f.e >= 0) f.e = (int)(((u64)f.e) & ~tailMaskV);
                sort(ns.begin(), ns.end()); nxt[ns] += c; }
        }
        cur.swap(nxt);
    }
    // final: exactly one fragment, with global start and global end (or n == 0)
    u64 h = 0; for (auto& kv : cur) { const State& st = kv.first; if (st.size() == 1 && st[0].s == -1 && st[0].e == -2) h += kv.second; }
    return h; }

static int cutwidthOf(int n, const vector<pair<int,int>>& D) { int w = 0; for (int t = 0; t + 1 < n; t++) { int c = 0; for (auto& a : D) if (a.second <= t && a.first > t) c++; w = max(w, c); } return w; }

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "modes: check n trials | chain s1 s2 ... | census w nmax X\n"); return 1; }
    string mode = argv[1];
    if (mode == "check") { int n = atoi(argv[2]); long trials = atol(argv[3]); mt19937_64 rng(20261104); long bad = 0, done = 0; int maxw = 0;
        if (n <= 7) { int P = n * (n - 1) / 2; vector<pair<int,int>> pairs; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pairs.push_back({j, i});
            for (long code = 0; code < (1L << P); code++) { vector<pair<int,int>> D; for (int k = 0; k < P; k++) if (code >> k & 1) D.push_back(pairs[k]);
                u64 a = hamiltonSubsetDP(n, spineOut(n, D)), b = hamiltonCutwidth(n, D); done++; if (a != b) { bad++; if (bad < 5) printf("MISMATCH n=%d code=%ld subset=%llu transfer=%llu\n", n, code, a, b); } maxw = max(maxw, cutwidthOf(n, D)); }
            printf("exhaustive n = %d: %ld tournaments, %ld mismatches, max cutwidth %d\n", n, done, bad, maxw); }
        else { for (long t = 0; t < trials; t++) { vector<pair<int,int>> D; double p = (double)(rng() % 1000) / 1000.0 * 0.5;
                for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if ((double)(rng() % 1000000) / 1000000.0 < p) D.push_back({j, i});
                if (cutwidthOf(n, D) > 12) continue; u64 a = hamiltonSubsetDP(n, spineOut(n, D)), b = hamiltonCutwidth(n, D); done++; if (a != b) { bad++; if (bad < 5) printf("MISMATCH n=%d |D|=%zu subset=%llu transfer=%llu\n", n, D.size(), a, b); } }
            printf("random n = %d: %ld tournaments compared, %ld mismatches\n", n, done, bad); }
        return 0; }
    if (mode == "chain") { vector<int> spans; for (int i = 2; i < argc; i++) spans.push_back(atoi(argv[i])); int n = 1; for (int s : spans) n += s;
        vector<pair<int,int>> D; int pos = 0; for (int s : spans) { D.push_back({pos + s, pos}); pos += s; }
        printf("chain spans"); for (int s : spans) printf(" %d", s); printf(": n = %d, cutwidth %d, H = %llu\n", n, cutwidthOf(n, D), hamiltonCutwidth(n, D)); return 0; }
    if (mode == "census") { int w = atoi(argv[2]), nmax = atoi(argv[3]); u64 X = strtoull(argv[4], 0, 10);
        // depth-first over spine vertices: at vertex v choose the set of backward arcs (v -> i), i < v, keeping cutwidth <= w; prune by H(prefix) <= X.
        // Configurations are canonical up to: a leading vertex with no backward arcs touching it is a source (H unchanged) -- we still enumerate them (cost only).
        vector<char> att(X + 1, 0); u64 configs = 0, pruned = 0; vector<pair<int,int>> D;
        // to keep the enumeration finite and useful: require every vertex except the first to be the tail of at least one backward arc or the head of one
        // (otherwise it is a source/sink relative to its neighbours and only pads); we keep it simple: enumerate all, prune by H.
        struct Rec { static void go(int v, int n, int nmax, int w, u64 X, vector<pair<int,int>>& D, vector<char>& att, u64& configs, u64& pruned) {
            // evaluate current prefix of n vertices (after placing vertex n-1)
            if (n >= 1) { u64 h = hamiltonCutwidth(n, D); configs++; if (h > X) { pruned++; return; } if (h <= X) att[h] = 1; }
            if (n == nmax) return;
            int nv = n; // new vertex index
            // choose subset of heads i < nv for backward arcs (nv -> i), subject to cutwidth <= w on every cut
            // crossing count at cut t after adding arcs: existing crossings + new arcs with i <= t
            vector<int> cross(nv, 0); for (auto& a : D) for (int t = a.second; t < a.first; t++) cross[t]++;
            // enumerate subsets of {0..nv-1} as heads; feasibility: for each t, cross[t] + #(chosen i <= t) <= w
            int m = nv; for (long sub = 0; sub < (1L << m); sub++) { bool ok = true; int acc = 0;
                for (int t = 0; t < m && ok; t++) { if (sub >> t & 1) acc++; if (cross[t] + acc > w) ok = false; }
                if (!ok) continue; size_t before = D.size(); for (int i = 0; i < m; i++) if (sub >> i & 1) D.push_back({nv, i});
                go(v + 1, nv + 1, nmax, w, X, D, att, configs, pruned); D.resize(before); } } };
        Rec::go(0, 0, nmax, w, X, D, att, configs, pruned);
        u64 cnt = 0; vector<u64> miss; for (u64 v = 1; v <= X; v += 2) { if (att[v]) cnt++; else miss.push_back(v); }
        printf("census cutwidth <= %d, n <= %d, X = %llu: %llu configurations evaluated (%llu pruned), %llu odd values attained of %llu; first misses:", w, nmax, X, configs, pruned, cnt, (X + 1) / 2);
        for (size_t i = 0; i < miss.size() && i < 40; i++) printf(" %llu", miss[i]); printf("%s\n", miss.size() > 40 ? " ..." : ""); return 0; }
    return 0;
}
