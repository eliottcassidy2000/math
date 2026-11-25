// Read nauty gentourng's default upper-triangle ASCII format on stdin.
// Enumerate strong sources with three marked transitive inflations and one
// actual edge flip having two unmarked endpoints. Output distinct paired
// coefficient tensors (with one concrete witness per pair) at zero mean constant.
#include <algorithm>
#include <array>
#include <iostream>
#include <set>
#include <vector>
#include <string>
using namespace std;
using Mask = unsigned;
using Polynomial = array<int,8>;
array<Mask,12> rows;
array<int,3> ports;
Mask all_bits, ordinary;
int n;
struct WordCounter {
    Polynomial totals{};
    void visit(Mask remaining, Mask unseen, Mask repeatable, Mask candidates) {
        if (!remaining && !unseen) {
            int tensor_index = 0;
            for (int k=0;k<3;++k) {
                if (!(repeatable & (1u<<ports[k]))) tensor_index |= 1<<k;
            }
            ++totals[tensor_index];
        }
        Mask options = candidates & (remaining|repeatable);
        while (options) {
            Mask bit = options & -options;
            options ^= bit;
            int v = __builtin_ctz(bit);
            Mask again = repeatable;
            if (!(ordinary&bit) && !(unseen&bit)) again ^= bit;
            visit(remaining&~bit, unseen&~bit, again, rows[v]);
        }
    }
    Polynomial evaluate() {
        Mask marked = all_bits ^ ordinary;
        visit(ordinary, marked, marked, all_bits);
        Polynomial ans{};
        for (int ix=0;ix<8;++ix) {
            for (int q=ix;;q=(q-1)&ix) {
                int delta = __builtin_popcount(ix^q);
                ans[q] += totals[ix]*(delta%2 ? -1 : 1)*(1<<delta);
                if (!q) break;
            }
        }
        return ans;
    }
};
bool two_run_caps() {
    for (int p:ports) {
        Mask incoming=all_bits^rows[p]^(1u<<p);
        vector<Mask> supports;
        for (Mask u=rows[p];u;u&=u-1) {
            int a=__builtin_ctz(u);
            for (Mask v=rows[a]&incoming;v;v&=v-1) {
                int b=__builtin_ctz(v);
                Mask support=((1u<<a)|(1u<<b))&ordinary;
                if (!support) return false;
                for (Mask prev:supports) if (!(prev&support)) return false;
                supports.push_back(support);
            }
        }
    }
    return true;
}
bool small_width_needed(const Polynomial& f, const Polynomial& g) {
    int c=f[0]+g[0];
    if (!c) return false; // handled separately by the zero-constant scan
    for (int t=0;t<8;++t) if (abs(f[t]+g[t]) >= (1<<14)) return true;
    for (int t:{1,2,4}) {
        int v=f[t]+g[t];
        if (v && -c%v==0) {
            int ratio=-c/v;
            if (ratio>=2 && !(ratio&(ratio-1))) return true;
        }
    }
    return false;
}
int main(int argc,char**argv) {
    n=stoi(argv[1]);
    all_bits=(1u<<n)-1;
    string s;
    set<pair<Polynomial,Polynomial>> emitted;
    long long sources=0,capped=0,reversals=0,zero=0;
    bool all_constants = argc>2 && string(argv[2])=="all";
    bool small_only = argc>2 && string(argv[2])=="small";
    while (cin>>s) {
        if (int(s.size())!=n*(n-1)/2) return 2;
        ++sources;
        rows.fill(0);
        int t=0;
        for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) {
            if(s[t++]=='1') rows[i]|=1u<<j; else rows[j]|=1u<<i;
        }
        for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) for(int k=j+1;k<n;++k) {
            ports={i,j,k};
            ordinary=all_bits^(1u<<i)^(1u<<j)^(1u<<k);
            if (!two_run_caps()) continue;
            ++capped;
            auto f=WordCounter().evaluate();
            for (Mask outer=ordinary;outer;outer&=outer-1) {
                int a=__builtin_ctz(outer);
                for (Mask inner=outer&(outer-1);inner;inner&=inner-1) {
                    int b=__builtin_ctz(inner);
                    rows[a]^=1u<<b;rows[b]^=1u<<a;
                    if (two_run_caps()) {
                        ++reversals;
                        auto g=WordCounter().evaluate();
                        if (f[0]+g[0]==0) ++zero;
                        if (all_constants || (small_only ? small_width_needed(f,g) : f[0]+g[0]==0)) {
                            if (emitted.insert({f,g}).second) {
                                cout<<s<<" "<<i<<" "<<j<<" "<<k<<" "<<a<<" "<<b;
                                for (int x:f) cout<<" "<<x;
                                for (int x:g) cout<<" "<<x;
                                cout<<"\n";
                            }
                        }
                    }
                    rows[a]^=1u<<b;rows[b]^=1u<<a;
                }
            }
        }
    }
    cerr<<n<<" sources "<<sources<<" marked capped "<<capped
        <<" capped reversals "<<reversals<<" zero constants "<<zero
        <<" distinct pairs "<<emitted.size()<<"\n";
}
