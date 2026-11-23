#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using namespace std;

using Rows=array<int,9>;
int height(const Rows& r) {
    int f[512][9]{};
    for(int i=0;i<9;i++) f[1<<i][i]=1;
    for(int s=1;s<512;s++) for(int l=s;l;l&=l-1) {
        int i=__builtin_ctz(static_cast<unsigned>(l));
        for(int o=r[i]&~s;o;o&=o-1)
            f[s|(o&-o)][__builtin_ctz(static_cast<unsigned>(o))]+=f[s][i];
    }
    return accumulate(f[511],f[511]+9,0);
}
int main() {
    vector<array<int,3>> lines;
    for(int d:{1,3,4,7}) {
        set<array<int,3>> parallel;
        for(int i=0;i<9;i++) {
            array<int,3> l{};
            for(int k=0;k<3;k++)
                l[k]=((i/3+k*(d/3))%3)*3+(i%3+k*(d%3))%3;
            sort(l.begin(),l.end());
            parallel.insert(l);
        }
        assert(parallel.size()==3);
        lines.insert(lines.end(),parallel.begin(),parallel.end());
    }
    map<pair<int,int>,int> incidences;
    for(const auto& l:lines) for(int i=0;i<3;i++) for(int j=i+1;j<3;j++)
        incidences[{l[i],l[j]}]++;
    assert(incidences.size()==36);
    for(const auto& entry:incidences) assert(entry.second==1);
    map<int,int> bins;
    for(int m=0;m<4096;m++) {
        Rows r{};
        for(int j=0;j<12;j++) for(int k=0;k<3;k++) {
            const auto& l=lines[j];
            r[l[k]]|=1<<l[(k+(((m>>j)&1)?1:2))%3];
        }
        for(int row:r) assert(__builtin_popcount(static_cast<unsigned>(row))==4);
        int h=height(r);
        bins[h]++;
        for(int i=0;i<9;i++) for(int j=i+1;j<9;j++) {
            r[i]^=1<<j; r[j]^=1<<i;
            int delta=height(r)-h;
            assert(delta!=2 && delta!=-2);
            r[i]^=1<<j; r[j]^=1<<i;
        }
    }
    const map<int,int> expected{
        {3159,64},{3243,288},{3255,576},{3275,864},
        {3281,864},{3299,864},{3303,288},{3333,288}};
    assert(bins==expected);
    cout<<"affine-nine exact census and all 147456 exits passed\n";
}
