#include <array>
#include <cstring>
#include <iostream>
#include <vector>
using Rows = std::array<int,9>;
using Runs = std::array<long long,5>;

Runs root_runs(const Rows& rows) {
    static long long dp[256][9][5];
    std::memset(dp,0,sizeof(dp));
    dp[0][8][1]=1;
    for (int i=0;i<8;i++) dp[1<<i][i][0]=1;
    for (int mask=0;mask<256;mask++) {
        for (int k=0;k<5;k++) {
            for (int i=0;i<8;i++) {
                long long count=dp[mask][i][k];
                if (!count) continue;
                int next=rows[i]&255&~mask;
                while (next) {
                    int j=__builtin_ctz(next); next&=next-1;
                    dp[mask|(1<<j)][j][k]+=count;
                }
                if ((rows[i]&(1<<8)) && k<4) dp[mask][8][k+1]+=count;
            }
            long long count=dp[mask][8][k];
            int next=rows[8]&255&~mask;
            while (count && next) {
                int j=__builtin_ctz(next); next&=next-1;
                dp[mask|(1<<j)][j][k]+=count;
            }
        }
    }
    Runs result{};
    for (int i=0;i<9;i++) for (int k=0;k<5;k++) result[k]+=dp[255][i][k];
    return result;
}

int main() {
    int cases; if (!(std::cin>>cases)) return 1;
    for (int case_index=0;case_index<cases;case_index++) {
        int r,c,number; std::cin>>r>>c>>number;
        if (r+c!=8 || r<2 || r>4) return 2;
        std::vector<std::vector<int>> patterns(number,std::vector<int>(r));
        for (auto& p:patterns) for (auto& row:p) std::cin>>row;
        std::vector<std::array<int,2>> inner;
        for (int i=0;i<8;i++) for (int j=i+1;j<8;j++)
            if ((i<r)==(j<r)) inner.push_back({i,j});
        long long tested=0,hits=0;
        for (int pattern=0;pattern<number;pattern++) {
            for (int mask=0;mask<(1<<inner.size());mask++) {
                Rows rows{}; rows[8]=(1<<r)-1;
                for (int i=r;i<8;i++) rows[i]|=1<<8;
                for (int i=0;i<r;i++) for (int j=r;j<8;j++) {
                    if (patterns[pattern][i]&(1<<(j-r))) rows[i]|=1<<j;
                    else rows[j]|=1<<i;
                }
                for (int b=0;b<(int)inner.size();b++) {
                    int i=inner[b][0],j=inner[b][1];
                    if (mask&(1<<b)) rows[i]|=1<<j;
                    else rows[j]|=1<<i;
                }
                Runs a=root_runs(rows); tested++;
                if (a[4]) return 3;
                long long amplitude=a[3],beta=a[1]-3*amplitude;
                if (amplitude>0 && a[2]==3*amplitude && (beta==1 || beta==-1)) {
                    hits++;
                    std::cout<<"UNIT "<<r<<' '<<pattern<<' '<<mask<<' '<<amplitude<<' '<<beta;
                    for (int row:rows) std::cout<<' '<<row;
                    std::cout<<'\n';
                }
            }
        }
        std::cout<<"CASE "<<r<<' '<<c<<' '<<number<<' '<<tested<<' '<<hits<<'\n';
    }
}
