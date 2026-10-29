#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Exact permutation/cut enumeration, pruned only after three forbidden joins.
// Ports are lexicographic ordered triples. Output has a fixed little-endian format.
struct Census {
    int n=12,d=1320;
    int a[12]{}, ports[12][12][12]{};
    unsigned rows[12]{};
    std::vector<std::uint64_t> matrix;
    std::uint64_t leaves=0,total=0;
    Census(bool flip):matrix(d*d,0) {
        int index=0;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(j!=i)
            for(int k=0;k<n;k++)if(k!=i&&k!=j)ports[i][j][k]=index++;
        for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){
            int b=u/4,c=v/4,q=u%4,r=v%4;
            bool arc=b!=c ? (c-b)%3==1 : q==0 || (q>0 && (r-q+3)%3==1);
            if(b==2&&c==2&&flip)arc=!arc;
            if(arc)rows[u]|=1u<<v;else rows[v]|=1u<<u;
        }
    }
    void add(int b,int c){
        int s=ports[a[0]][a[b]][a[c]];
        int t=ports[a[b-1]][a[c-1]][a[n-1]];
        matrix[s*d+t]++;total++;
    }
    void visit(int depth,unsigned used,int bad,int first,int second){
        if(depth==n){
            leaves++;
            if(bad==2)add(first,second);
            else if(bad==1){
                for(int cut=1;cut<n;cut++)if(cut!=first)
                    add(cut<first?cut:first,cut<first?first:cut);
            }else for(int b=1;b<n;b++)for(int c=b+1;c<n;c++)add(b,c);
            return;
        }
        for(int v=0;v<n;v++)if(!(used>>v&1u)){
            bool forbidden=depth>0 && !(rows[a[depth-1]]>>v&1u);
            if(forbidden&&bad==2)continue;
            a[depth]=v;
            visit(depth+1,used|(1u<<v),bad+forbidden,
                  forbidden&&bad==0?depth:first,forbidden&&bad==1?depth:second);
        }
    }
};
int main(int argc,char**argv){
    if(argc!=3){std::cerr<<"usage: three_run_hostile flip output\n";return 2;}
    Census c(std::string(argv[1])=="1");c.visit(0,0,0,0,0);
    std::ofstream out(argv[2],std::ios::binary);
    for(auto x:c.matrix)for(int byte=0;byte<8;byte++)out.put(char((x>>(8*byte))&255));
    if(!out){std::cerr<<"output failure\n";return 3;}
    std::cout<<"covers "<<c.total<<" permutations "<<c.leaves<<"\n";
    std::cout<<"rows";for(auto row:c.rows)std::cout<<" "<<row;std::cout<<"\n";
}
