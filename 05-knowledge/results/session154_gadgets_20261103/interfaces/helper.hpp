#pragma once
#include <bits/stdc++.h>
using namespace std;
int hp(const vector<int>& a){int n=a.size(),N=1<<n; vector<array<int,11>> d(N); for(int i=0;i<n;i++)d[1<<i][i]=1;for(int s=1;s<N;s++)for(int i=0;i<n;i++)if(s>>i&1){int p=s^(1<<i);if(p)for(int j=0;j<n;j++)if((p>>j&1)&&(a[j]>>i&1))d[s][i]+=d[p][j];}int z=0;for(int i=0;i<n;i++)z+=d[N-1][i];return z;}
vector<int> ear(vector<int>a,int m){int n=a.size();a.push_back(m);for(int i=0;i<n;i++)if(!(m>>i&1))a[i]|=1<<n;return a;}
