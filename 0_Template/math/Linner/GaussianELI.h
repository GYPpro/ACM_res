using ld = long double;
// A:[NxN] b[1xN]
//返回：{1,[]}代表有唯一解，{-1,[]}代表无解或无穷解
pair<int,vector<ld>> GaussianELI(vector<vector<ld>> a,vector<ld> b) {
    assert(a.size()),assert(a.size() == a[0].size()),assert(a.size() == b.size());
    int n = a.size();
    vector<int> p(n);
    iota(all(p),0);
    for(int i = 0;i < n;i ++){
        sort(p.begin()+i,p.end(),[&](int x,int y) -> bool {return abs(a[x][i]) > abs(a[y][i]); });
        if(a[p[i]][i] == 0) continue;
        for(int j = i+1;j < n;j ++)
        {
            ld k = a[p[j]][i]/a[p[i]][i] ;
            for(int f = i;f < n;f ++) a[p[j]][f] -= a[p[i]][f] * k;
            b[p[j]] -= b[p[i]] * k;
        }
    }
    vector<ld> res(n);
    for(int i = n-1;i >= 0;i --)
    {
        ld dev = 0;
        for(int f = i+1;f < n;f ++) dev += a[p[i]][f] * res[f];
        if(fabs(a[p[i]][i]) <= 1e-8){ return {-1,{}};} else res[i] = (b[p[i]] - dev)/a[p[i]][i];
    }
    return {f,res};
}