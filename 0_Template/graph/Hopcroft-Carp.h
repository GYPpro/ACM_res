/**************************
* 二分图匹配(Hopcroft-Carp算法)
* 复杂度 O (sqrt(n) * E) 
* 邻接表存图 
* Nx 为左端的端点数,使用前先赋值
**************************/ 
const int N = 3005, INF = 0x3f3f3f3f;
vector<int> G[N];  //邻接表
int Nx,Ny,k;  //x、y部的点数、k边数
int Mx[N],My[N];  //记录匹配的点的序号
int dx[N],dy[N];  //标记数组
int dis,u,v;
bool used[N];
bool bfs()
{
    queue<int> Q;
    dis = INF;  //初始化增广路的最短长度
    memset(dx,-1,sizeof(dx)); //初始化标记
    memset(dy,-1,sizeof(dy));
    for(int i=1;i<=Nx;++i)
    {
        if(Mx[i] == -1)
        {
            Q.push(i), dx[i] = 0;//将x部未匹配的点加入到队列中
        }
    }
    while(!Q.empty())
    {
        int u = Q.front();Q.pop();
        if(dx[u] > dis) break;
        for(auto v:G[u])
        {
            if(dy[v] == -1)   //取未标记过的y部的点
            {
                dy[v] = dx[u] + 1;  //做上标记
                if(My[v] == -1) dis = dy[v];  //若点v未匹配，此即为最短增广路
                else dx[My[v]] = dy[v] + 1, Q.push(My[v]);  //若已经匹配，则将v的匹配点纳入到增广路，并入队
            }
        }
    }
    return dis != INF;  //若dis==inf，说明无增广路
}
bool DFS(int u)
{
    for(auto v:G[u])
    {
        if(!used[v] && dy[v] == dx[u] + 1)//点u只能和增广路上的下一个点匹配
        {//used[v]表示这个点已被查询过，它要么已经匹配了，要么不可匹配
            used[v] = true;
            if(My[v] != -1 && dy[v] == dis) continue; //若v匹配过，但是是增广路的终点，则忽略它
            if(My[v] == -1 || DFS(My[v]))//若点v未匹配，或者点v的匹配点可以找到新的匹配点，则匹配u和v
            {
                My[v] = u;
                Mx[u] = v;
                return true;  //匹配成功就返回
            }
        } 
    }
    return false;
}
int MaxMatch()
{
    int res = 0;
    memset(Mx,-1,sizeof(Mx));
    memset(My,-1,sizeof(My));
    while(bfs())
    {
        memset(used,false,sizeof(used));
        for(int i = 1;i <=Nx;++i)
            if(Mx[i] == -1 && DFS(i)) //取x部未匹配的点进行深搜
                ++res;   //若成功匹配则做出贡献
    }
    return res;
}
void solve()
{
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    cin>>Nx>>Ny>>k;
    while(k--)
    {
        cin>>u>>v;
        G[u].push_back(v);
    }
    cout<<MaxMatch()<<"\n";
}