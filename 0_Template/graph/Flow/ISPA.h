
int n,m;
struct p3376
{
    int to;
    ll w=0;
    int nex;
}a[540005];
map<pair<int,int>,int> key;
int head[1580]={0},cnt=1,dept[1580]={0};
void add(int u,int v,ll w)
{
    if(key[{u,v}]!=0)     //去除重边
    {
        a[key[{u,v}]].w+=w;
        return;
    }
    a[++cnt].to=v;
    a[cnt].w=w;
    a[cnt].nex=head[u];
    key[{u,v}]=cnt;
    head[u]=cnt;
    a[++cnt].to=u;
    a[cnt].w=0;
    a[cnt].nex=head[v];
    key[{v,u}]=cnt;
    head[v]=cnt;
}
//..........................//
//Dinic算法会调用太多次bfs，考虑优化
//于是就有了ISPA算法
//从汇点开始bfs，标记深度
//再从源点开始dfs，对经过的点进行深度修改，当出现断层时，算法结束
ll ans=0;
int g[1580]={0},maxn=0;

void bfs(int t)  //从汇点bfs
{
    queue<int> q;
    bool vis[1580];
    memset(vis,0,sizeof(vis));
    q.push(t);
    vis[t]=true;
    //g[0]++;
    while(!q.empty())
    {
        int temp=q.front();
        q.pop();
        g[dept[temp]]++;
        maxn=max(dept[temp],maxn);
        //vis[temp]=true;
        for(int i=head[temp];i;i=a[i].nex)
        {
            int v=a[i].to;
            if(!vis[v])
            {
                vis[v]=true;
                dept[v]=dept[temp]+1;
                q.push(v);
            }
        }
    }
}
ll dfs(int p,int s,int t,ll tot)
{
    int i;
    if(p==t)  //到达汇点
    {
        ans+=tot;  //更新答案
        return tot; //返回
    }
    ll k=0,res=0;
    for(i=head[p];i&&tot;i=a[i].nex)
    {
        int v=a[i].to;
        if(dept[v]!=dept[p]-1||a[i].w<=0)continue;  //残量为0，或非递减的深度，跳过
        k=dfs(v,s,t,min(a[i].w,tot)); //取p->v的流量  

        a[i].w-=k; //更新残量
        a[i^1].w+=k;

        res+=k;  //增加总流量

        tot-=k;  //残量减少
    }
    if(tot>0)   //流过来的流量还有余，使深度增加
    {
        dept[p]++;
        g[dept[p]-1]--;
        g[dept[p]]++;
        if(g[dept[p]-1]<=0)dept[s]=n+1;    //出现断层，对s的深度进行标记
    }
    return res;
}
void isap(int s,int t)
{
    while(dept[s]<=n)dfs(s,s,t,100000000000);  //未出现断层就反复dfs
}
void solve()
{
    int s,t;
    cin>>n>>m>>s>>t;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        ll w;
        cin>>u>>v>>w;
        add(u,v,w);
    }
    bfs(t);
    isap(s,t);
    cout<<ans<<"\n";
}