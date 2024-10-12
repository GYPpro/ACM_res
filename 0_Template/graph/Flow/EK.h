
int n,m,s,t;
ll max_flow=0,min_cost=0;
struct p3381
{
    int to;
    int nex;
    ll w;
    ll c;
}a[200005];
int head[5005]={0},cnt=1,x[5005][2]={0};
void add(int u,int v,ll w,ll c)
{
    a[++cnt].to=v;
    a[cnt].nex=head[u];
    a[cnt].w=w;
    a[cnt].c=c;
    head[u]=cnt;

    a[++cnt].to=u;
    a[cnt].nex=head[v];
    a[cnt].w=0;
    a[cnt].c=-c;
    head[v]=cnt;
}
ll dist[5005]={0};   //最小费用数组
/////////////////////////////////////
//基于EK算法的最小费用最大流
//用SPFA对费用求增广路
//用EK算法中的更新操作求解
bool SPFA() //SPFA找增广路
{
    int vis[5005]={0};
    queue<int> q;
    q.push(s);
    for(int i=1;i<=n;i++)dist[i]=inf;
    dist[s]=0;
    while(!q.empty())
    {
        int temp=q.front();
        q.pop();
        vis[temp]=0;
        for(int i=head[temp];i;i=a[i].nex)
        {
            int v=a[i].to;
            if(dist[v]>dist[temp]+a[i].c&&a[i].w>0)
            {
                dist[v]=dist[temp]+a[i].c;
                if(!vis[v])
                {
                    vis[v]=1;
                    q.push(v);
                }
                x[v][0]=temp;
                x[v][1]=i;
            }
        }
    }
    if(dist[t]==inf)return false;
    return true;
}
void EK()//EK算法求解
{
    ll cost=0,minn=inf;
    for(int p=t;p!=s;p=x[p][0])
    {
        minn=min(minn,a[x[p][1]].w);
    }
    for(int p=t;p!=s;p=x[p][0])
    {
        a[x[p][1]].w-=minn;
        a[x[p][1]^1].w+=minn;
    }
    cost=dist[t]*minn;
    min_cost+=cost;
    max_flow+=minn;
}
void solve()
{
    cin>>n>>m>>s>>t;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        ll w,c;
        cin>>u>>v>>w>>c;
        add(u,v,w,c);
    }
    while(SPFA())
    {
        EK();
    }
    cout<<max_flow<<" "<<min_cost<<"\n";
}