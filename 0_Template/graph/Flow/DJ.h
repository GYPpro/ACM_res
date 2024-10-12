
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
ll h[5005]={0};    //势数组，dj不能跑负边，所以要操作
struct node
{
    int id;
    bool operator < (const node &x) const {return dist[x.id]<dist[id];}
};
/////////////////////////////////////
bool DJ()
{
    int vis[5005]={0};
    priority_queue<node> q;
    q.push({s});
    for(int i=1;i<=n;i++)dist[i]=inf;
    dist[s]=0;
    while(!q.empty())
    {
        node temp=q.top();
        q.pop();
        vis[temp.id]=0;
        for(int i=head[temp.id];i;i=a[i].nex)
        {
            int v=a[i].to;
            if(dist[v]>dist[temp.id]+a[i].c+h[temp.id]-h[v]&&a[i].w>0)
            {
                dist[v]=dist[temp.id]+a[i].c+h[temp.id]-h[v];
                if(!vis[v])
                {
                    vis[v]=1;
                    q.push({v});
                }
                x[v][0]=temp.id;
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
    cost=(dist[t]-(h[s]-h[t]))*minn;
    min_cost+=cost;
    max_flow+=minn;
    for(int i=1;i<=n;i++)h[i]+=dist[i];
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
    while(DJ())
    {
        EK();
    }
    cout<<max_flow<<" "<<min_cost<<"\n";
}