#include <template_overAll.h>

vector<ll> fath ;//: fath[i]为i结点的父节点
vector<ll> deep ;//: deep[i]为i结点的深度
vector<vector<ll>> fa ;//: fa[i][j]为i的第2^i个祖先结点
vector<vector<ll>> g ;//: 邻接表

// 朴素算法
ll lca_p(vector<ll> fath, vector<ll> deep, ll u, ll v)
{
    while (u != v)
    {
        if (deep[u] == deep[v])
        {
            u = fath[u];
            v = fath[v];
        }
        else if (deep[u] > deep[v])
        {
            u = fath[u];
        }
        else if (deep[u] < deep[v])
        {
            v = fath[v];
        }
    }
    return u;
}

// 倍增算法
/*
root: 起始节点
fno: root的父亲节点
*/
// 初始化
void init(){
deep[0] = 0;
for (ll i = 0; i < 31; i++)
    fa[0].push_back(0);
}
// 记录每个节点的祖先和深度
void dfs(vector<vector<ll>> &g, vector<vector<ll>> &fa, vector<ll> &deep, ll fno, ll root) // root为当前结点，fno为root的父节点
{
    fa[root].push_back(fno); // 初始化:root的第一个祖先结点是父节点
    deep[root] = deep[fa[root][0]] + 1;
    for (ll i = 1; i < 31; i++)
    {
        fa[root].push_back(fa[fa[root][i - 1]][i - 1]);
    }
    ll sz = g[root].size();
    for (ll i = 0; i < sz; i++) // 遍历子节点
    {
        if (g[root][i] == fno)
            continue;
        dfs(g, fa, deep, root, g[root][i]);
    }
    return;
}

// fa:lca倍增值 deep:深度 tr:树
auto dfs = [&](int fno,int root) -> void{
    fa[root].push_back(fno);
    deep[root] = deep[fa[root][0]] + 1;
    for (ll i = 1; i < 31; i++)
    {
        fa[root].push_back(fa[fa[root][i - 1]][i - 1]);
    }
    ll sz = tr[root].size();
    for (ll i = 0; i < sz; i++) // 遍历子节点
    {
        if (tr[root][i] == fno)
            continue;
        dfs(root, tr[root][i]);
    }
    return;
};

auto lca = [&](int u,int v) -> int{
    if (deep[u] > deep[v])
    {
        ll swap = u;
        u = v;
        v = swap;
    }
    ll tmp = deep[v] - deep[u];
    for (ll j = 0; tmp; j++, tmp >>= 1)
    {
        if (tmp & 1)
            v = fa[v][j];
    }
    if (u == v)
        return v;
    for (ll j = 30; j >= 0 && v != u; j--)
    {
        if (fa[u][j] != fa[v][j])
        {
            u = fa[u][j];
            v = fa[v][j];
        }
    }
    return fa[u][0]; // 返回最近公共祖先
};

// 寻找u和v的最近公共祖先
ll lca_b(vector<vector<ll>> &fa, vector<ll> &deep, ll u, ll v)
{
    if (deep[u] > deep[v])
    {
        ll swap = u;
        u = v;
        v = swap;
    }
    ll tmp = deep[v] - deep[u];
    for (ll j = 0; tmp; j++, tmp >>= 1)
    {
        if (tmp & 1)
            v = fa[v][j];
    }
    if (u == v)
        return v;
    for (ll j = 30; j >= 0 && v != u; j--)
    {
        if (fa[u][j] != fa[v][j])
        {
            u = fa[u][j];
            v = fa[v][j];
        }
    }
    return fa[u][0]; // 返回最近公共祖先
}

// Tarjan算法