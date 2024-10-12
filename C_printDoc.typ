#import "@preview/codelst:2.0.1": sourcecode
// Display inline code in a box
#set text(font:("Times New Roman","Source Han Serif SC"))
#show raw.where(block: false): box.with(
  fill: luma(230),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
)
#show raw.where(block: true): block.with(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
)
#show raw: set text(
    font: ("consolas", "Source Han Serif SC")
  )
#set page(
  // flipped: true,
  // background: [#image("background.png")]
  paper: "a4",
)
#set text(
    font:("Times New Roman","Source Han Serif SC"),
    style:"normal",
    weight: "regular",
    size: 13pt,
)
#show math.equation:set text(font:("New Computer Modern Math","Source Han Serif SC"))
#let nxtIdx(name) = box[ #counter(name).step()#counter(name).display()]
#set math.equation(numbering: "(1)")

#outline(depth: 2)

#set page(
  paper:"a4",
  number-align: right,
  margin: (x:2cm,y:2.5cm),
  header: [
    #box(baseline:5pt)[#set text(
      size: 11pt,
    )
    #align(
      left+bottom,
      [
        #smallcaps[Templetes]
        #h(1fr)#text("Github GYPpro/Acm_res",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)

= #smallcaps[Basic]

== `bitset.h`


 #sourcecode[```cpp
#include<bitset>
int n;
auto bin = bitset<32>(n);

cout << bin;
//输出二进制位

cout << bin.to_ullong();
//输出十进制

bin[1] = 1
//随机访问

bin = !bin ^ (bin & bin | bitset<32>(1))
//位运算

bin != bin
//比较运算符

bin.count()
//1的数量

bin.test(i)
//随机访问，类似std::vector::pos()

bin.any()
//有一位1就true

bin.none()
//全0就返回true

bin.all()
//全1就返回true

bin.flip()
//翻转全部

bin.flip(i)
//a[i] = !a[i]

bin._Find_first()
//第一个1的下标

bin._Find_next(i)
//从下标n往后第一个1的下标
```]
 #pagebreak() 
== `pbds.h`


 #sourcecode[```cpp
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
using ord_set = tree<int, null_type, less<int>, rb_tree_tag,
tree_order_statistics_node_update>;
using ord_mset =  tree<int, null_type, less_equal<int>, rb_tree_tag,
tree_order_statistics_node_update>;
//find_by_order
//order_of_key
```]
 #pagebreak() 
= #smallcaps[Ds]

== `bst.h`


 #sourcecode[```cpp

const int N=1e5+100;
struct node{
	int s[2];
	int v,p,cnt,sz;
	void init(int p1,int v1){
		p=p1;v=v1;
		cnt=sz=1;
	}
}tr[N];
int root=0,idx=0;
void pushup(int x){
	tr[x].sz=tr[x].cnt+tr[tr[x].s[1]].sz+tr[tr[x].s[0]].sz;
}
void rotate(int x){
	int y=tr[x].p;
	int z=tr[y].p;
	int k=tr[y].s[1]==x;
	tr[y].s[k]=tr[x].s[k^1];
	tr[tr[x].s[k^1]].p=y;
	tr[z].s[tr[z].s[1]==y]=x;
	tr[x].p=z;
	tr[y].p=x;
	tr[x].s[k^1]=y;
	pushup(y);pushup(x);
}
void splay(int x,int k){s
	while(tr[x].p!=k){
		int y=tr[x].p;
		int z=tr[y].p;
		if(z!=k) (tr[y].s[0]==x)^(tr[z].s[0]==y) ? rotate(x) : rotate(y);
		rotate(x);
	}
	if(k==0) root=x;
}
void find(int v){
	int x=root;
	while(tr[x].v!=v && tr[x].s[v>tr[x].v] ) x=tr[x].s[v>tr[x].v];
	splay(x,0);
}
int get_pre(int v){
	find(v);
	int x=root;
	if(tr[x].v<v) return x;
	x=tr[x].s[0];
	while(tr[x].s[1]) x=tr[x].s[1];
	splay(x,0);
	return x;
}
int get_suc(int v){
	find(v);
	int x=root;
	if(tr[x].v>v) return x;
	x=tr[x].s[1];
	while(tr[x].s[0]) x=tr[x].s[0];
	splay(x,0);
	return x;
}
void del(int v){
	int pre=get_pre(v);
	int suc=get_suc(v);
	splay(pre,0);splay(suc,pre);
	int d=tr[suc].s[0];
	if(tr[d].cnt>1){
		tr[d].cnt--;splay(d,0);
	}
	else{
		tr[suc].s[0]=0;splay(suc,0);
	}
}
void insert(int v){
	int x=root;
	int p=0;
	while(x && tr[x].v!=v){
		p=x;x=tr[x].s[v>tr[x].v];
	}
	if(x) tr[x].cnt++;
	else{
		x=++idx;
		tr[p].s[v>tr[p].v]=x;
		tr[x].init(p,v);
	}
	splay(x,0);
}
int get_rank(int v){
	insert(v);
	int res=tr[tr[root].s[0]].sz;
	del(v);
	return res;
}
int get_val(int k){
	int x=root;
	while(1){
		int y=tr[x].s[0];
		if(tr[x].cnt+tr[y].sz<k){
			k-=tr[y].sz+tr[x].cnt;
			x=tr[x].s[1];
		}
		else{
			if(tr[y].sz>=k) x=tr[x].s[0];
			else break;
		}
	}
	splay(x,0);
	return tr[x].v;
}
```]
 #pagebreak() 
== `dsu.h`


 #sourcecode[```cpp

class DSU {
    std::vector<int> f, siz;
public:
    DSU() {}
    DSU(int n) {
        init(n);
    }
    
    void init(int n) {
        f.resize(n);
        for(int i = 0;i < n;i ++) f[i] = i;
        siz.assign(n, 1);
    }
    
    int find(int x) {
        while (x != f[x]) {
            x = f[x] = f[f[x]];
        }
        return x;
    }
    
    bool same(int x, int y) {
        return find(x) == find(y);
    }
    
    bool merge(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) {
            return false;
        }
        siz[x] += siz[y];
        f[y] = x;
        return true;
    }
    
    int size(int x) {
        return siz[find(x)];
    }
};


```]
 #pagebreak() 
== `segTree_add.h`


 #sourcecode[```cpp

// AC 带懒惰标记线段树 
template <class TYPE_NAME>
class lazyTree
{
    /*
     *  TYPE_NAME需要支持：+ += != 和自定义的合并运算符
     *  实现了懒惰标记，仅支持区间批量增加
     *  区间批量减需要TYPE_NAME支持-，且有-a = 0 - a
     *  额外处理了一个单点修改为对应值的函数，非原生实现，其复杂度为 4logn
     *  不提供在线
     *  不提供持久化
     */
private:
    vector<TYPE_NAME> d;
    vector<TYPE_NAME> a;
    vector<TYPE_NAME> b;
    int n;
    const TYPE_NAME INI = 0; // 不会影响合并运算的初始值，如max取INF，min取0，mti取1

    void subbuild(int s, int t, int p)
    {
        if (s == t)
        {
            d[p] = a[s];
            return;
        }
        int m = s + ((t - s) >> 1); //  (s+t)/2
        subbuild(s, m, p * 2);
        subbuild(m + 1, t, p * 2 + 1);
        d[p] = d[p * 2] + d[p * 2 + 1];
        //    合并运算符 ↑
    }

    TYPE_NAME subGetSum(int l, int r, int s, int t, int p)
    {
        if (l <= s && t <= r)
            return d[p];
        int m = s + ((t - s) >> 1);
        if (b[p] != 0)
        {
            d[p * 2] += b[p] * (m - s + 1); // 合并运算符的高阶运算 此处运算为应用懒惰标记
            d[p * 2 + 1] += b[p] * (t - m); // 合并运算符的高阶运算 此处运算为应用懒惰标记
            b[p * 2] += b[p];               // 下传标记，无需修改
            b[p * 2 + 1] += b[p];           // 下传标记，无需修改
            b[p] = 0;
        }
        TYPE_NAME ansl = INI;
        TYPE_NAME ansr = INI;
        if (l <= m)
            ansl = subGetSum(l, r, s, m, p * 2);
        if (r > m)
            ansr = subGetSum(l, r, m + 1, t, p * 2 + 1);
        return ansl + ansr;
        // 合并运算符↑
    }

    void subUpdate(int l, int r, TYPE_NAME c, int s, int t, int p)
    {
        if (l <= s && t <= r)
        {
            d[p] += (t - s + 1) * c; // 合并运算符的高阶运算 此处运算为修改整匹配区间值
            b[p] += c;               // 记录懒惰标记，无需修改
            return;
        }
        int m = s + ((t - s) >> 1);
        if (b[p] != 0 && s != t)
        {
            d[p * 2] += b[p] * (m - s + 1); // 合并运算符的高阶运算 此处运算为应用懒惰标记
            d[p * 2 + 1] += b[p] * (t - m); // 合并运算符的高阶运算 此处运算为应用懒惰标记
            b[p * 2] += b[p];               // 下传标记，无需修改
            b[p * 2 + 1] += b[p];           // 下传标记，无需修改
            b[p] = 0;
        }
        if (l <= m)
            subUpdate(l, r, c, s, m, p * 2);
        if (r > m)
            subUpdate(l, r, c, m + 1, t, p * 2 + 1);
        d[p] = d[p * 2] + d[p * 2 + 1];
        //    合并运算符 ↑
    }

public:
    lazyTree(TYPE_NAME _n)
    {
        n = _n;
        d.resize(4 * n + 5);
        a.resize(4 * n + 5);
        b.resize(4 * n + 5);
    }

    void build(vector<TYPE_NAME> _a)
    {
        a = _a;
        subbuild(1, n, 1);
    }

    TYPE_NAME getsum(int l, int r)
    {
        return subGetSum(l, r, 1, n, 1);
    }

    void update(int l, int r, TYPE_NAME c) // 修改区间
    {
        subUpdate(l, r, c, 1, n, 1);
    }

    void update(int idx, TYPE_NAME tar)
    { // 修改单点，未测试
        TYPE_NAME tmp = getsum(idx, idx);
        tar -= tmp;
        subUpdate(idx, idx, tar, 1, n, 1);
    }
};
```]
 #pagebreak() 
== `segTree_add_setto.h`


 #sourcecode[```cpp
template <typename T>
class SegTreeLazyRangeSet {
  vector<T> tree, lazy;
  vector<T> *arr;
  vector<bool> ifLazy;
  int n, root, n4, end;

  void maintain(int cl, int cr, int p) {
    int cm = cl + (cr - cl) / 2;
    if (cl != cr && ifLazy[p]) {
      lazy[p * 2] = lazy[p],ifLazy[p*2] = 1;
      lazy[p * 2 + 1] = lazy[p],ifLazy[p*2+1] = 1;
      tree[p * 2] = lazy[p] * (cm - cl + 1);
      tree[p * 2 + 1] = lazy[p] * (cr - cm);
      lazy[p] = 0;
      ifLazy[p] = 0;
    }
  }

  T range_sum(int l, int r, int cl, int cr, int p) {
    if (l <= cl && cr <= r) return tree[p];
    int m = cl + (cr - cl) / 2;
    T sum = 0;
    maintain(cl, cr, p);
    if (l <= m) sum += range_sum(l, r, cl, m, p * 2);
    if (r > m) sum += range_sum(l, r, m + 1, cr, p * 2 + 1);
    return sum;
  }

  void range_set(int l, int r, T val, int cl, int cr, int p) {
    if (l <= cl && cr <= r) {
      lazy[p] = val;
      ifLazy[p] = 1;
      tree[p] = (cr - cl + 1) * val;
      return;
    }
    int m = cl + (cr - cl) / 2;
    maintain(cl, cr, p);
    if (l <= m) range_set(l, r, val, cl, m, p * 2);
    if (r > m) range_set(l, r, val, m + 1, cr, p * 2 + 1);
    tree[p] = tree[p * 2] + tree[p * 2 + 1];
  }

  void build(int s, int t, int p) {
    if (s == t) {
      tree[p] = (*arr)[s];
      return;
    }
    int m = s + (t - s) / 2;
    build(s, m, p * 2);
    build(m + 1, t, p * 2 + 1);
    tree[p] = tree[p * 2] + tree[p * 2 + 1];
  }

 public:
  explicit SegTreeLazyRangeSet<T>(vector<T> v) {
    n = v.size();
    n4 = n * 4;
    tree = vector<T>(n4, 0);
    lazy = vector<T>(n4, 0);
    ifLazy = vector<bool>(n4,0);
    arr = &v;
    end = n - 1;
    root = 1;
    build(0, end, 1);
    arr = nullptr;
  }

  void show(int p, int depth = 0) {
    if (p > n4 || tree[p] == 0) return;
    show(p * 2, depth + 1);
    for (int i = 0; i < depth; ++i) putchar('\t');
    printf("%d:%d\n", tree[p], lazy[p]);
    show(p * 2 + 1, depth + 1);
  }

  T range_sum(int l, int r) { return range_sum(l, r, 0, end, root); }

  void range_set(int l, int r, T val) { range_set(l, r, val, 0, end, root); }
};  
```]
 #pagebreak() 
== `segTree_mul_add.h`


 #sourcecode[```cpp
/*
三种操作：（模数为mod） 
1.区间乘val 
2.区间加val
3.区间和 
*/
#define int long long 
#define lc p<<1
#define rc p<<1|1
const int N=1e5+10;
int mod=1e9+7;
struct node{
	int l,r,sum,mul,add;
}tr[4*N];
vector<int> a(N);
void pushup(int p){
	tr[p].sum=tr[lc].sum+tr[rc].sum;
	tr[p].sum%=mod;
}
void build(int p,int l,int r){
	tr[p].l=l;tr[p].r=r; 
	tr[p].mul=1;tr[p].add=0;
	if(l==r){
		tr[p].sum=a[l];
		return;
	}
	int m=l+r>>1;
	build(lc,l,m);
	build(rc,m+1,r);
	pushup(p);
}
void pushdown(int p){
	tr[lc].sum=(tr[lc].sum*tr[p].mul+tr[p].add*(tr[lc].r-tr[lc].l+1))%mod;
	tr[rc].sum=(tr[rc].sum*tr[p].mul+tr[p].add*(tr[rc].r-tr[rc].l+1))%mod;
	tr[lc].mul=tr[lc].mul*tr[p].mul%mod;
	tr[rc].mul=tr[rc].mul*tr[p].mul%mod;
	tr[lc].add=(tr[lc].add*tr[p].mul+tr[p].add)%mod;
	tr[rc].add=(tr[rc].add*tr[p].mul+tr[p].add)%mod;
	tr[p].mul=1;tr[p].add=0;
}
void update_mul(int p,int l,int r,int val){
	if(l<=tr[p].l&&tr[p].r<=r){
		tr[p].sum=tr[p].sum*val%mod;
		tr[p].mul=tr[p].mul*val%mod;
		tr[p].add=tr[p].add*val%mod;
		return;
	}
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	if(l<=m) update_mul(lc,l,r,val);
	if(m<r) update_mul(rc,l,r,val);
	pushup(p);
	return;
}
void update_add(int p,int l,int r,int val){
	if(l<=tr[p].l&&tr[p].r<=r){
		tr[p].add=(tr[p].add+val)%mod;
		tr[p].sum=(tr[p].sum+val*(tr[p].r-tr[p].l+1))%mod;
		return;
	}
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	if(l<=m) update_add(lc,l,r,val);
	if(m<r) update_add(rc,l,r,val);
	pushup(p);
	return;
}
int query(int p,int l,int r){
	if(l<=tr[p].l&&tr[p].r<=r) return tr[p].sum;
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	int sum=0;
	if(l<=m) sum+=query(lc,l,r);
	if(m<r) sum+=query(rc,l,r);
	return sum%mod;
}
void solve(){
	int n,m;
	cin>>n>>m>>mod;
	for(int i=1;i<=n;i++) cin>>a[i];
	build(1,1,n);
	while(m--){
		int op;
		cin>>op;
		if(op==1){
			int x,y,val;
			cin>>x>>y>>val;
			update_mul(1,x,y,val);
		}
		else if(op==2){
			int x,y,val;
			cin>>x>>y>>val;
			update_add(1,x,y,val);
		}
		else{
			int x,y;
			cin>>x>>y;
			cout<<query(1,x,y)<<"\n";
		}
	}
	return;
}

```]
 #pagebreak() 
== `segTree_MX_MI.h`


 #sourcecode[```cpp

//AC MJ的MIN/MAX树
template<class Info>
struct SegmentTree {
    int n;
    std::vector<Info> info;
    SegmentTree() : n(0) {}
    SegmentTree(int n_, Info v_ = Info()) {
        init(n_, v_);
    }
    template<class T>
    SegmentTree(std::vector<T> init_) {
        init(init_);
    }
    void init(int n_, Info v_ = Info()) {
        init(std::vector(n_, v_));
    }
    template<class T>
    void init(std::vector<T> init_) {
        n = init_.size();
        info.assign(4 << std::__lg(n), Info());
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init_[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = info[2 * p] + info[2 * p + 1];
    }
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
    template<class F>
    int findFirst(int p, int l, int r, int x, int y, F &&pred) {
        if (l >= y || r <= x) {
            return -1;
        }
        if (l >= x && r <= y && !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        int res = findFirst(2 * p, l, m, x, y, pred);
        if (res == -1) {
            res = findFirst(2 * p + 1, m, r, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findFirst(int l, int r, F &&pred) {
        return findFirst(1, 0, n, l, r, pred);
    }
    template<class F>
    int findLast(int p, int l, int r, int x, int y, F &&pred) {
        if (l >= y || r <= x) {
            return -1;
        }
        if (l >= x && r <= y && !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        int res = findLast(2 * p + 1, m, r, x, y, pred);
        if (res == -1) {
            res = findLast(2 * p, l, m, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findLast(int l, int r, F &&pred) {
        return findLast(1, 0, n, l, r, pred);
    }
};
const int inf = 1E9;
struct Info
{
    int mn {inf}, mnId, mx {-inf}, mxId;
} ;
Info operator+(Info a, Info b)  
{
    if (a.mn > b.mn)
        a.mn = b.mn, a.mnId = b.mnId;
    if (a.mx < b.mx)
        a.mx = b.mx, a.mxId = b.mxId;
    return a;
}


```]
 #pagebreak() 
== `ST.h`


 #sourcecode[```cpp

class SparseTable
{
    using func_type = function<int(const int &, const int &)>;
    
    vector<vector<int>> ST;
    int len;
    vector<int> preLog;
    func_type op;
    static int default_func(const int &t1, const int &t2) { return max(t1, t2); }

public:
    void build(const vector<int> &v, func_type _func = default_func)
    {
        op = _func;
        len = v.size();
        int l1 = ceil(log2(len)) + 1;
        ST.assign(len, vector<int>(l1, 0));
        for (int i = 0; i < len; ++i)
        {
            ST[i][0] = v[i];
        }
        for (int j = 1; j < l1; ++j)
        {
            int pj = (1 << (j - 1));
            for (int i = 0; i + pj < len; ++i)
            {
                ST[i][j] = op(ST[i][j - 1], ST[i + (1 << (j - 1))][j - 1]);
            }
        }
        preLog.resize(len + 1);
        lop(i, 1, len + 1) preLog[i] = floor(log2(i));
    }

    int getsum(int l, int r)
    {
        if (r < l)
            return 0;
        l = max(0, l), r = min(len - 1, r);
        if (r == 0)
            return 0;
        int lt = r - l + 1;
        int q = preLog[lt];
        return op(ST[l][q], ST[r - (1 << q) + 1][q]);
    }
};

```]
 #pagebreak() 
== `twoDimPrfxSum.h`


 #sourcecode[```cpp

class prfx_2{
public:
    vector<vector<int>> mtx;
    int n,m;
    public:
    prfx_2(vector<vector<int>> _mtx){init(_mtx);};
    prfx_2() { };
    void init(vector<vector<int>> _mtx)
    {
        n = _mtx.size();
        m = _mtx[0].size();
        mtx.resize(n+1);
        for(auto &x:mtx) x.resize(m+1);
        lop(i,1,n+1)
            lop(j,1,m+1)
                mtx[i][j] = mtx[i-1][j] + mtx[i][j-1] - mtx[i-1][j-1] + _mtx[i-1][j-1];
    }

    int getsum(int x1,int y1,int x2,int y2)
    {
        x1 ++,x2 ++,y1 ++,y2 ++;
        return mtx[x2][y2] - mtx[x1-1][y2] - mtx[x2][y1-1] + mtx[x1-1][y1-1];
    }

    int getsum(pii d1,pii d2)
    {
        auto [x1,y1] = d1;
        auto [x2,y2] = d2;
        x1 ++,x2 ++,y1 ++,y2 ++;
        return mtx[x2][y2] - mtx[x1-1][y2] - mtx[x2][y1-1] + mtx[x1-1][y1-1];
    }

    vector<int> & operator[](std::size_t i)
    {
        return mtx[i];
    }

};

class conj_diff_2{
    vector<vector<int>> mtx;
    vector<vector<int>> prf;
    int n,m;
    public:

    conj_diff_2(int _n,int _m)
    {
        n = _n+1,m = _m+1;
        vector<vector<int>> initmp(n,vector<int> (m,0));
        init(initmp);
    }

    conj_diff_2(vector<vector<int>> _mtx)
    {
        this->init(_mtx);
    }

    conj_diff_2(){ }

    void init(vector<vector<int>> _mtx)
    {
        n = _mtx.size();
        m = _mtx[0].size();
        mtx.resize(n+2);
        for(auto &x:mtx) x.resize(m+2);
        prf.resize(n+2);
        for(auto &x:prf) x.resize(m+2);
        lop(i,1,n+1)
            lop(j,1,m+1)
                prf[i][j] = _mtx[i-1][j-1];
    }

    void modify(int x1,int y1,int x2,int y2,int k)
    {
        x1 ++,x2 ++,y1 ++,y2 ++;
        mtx[x1][y1] += k;
        mtx[x2+1][y1] -= k,mtx[x1][y2+1] -= k;
        mtx[x2+1][y2+1] += k;
    }

    void modify(pii d1,pii d2,int k)
    {
        this->modify(d1.fi,d1.se,d2.fi,d2.se,k);
    }

    vector<vector<int>> cacu()
    {
        auto res = prfx_2(mtx);
        vector<vector<int>> rst(n,vector<int>(m));
        lop(i,1,n+1)
            lop(j,1,m+1)
                rst[i-1][j-1] = prf[i][j] +  res[i+1][j+1];
        return rst;
    }

    vector<int> & operator[](std::size_t i)
    {
        return mtx[i];
    }
};


```]
 #pagebreak() 
= #smallcaps[Geo]

== `3_Dim.h`


 #sourcecode[```cpp
using ld = long double;

struct Point3
{
    ld x, y, z;
    Point3(ld x_ = 0, ld y_ = 0, ld z_ = 0) : x(x_), y(y_), z(z_) {}
    Point3 &operator+=(Point3 p) &
    {
        return x += p.x, y += p.y, z += p.z, *this;
    }
    Point3 &operator-=(Point3 p) &
    {
        return x -= p.x, y -= p.y, z -= p.z, *this;
    }
    Point3 &operator*=(Point3 p) &
    {
        return x *= p.x, y *= p.y, z *= p.z, *this;
    }
    Point3 &operator*=(ld t) &
    {
        return x *= t, y *= t, z *= t, *this;
    }
    Point3 &operator/=(ld t) &
    {
        return x /= t, y /= t, z /= t, *this;
    }
    friend Point3 operator+(Point3 a, Point3 b) { return a += b; }
    friend Point3 operator-(Point3 a, Point3 b) { return a -= b; }
    friend Point3 operator*(Point3 a, Point3 b) { return a *= b; }
    friend Point3 operator*(Point3 a, ld b) { return a *= b; }
    friend Point3 operator*(ld a, Point3 b) { return b *= a; }
    friend Point3 operator/(Point3 a, ld b) { return a /= b; }
    friend auto &operator>>(istream &is, Point3 &p)
    {
        return is >> p.x >> p.y >> p.z;
    }
    friend auto &operator<<(ostream &os, Point3 p)
    {
        return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
    }
};
using P3 = Point3;
struct Line3
{
    Point3 a, b;
};
using L3 = Line3;
struct Plane
{
    Point3 u, v, w;
};

ld len(P3 p)
{ // 原点到当前点的距离计算
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}
P3 crossEx(P3 a, P3 b)
{ // 叉乘
    P3 ans;
    ans.x = a.y * b.z - a.z * b.y;
    ans.y = a.z * b.x - a.x * b.z;
    ans.z = a.x * b.y - a.y * b.x;
    return ans;
}
ld cross(P3 a, P3 b)
{
    return len(crossEx(a, b));
}
ld dot(P3 a, P3 b)
{ // 点乘
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
P3 getVec(Plane s)
{ // 获取平面法向量
    return crossEx(s.u - s.v, s.v - s.w);
}
ld dis(P3 a, P3 b)
{ // 三维欧几里得距离公式
    ld val = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
    return sqrt(val);
}
P3 standardize(P3 vec)
{ // 将三维向量转换为单位向量
    return vec / len(vec);
}

bool onLine(P3 p1, P3 p2, P3 p3)
{ // 三点是否共线
    return sign(cross(p1 - p2, p3 - p2)) == 0;
}
bool onLine(Plane s)
{
    return onLine(s.u, s.v, s.w);
}
bool onPlane(P3 p1, P3 p2, P3 p3, P3 p4)
{ // 四点是否共面
    ld val = dot(getVec({p1, p2, p3}), p4 - p1);
    return sign(val) == 0;
}
bool pointOnSegment(P3 p, L3 l)
{
    return sign(cross(p - l.a, p - l.b)) == 0 && min(l.a.x, l.b.x) <= p.x &&
           p.x <= max(l.a.x, l.b.x) && min(l.a.y, l.b.y) <= p.y && p.y <= max(l.a.y, l.b.y) &&
           min(l.a.z, l.b.z) <= p.z && p.z <= max(l.a.z, l.b.z);
}
bool pointOnSegmentEx(P3 p, L3 l)
{ // pointOnSegment去除端点版
    return sign(cross(p - l.a, p - l.b)) == 0 && min(l.a.x, l.b.x) < p.x &&
           p.x < max(l.a.x, l.b.x) && min(l.a.y, l.b.y) < p.y && p.y < max(l.a.y, l.b.y) &&
           min(l.a.z, l.b.z) < p.z && p.z < max(l.a.z, l.b.z);
}
bool pointOnSegmentSide(P3 p1, P3 p2, L3 l)
{
    if (!onPlane(p1, p2, l.a, l.b))
    { // 特判不共面
        return 0;
    }
    ld val = dot(crossEx(l.a - l.b, p1 - l.b), crossEx(l.a - l.b, p2 - l.b));
    return sign(val) == 1;
}
bool pointOnPlaneSide(P3 p1, P3 p2, Plane s)
{
    ld val = dot(getVec(s), p1 - s.u) * dot(getVec(s), p2 - s.u);
    return sign(val) == 1;
}
bool lineParallel(L3 l1, L3 l2)
{ // 平行
    return sign(cross(l1.a - l1.b, l2.a - l2.b)) == 0;
}
bool lineVertical(L3 l1, L3 l2)
{ // 垂直
    return sign(dot(l1.a - l1.b, l2.a - l2.b)) == 0;
}

```]
 #pagebreak() 
== `Cood_Convert.h`


 #sourcecode[```cpp
using p = pair<int,int>;
```]
 #pagebreak() 
== `Distance.typ`


 
=== 曼哈顿距离
$ d(A,B) = |x_1 - x_2| + |y_1 - y_2| $

=== 欧几里得距离
$ d(A,B) = sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2) $

=== 切比雪夫距离
$ d(A,B) = max(|x_1 - x_2|, |y_1 - y_2|) $

=== 闵可夫斯基距离
$ d(A,B) = (|x_1 - x_2|^p + |y_1 - y_2|^p)^{1/p} $

=== 曼哈顿转切比雪夫

对于直角坐标中的$A(x_1,y_1),B(x_2,y_2)$

其曼哈顿距离
$ d(A,B) = max(|(x_1+y_1) - (x_2+y_2)|,|(x_1-y_1)-(x_2-y_2|)) $
即为点$A'(x_1+y_1,x_1-y_1),B'(x_2+y_2,x_2-y_2)$的切比雪夫距离。

同理，其切比雪夫距离
$ d(A,B) = max(|(x_1+y_1)/2-(x_2+y_2)/2| + |(x_1-y_1)/2-(x_2-y_2)/2|) $
即为点$A'((x_1+y_1)/2,(x_1-y_1)/2),B'((x_2+y_2)/2, (x_2-y_2)/2)$的曼哈顿距离。

综上：

$
"曼哈顿距离" & =>"切比雪夫距离：" \

(x,y) & => (x+y,x-y) \

"切比雪夫距离"&=>"曼哈顿距离："\

(x,y) &=> ((x+y)/2,(x-y)/2) $
 #pagebreak() 
== `Rotating_Calipers.h`


 #sourcecode[```cpp

//Rotating_Calipers
template<typename VALUE_TYPE>
class Rotating_Calipers
{
public:
    using pv = pair<VALUE_TYPE, VALUE_TYPE>;
    using vec_pv = vector<pair<VALUE_TYPE, VALUE_TYPE>>;
    vec_pv p;

    static VALUE_TYPE cross(pv p1, pv p2, pv p0)
    {
        pv t1 = {p1.fi - p0.fi, p1.se - p0.se},
           t2 = {p2.fi - p0.fi, p2.se - p0.se};
        return t1.fi * t2.se - t1.se * t2.fi;
    }

    static VALUE_TYPE dis(const pv &p1,const pv &p2){
        return (p1.fi - p2.fi) * (p1.fi - p2.fi) + (p1.se - p2.se) * (p1.se - p2.se);
    };

public:
    
    Rotating_Calipers() {}

    Rotating_Calipers(vec_pv _A) {
        build(_A);
    }

    void build(const vec_pv & _A) {
        p = ConvexHull(_A);
    }

    static vec_pv ConvexHull(vec_pv A, VALUE_TYPE flag = 1)
    {
        int n = A.size();
        if (n <= 2) return A; 
        vec_pv ans(n * 2);
        sort(A.begin(), A.end(),
        [](pv a,pv b) -> bool {
            if(fabs(a.fi - b.fi) < 1e-10)
                return a.se < b.se;
            else return a.fi < b.fi;}    );
        int now = -1;
        for (int i = 0; i < n; i++)
        { // 维护下凸包
            while (now > 0 && cross(A[i], ans[now], ans[now - 1]) < flag) now--;
            ans[++now] = A[i];
        }
        int pre = now;
        for (int i = n - 2; i >= 0; i--)
        { // 维护上凸包
            while (now > pre && cross(A[i], ans[now], ans[now - 1]) < flag) now--;
            ans[++now] = A[i];
        }
        ans.resize(now);
        return ans;
    }

    VALUE_TYPE getDiameter() {
        int j = 1;
        VALUE_TYPE ans = 0;
        int m = p.size();
        p.push_back(p[0]);
        for(int i = 0;i < m;i ++)
        {
            while( cross(p[i+1],p[j],p[i]) > cross(p[i+1],p[j+1],p[i]) ) j = (j+1)%m;
            ans = max(ans, max( dis(p[i],p[j]) , dis(p[i+1],p[j]) ) );
        }
        p.pop_back();
        return ans;
    }

    VALUE_TYPE getPerimeter() {
        VALUE_TYPE sum = 0;
        p.pb(p[0]);
        for(int i = 0;i < (int)p.size() - 1;i ++)
        {
            sum += sqrtl(dis(p[i],p[i+1]));
        }
        p.pop_back();
        return sum;
    }

};


```]
 #pagebreak() 
== `segIntersection.h`


 #sourcecode[```cpp
using pii = pair<int, int>

#define fi first
#define se second
    const long double EPS = 1e-9;

template <class T>
int sign(T x)
{
    if (-EPS < x && x < EPS)
        return 0;
    return x < 0 ? -1 : 1;
}

// 叉乘
template <class T>
T cross(pair<T, T> a, pair<T, T> b)
{
    return a.fi * b.se - a.se * b.fi;
}

// 二维快速跨立实验
template <class T>
bool segIntersection(pair<T, T> l1, pair<T, T> l2)
{
    auto [s1, e1] = l1;
    auto [s2, e2] = l2;
    auto A = max(s1.fi, e1.fi), AA = min(s1.fi, e1.fi);
    auto B = max(s1.se, e1.se), BB = min(s1.se, e1.se);
    auto C = max(s2.fi, e2.fi), CC = min(s2.fi, e2.fi);
    auto D = max(s2.se, e2.se), DD = min(s2.se, e2.se);
    return A >= CC && B >= DD && C >= AA && D >= BB &&
           sign(cross(s1, s2, e1) * cross(s1, e1, e2)) == 1 &&
           sign(cross(s2, s1, e2) * cross(s2, e2, e1)) == 1;
}

//三维线段交点，需要P3封装，不相交返回{0,{}}
pair<bool, P3> lineIntersection(L3 l1, L3 l2)
{
    if (!onPlane(l1.a, l1.b, l2.a, l2.b) || lineParallel(l1, l2))
    {
        return {0, {}};
    }
    auto [s1, e1] = l1;
    auto [s2, e2] = l2;
    ld val = 0;
    if (!onPlane(l1.a, l1.b, {0, 0, 0}, {0, 0, 1}))
    {
        val = ((s1.x - s2.x) * (s2.y - e2.y) - (s1.y - s2.y) * (s2.x - e2.x)) /
              ((s1.x - e1.x) * (s2.y - e2.y) - (s1.y - e1.y) * (s2.x - e2.x));
    }
    else if (!onPlane(l1.a, l1.b, {0, 0, 0}, {0, 1, 0}))
    {
        val = ((s1.x - s2.x) * (s2.z - e2.z) - (s1.z - s2.z) * (s2.x - e2.x)) /
              ((s1.x - e1.x) * (s2.z - e2.z) - (s1.z - e1.z) * (s2.x - e2.x));
    }
    else
    {
        val = ((s1.y - s2.y) * (s2.z - e2.z) - (s1.z - s2.z) * (s2.y - e2.y)) /
              ((s1.y - e1.y) * (s2.z - e2.z) - (s1.z - e1.z) * (s2.y - e2.y));
    }
    return {1, s1 + (e1 - s1) * val};
}

```]
 #pagebreak() 
= #smallcaps[Graph]

== #smallcaps[Flow]

=== `max_Flow_print.h`


 #sourcecode[```cpp

class maxFlow//AC
{
private:
    class edge
    {
    public:
        ll int nxt, cap, flow;                  
        pair<int, int> revNodeIdx; 
    public:
        edge(int _nxt, int _cap)
        {
            nxt = _nxt,cap = _cap,flow = 0;
        }
        void setRevIdx(int _i, int _j) {  revNodeIdx = {_i,_j}; }
    };
    vector<vector<edge>> edgeList; 
    vector<int> dep,fir;
    ll int maxFlowAns;
    int T, S;
public:
    maxFlow(int _n)
    {
        maxFlowAns = 0;
        S = 1;
        T = _n;
        edgeList.resize(_n + 1);
        dep.resize(_n + 1);
        fir.resize(_n+1);
    }
    void resetTS(int _T, int _S) { T = _T,S = _S; }

    void addedge(int _u, int _v, int _w)
    {
        edgeList[_u].push_back(edge(_v, _w));
        edgeList[_v].push_back(edge(_u, 0)); 
        edgeList[_u][edgeList[_u].size() - 1].setRevIdx(_v, edgeList[_v].size() - 1);
        edgeList[_v][edgeList[_v].size() - 1].setRevIdx(_u, edgeList[_u].size() - 1);
    }

    bool bfs()
    {
        queue<int> que;
        for (auto x = dep.begin(); x != dep.end(); x++) (*x) = 0; 
        dep[S] = 1;
        que.push(S);
        while (que.size())
        {
            ll int at = que.front();
            que.pop();
            for (int i = 0; i < edgeList[at].size(); i++)
            {
                auto tar = edgeList[at][i];
                if ((!dep[tar.nxt]) && (tar.flow < tar.cap))
                {
                    dep[tar.nxt] = dep[at] + 1;
                    que.push(tar.nxt);
                }
            }
        }
        return dep[T];
    }

    ll int dfs(int at, ll int flow)
    {
        if ((at == T) || (!flow))
            return flow; // 到了或者没了
        ll int ret = 0;  // 本节点最大流
        for (int &i = fir[at]; i < edgeList[at].size(); i++)
        {
            auto tar = edgeList[at][i];   // 目前遍历的边
            int tlFlow = 0;   // 目前边的最大流
            if (dep[at] == dep[tar.nxt] - 1) // 遍历到的边为合法边
            {
                tlFlow = dfs(tar.nxt, min((ll)tar.cap - tar.flow, flow - ret));
                if (!tlFlow)
                    continue;   // 若最大流为0，什么都不做
                ret += tlFlow;   // 本节点最大流累加
                edgeList[at][i].flow += tlFlow;  // 本节点实时流量
                edgeList[tar.revNodeIdx.first][tar.revNodeIdx.second].flow -= tlFlow; // 反向边流量
                if (ret == flow)
                    return ret; // 充满了就不继续扫了
            }
        }
        return ret;
    }

    ll int dinic()
    {
        if (maxFlowAns)
            return maxFlowAns;
        while (bfs())
        {
            for(auto x = fir.begin();x != fir.end();x ++) (*x) = 0;
            maxFlowAns += dfs(S, INT_MAX);
        }
        return maxFlowAns;
    }
};

```]
 #pagebreak() 
=== `min_Cost.h`


 #sourcecode[```cpp

const int INF = 0x3f3f3f3f

class PD//AC
{
public:
    class edge
    {
    public:
        int v, f, c, next;
        edge(int _v,int _f,int _c,int _next)
        {
            v = _v,f = _f,c = _c,next = _next;
        }
        edge() { }
    } ;

    void vecset(int value,vector<int> &arr)
    {
        for(int i = 0;i < arr.size();i ++) arr[i] = value;
        return;
    }

    class node
    {
    public:
        int v, e;
    } ;

    class mypair
    {
    public:
        int dis, id;

        bool operator<(const mypair &a) const { return dis > a.dis; }

        mypair(int d, int x)
        {
            dis = d;
            id = x;
        }
    };

    vector<int> head,dis,vis,h;
    vector<edge> e;
    vector<node> p;
    int n, m, s, t, cnt = 1, maxf, minc;

    PD(int _n,int _m,int _s,int _t)
    {
        n = _n, m = _m,s = _s,t = _t;
        maxf = 0, minc = 0;
        head.resize(n+2), dis.resize(n+2),vis.resize(n+2);
        e.resize(2);
        h.resize(n+2), p.resize(m+2);
    }

    void addedge(int u, int v, int f, int c)
    {
        e.push_back(edge(v,f,c,head[u]));
        head[u] = e.size()-1;
        e.push_back(edge(u,0,-c,head[v]));
        head[v] = e.size()-1;
    }

    bool dijkstra()
    {
        priority_queue<mypair> q;
        vecset(INF,dis);
        vecset(0,vis);
        dis[s] = 0;
        q.push(mypair(0, s));
        while (!q.empty())
        {
            int u = q.top().id;
            q.pop();
            if (vis[u])
                continue;
            vis[u] = 1;
            for (int i = head[u]; i; i = e[i].next)
            {
                int v = e[i].v, nc = e[i].c + h[u] - h[v];
                if (e[i].f && dis[v] > dis[u] + nc)
                {
                    dis[v] = dis[u] + nc;
                    p[v].v = u;
                    p[v].e = i;
                    if (!vis[v])
                        q.push(mypair(dis[v], v));
                }
            }
        }
        return dis[t] != INF;
    }

    void spfa()
    {
        queue<int> q;
        vecset(63,h);
        h[s] = 0, vis[s] = 1;
        q.push(s);
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            vis[u] = 0;
            for (int i = head[u]; i; i = e[i].next)
            {
                int v = e[i].v;
                if (e[i].f && h[v] > h[u] + e[i].c)
                {
                    h[v] = h[u] + e[i].c;
                    if (!vis[v])
                    {
                        vis[v] = 1;
                        q.push(v);
                    }
                }
            }
        }
    }

    int pd()
    {
        spfa();
        while (dijkstra())
        {
            int minf = INF;
            for (int i = 1; i <= n; i++)
                h[i] += dis[i];
            for (int i = t; i != s; i = p[i].v)
                minf = min(minf, e[p[i].e].f);
            for (int i = t; i != s; i = p[i].v)
            {
                e[p[i].e].f -= minf;
                e[p[i].e ^ 1].f += minf;
            }
            maxf += minf;
            minc += minf * h[t];
        }
        return 0;
    }

    pair<int,int> get()
    {
        return {maxf,minc};
    }
};

```]
 #pagebreak() 
== #smallcaps[Path]

=== `SCSP.h`


 #sourcecode[```cpp
// Sides Change Shortest Path
#include<bits/stdc++.h>
using namespace std;
#define ll long long
const ll inf=1e18;
#define N 200505
#define mod 1000000007
//换边最短路
int n,m,qu;
struct Original_edge{int u,v;ll len;}OE[N];
struct tabl{int to,id;ll len;};
struct node{int id; ll len; bool operator < (const node &x) const {return x.len<len;}};
vector<tabl> edge[N<<1];
priority_queue<node> q;
int lstvis[N]={0},l[N]={0},r[N]={0},ind[N]={0},vis[N];
ll t[N<<4]={0};
ll disT[N],disN[N];
bool on_path[N];
int path_cnt=0;

void dj(int p,ll dis[],int f)
{
    for(int i=1;i<=n;i++)
    {
        vis[i]=0;
        dis[i]=inf;
    }
    dis[p]=0;
    q.push((node){p,0});
    while(!q.empty())
    {
        node temp=q.top();
        q.pop();
        int u=temp.id;
        ll w=temp.len;
        if(vis[u])continue;
        vis[u]++;
        dis[u]=w;
        for(auto re:edge[u])
        {
            int v=re.to;
            int id=re.id;
            ll tw=re.len;

            if(dis[v]<=w+tw)continue;
            
            dis[v]=w+tw;

            lstvis[v]=id;  //记前驱

            q.push((node){v,w+tw});

            if(f==1&&!on_path[v])l[v]=l[u];  //记住前后缀
            if(f==2&&!on_path[v])r[v]=r[u];
        }
    }
}
void trace()
{
    int u=1;
    on_path[u]=true;
    l[u]=r[u]=0;

    for(int i=1;u!=n;i++)
    {
        int e_id=lstvis[u];  //取前驱
        ind[e_id]=i;      //做标记

        u^=OE[e_id].u^OE[e_id].v;  //取前点
        on_path[u]=true;       //做标记
        l[u]=r[u]=i;      //做标记
        path_cnt++;      //路长
    }
}
void build(int le,int ri,int p)
{
    t[p]=inf;
    if(le==ri)return;
    int mid=(le+ri)>>1,lef=(p<<1),rig=lef|1;
    build(le,mid,lef);
    build(mid+1,ri,rig);
}
void update(int L,int R,int x,int y,int p,ll k)
{
    if(x>y)return;
    if(x<=L&&R<=y)
    {
        t[p]=min(t[p],k);
        return;
    }
    int mid=(L+R)>>1,lef=(p<<1),rig=lef|1;
    if(x<=mid)update(L,mid,x,y,lef,k);
    if(y>mid)update(mid+1,R,x,y,rig,k);
}
ll query(int L,int R,int x,int y,int p)
{
    ll ans=t[p];
    if(L==R)
    {
        return ans;
    }
    int mid=(L+R)>>1,lef=(p<<1),rig=lef|1;
    if(y<=mid)ans=min(ans,query(L,mid,x,y,lef));
    else ans=min(ans,query(mid+1,R,x,y,rig));
    return ans;
}
void solve()
{
    memset(on_path,false,sizeof(on_path));
    cin>>n>>m>>qu;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        ll w;
        cin>>u>>v>>w;

        OE[i].u=u;
        OE[i].v=v;
        OE[i].len=w;

        edge[u].push_back({v,i,w});
        edge[v].push_back({u,i,w});
    }   

    dj(n,disN,0);
    trace();
    dj(1,disT,1);
    dj(n,disN,2);

    build(1,path_cnt,1);

    for(int i=1;i<=m;i++)
    {
        int u=OE[i].u,v=OE[i].v;
        ll w=OE[i].len;
        if(ind[i])continue;
        update(1,path_cnt,l[u]+1,r[v],1,disT[u]+w+disN[v]);
        update(1,path_cnt,l[v]+1,r[u],1,disT[v]+w+disN[u]);
    }

    while(qu--)
    {
        int ti;
        ll ch;
        ll ans;
        cin>>ti>>ch;
        if(ind[ti])
        {
            ans=disT[n]-OE[ti].len+ch;
            if(ch>OE[ti].len)
            {
                ans=min(ans,query(1,path_cnt,ind[ti],ind[ti],1));
            }
        }
        else 
        {
            ans=disT[n];
            if(OE[ti].len>ch)
            {
                int u=OE[ti].u,v=OE[ti].v;
                ans=min(ans,min(disT[u]+ch+disN[v],disT[v]+ch+disN[u]));
            }
        }
        cout<<ans<<"\n";
    }
}
```]
 #pagebreak() 
=== `SPFA+SLE.h`


 #sourcecode[```cpp
#include<bits/stdc++.h>
using namespace std;
#define ll long long 
const int inf =1e9+7;
////////////////////
struct p3008  //链式前向星
{
    int to;
    int nex;
    int v;
}a[300000];
int head[25006]={0},cnt=0;

int dis[25006]={0};  //距离

void add(int u,int v,int w)
{
    a[++cnt].nex=head[u];
    a[cnt].to=v;
    a[cnt].v=w;
    head[u]=cnt;
}
/////////////////////////////////
//SPFA，但是双端队列优化
//小的放队头，大的放队尾
void spfa(int t,int s)
{
    for(int i=1;i<=t;i++)dis[i]=inf;  //初始化距离
    deque<int> q;
    vector<int> vis(t+1),cou(t+1);

    q.push_front(s);  //起点入队

    dis[s]=0;  //起点初始化
    vis[s]++;

    while(!q.empty())
    {
        int u=q.front();  //取出队头
        q.pop_front();

        cou[u]++;

        if(cou[u]>t)return;  //遍历次数过大，出现负环

        for(int i=head[u];i;i=a[i].nex)  //遍历儿子
        {
            int v=a[i].to;
            int w=a[i].v;
            if(dis[v]>dis[u]+w)  //可松弛
            {
                dis[v]=dis[u]+w;  //松弛操作

                if(!vis[v])
                {
                    if(q.empty())q.push_back(v);       //特判！队空则随便入队
                    else if(dis[q.front()]>=dis[v])q.push_front(v);  //小的放队头
                    else q.push_back(v);  //大的放队尾
                    vis[v]++;
                }
            }
        }
        vis[u]--;  //取消入队标记
    }
}
void solve()
{
    int t,r,p,s;
    cin>>t>>r>>p>>s;
    for(int i=1;i<=r;i++)
    {
        int u,v;
        int w;
        cin>>u>>v>>w;
        add(u,v,w);
        add(v,u,w);
    }
    for(int i=1;i<=p;i++)
    {
        int u,v;
        int w;
        cin>>u>>v>>w;
        add(u,v,w);
    }

    spfa(t,s);  //建好边调用即可，t是点数，s是起点
    
    for(int i=1;i<=t;i++)
    {
        if(dis[i]>=inf)cout<<"NO PATH\n";
        else cout<<dis[i]<<"\n";
    }
}
```]
 #pagebreak() 
=== `SPFA.h`


 #sourcecode[```cpp
#include<bits/stdc++.h>
using namespace std;
#define ll long long
/////////////////////////////////////
struct p4779SPFA  //链式前向星
{
    int to;
    int nex;
    ll w;
}star[550000];
int head[100005]={0},cnt=0;
void add(int u,int v,ll w)
{
    star[++cnt].to=v;
    star[cnt].nex=head[u];
    star[cnt].w=w;
    head[u]=cnt;
}
///////////////////////////////

//SPFA,一款暴力的最短路算法
//常用于网络流、判负环、带负边的最短路
//暴力版的dj，每次松弛都让未入队的点入队，直到无法松弛为止
int n,m,s;
ll dis[100005]={0};   //距离
int times[100005]={0}; //遍历次数
bool vis[100005];  //入队情况
void SPFA()
{
    memset(vis, 0 ,sizeof(vis)); //清空vis
    for(int i=1;i<=n;i++)dis[i]=2147483647;  //距离初始化为无穷大
    dis[s]=0;  //起点归零
    queue<int> q;  //队列
    q.push(s);
    while(!q.empty())
    {
        int temp=q.front();
        q.pop();
        times[temp]++;
        vis[temp]=false;  //出队判定

        if(times[temp]>2+n)  //判负环  ，遍历次数过多说明出现了负环
        {
            times[s]=n+10;  //做标记
            break;
        }

        for(int i=head[temp];i;i=star[i].nex)  //遍历儿子
        {
            int v=star[i].to;
            if(dis[v]>dis[temp]+star[i].w)  //可松弛
            {
                dis[v]=dis[temp]+star[i].w;//松弛
                if(!vis[v])  //未入队则入队
                {
                    q.push(v);
                    vis[v]=true;
                }
            }
        }
    }

}
void solve()
{
    //spfa用于求最短路
    cin>>n>>m>>s;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        ll w;
        cin>>u>>v>>w;
        add(u,v,w);
    }

    SPFA();//建边后调用即可

    if(times[s]>n)  //找负环
    {
        cout<<"minus circle!\n";
        return;
    }
    for(int i=1;i<=n;i++)cout<<dis[i]<<" ";
    cout<<"\n";
}
```]
 #pagebreak() 
== #smallcaps[Tree]

=== `hld.h`


 #sourcecode[```cpp
void dfs1(int o) {
  son[o] = -1;
  siz[o] = 1;
  for (int j = h[o]; j; j = nxt[j])
    if (!dep[p[j]]) {
      dep[p[j]] = dep[o] + 1;
      fa[p[j]] = o;
      dfs1(p[j]);
      siz[o] += siz[p[j]];
      if (son[o] == -1 || siz[p[j]] > siz[son[o]]) son[o] = p[j];
    }
}

void dfs2(int o, int t) {
  top[o] = t;
  cnt++;
  dfn[o] = cnt;
  rnk[cnt] = o;
  if (son[o] == -1) return;
  dfs2(son[o], t);  // 优先对重儿子进行 DFS，可以保证同一条重链上的点 DFS 序连续
  for (int j = h[o]; j; j = nxt[j])
    if (p[j] != son[o] && p[j] != fa[o]) dfs2(p[j], p[j]);
}

int lca(int u, int v) {
  while (top[u] != top[v]) {
    if (dep[top[u]] > dep[top[v]])
      u = fa[top[u]];
    else
      v = fa[top[v]];
  }
  return dep[u] > dep[v] ? v : u;
} 

// st 是线段树结构体
int querymax(int x, int y) {
  int ret = -inf, fx = top[x], fy = top[y];
  while (fx != fy) {
    if (dep[fx] >= dep[fy])
      ret = max(ret, st.query1(1, 1, n, dfn[fx], dfn[x])), x = fa[fx];
    else
      ret = max(ret, st.query1(1, 1, n, dfn[fy], dfn[y])), y = fa[fy];
    fx = top[x];
    fy = top[y];
  }
  if (dfn[x] < dfn[y])
    ret = max(ret, st.query1(1, 1, n, dfn[x], dfn[y]));
  else
    ret = max(ret, st.query1(1, 1, n, dfn[y], dfn[x]));
  return ret;
}
```]
 #pagebreak() 
=== `lca.h`


 #sourcecode[```cpp

class LCA{ 
public:
    vector<vector<pii>> cnj;
    vector<int> lg,dep;
    vector<array<int,32>> fa,wei;
    int n;

    LCA(int _n) {
        n = _n;
        cnj.resize(n+1);
        lg.resize(n+1),fa.resize(n+1),dep.resize(n+1),wei.resize(n+1);
        for(int i = 1; i <= n; i ++)
            lg[i] = lg[i-1] + (1 << lg[i-1] == i);
    }

    void addEdge(int u,int v,int w) {
        cnj[u].push_back({v,w});
        cnj[v].push_back({u,w});
    }

    void build(int rt = 1) {
        using itf = function<void(int,int)>;
        itf dfs = [&](int p,int f) -> void {
            fa[p][0] = f,dep[p] = dep[f] + 1;
            // wei[p][0] = 0;
            for(int i = 1;i <= lg[dep[p]];i ++) fa[p][i] = fa[fa[p][i-1]][i-1];
            for(int i = 1;i <= lg[dep[p]];i ++) wei[p][i] = max(wei[p][i-1],wei[fa[p][i-1]][i-1]);
            for(auto [x,w]:cnj[p]) if(x == f) continue;
            else wei[x][0] = w,dfs(x,p);
        };
        dfs(rt,0);
    }

    int get(int x,int y) {
        if(dep[x] < dep[y]) swap(x,y);
        while(dep[x] > dep[y]) x = fa[x][lg[dep[x] - dep[y]] - 1];
        if(x == y) return x;
        for(int k = lg[dep[x]]-1;k >= 0;k --) if(fa[x][k] != fa[y][k]) x = fa[x][k],y = fa[y][k];
        return fa[x][0];
    }

    int getmaxw(int x,int y) {
        int curmx = 0;
        if(dep[x] < dep[y]) swap(x,y);
        while(dep[x] > dep[y]) curmx = max(curmx,wei[x][lg[dep[x] - dep[y]] - 1]), x = fa[x][lg[dep[x] - dep[y]] - 1];
        if(x == y) return curmx;
        for(int k = lg[dep[x]]-1;k >= 0;k --) 
            if(fa[x][k] != fa[y][k]) 
                curmx = max(curmx,wei[x][k]),x = fa[x][k],
                curmx = max(curmx,wei[y][k]),y = fa[y][k];
        return max({curmx,wei[x][0],wei[y][0]});
    } 
};


```]
 #pagebreak() 
= #smallcaps[Math]

== #smallcaps[Number_theory]

=== `basic.h`


 #sourcecode[```cpp
__builtin_ffsll(x)
// 返回 x 的二进制末尾最后一个 1 的位置

__builtin_clzll(x)
// 返回 x 的二进制的前导 0 的个数。

__builtin_ctzll(x)
// 返回 x 的二进制末尾连续 0 的个数。

__builtin_clrsbll(x)
// 当 x 的符号位为 0 时返回 x 的二进制的前导 0 的个数减一，否则返回 x 的二进制的前导 1 的个数减一。

__builtin_popcountll(x)
// 返回 x 的二进制中 1 的个数。

__builtin_parity(x)
// 判断 x 的二进制中 1 的个数的奇偶性。

int binpow(int x, int y)
{
    int res = 1;
    while (y > 0)
    {
        if (y & 1)
            res = res * x % mod;
        x = x * x % mod;
        y >>= 1;
    }
    return res;
}

void exgcd(int a, int b, int& x, int& y) {
  if (b == 0) {
    x = 1, y = 0;
    return;
  }
  exgcd(b, a % b, y, x);
  y -= a / b * x;
}

binpow(x,mod-2)
```]
 #pagebreak() 
=== `Comb.h`


 #sourcecode[```cpp

const int N = 1e6;
const int mod = 1e9+7;

int binpow(int x, int y)
{
    int ans = 1;
    while (y)
    {
        if (y & 1) ans = ans * x % mod;
        x = x * x % mod;
        y >>= 1;
    }
    return ans;
}

vector<int> fac(N), inv(N);

void init()
{
    fac[0] = inv[0] = 1;
    for (int i = 1; i < N; i++) fac[i] = fac[i - 1] * i % mod;
    inv[N - 1] = binpow(fac[N - 1], mod - 2);
    for (int i = N - 2; i >= 1; i--)
    {
        inv[i] = inv[i + 1] * (i + 1) % mod;
    }
}

auto C = [&](int x, int y) -> int
{
    return (fac[x] * inv[y] % mod) * inv[x - y] % mod;
};
```]
 #pagebreak() 
=== `CRT.h`


 #sourcecode[```cpp

int CRT(vector<int> &r, vector<int> &a)
{ // % r === a
    int n = a.size();
    __int128 k = 1, ans = 0;
    for (int i = 0; i < n; i++) k *= r[i];
    for (int i = 0; i < n; i++)
    {
        __int128 m = k / r[i];
        int b, y;
        exgcd(m, r[i], b, y); // b * m mod r[i] = 1
        ans = (ans + a[i] * m * b % k) % k;
    }
    return (ans % k + k) % k;
}



int mul(int a, int b, int m) {
    return (__int128)a * b % m;
}

int exgcd (int a,int b,int &x,int &y) {
    if (b == 0) { x = 1, y = 0; return a; }
    int g = exgcd(b, a % b, x, y), tp = x;
    x = y, y = tp - a / b * y;
    return g;
};

int EXCRT(vector<int> &a,vector<int> &r) { // % r == a 
    int x, y, k;
    int n = r.size();
    int M = a[0], ans = r[0];
    for (int i = 1; i < n; ++ i) {
        int ca = M, cb = a[i], cc = (r[i] - ans % cb + cb) % cb;
        int gcd = exgcd(ca, cb, x, y), bg = cb / gcd;
        if (cc % gcd != 0) return -1;
        x = mul(x, cc / gcd, bg);
        ans += x * M;
        M *= bg;
        ans = (ans % M + M) % M;
    }
    return (ans % M + M) % M;
}
```]
 #pagebreak() 
=== `Eular_phi.h`


 #sourcecode[```cpp
int euler_phi(int n) {
  int ans = n;
  for (int i = 2; i * i <= n; i++)
    if (n % i == 0) {
      ans = ans / i * (i - 1);
      while (n % i == 0) n /= i;
    }
  if (n > 1) ans = ans / n * (n - 1);
  return ans;
}
```]
 #pagebreak() 
=== `Eular_sieve.h`


 #sourcecode[```cpp


vector<int> init(int n)
{
    vector<int> pri;
    vector<bool> vis(n, 0); 
    for (int i = 2; i <= n; i++)
    {
        if (!vis[i])
            pri.push_back(i);
        for (int j = 0; j < pri.size(); j++)
        {
            if (i * pri[j] > n)
                break;
            vis[pri[j] * i] = 1;
            if (i % pri[j] == 0)
                break;
        }
    }
    return pri;
}

```]
 #pagebreak() 
=== `factor_pr.h`


 #sourcecode[```cpp

#define int long long
#define pii pair<int, int>
const int INF = 1145141919810LL;
using namespace std;

class Pollard_Rho
{
private:

    vector<int> B;

    int mul(int a, int b, int m)
    {
        int r = a * b - m * (int)(1.L / m * a * b);
        return r - m * (r >= m) + m * (r < 0);
    }

    int mypow(int a, int b, int m)
    {
        int res = 1 % m;
        for (; b; b >>= 1, a = mul(a, a, m))
        {
            if (b & 1)
            {
                res = mul(res, a, m);
            }
        }
        return res;
    }

    bool MR(int n)
    {
        if (n <= 1)
            return 0;
        for (int p : B)
        {
            if (n == p)
                return 1;
            if (n % p == 0)
                return 0;
        }
        int m = (n - 1) >> __builtin_ctz(n - 1);
        for (int p : B)
        {
            int t = m, a = mypow(p, m, n);
            while (t != n - 1 && a != 1 && a != n - 1)
            {
                a = mul(a, a, n);
                t *= 2;
            }
            if (a != n - 1 && t % 2 == 0)
                return 0;
        }
        return 1;
    }

    inline const int getfecsum(int _n)
    {
        int sum = 0;
        while (_n)
        {
            sum += _n % 10;
            _n /= 10;
        }
        return sum;
    };

    int PR(int n)
    {
        for (int p : B)
        {
            if (n % p == 0)
                return p;
        }
        auto f = [&](int x) -> int
        {
            x = mul(x, x, n) + 1;
            return x >= n ? x - n : x;
        };
        int x = 0, y = 0, tot = 0, p = 1, q, g;
        for (int i = 0; (i & 255) || (g = gcd(p, n)) == 1; i++, x = f(x), y = f(f(y)))
        {
            if (x == y)
            {
                x = tot++;
                y = f(x);
            }
            q = mul(p, abs(x - y), n);
            if (q)
                p = q;
        }
        return g;
    }

    vector<int> fac(int n)
    {
        // if(n == 0)
        // #define pb emplace_back
        if (n <= 1)
            return {};
        if (MR(n))
            return {n};
        int d = PR(n);
        auto v1 = fac(d), v2 = fac(n / d);
        auto i1 = v1.begin(), i2 = v2.begin();
        vector<int> ans;
        while (i1 != v1.end() || i2 != v2.end())
        {
            if (i1 == v1.end())
            {
                ans.pb(*i2++);
            }
            else if (i2 == v2.end())
            {
                ans.pb(*i1++);
            }
            else
            {
                if (*i1 < *i2)
                {
                    ans.pb(*i1++);
                }
                else
                {
                    ans.pb(*i2++);
                }
            }
        }
        return ans;
    }

public:

    Pollard_Rho(){
        B = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    }

    vector<pii> fac_Comp(int n)
    {
        auto srt = fac(n);
        map<int, int> cnt;
        for (auto x : srt)
            cnt[x]++;
        vector<pii> rt;
        for (auto x : cnt)
            rt.push_back(x);
        return rt;
    }

    vector<int> fac_pri(int n)
    {
        return fac(n);
    }
    
    vector<int> fac_all(int n)
    {
        vector<pii> rt = fac_Comp(n);
        vector<int> v;
        function<void(int, int)> dfs = [&](int id, int x)
        {
            if (id == rt.size())
            {
                v.push_back(x);
                return;
            }
            for(int i = 0;i <= rt[id].se;i ++)
            {
                dfs(id + 1, x * (mypow(rt[id].fi, i, INF)));
            }
        };
        dfs(0, 1);
        return v;
    }
};
 
```]
 #pagebreak() 
== #smallcaps[Other]

=== `Frac.h`


 #sourcecode[```cpp
template<class T>
struct Frac {
    T num;
    T den;
    Frac(T num_, T den_) : num(num_), den(den_) {
        if (den < 0) {
            den = -den;
            num = -num;
        }
    }
    Frac() : Frac(0, 1) {}
    Frac(T num_) : Frac(num_, 1) {}
    explicit operator double() const {
        return 1. * num / den;
    }
    Frac &operator+=(const Frac &rhs) {
        num = num * rhs.den + rhs.num * den;
        den *= rhs.den;
        return *this;
    }
    Frac &operator-=(const Frac &rhs) {
        num = num * rhs.den - rhs.num * den;
        den *= rhs.den;
        return *this;
    }
    Frac &operator*=(const Frac &rhs) {
        num *= rhs.num;
        den *= rhs.den;
        return *this;
    }
    Frac &operator/=(const Frac &rhs) {
        num *= rhs.den;
        den *= rhs.num;
        if (den < 0) {
            num = -num;
            den = -den;
        }
        return *this;
    }
    friend Frac operator+(Frac lhs, const Frac &rhs) {
        return lhs += rhs;
    }
    friend Frac operator-(Frac lhs, const Frac &rhs) {
        return lhs -= rhs;
    }
    friend Frac operator*(Frac lhs, const Frac &rhs) {
        return lhs *= rhs;
    }
    friend Frac operator/(Frac lhs, const Frac &rhs) {
        return lhs /= rhs;
    }
    friend Frac operator-(const Frac &a) {
        return Frac(-a.num, a.den);
    }
    friend bool operator==(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den == rhs.num * lhs.den;
    }
    friend bool operator!=(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den != rhs.num * lhs.den;
    }
    friend bool operator<(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den < rhs.num * lhs.den;
    }
    friend bool operator>(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den > rhs.num * lhs.den;
    }
    friend bool operator<=(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den <= rhs.num * lhs.den;
    }
    friend bool operator>=(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den >= rhs.num * lhs.den;
    }
    friend std::ostream &operator<<(std::ostream &os, Frac x) {
        T g = std::gcd(x.num, x.den);
        if (x.den == g) {
            return os << x.num / g;
        } else {
            return os << x.num / g << "/" << x.den / g;
        }
    }
};
 
using F = Frac<int>;
```]
 #pagebreak() 
= #smallcaps[String]

== `AC_automaton.h`


 #sourcecode[```cpp
struct ACAutomaton
{
    static constexpr int N = 1e6 + 10;
    int ch[N][26], fail[N], cntNodes;
    int cnt[N];
    ACAutomaton()
    {
        cntNodes = 1;
    }
    void insert(string s)
    {
        int u = 1;
        for (auto c : s)
        {
            int &v = ch[u][c - 'a'];
            if (!v)
                v = ++cntNodes;
            u = v;
        }
        cnt[u]++;
    }
    void build()
    {
        fill(ch[0], ch[0] + 26, 1);
        queue<int> q;
        q.push(1);
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            for (int i = 0; i < 26; i++)
            {
                int &v = ch[u][i];
                if (!v)
                    v = ch[fail[u]][i];
                else
                {
                    fail[v] = ch[fail[u]][i];
                    q.push(v);
                }
            }
        }
    }
    LL query(string t)
    {
        LL ans = 0;
        int u = 1;
        for (auto c : t)
        {
            u = ch[u][c - 'a'];
            for (int v = u; v && ~cnt[v]; v = fail[v])
            {
                ans += cnt[v];
                cnt[v] = -1;
            }
        }
        return ans;
    }
};
```]
 #pagebreak() 
== `compress_print.h`


 #sourcecode[```cpp

const int N = 1 << 21;
static const int mod1 = 1E9 + 7, base1 = 127;
static const int mod2 = 1E9 + 9, base2 = 131;
vector<int> val1;
vector<int> val2;
void init(int n = N)
{
    val1.resize(n + 1), val2.resize(n + 2);
    val1[0] = 1, val2[0] = 1;
    for (int i = 1; i <= n; i++)
    {
        val1[i] = val1[i - 1] * base1 % mod1;
        val2[i] = val2[i - 1] * base2 % mod2;
    }
}

string compress(vector<string> in)
{ // 前后缀压缩
    vector<int> h1{1};
    vector<int> h2{1};
    string ans = "#";
    for (auto s : in)
    {
        s = "#" + s;
        int st = 0;
        int chk1 = 0;
        int chk2 = 0;
        for (int j = 1; j < s.size() && j < ans.size(); j++)
        {
            chk1 = (chk1 * base1 % mod1 + s[j]) % mod1;
            chk2 = (chk2 * base2 % mod2 + s[j]) % mod2;
            if ((h1.back() == (h1[ans.size() - 1 - j] * val1[j] % mod1+ chk1) % mod1) &&
                (h2.back() == (h2[ans.size() - 1 - j] * val2[j] % mod2+ chk2) % mod2)    ) 
                st = j;
        }
        for (int j = st + 1; j < s.size(); j++)
        {
            ans += s[j];
            h1.push_back((h1.back() * base1 % mod1 + s[j]) % mod1);
            h2.push_back((h2.back() * base2 % mod2 + s[j]) % mod2);
        }
    }
    return ans.substr(1);
}


```]
 #pagebreak() 
== `get_occr.h`


 #sourcecode[```cpp
#include <template_overAll.h>

/*
 * 找到某一堆短字符串在长字符串中的出现位置
 * dira=1为最早出现的后端点下标 dira=0为最晚出现的前端点下标
 * 源字符串s长度为|s|，查找字符串列表中所有字符串长度和为|_s|
 * 则时间复杂度为O(max(|_s|log(|_s|),|s|))
 */
class get_occr
{
private:
    string s;
public:
    get_occr(string _s) { s = _s; }
    vector<int> locate(vector<string> _s,bool dira = 1)
    {
        int n = _s.size();
        vector<int> occr(n,-1);
        map<char,vector<pair<int,int>>> gncing;
        if(dira == 1)
        {
            for(int i = 0;i < n;i++)
                gncing[_s[i][0]].push_back({i,0});
            for(int i = 0;i < s.size();i ++)
            {
                vector<pair<int,int>> gnctmp = gncing[s[i]];
                gncing[s[i]].clear();
                for(int j = 0;j < gnctmp.size();j ++)
                {
                    if(gnctmp[j].se+1 < _s[gnctmp[j].fi].size())
                            gncing[_s[gnctmp[j].fi][gnctmp[j].se+1]].push_back({gnctmp[j].fi,gnctmp[j].se+1});
                    else occr[gnctmp[j].fi] = i;
                }
            }
        } else {
            for(int i = 0;i < n;i++) gncing[_s[i][_s[i].size()-1]].push_back({i,_s[i].size()-1});
            for(int i= s.size()-1;i >=0;i --)
            {
                vector<pair<int,int>> gnctmp = gncing[s[i]];
                gncing[s[i]].clear();
                for(int j = 0;j < gnctmp.size();j ++)
                {
                    if(gnctmp[j].se -1 >= 0)
                            gncing[_s[gnctmp[j].fi][gnctmp[j].se-1]].push_back({gnctmp[j].fi,gnctmp[j].se-1});
                    else occr[gnctmp[j].fi] = i;
                }
            }
        }
        return occr;
    }
};

```]
 #pagebreak() 
== `hash_print.h`


 #sourcecode[```cpp
#define int long long
const int N = 1 << 21;
static const int mod1 = 1E9 + 7, base1 = 127;
static const int mod2 = 1E9 + 9, base2 = 131;
vector<int> val1;
vector<int> val2;
using puv = pair<int,int>;
void init(int n = N)
{
    val1.resize(n + 1), val2.resize(n + 2);
    val1[0] = 1, val2[0] = 1;
    for (int i = 1; i <= n; i++)
    {
        val1[i] = val1[i - 1] * base1 % mod1;
        val2[i] = val2[i - 1] * base2 % mod2;
    }
}
class hstring
{
public:
    vector<int> h1;
    vector<int> h2;
    string s;

    hstring(string s_) : s(s_), h1{0}, h2{0}
    {
        build();
    }

    void build()
    {
        for (auto it : s)
        {
            h1.push_back((h1.back() * base1 % mod1 + it) % mod1);
            h2.push_back((h2.back() * base2 % mod2 + it) % mod2);
        }
    }

    puv get()
    { // 输出整串的哈希值
        return {h1.back(), h2.back()};
    }

    puv substring(int l, int r)
    { // 输出子串的哈希值
        if (l > r) swap(l, r);
        int ans1 = (mod1 + h1[r + 1] - h1[l] * val1[r - l + 1] % mod1) % mod1;
        int ans2 = (mod2 + h2[r + 1] - h2[l] * val2[r - l + 1] % mod2) % mod2;
        return {ans1, ans2};
    }

    puv modify(int idx, char x) 
    { //修改 idx 位为 x
        int n = s.size() - 1;
        int ans1 = (h1.back() + val1[n - idx] * (x - s[idx]) % mod1) % mod1;
        int ans2 = (h2.back() + val2[n - idx] * (x - s[idx]) % mod2) % mod2;
        return {ans1, ans2};
    }
};


```]
 #pagebreak() 
== `KMP.h`


 #sourcecode[```cpp
#include <template_overAll.h>

class KMP
{
private:
    string s;
    string inis;
public:
    vector<int> pi;
    KMP(string _s)
    {
        s = _s;
        inis = s;
    }
    void prefix_function()
    {
        pi.clear();
        int n = (int)s.length();
        pi.resize(n);
        for (int i = 1; i < n; i++)
        {
            int j = pi[i - 1];
            while (j > 0 && s[i] != s[j])
                j = pi[j - 1];
            if (s[i] == s[j])
                j++;
            pi[i] = j;
        }
        return;
    }
    vector<int> find_occr(string p)
    {
        s = inis;
        s = p + "#" + s;
        prefix_function();
        vector<int> v;
        for (int i = p.size() + 1; i < s.size(); i++)
            if (pi[i] == p.size())
                v.push_back(i - 2 * p.size());
        return v;
    }
};
 
```]
 #pagebreak() 
== `trie_Tree.h`


 #sourcecode[```cpp
#include <template_overAll.h>

class Trie//AC
{
public:
    vector<map<char, int>> t;
    int root = 0;
    Trie()
    {
        t.resize(1);
    }
    void addedge(string _s)
    {
        int pvidx = root;
        _s.push_back('-');
        for (int i = 0; i < _s.size(); i++)
        {
            if (t[pvidx].find(_s[i]) != t[pvidx].end())
            {
                pvidx = t[pvidx][_s[i]];
            }
            else
            {
                t[pvidx][_s[i]] = t.size();
                t.push_back(map<char, int>());
                pvidx = t[pvidx][_s[i]];
            }
        }
    }
    bool ifcmp(string &s)
    {
        int pvidx = root;
        for(int i = 0;i < s.size();i ++)
        {
            if(t[pvidx].find(s[i]) != t[pvidx].end()) pvidx = t[pvidx][s[i]];
            else return 0;
        }
        return t[pvidx].find('-') != t[pvidx].end();
    }
};
```]
 #pagebreak() 
