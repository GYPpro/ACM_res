class PA {
 private:
  static const int N = 100010;
  struct Node {
    int len;
    int ptr[26], fail;
    Node(int len = 0) : len(len), fail(0) { memset(ptr, 0, sizeof(ptr)); }
  } nd[N];

  int size, cnt;  // size为字符串长度，cnt为节点个数
  int cur;  //当前指针停留的位置，即最后插入字符所对应的节点
  char s[N];

  int getfail(int x)  //沿着fail指针找到第一个回文后缀
  {
    while (s[size - nd[x].len - 1] != s[size]) {
      x = nd[x].fail;
    }
    return x;
  }

 public:
  PA() : size(0), cnt(0), cur(0) {
    nd[cnt] = Node(0);
    nd[cnt].fail = 1;
    nd[++cnt] = Node(-1);
    nd[cnt].fail = 0;
    s[0] = '$';
  }

  void extend(char c) {
    s[++size] = c;
    int now = getfail(cur);  //找到插入的位置
    if (!nd[now].ptr[c - 'a'])  //若没有这个节点，则新建并求出它的fail指针
    {
      int tmp = ++cnt;
      nd[tmp] = Node(nd[now].len + 2);
      nd[tmp].fail = nd[getfail(nd[now].fail)].ptr[c - 'a'];
      nd[now].ptr[c - 'a'] = tmp;
    }
    cur = nd[now].ptr[c - 'a'];
  }

  int qlen() { return nd[cur].len; }
}