
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
