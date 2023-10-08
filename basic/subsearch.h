#include <template_overAll.h>

// ac
template <class T>
int subsearch_1(T tar, vector<T> &arr) // 第一类二分："非下降"序列，"找到"就返回元素下标，相等元素返回最"左侧"下标，否则返回-1
{
    int l = 0;
    int r = arr.size();
    int ans = -1;
    while (r >= l)
    {
        int m = l + (r - l) >> 1;
        if (arr[m] >= tar)
        {
            if (arr[m] == tar) ans = m;
            r = m - 1;
        }
        else l = m + 1;
    }
    return ans;
}

template <class T>
int subsearch_2(T tar, vector<T> &arr) // 第二类二分：取小于或等于tar的最大的元素
{
    int l = 0;
    int r = arr.size();
    int m = 0;
    while (r > l + 1)
    {
        // m = l + (r - l) >> 1;
        m = (l + r) >> 1;
        if (arr[m] > tar) r = m;
        else l = m;
    }
    return l;
}

template <class T>
int subsearch_3(T tar, vector<T> &arr)//第三类二分：取大于或等于tar的最小的元素
{
    int l = 0;
    int r = arr.size();
    int m = 0;
    while (r > l + 1)
    {
        // m = l + (r - l) >> 1;
        m = (l + r) >> 1;
        if (arr[m] > tar) r = m;
        else l = m;
    }
    return r;
}
