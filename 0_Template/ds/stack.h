// ##IGNORE##
template <class TYPE_NAME>
class stackElememt{
public:
    stackElememt *top;
    stackElememt *base;
    TYPE_NAME data;
};

template <class TYPE_NAME>
class mystack{
private:
    // stackElememt<T> *bedrock;
    stackElememt<TYPE_NAME> *upper;
public:
    mystack()
    {
        stackElememt<TYPE_NAME> *st = new stackElememt<TYPE_NAME>;
        // bedrock = st;
        // bedrock->base = bedrock;
        // bedrock->top = bedrock;
        upper = st;
        upper->base = nullptr;
        upper->top = upper;
    }
    void push(TYPE_NAME _t)
    {
        stackElememt<TYPE_NAME> *st = new stackElememt<TYPE_NAME>;
        st->data = _t;
        upper->top = st;
        st->top = st;
        st->base = upper;
        upper = upper->top;
    }
    TYPE_NAME pop()
    {
        TYPE_NAME rt = upper->data;
        stackElememt<TYPE_NAME> *st = upper;
        upper->base->top = upper->base;
        upper = upper->base;
        free(st);
        return rt;
    }
    ~mystack()
    {
        for(;;)
        {
            if(upper == nullptr) break;
            stackElememt<TYPE_NAME> *st = upper;
            upper = upper->base;
            free(st);
        }
    }
};
