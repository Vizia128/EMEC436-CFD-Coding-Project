using Latexify

struct ttest
    aa
    bb
    cc
    dd
end

function ttest(;aa=1, bb=2, cc=3, dd=4)
    return ttest(aa, bb, cc, dd)
end

function foo(TT::ttest)
    (;aa, bb, dd) = TT
    println(dd)
    return bb
end

TT = ttest()

ff = foo(TT)