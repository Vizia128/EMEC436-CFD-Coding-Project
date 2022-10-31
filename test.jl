using Latexify, Interpolations

xx = [i+j+k for i in 0:32, j in 0:32, k in 0:32]
yy = [i^2 + j^2 + k^2 for i in -16:16, j in -16:16, k in -16:16]

xxi = interpolate(xx, BSpline(Linear()))
yyi = interpolate(yy, BSpline(Linear()))

foo(x,y) = Point(xxi[x,y,8], yyi[x,y,8])
streamplot(foo, 1:32, 1:32)

tt = Observable(1)

xxit = @lift()

function getPoint(x,y,t)
    p = Point()
    function gp(x,y)

fii(x,y) = Point(xxi[x,y,tt], yyi[x,y,tt])
streamplot(fii, 1:32, 1:32)