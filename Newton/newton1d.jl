#=
f(x) = x^3 + 2x^2 -5x +6
-5<= x <= 3 で図示せよ
=#

# 開始値:ステップ:終了値
x = -5:0.1:3.0
# 関数
f(x) = x^3+2x^2-5x+6
y = f

using Plots
plot(x,y)
savefig("newton1d.png")