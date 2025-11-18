using Plots
using DelimitedFiles

dT = 5.2
Tsat = 373.15
alphaL = 0.1689e-6
rhol = 958.4
rhov = 0.597
hlv = 2.26e6
cp = 4216
Ja = (cp * dT * rhol)/(hlv * rhov)
B = Ja * (12 * alphaL / pi)^0.5

gr()
# plotlyjs()

inlist = [
    "yourOutputFile.csv",
]

diameter(v) = (6*v/pi)^(1/3)
radius(v) = (3*v/(4*pi))^(1/3)

maity = [readdlm("maity_bubble1.csv", ',', Float64),
         readdlm("maity_bubble2.csv", ',', Float64),
         readdlm("maity_bubble3.csv", ',', Float64)]

maityTime = [i[:,1] for i in maity]
maityDiameter = [i[:,2] for i in maity]

t0 = maityTime[1][1]

maityShapes = [:xcross :diamond :utriangle]
maityLegend = ["Bubble 1" "Bubble 2" "Bubble 3"]

scatter(maityTime.*1e3, maityDiameter.*1e3, markershape=maityShapes, label=maityLegend, grid=false, ylabel="Bubble diameter [mm]",xlabel="Time [ms]", legend=:bottomright, markersize=6,framestyle=:box,fontfamily="Computer Modern",xtickfontsize=10,ytickfontsize=10,legendfontsize=10,guidefontsize=10, xlim=(0,:auto))

simLabels = "Numerical"

sim = [readdlm(i, '\t', Float64) for i in inlist]

simTime = [i[:,1] .+ t0 for i in sim]
simDiameter = [radius.(i[:,2].*2).*2 for i in sim]
plot!(simTime.*1e3, simDiameter.*1e3, label=simLabels, lw=3)
