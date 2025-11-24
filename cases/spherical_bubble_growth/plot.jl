using Plots
Plots.default(show = true)
using DelimitedFiles
using SpecialFunctions

gr()
# plotlyjs()

# Set constants

dT = 5
Tsat = 373.15
alphaL = 0.1689e-6
rhol = 958.4
rhov = 0.597
hlv = 2.26e6
cp = 4216
R0 = 9.69812679359357e-05
Ja = (cp * dT * rhol)/(hlv * rhov)
B = Ja * (12 * alphaL / pi)^0.5
t0 = R0^2 / B^2

# Define functions

z(r) = 0.5*(3^0.5)*r/(alphaL * t0)^0.5

R_from_V(v) = ((3/4)*v/pi)^(1/3)   # calculates radius of bubble from wedge volume
R(t) = B*(t^0.5)       # analytical radius
V(r) = (r^3)*pi/54              # calculates volume of bubble wedge from radius

offset_time(t) = t.+t0
offset_volume(v) = 2*(360/5)*v

get_time(x) = offset_time(x[:,1])
get_vol(x) = offset_volume(x[:,2])

R_array(vol) = [R_from_V(i) for i in vol]

function T(R)
    if R < R0
        return Tsat
    else
        return Tsat + dT*erf(z(R-R0))
    end
end

# Import data

data = [
    readdlm("sampled_05micron.csv",'\t',Float64),
    readdlm("sampled_1micron.csv",'\t',Float64),
    readdlm("sampled_2micron.csv",'\t',Float64),
]


# Separate time and volume data

ana_time = 0:1e-6:0.6e-3       # time range for plotting analytical sol
time = get_time.(data)    # numerical

volume = get_vol.(data)

# Calculate radius for each

R_analytical = [R(i) for i in ana_time]
Rad = R_array.(volume)

# Plot data

lls = ["Numerical, Δx = 0.5μm" "Numerical, Δx = 1μm" "Numerical, Δx = 2μm"]

plot(ana_time*10^(3), R_analytical/R0, label="Analytical solution", xlabel="Time (ms)", ylabel="R(t)/R0", grid=false, lw=2,legendposition=:bottomright,framestyle=:box,xlim=(0,0.6))
scatter!(time.*10^(3), (Rad)./R0, label=lls, markersize=5, markershape=[:circle :dtriangle :square])
