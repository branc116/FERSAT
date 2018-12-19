using Makie
scene = Scene()
ti = Node(0.1)
f(v, t) = sin.(v)

scene = lines!(scene,lift(t -> f.((0:0.01:10)*t, t), ti))


push!(ti, 10)

for i=0.1:0.1:50
    push!(ti, i)
    text!(scene, "t=$i")
    force_update!();
    sleep(1/24)
end
