using Makie;
using Statistics;
using FFTW;
include("libs.jl")

sampleTime = 1800000;
a=readWhole("RawFiles\\gqrx_20181207_153907_437485000_1800000_fc.raw", UInt16(0));
s = Scene()
ti = Node(1);
push!(ti, 1);
s = lines!(s,lift(t -> powerSpectarElement(a, sampleTime, t)[10:end] .|> log, ti))
for i=1:1:600
    push!(ti, i)
    force_update!();
    sleep(1/60)
end
displayRealTime!(s, ti, specs)
arr=specs
