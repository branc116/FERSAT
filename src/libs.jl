using Makie;
using Statistics;
using FFTW;
using AbstractPlotting;
using Observables;

include("GuiLibs.jl");
include("commonLibs.jl");

function analizeComplex(filePath, sampleRate=1800000)
    s1000, a1000 = AbstractPlotting.textslider(0:10, "1000t", start = 0);
    s100, a100 = AbstractPlotting.textslider(0:10, "100t", start = 0);
    s10, a10 = AbstractPlotting.textslider(0:10, "10t", start = 0);
    s1, a1 = AbstractPlotting.textslider(1:10, "t", start = 1);
    sSegSize, segSize = AbstractPlotting.textslider(1000:1000:50000, "Segment Size", start = 1000);
    sDecimate, decim = AbstractPlotting.textslider(0:10, "Decimate", start = 0);

    t = lift((t1, t2, t3, t4) -> t1*1000 + t2*100 + t3*10 + t4, a1000, a100, a10, a1);
    tSeconds = lift((t, sampleRate, segmentSize) -> "$(round(t*segmentSize/sampleRate, digits=3))s", t, sampleRate, segSize);

    _curSegment = lift(readFileSegment, Node(Complex{Float32}), filePath, t, segSize);
    _fft = lift(myFft, _curSegment, decim);
    _shiftF = lift(getShiftFreq, _fft, sampleRate);
    
    _shifted = lift(shift, _curSegment, t, _shiftF, sampleRate);
    _lopas = lift(lopass, _shifted);
    _unwrap = lift(unwrapCpl, _lopas);
    _unwrapDeriv = lift(diffNaiv, _unwrap);
    _sign = lift(sign, _unwrapDeriv);
    
    _text = FerSatGui.createNormalText(tSeconds, (10, 10));
    _shiftFreq = FerSatGui.createNormalText(lift((sf) -> "Shift Frequency: $sf", _shiftF), (10, 10));

    FerSatGui.join(
            FerSatGui.createMultiPlot([
                FerSatGui.NamedObsevableArray(_curSegment, "Current segment"),
                FerSatGui.NamedObsevableArray(_fft, "Fft"),
                FerSatGui.NamedObsevableArray(_shifted, "Shifted"),
                FerSatGui.NamedObsevableArray(_lopas, "Lo pass"),
                FerSatGui.NamedObsevableArray(_unwrap, "Unwrapped"),
                FerSatGui.NamedObsevableArray(_unwrapDeriv, "Unwrapped derivation"),
                FerSatGui.NamedObsevableArray(_sign, "Sign")
                 ], 2

        ), [s1, s10, s100, s1000, sSegSize, sDecimate, _shiftFreq, _text])
    # sFps, fps = AbstractPlotting.textslider(10:500, "fps", start = 30);
end
