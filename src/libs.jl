using Makie;
using Statistics;
using FFTW;
using AbstractPlotting;
using Observables;

include("GuiLibs.jl");
include("commonLibs.jl");


function analizeFile(filePath, sampleRate=1800000)
        
        s1, a = AbstractPlotting.textslider(0:100, "1000t", start = 0);
        smin, amin = AbstractPlotting.textslider(0:100, "10tmin", start = 0);
        smicro, amicro = AbstractPlotting.textslider(1:10, "t", start = 1);
        sFps, fps = AbstractPlotting.textslider(10:500, "fps", start = 30);
        ar=readWhole(filePath, UInt16(0));
        treal = lift((t1, t2, t3) -> t1*1000 + t2*10 + t3, a, amin, amicro);
        timeLift = lift((treal, fps) -> "$(round(treal/fps, digits=3))s", treal, fps);
        myLift = lift(powerSpectarElement, treal, ar, sampleRate, fps);
        myLift2 = lift(bpskFistStage, treal, ar, sampleRate, fps);
        integrator_lift = lift(integration, myLift2, sampleRate);
        seconds_label = text(timeLift)


        FerSatGui.join(
            FerSatGui.createMultiPlot([
                FerSatGui.NamedObsevableArray(integrator_lift, "Integrate"),
                FerSatGui.NamedObsevableArray(myLift2, "Time domain"),
                FerSatGui.NamedObsevableArray(myLift, "FFt")], 3
        ), [s1, smin, smicro, sFps, seconds_label])
        # vbox(hbox(s1, smin, smicro, sFps, seconds_label), fftElements, integrator_plot, parent=Scene())
end
function analizeComplex(filePath, sampleRate=1800000)
    s1000, a1000 = AbstractPlotting.textslider(0:10, "1000t", start = 0);
    s100, a100 = AbstractPlotting.textslider(0:10, "100t", start = 0);
    s10, a10 = AbstractPlotting.textslider(0:10, "10t", start = 0);
    s1, a1 = AbstractPlotting.textslider(1:10, "t", start = 1);
    sSegSize, segSize = AbstractPlotting.textslider(1000:1000:50000, "Segment Size", start = 1000);
    sShiftFreq, _shiftFreq = AbstractPlotting.textslider(20000:10:25000, "Shift Frequency", start = 23000);
    # sTravleMean, _travleMean = AbstractPlotting.textslider(1:1:100, "Travle mean", start = 10);

    t = lift((t1, t2, t3, t4) -> t1*1000 + t2*100 + t3*10 + t4, a1000, a100, a10, a1);
    tSeconds = lift((t, sampleRate, segmentSize) -> "$(round(t*segmentSize/sampleRate, digits=3))s", t, sampleRate, segSize);
    # curInterval = lift((t, segmentSize) -> ((t)*segmentSize):(((t + 1)*segmentSize) - 1), t, segSize);

    _curSegment = lift(readFileSegment, Node(Complex{Float32}), filePath, t, segSize);
    # curSegmentTuple = lift(toTuple, _segment);

    _shifted = lift(shift, _curSegment, t, _shiftFreq, sampleRate);
    # shiftedTuple = lift(toTuple, shifted);

    _lopas = lift(lopass, _shifted);
    # _lopasTuple = lift(toTuple, _lopas);

    _unwrap = lift(unwrapCpl, _lopas);
    # _unwrapTuple = lift(toTuple, _unwrap);

    _unwrapDeriv = lift(diffNaiv, _unwrap);
    # _unwrapDerivTuple = lift(toTuple, _unwrapDeriv);

    _sign = lift(sign, _unwrapDeriv);
    # _signTuple = lift(toTuple, _sign);
    _signTuple = lift(toTuple, _sign);
    _fft = lift(myFft, _lopas);
    _fftTuple = lift(toTuple, _fft);
    _text = FerSatGui.createNormalText(tSeconds, (10, 10));
    FerSatGui.join(
            FerSatGui.createMultiPlot([
                FerSatGui.NamedObsevableArray(_curSegment, "Current segment"),
                FerSatGui.NamedObsevableArray(_shifted, "Shifted"),
                FerSatGui.NamedObsevableArray(_lopas, "Lo pass"),
                FerSatGui.NamedObsevableArray(_unwrap, "Unwrapped"),
                FerSatGui.NamedObsevableArray(_unwrapDeriv, "Unwrapped derivation"),
                FerSatGui.NamedObsevableArray(_sign, "Sign"), ], 2

        ), [s1, s10, s100, s1000, sSegSize, sShiftFreq, _text])
    # sFps, fps = AbstractPlotting.textslider(10:500, "fps", start = 30);
end
