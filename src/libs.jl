using Makie;
using Statistics;
using FFTW;
using AbstractPlotting;
using Observables;

function readWhole(fileName, tip::T) where {T <: Number}
        buff = (fileName |> stat).size |> i -> zeros(T, i / sizeof(T) |> Int64);
        stream = open(fileName);
        buff = read!(stream, buff);
        close(stream);
        return buff;
end

function rdc(a, n)
        sz = size(a)[1]
        [Statistics.mean([abs(a[i]) for i=j-n:j]) for j=(1+n:n:sz)]
end
function rdc(a::AbstractArray{T, 1}, n) where {T<:Unsigned}
        sz = size(a)[1]
        [(Statistics.mean([abs(a[i]) for i=j-n:j]) |> floor |> T) for j=(1+n:n:sz)]
end
function removeNoise!(arr::Array{Complex{Float64}, 1}, fact=10)
        mi = [abs(i) for i=arr[3:end]] |> ii -> Statistics.mean(ii)*fact;
        for i in 1:length(arr)
                if abs(arr[i]) < mi
                        setindex!(arr, complex(Float64(0)), i);
                end
        end
end
function removeNoiseBPSH!(arr::Array{Complex{Float64}, 1}, takeAround=1000000)
        maxim = [abs(i) for i=arr[3:end]] |> findmax;
        for i in 1:(maxim[2] - takeAround)
                setindex!(arr, complex(Float64(0)), i);
        end
        len = length(arr);
        for i in (maxim[2] + takeAround):len
                setindex!(arr, complex(Float64(0)), i);
        end
end
function removeIm!(arr::Array{Complex{Float64}, 1})
        aa = Array{Float64, 1}()
        for i in 1:length(arr)
                setindex!(arr, abs(arr[i]), i);
        end
        Array{Float64, 1}(arr)
end
function plotarray(arr, everyEnth=1, sampleTime::Real = 1, odPosto=0, doPosto=1)
        aView = a[(end*odPosto |> floor |> Int64):(end*doPosto |> floor |> Int64)];
        startTime = length(a)*odPosto;
        plot(i->(i + startTime)/sampleTime, i-> aView[i], 1:everyEnth:length(aView), xlabel=:time, ylabel=:amplituda)
end
function fftss(arr, odPosto, doPosto)
        aView = a[(end*odPosto |> floor |> Int64):(end*doPosto |> floor |> Int64)];
        FFTW.fft(aView)
end
function powerSpectar(arr, sampTime, fps=60, seconds=1)
        dataPFrame = sampTime/fps |> floor |> Int64;
        samples = seconds*fps;
        [fft(arr[(i*dataPFrame+1):((i+1)*dataPFrame+1)]) for i=0:samples] |> tm -> [removeIm!(i) for i=tm] |> tm -> [i[(end/2 |> floor |> Int64):end] for i=tm]
end
function powerSpectar(arr, sampTime, ti, fps=60, seconds=1)
        dataPFrame = sampTime/fps |> floor |> Int64;
        samples = seconds*fps;
        [i for i=0:samples]
end

function powerSpectarElement(i, arr, sampTime, fps=60)
        dataPFrame = sampTime/fps |> floor |> Int64;
        fftsplice(arr[i*dataPFrame+1:(i+1)*dataPFrame+1])
end
function timeDomainElement(arr, sampTime, i, fps=60)
        dataPFrame = sampTime/fps |> floor |> Int64;
        arr[i*dataPFrame+1:(i+1)*dataPFrame+1]
end
function displayRealTime!(s, ti, arr)
        for i=1:1:600
            push!(ti, i%60 + 1)
            force_update!();
            sleep(1/60)
        end
        return s;
end
function fftsplice(arr)
        arr[1:end+(end%2-1)] |>
        fft .|>
        abs |>
        i -> i[1:(end/2 |> ceil |> Int64)] |> #+ i[(end/2 |> ceil |> Int64):-1:1] |>
        i -> i[3:end-1]
end
function bpskFistStage(t, ar, sampleTime = 1800000, fps=60)
    raw = timeDomainElement(ar, sampleTime, t, fps)
    signalFreq = sampleTime*(7431/30000+0.5)
    steper = (1+length(raw)*t:length(raw)*(t+1))/sampleTime
    [cos(2*pi*signalFreq*i) for i=steper]  .* raw
end


function integration(arr::Array{T, 1}, sampTime) where {T <: Number}
        retArr = zeros(Float64, length(arr));
        retArr[1] = arr[1]/sampTime;
        for i=2:length(arr)
                retArr[i] = retArr[i - 1] + arr[i]/sampTime;
        end
        retArr
end

function analizeFile(filePath, sampleRate=1800000)
        ar=readWhole(filePath, UInt16(0));
        s1, a = AbstractPlotting.textslider(0:100, "1000t", start = 0);
        smin, amin = AbstractPlotting.textslider(0:100, "10tmin", start = 0);
        smicro, amicro = AbstractPlotting.textslider(1:10, "t", start = 1);
        AbstractPlotting.button()
        treal = lift((t1, t2, t3) -> t1*1000 + t2*10 + t3, a, amin, amicro);
        sFps, fps = AbstractPlotting.textslider(10:500, "fps", start = 30);
        timeLift = lift((treal, fps) -> "$(round(treal/fps, digits=3))s", treal, fps);

        myLift = lift(powerSpectarElement, treal, ar, sampleRate, fps);
        fftElements = lines(myLift);
        myLift2 = lift(bpskFistStage, treal, ar, sampleRate, fps);
        timeElements = lines(myLift2);
        integrator_lift = lift(integration, myLift2, sampleRate);
        integrator_plot = lines(integrator_lift);
        integrator_plot = scatter!(integrator_plot, [(0, 5), (0, -5)]);
        seconds_label = text(timeLift)

        vbox(hbox(s1, smin, smicro, sFps, seconds_label), fftElements, integrator_plot, parent=Scene())
end
