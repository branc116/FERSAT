
function readWhole(fileName, tip::T) where {T <: Integer}
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

function powerSpectarElement(arr, sampTime, i, fps=60)
        dataPFrame = sampTime/fps |> floor |> Int64;
        fftsplice(arr[i*dataPFrame+1:(i+1)*dataPFrame+1])
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
        i -> i[(end/2 |> ceil |> Int64):end] + i[(end/2 |> ceil |> Int64):-1:1] |>
        i -> i[3:end-1]
end
