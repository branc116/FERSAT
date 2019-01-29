# using Plots;
using FFTW;
using DSP;
using Statistics;
# read I/Q samples

function shift(data::T, segment::Number, shiftFreq::Number=23000, sampleTime=1800000) where {T<:AbstractArray{<:Number}} 
    len = length(data);
    constant = (-im*2*pi*shiftFreq/sampleTime);
    return exp.(([(len*(segment - 1)):(len*segment - 1)] .* constant)...) .* data;
end
function lopass(data::T, sampleTime=1800000) where {T<:AbstractArray{<:Number}}
    responsetype = Lowpass(1500; fs=sampleTime);
    prototype = Butterworth(40);
    return filt(digitalfilter(responsetype, prototype), data)
end
function unwrapCpl(data::T) where {T<:AbstractArray{<:Complex}}
    return angle.(data) |> unwrap
end
function diffNaiv(x::T, i::Integer) where {T1<:Real, T<:AbstractArray{T1}}
     return i > 1 && i <= length(x) ? x[i] - x[i-1] : T1(0)
end
function diffNaiv(x::T) where {T1<:Real, T<:AbstractArray{T1}}
    return [diffNaiv(x, i) for i=1:length(x)]
end
function sign(data::T) where {T<:AbstractArray{<:Real}}
    return Base.sign.(data);
end


function readWhole(fileName, tip::T) where {T <: Number}
    buff = (fileName |> stat).size |> i -> zeros(T, i / sizeof(T) |> Int64);
    stream = open(fileName, 'r');
    buff = read!(stream, buff);
    close(stream);
    return buff;
end
function readFileSegment(::Type{G}, fileName::String, segment::UnitRange{T}, ) where {T <: Integer, G <: Number}
    buff = zeros(G, length(segment));
    stream = open(fileName);
    seek(stream, sizeof(G)*(segment[1] - 1));
    read!(stream, buff);
    close(stream);
    buff
end
function readFileSegment(::Type{G}, fileName::String, segment::T2, segmentSize::T2=10000)  where {T <: Integer, G <: Number, T1 <: Integer, T2 <:Integer} 
    return readFileSegment(G, fileName, ((segment)*segmentSize):(((segment + 1)*segmentSize) - 1));
end

function rdc(a::T, n::Integer) where {T <: AbstractArray{<:Real}}
    sz = length(a)
    [Statistics.mean(a[(j-n):j]) for j=((1+n):n:sz)]
end
function travelingMean(a::T, n::Integer) where {T <: AbstractArray{<:Real}}
    sz = length(a)
    [Statistics.median(a[(j-n):j]) for j=((1+n):1:sz)]
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
function fftss(arr, odPosto, doPosto)
        aView = arr[(end*odPosto |> floor |> Int64):(end*doPosto |> floor |> Int64)];
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

function timeDomainElement(arr, sampTime, i, fps=60)
        dataPFrame = sampTime/fps |> floor |> Int64;
        arr[i*dataPFrame+1:(i+1)*dataPFrame+1]
end
function bpskFistStage(t, ar, sampleTime = 1800000, fps=60)
    raw = timeDomainElement(ar, sampleTime, t, fps)
    signalFreq = sampleTime*(7431/30000+0.5)
    steper = (1+length(raw)*t:length(raw)*(t+1))/sampleTime
    [cos(2*pi*signalFreq*i) for i=steper]  .* raw
end

function toTuple(data::T) where {T<:AbstractArray{<:Complex}}
    return [(Float64(i.re), Float64(i.im)) for i=data];
end
function toTuple(data::T) where {T<:AbstractArray{<:Real}}
    return [(Float64(i), Float64(data[i])) for i=1:length(data)];
end
function toTuple(data::T) where {T<:AbstractArray{Tuple{Float64, Float64}}}
    return data;
end

function integration(arr::Array{T, 1}, sampTime) where {T <: Number}
        retArr = zeros(Float64, length(arr));
        retArr[1] = arr[1]/sampTime;
        for i=2:length(arr)
                retArr[i] = retArr[i - 1] + arr[i]/sampTime;
        end
        retArr
end

# function gimme()
#     data = reinterpret(Complex{Float32}, read("/home/domagoj/gqrx_20181207_153804_437485000_1800000_fc.raw"));
#     dataLength = floor(Int, size(data,1));

#     # multiply by complex exponential -> frequency shift; desired carrier frequency is shifted to origin
#     T = 1/1800000;
#     shiftFreq = 23000;
#     shiftExp = [exp(-im*2*pi*shiftFreq*n*T) for n in 1:dataLength];
#     data = data.*shiftExp;

#     # create Butterworth lowpass filter
#     responsetype = Lowpass(1500; fs=1800000);
#     prototype = Butterworth(40);

#     # apply created filter to sample array
#     db = filt(digitalfilter(responsetype, prototype), data);
# end
# try
#     db;
# catch
#     db = gimme();
# end

# angularv = (x, i) -> i> 1 && i <= length(x) ? angle(x[i]/x[i-1]) : 0
# dbAng = angle.(db) |> unwrap;
# #dbd = [diffNaiv(dbAng, i) for i=1:length(dbAng)];
# #dbav = [angularv(dbAng, i) for i=1:length(dbAng)];
# #dbdAvg = [Statistics.median((dbd[i-30:i])) for i=31:30:length(dbd)];
# # Q = imag.(data) - real.(data) .|> sign ;
# Qd = [diffNaiv(dbAng, i) for i=2:length(dbAng)];
# # dbd = sign.(dbd);
# # plot d_phi / dt
# windowSize = 1024;
# n = floor(Int, size(Q,1)/windowSize);
# b = true;
# function a(arrr, amp, frame, ws)
#     for i in (18 * ws):(20 * ws)
#         currr = arrr[(i-1)*windowSize+1:i*windowSize];
# #		if (currr |> findmin)[1] >  | true
#             plot(currr)
#             if (amp > 0)
#                 ylims!(-amp,amp)
#             end
#             #xlims!(-.01,.01)
#             gui()
#             sleep(frame)
#             if !b
#                 break;
#             end
#             println(i);
# #		end
#     end
#     println("Stopped");
# end
# @async a(Qd, 0.01, 0.04, windowSize);
#=
# plot "live" constellation diagram
windowSize = 10240;
n = floor(Int, size(db,1)/windowSize);
for i in 1:n
    plot(db[(i-1)*windowSize+1:i*windowSize])
    ylims!(-.01u,.01)
    xlims!(-.01,.01)
    gui();
    sleep(0.04)
end
=#

# EOF
