# using Plots;
using FFTW;
using DSP;
# read I/Q samples

function shift(data::T, segment::Number, shiftFreq::Number, sampleTime=1800000) where {T<:AbstractArray{<:Number}} 
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
function myFft(data)
    fft(data) |> fftshift .|> abs .|> log10
end
function isData(fftData, sampleRate=1800000)
    a = fftData |> findmax
    if (a[1] > 1)
        return true
    end
    return false
end
function getShiftFreq(fftData, sampleRate=1800000)
    a = fftData |> findmax
    len=length(fftData);
    if (a[1] > 1)
        bin = a[2];
        return bin*sampleRate/len;
    end
end
function bitrate(fftdata, )
    
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
