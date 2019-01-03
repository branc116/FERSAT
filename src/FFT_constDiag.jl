using Plots;
using FFTW;
using DSP;

# read I/Q samples
data = reinterpret(Complex{Float32}, read("/home/domagoj/gqrx_20190101_155539_97600000_1800000_fc.raw"));
dataLength = floor(Int, size(data,1));

# signal is complex, doesn't have symmetrical FFT spectrum (with respect to N/2 bin)
# problem -> any real filter has symmetrical characteristic -> undesired frequencies remain unfiltered
# so we need to multiply by complex exponential -> frequency shift
# desired carrier frequency is shifted to origin, therefore we can apply any low-pass filter
T = 1/1800000;
shiftFreq = -57000;
shiftExp = [exp(-im*2*pi*shiftFreq*n*T) for n in 1:dataLength];
data = data.*shiftExp; 

# create Butterworth lowpass filter
responsetype = Lowpass(10000; fs=1800000);
prototype = Butterworth(40);

# apply created filter to data array
db = filt(digitalfilter(responsetype, prototype), data);

# generate and plot 32768-FFT results
windowSize = 32768;
n = floor(Int, size(db,1)/windowSize);
for i in 1:n
    F = fft(db[(i-1)*windowSize+1:i*windowSize]);
    F = fftshift(F);
    Fr = [20*log10(abs(i)) for i in F];
    plot(Fr)
    ylims!(-30, 90);
    gui();
    sleep(0.04)
end

# downsampling (don't go below 2*filter_corner_frequency)
decimat = 30;
db = db[1:decimat:end]; 

# plot "live" constellation diagram
windowSize = 1024;
n = floor(Int, size(db,1)/windowSize);
for i in 1:n
    plot(db[(i-1)*windowSize+1:i*windowSize])
    ylims!(-2,2)
    xlims!(-2,2)
    gui();
    sleep(0.04)
end

# TODO
# -demodulation
# -interactivity

# EOF

