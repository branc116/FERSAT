using Plots;
using FFTW;
using DSP;

# read I/Q samples
function gimme()
	data = reinterpret(Complex{Float32}, read("/home/domagoj/gqrx_20181207_153804_437485000_1800000_fc.raw"));
	dataLength = floor(Int, size(data,1));

	# multiply by complex exponential -> frequency shift; desired carrier frequency is shifted to origin
	T = 1/1800000;
	shiftFreq = 23000;
	shiftExp = [exp(-im*2*pi*shiftFreq*n*T) for n in 1:dataLength];
	data = data.*shiftExp;

	# create Butterworth lowpass filter
	responsetype = Lowpass(1500; fs=1800000);
	prototype = Butterworth(40);

	# apply created filter to sample array
	db = filt(digitalfilter(responsetype, prototype), data);
end
try
	db;
catch
	db = gimme();
end

angularv = (x, i) -> i> 1 && i <= length(x) ? angle(x[i]/x[i-1]) : 0
dbAng = angle.(db) |> unwrap;
#dbd = [diffNaiv(dbAng, i) for i=1:length(dbAng)];
#dbav = [angularv(dbAng, i) for i=1:length(dbAng)];
#dbdAvg = [Statistics.median((dbd[i-30:i])) for i=31:30:length(dbd)];
diffNaiv = (x, i) -> i> 1 && i <= length(x) ? x[i] - x[i-1] : typeof(x[1])(0);
Q = imag.(data) - real.(data) .|> sign ;
Qd = [diffNaiv(dbAng, i) for i=2:length(dbAng)];
dbd = sign.(dbd);
# plot d_phi / dt
windowSize = 1024;
n = floor(Int, size(Q,1)/windowSize);
b = true;
function a(arrr, amp, frame, ws)
	for i in (18 * ws):(20 * ws)
		currr = arrr[(i-1)*windowSize+1:i*windowSize];
#		if (currr |> findmin)[1] >  | true
			plot(currr)
			if (amp > 0)
				ylims!(-amp,amp)
			end
			#xlims!(-.01,.01)
			gui()
			sleep(frame)
			if !b
				break;
			end
			println(i);
#		end
	endcd
	println("Stopped");
end
@async a(Qd, 0.01, 0.04, windowSize)
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
