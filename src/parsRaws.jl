using FerSat
sampleTime = 1800000;
analizeFile("RawFiles\\gqrx_20181207_153907_437485000_1800000_fc.raw", sampleTime)
ar=readWhole("RawFiles\\gqrx_20181207_153804_437485000_1800000_fc.raw", Float32(1)im);
using Plots;plotly();
using FFTW;
length(ar)
interest = (ar .|> abs |>  findmax)[2]
unique(ar)
ar |> unique
ar[interest - 10000:interest + 10000] |> fft |> fftshift .|> abs .|> log10  |> plot
