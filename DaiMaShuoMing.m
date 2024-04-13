% The following code snippet comprises publicly available code for simulating experimental procedures, with accompanying descriptions as follows:

% wavelet_denoising: Implements wavelet denoising.
% VMD_denoising: Performs denoising using Variational Mode Decomposition (VMD).
% DDTF_denoising: Implements denoising using Discrete Dyadic Transform Filtering (DDTF).
% waveletDDTF_denoising: Applies joint sparse representation with wavelet and DDTF for denoising.
% col2imstep: Facilitates data reordering.
% countcover: Counts the number of times each data point is covered by blocks during data reordering.
% displayDictionaryElementsAsImage: Visualizes dictionary elements.
% learnt_dict1=filter_learning(): Executes threshold processing.
% OMP: Utilizes Orthogonal Matching Pursuit (OMP) for sparse coding of columns using learned dictionaries.
% Parameters for simulated signals:

% fs = 30e3: Sampling frequency.
% fn = 2e3/1: Natural frequency.
% y0 = 10: Displacement constant.
% g = 0.1: Damping coefficient.
% T = 0.005*2: Repetition period.
% N = 4096: Number of sample points.
% NT = round(fs*T): Number of sample points per single period.
% t = 0:1/fs:(N-1)/fs: Sampling instances.
% t0 = 0:1/fs:(NT-1)/fs: Sampling instances per single period.
% K = ceil(N/NT) + 1: Repetition count.