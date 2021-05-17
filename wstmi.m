function d = wstmi(clean_speech, degraded_speech, Fs)
%     d = wstmi(clean_speech, degraded_speech, Fs) returns the output of the weighted
%     spectro-temporal modulation index (wSTMI) predictor.
% 
%     Implementation of the weighted spectro-temporal modulation index (wSTMI) predictor, 
%     described in [1].
% 
%     Inputs:
%            clean_speech:        clean reference time domain signal
%            degraded_speech:     noisy/processed time domain signal
%            Fs:                  sampling rate [Hz]
% 
%     Output:
%            d: intelligibility index

%     References:
%     [1] Edraki, A., Chan, W. Y., Jensen, J., & Fogerty, D. (2020). 
%     Speech Intelligibility Prediction Using Spectro-Temporal Modulation
%     Analysis. IEEE/ACM Transactions on Audio, Speech, and Language 
%     Processing, 29, 210-225.


%     Copyright 2021, Amin Edraki, a.edraki@queensu.ca, Queen's University
%     at Kingston, Multimedia Coding and Communications Laboratory.

%     wstmi is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     wstmi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with wstmi.  If not, see <https://www.gnu.org/licenses/>.

if length(clean_speech) ~= length(degraded_speech)
    error('clean and degraded signals should have the same length!');
end

clean_speech    = clean_speech(:);
degraded_speech = degraded_speech(:);


% resampling the input signals
if Fs ~= 10000
    clean_speech    = resample(clean_speech, 10000, Fs);
    degraded_speech = resample(degraded_speech, 10000, Fs);
    Fs              = 10000;
end


% Parameters of the log Mel-spectrogram
win_length   = 25.6;
win_shift    = win_length/2;
freq_range   = [64 min(floor(Fs./2), 12000)];
num_bands    = 130;
band_factor  = [];

% spectro-temporal modulation channel weights
w   = [0.000,0.031,0.140;
    0.013,0.041,0.055;
    0.459,0.528,0.000;
    0.151,0.000,0.000];

% linear model offset
b   = 0.16;


% voice activity detection
dyn_range = 40;
N_frame   = 256;
[clean_speech, degraded_speech] = removeSilentFrames(clean_speech, degraded_speech, dyn_range, N_frame, N_frame/2);


% spectro-temporal modulation envelope extraction
[X_spec, ~] = log_mel_spectrogram(clean_speech, Fs, win_shift, win_length, freq_range, num_bands, band_factor);
[Y_spec, ~] = log_mel_spectrogram(degraded_speech, Fs, win_shift, win_length, freq_range, num_bands, band_factor);
X           = sgbfb_4D(X_spec);
Y           = sgbfb_4D(Y_spec);


% row normalization
X   = X - mean(X, 4);
Y   = Y - mean(Y, 4);
X   = X ./ sqrt(sum(X.*X, 4));
Y   = Y ./ sqrt(sum(Y.*Y, 4));

% cross correlation
rho = nanmean(sum(X.*Y, 4), 3);

% intelligibility prediction
d   = sum(w(:).*rho(:)) + b;

end






function [x_sil, y_sil] = removeSilentFrames(x, y, range, N, K)
%     [X_SIL Y_SIL] = REMOVESILENTFRAMES(X, Y, RANGE, N, K) X and Y
%     are segmented with frame-length N and overlap K, where the maximum energy
%     of all frames of X is determined, say X_MAX. X_SIL and Y_SIL are the
%     reconstructed signals, excluding the frames, where the energy of a frame
%     of X is smaller than X_MAX-RANGE

%     Copyright 2009: Delft University of Technology, Signal & Information
%     Processing Lab. The software is free for non-commercial use. This program
%     comes WITHOUT ANY WARRANTY.

x       = x(:);
y       = y(:);

frames  = 1:K:(length(x)-N);
w       = hanning(N);
msk     = zeros(size(frames));

for j = 1:length(frames)
    jj      = frames(j):(frames(j)+N-1);
    msk(j) 	= 20*log10(norm(x(jj).*w)./sqrt(N));
end

msk     = (msk-max(msk)+range)>0;
count   = 1;

x_sil   = zeros(size(x));
y_sil   = zeros(size(y));

for j = 1:length(frames)
    if msk(j)
        jj_i            = frames(j):(frames(j)+N-1);
        jj_o            = frames(count):(frames(count)+N-1);
        x_sil(jj_o)     = x_sil(jj_o) + x(jj_i).*w;
        y_sil(jj_o)  	= y_sil(jj_o) + y(jj_i).*w;
        count           = count+1;
    end
end

x_sil = x_sil(1:jj_o(end));
y_sil = y_sil(1:jj_o(end));
end


