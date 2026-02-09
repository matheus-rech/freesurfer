function [fftaxis, deltafreq, indaxis, indaxisr] = fast_fftaxis(Ntp,TR)
% [fftaxis deltafreq indaxis indaxisr] = fast_fftaxis(Ntp,TR)
%
% Returns the frequencies at which the fft is computed, 
% from DC to the nyquist frequency. There will be Ntp/2+1
% frequencies. The last (Ntp/2-1) will be a replication of
% the 2:(Ntp/2) of the first section, but reversed.
%
% Example:
% fftaxis_pos = fast_fftaxis(Nf,TR);
%
% To get the negative frequencies:
% fftaxis_neg = conj(fliplr(fftaxis_pos(2:end-1)));
%
% EG, 
% h = randn(1,128);
% hfft = fft(h);
% [fftaxis deltafreq indaxis indaxisr] = fast_fftaxis(128,1);
% plot(fftaxis,abs(hfft(indaxis)))
% 
% If you only had the pos freq, eg,
% a = hfft(indaxis);
% b = conj(a(indaxisr))
% c = [a b];
% c will now be equal to hfft


%
% fast_fftaxis.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(nargin ~= 2) 
  msg = 'USAGE: [fftaxis deltafreq indaxis] = fast_fftaxis(Ntp,TR)';
  qoe(msg); error(msg);
end

Ntpdiv2 = round(Ntp/2);
nn = 0:Ntpdiv2;
freqmax = (1/TR)/2;         % Nyquist
deltafreq = freqmax/(Ntp/2); % Measured from 0 to Nyquist
fftaxis = deltafreq*nn;
indaxis = nn+1;
indaxisr = [Ntpdiv2:-1:2];
return;
