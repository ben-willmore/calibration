function spect=golAnal(golA, golB, golayA, golayB);
% function spect=golAnal(golA, golB, golayA, golayB)
%	Golay code frequency analysis.
% 	Calculates the complex frequency spectrum from 
%	a pair of golay code signals (golA and golB) 
%	using the original codes golayA and golayB. 
%	If golA and golB are longer than golayA and golayB the
%	original golay codes will be zero padded.
%
%	see Zhou, Green & Middlebrooks (1992) J.Acoust.Soc.Am. 92:1169-71 
%	for theory
%  
%	see also function [a b]=golay(X) for generating golay codes.

if sum(size(golA)==size(golayA'))==2 
    golA=golA';
end;

if sum(size(golB)==size(golayB'))==2 
    golB=golB';
end;

ldiff = length(golA) - length(golB); 
if ldiff ~= 0
  error(['codes to analyse must be of equal length ' num2str(ldiff)]);
end;

L=length(golayA);

% make sure that input codes and complementary codes are of equal length
ldiff = length(golA) - L; 
if ldiff<0
  golA=[golA zeros(1,-ldiff)];
  golB=[golB zeros(1,-ldiff)];
end;

if ldiff>0
  golayA=[golayA zeros(1,ldiff)];
  golayB=[golayB zeros(1,ldiff)];
end;

% calculate spectrum
spect = ((fft(golA).*conj(fft(golayA)))+(fft(golB).*conj(fft(golayB))))/(2*L);

% get rid of mirrored part of spectrum
% spect = spect(1:length(spect)/2);