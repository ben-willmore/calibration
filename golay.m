function [ga, gb] = golay(x);
% function [ga, gb] = golay(x);
%	generates a pair of golay codes [ga, gb]
%	of length 2^x
%
%	see Zhou, Green & Middlebrooks (1992) J.Acoust.Soc.Am. 92:1169-71 
%	for theory
%

ga = [1 1];
gb = [1 -1];

for idx=1:x-1,
  ha=[ga gb];
  hb=[ga -gb];
  ga=ha;
  gb=hb;
end;

