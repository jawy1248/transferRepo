function [xt,ts] = ifftp(Xfft,ws,varargin)
% [xt,ts] = ifftp(Xfft,ws,DIM);
% 
% This m-file is the inverse of fftp.  It takes the inverse FFT of the
% complex FFT "Xfft" (produced by fftp) and returns the time record and
% corresponding time vector.  This is accomplished by mirroring the complex
% conjugates of the coefficients in Xfft and then using Matlab's 'ifft.m'.
%   Xfft - FFT or array of FFTs.  If Xfft is a matrix,
%        the FFT operates on the larger dimension of Xfft, unless the
%        optional argument DIM is given.
%   DIM - (optional) the dimension over which the FFT acts
%   ws - frequency vector
% 
% If no output arguments are specified, a plot of the time history(ies) is
% created in a new figure window.
%
% This routine assumes that Xfft has not been scaled after applying
% 'fftp'.
%
% Matt Allen Fall 2018
% msallen@engr.wisc.edu
%

w1 = ws(2)-ws(1);

if nargin > 2
    DIM = varargin{1};
    if DIM == 2;
        Xfft = Xfft.';
    end
else
    [a,b] = size(Xfft);
    if b > a;
        Xfft = Xfft.';
    end
end

% The time sequences are now the columns of xt.
L = size(Xfft,1);
if length(ws) ~= L;
    warning('Frequency Vector should be the same length as FFT sequence');
end

N = 2*(L-1);
xt = ifft([Xfft; conj(Xfft([(end-1):-1:2],:))],[],1);
    %norm(imag(xt))/norm(real(xt))
xt = real(xt);

T = 2*pi/w1;
ts = [0:T/N:T-T/N];

if nargout < 1;
    figure
    plot(ts,xt); grid on;
    title('IFFT(x)');
    xlabel('Time (s)'); ylabel('x(t)');
end
    