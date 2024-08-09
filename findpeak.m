function [y,x,xindex] = findpeak(Y,X,varargin)
% 
% [y,x,xindex] = findpeak(Y,X)
% find peaks in selected interval in current plot, and
% return the max y value(s), corresponding x value(s).
% X must be a vector (horizontal axis). 
% If Y is a 2-D matrix, the function runs on the dimension with the same 
% length as X. 
% If Y is a complex vector (matrix), the function finds the maximum
% absolute value(s).
%
%  Mouse Operation:
%     select interval - click left mouse button above and below the region
%                       of interest
%     cancel input    - click middle button
%     done picking peaks - click right button
%
% Created by Shifei Yang, 2010.07.15
% Modified by MSA to be compatible with prior versions, 2010.07.19
% University of Wisconsin-Madison, msallen@engr.wisc.edu

if size(X,1) ~= 1 && size(X,2) ~= 1
    error('X must be a vector')
end
if size(Y,1) ~= length(X) && size(Y,2) == length(X)
    Y = Y.';
elseif size(Y,1) ~= length(X) && size(Y,2) ~= length(X)
    error('X and Y must have the same length')
end

if nargin > 2;
    npeaks = varargin{1}; % find exactly this many peaks.
    for k = 1:npeaks
        [x12,junk,mbutton] = ginput(2);
        x1 = x12(1); x2 = x12(2);
        xg = sort([x1,x2]);
        xind1 = find(X <= xg(1),1,'last');
        xind2 = find(X >= xg(2),1,'first');
        xind = (xind1:xind2).';

        maxY = zeros(size(Y,2),1);
        maxind = zeros(size(Y,2),1);
        for ii = 1:size(Y,2)
            [maxY(ii),maxind(ii)] = max(Y(xind,ii));
            maxind(ii) = maxind(ii)+xind1-1;
        end; clear xind
        [junk,psind] = max(maxY);
        indx = maxind(psind);
        x(k) = X(indx);
        y(k,:) = Y(indx,:);
        xindex(k) = indx;
        line(x(k),y(k,:),'Marker','*','Color','k');
    end
else
    
    disp('To select a peak click the LEFT button twice');
    disp('To cancel selection for the current peak click the MIDDLE mouse button.');
    disp('To quit click the RIGHT mouse button.');

    jj = 1;

    while true
        [x1,junk,mbutton] = ginput(1);
        if mbutton == 1
            [x2,junk,mbutton] = ginput(1);
            if mbutton == 2
                continue
            elseif mbutton == 3
                break
            else
                xg = sort([x1,x2]);
                xind1 = find(X <= xg(1),1,'last');
                xind2 = find(X >= xg(2),1,'first');
                xind = (xind1:xind2).';

                maxY = zeros(size(Y,2),1);
                maxind = zeros(size(Y,2),1);
                for ii = 1:size(Y,2)
                    [maxY(ii),maxind(ii)] = max(Y(xind,ii));
                    maxind(ii) = maxind(ii)+xind1-1;
                end; clear xind
                [junk,psind] = max(maxY);
                indx = maxind(psind);
                x(jj) = X(indx);
                y(jj,:) = Y(indx,:);
                xindex(jj) = indx;
                line(x(jj),abs(y(jj,:)),'Marker','*','Color','k');
                jj = jj + 1;
            end
        elseif mbutton == 2
            continue
        else
            break
        end
    end
end
[x,sortind] = sort(x);
x = x.';
y = y(sortind,:).';
xindex = xindex(sortind);

