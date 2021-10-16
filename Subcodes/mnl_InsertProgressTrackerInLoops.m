function mnl_InsertProgressTrackerInLoops(i,n)
%Inputs
%   i=the current number
%   n=the total number of loops to do
persistent i_prev_prct;
PcValue=(i/n)*100;
rPcValue=round(PcValue,1);
if rPcValue==100 && i<n
    PcValue=99.9;
else
    PcValue=rPcValue;
end
msg=sprintf('%s\n','Percent Complete...');
%% Print message when counting starts.
if isempty(i_prev_prct) || PcValue < i_prev_prct
    i_prev_prct = 0;
    S_prev = getPrctStr(i_prev_prct);    
    fprintf('%s%s\n',msg, S_prev);
end
%% Print updated percentage.
if PcValue~=i_prev_prct
    S_prev = getPrctStr(i_prev_prct);
    fprintf(getBackspaceStr(numel(S_prev)));
    S = getPrctStr(PcValue);
    fprintf('%s', S);    
    i_prev_prct=PcValue;
end
%% Clear percentage variable.
if PcValue == 100
    fprintf(' Done.\n');
    clear i_prev_prct;
end
end
function S = getPrctStr(prct)
S = sprintf('%d%%  %s',prct,getDotStr(prct));
if prct < 10
    S = ['  ',S];
elseif prct < 100
    S = [' ',S];
end
end
function S = getDotStr(prct)
S = repmat(' ',1,10);
S(1:floor(prct/10)) = '.';
end
function S = getBackspaceStr(N)
S = repmat('\b',1,N);
end