function x = loc(mstring,sstring)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Returns number of the row of mstring that has the same 
% non-blanck characters as sstring. 
% *Note: Strings must be placed in between single quotation marks.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Info
[rm,cm] = size(mstring);
cs      = max(size(sstring));

% If necessary, add blanck columns to sstring so it will have the
%  same number of columns as mstring.
if cm > cs
    nblancks = cm - cs;
    for i = 1:nblancks
        sstring = [sstring,' '];
    end
end

% Adding spaces.
if(cm ~= max(size(sstring)))
    disp(['problem with padding ',sstring])
    disp('The character string might be longer than name list')
    mstring  
    pause
end

% Finding value.
x = [];
for r = 1:rm
    if(length(find(mstring(r,:)==sstring))==cm)
        x = r;
    end
end

% No match.
if x == 0
    if(~exist('switchmod'))
        disp(['Could not find ',sstring]); 
    end
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%