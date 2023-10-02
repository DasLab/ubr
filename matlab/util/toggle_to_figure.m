function toggle_to_figure( fignum );
% toggle_to_figure( fignum );
%
% Toggle to figure; if figure exists,
%  do *not* shift focus.
%
% (C) R. Das, Stanford & HHMI, 2023

if ~any(findobj(allchild(0), 'flat', 'type', 'figure')==fignum)
    figure(fignum);
end
set(groot,'CurrentFigure',fignum);


