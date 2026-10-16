function h = histogram(x, varargin)
%HISTOGRAM  Octave stand-in for MATLAB's histogram (ADDED 2026-09-17).
%   h = HISTOGRAM(x, ..., 'Normalization', 'probability') draws a bar chart of
%   the histogram of x(:) (NaN dropped) and returns the bar handle, so that
%   delete(h) works as in the MATLAB scripts. Supported options:
%   'Normalization' ('count' | 'probability' | 'pdf'), 'NumBins' (default: the
%   Freedman-Diaconis rule, at least 10, at most 60), 'BinEdges', 'FaceColor',
%   'FaceAlpha'. Anything else is ignored.
%   Lives in ReliabilityToolbox/octave_compat, added to the path ONLY under Octave.
x = x(:);  x = x(~isnan(x) & ~isinf(x));
if isempty(x)                          % nothing to draw (all-NaN column): invisible placeholder so delete(h) still works
    h = plot(NaN, NaN, 'HandleVisibility', 'off');
    return
end
norm_mode = 'count';  nbins = [];  edges = [];  face = [];  alpha = 0.6;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'normalization', norm_mode = lower(varargin{k+1});
        case 'numbins',       nbins = varargin{k+1};
        case 'binedges',      edges = varargin{k+1};
        case 'facecolor',     face  = varargin{k+1};
        case 'facealpha',     alpha = varargin{k+1};
    end
end
if isempty(edges)
    if isempty(nbins)
        n = numel(x);
        iqr_x = diff(quantile(x, [0.25 0.75]));
        if iqr_x > 0 && n > 1
            bw = 2 * iqr_x / n^(1/3);
            nbins = min(60, max(10, ceil((max(x) - min(x)) / bw)));
        else
            nbins = 10;
        end
    end
    lo = min(x); hi = max(x);
    if hi == lo, hi = lo + 1; end
    edges = linspace(lo, hi, nbins + 1);
end
counts  = histc(x, edges);  counts = counts(1:end-1);  counts(end) = counts(end) + sum(x == edges(end));
centers = (edges(1:end-1) + edges(2:end)) / 2;
widths  = diff(edges);
switch norm_mode
    case 'probability', vals = counts / max(sum(counts), 1);
    case 'pdf',         vals = counts ./ (max(sum(counts), 1) * widths(:));
    otherwise,          vals = counts;
end
h = bar(centers, vals(:)', 1.0);
if isempty(face)
    ax = gca; co = get(ax, 'ColorOrder'); idx = mod(numel(findobj(ax, 'type', 'hggroup')) + numel(findobj(ax, 'type', 'patch')) - 1, size(co, 1)) + 1;
    face = co(max(idx, 1), :);
end
set(h, 'FaceColor', face, 'EdgeColor', 'none');
try                                   % bar returns an hggroup without FaceAlpha; its patch children have it
    set(findobj(h, 'type', 'patch'), 'FaceAlpha', alpha);
catch
end
end
