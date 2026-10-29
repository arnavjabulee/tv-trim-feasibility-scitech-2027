function split_contour_levels(fname, prefix)
T = readtable(fname);

V = T{:,'Vt_mps'};
H = T{:,'Altitude_m'};
L = T{:,'Level'};

levels = unique(L(~isnan(L)));

for k = 1:numel(levels)
    lev = levels(k);
    idx = (L == lev) | isnan(L);

    Tk = T(idx,:);

    % trim leading/trailing NaN separator rows
    while ~isempty(Tk) && all(ismissing_or_nan_row(Tk(1,:)))
        Tk(1,:) = [];
    end
    while ~isempty(Tk) && all(ismissing_or_nan_row(Tk(end,:)))
        Tk(end,:) = [];
    end

    if isempty(Tk)
        continue
    end

    if lev < 0
        out = sprintf('%s_m%g.csv', prefix, abs(lev));
    else
        out = sprintf('%s_%g.csv', prefix, lev);
    end


    writetable(Tk, out);
    fprintf('Saved %s\n', out);
end
end

function tf = ismissing_or_nan_row(row)
vals = row{1,:};
tf = all(isnan(vals));
end