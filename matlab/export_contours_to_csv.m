function export_contours_to_csv(X, Y, Z, levels, stem)

C = contourc(X(1,:), Y(:,1), Z, levels);

rows = zeros(0,3);   % force 3 columns from the start
i = 1;

while i < size(C,2)
    level = C(1,i);
    npts  = C(2,i);

    if npts <= 0 || (i + npts) > size(C,2)
        break
    end

    pts = C(:, i+1:i+npts).';   % npts x 2

    if size(pts,2) ~= 2
        i = i + npts + 1;
        continue
    end

    seg = [pts, repmat(level, npts, 1)];  % npts x 3
    rows = [rows; seg; nan(1,3)]; %#ok<AGROW>

    i = i + npts + 1;
end

if isempty(rows)
    warning('No contours found for %s', stem);
    rows = nan(1,3);
end

Tout = table(rows(:,1), rows(:,2), rows(:,3), ...
    'VariableNames', {'Vt_mps','Altitude_m','Level'});

writetable(Tout, [stem '.csv']);
fprintf('Saved %s.csv\n', stem);

end