%% Gibbs simplex constraint
function [phia,phib,phic] = calc_gibbs_simplex(phia,phib,phic)
    [row,col] = find(phia<0 | phib<0 | phic<0);
    violations = row+size(phia,1)*(col-1);

    while length(violations)>1
        for k=violations'
            if phia(k)<0.0
                % Check all dual combinations
                if phib(k)<=0.0
                    phia(k) = 0.0;
                    phib(k) = 0.0;
                    phic(k) = 1.0;
                elseif phic(k)<=0.0
                    phia(k) = 0.0;
                    phib(k) = 1.0;
                    phic(k) = 0.0;
                else
                    phib(k) = phib(k) + 0.5*phia(k);
                    phic(k) = phic(k) + 0.5*phia(k);
                    phia(k) = 0.0;
                end
            elseif phib(k)<0
                % Check remaining combination with c
                if phic(k)<=0
                    phia(k) = 1.0;
                    phib(k) = 0.0;
                    phic(k) = 0.0;
                else
                    phia(k) = phia(k) + 0.5*phib(k);
                    phic(k) = phic(k) + 0.5*phib(k);
                    phib(k) = 0.0;
                end  
            elseif phic(k)<0
                phia(k) = phia(k) + 0.5*phic(k);
                phib(k) = phib(k) + 0.5*phic(k);
                phic(k) = 0.0;
            end
        end
        [row,col] = find(phia<0 | phib<0 | phic<0);
        violations = row+size(phia,1)*(col-1);
end