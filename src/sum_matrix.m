function [sum] = sum_matrix(M, m)
    [row, col] = size(M);
    %disp(row);
    %disp(col);
    sum = 0;
    for r=1:row
        for c=1:col
            if r==m ||c==m
                %do nothing                
            else
                sum = sum + M(r,c);
            end
        end
    end    
end