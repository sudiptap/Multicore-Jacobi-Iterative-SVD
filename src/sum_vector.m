function [sum] = sum_vector(V, m)
    [row, col] = size(V);
    %disp(row);
    %disp(col);
    sum = 0;    
    for c=1:row
        if c==m
            %do nothing            
        else
            sum = sum + V(c);
        end
    end        
end