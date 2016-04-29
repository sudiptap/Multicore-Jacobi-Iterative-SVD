function [T2] = mode1_mul()
    l = 5;
    T = rand(l,l,l);
    A = rand(l,l);
    for j1=1:l
        for i2=1:l
            for i3=1:l
                sum=0;
                for i1=1:l
                    sum = sum + T(i1,i2,i3)*A(j1,i1);
                end
                T2(j1,i2,i3)=sum;
            end
        end
    end       
end