function [Tans] = mode1_mul()
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
                T1(j1,i2,i3)=sum;
            end
        end
    end
    
    for j2=1:l
		for i1=1:l
			for i3=1:l
				sum = 0;		
				for i2=1:l
					sum = sum + T1(i1,i2,i3) * A(j2,i2);		
                end
				T2(i1,j2,i3) =sum;
            end
        end
    end
    
    for j3=1:l
		for i1=1:l
			for i2=1:l
				sum = 0;		
				for i3=1:l
					sum = sum + T2(i1,i2,i3) * A(j3,i3);			
                end
				T3(i1,i2,j3)=sum;
            end
        end
    end
    
    Tttm = ttm(tensor(T),{A,A,A});
    Tans = Tttm-T3;
    %disp(T3);   
end