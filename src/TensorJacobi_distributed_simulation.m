function [ACap,fnormarray,flag] = TensorJacobi_distributed_simulation(A,I,R)
%generate a random symmetric tensor
%I=10;
%M = rand(1,I);
%A = zeros(I,I,I);
%for i = 1:I
    %for j = 1:I
        %for k=1:I
            %A(i,j,k) = M(i)*M(j)*M(k);
        %end
    %end
%end

% Input part - put into argument of the function later
%R = 3;
Q = eye(I);
T = A;
%A = tensor(A);

%Step 1
lowerBound = 0; upperBound = (2/I);
epsilon=lowerBound+rand(1,1)*(upperBound-lowerBound);

%Step 2
maxiter = 14;
theta = 117;
fnormarray = zeros(1,100);
fprintf('iteration starts after this\n');

%this part is done for distributed implementation
G =  cell(1,maxiter);
%Q = cell(1,maxiter);
%T1 = cell(1,maxiter);
%U = cell(1,maxiter);

%for mi=1:maxiter
    %Q{mi} = eye(I);
    %T1{mi} = A;
%end

flag = 0;
for mmi = 1:300
    fprintf('i = %d\n',mmi);
    for mi=1:maxiter 
    %fprintf('i = %d\n',mi);
    %step 3
        G{mi} = eye(I); 
    %T = T1{mi};
    
    %for m=1:R
        %for n=(R+1):I
        m = randi([1 R]);
        n = randi([R+1 I]);
            %algorithm 2 to maximize theta
            T = double(T);
            sum_Tijm_Tijn = sum(sum(T(1:R-1,1:R-1,m) .* T(1:R-1,1:R-1,n)));
            sum_Tinn_Timn = sum(T(1:R-1,n,n) .* T(1:R-1,m,n));
            coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T(n,n,n)*T(m,n,n);
            %fprintf('%f\n',coeff_t6);
            
            sum_Tijm_Tijn = sum(sum(T(1:R-1,1:R-1,m) .* T(1:R-1,1:R-1,n)));            
            sum_Timm_Timn = sum(T(1:R-1,m,m) .* T(1:R-1,m,n));            
            coeff_t4 = 6*sum_Tijm_Tijn - 12*sum_Tijm_Tijn - 12*sum_Tinn_Timn - 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*T(m,n,n)*T(n,n,n) - 6*T(n,n,n)*T(m,m,m) - 18*T(m,m,n)*T(m,n,n) - 36*T(m,n,n)*T(m,m,n) + 18*T(m,n,n)*T(n,n,n);
            %fprintf('%f\n',coeff_t4);
            
            sum_Tinn_Tinn = sum(T(1:R-1,n,n) .* T(1:R-1,n,n));
            sum_Tinn_Timm = sum(T(1:R-1,n,n) .* T(1:R-1,m,m));
            sum_Timn_Timn = sum(T(1:R-1,m,n) .* T(1:R-1,m,n));
            coeff_t5 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm - 24*sum_Timn_Timn + 6*T(n,n,n)*T(n,n,n) -12*T(n,n,n)*T(m,m,n) - 18*T(m,n,n)*T(m,n,n);
            %fprintf('%f\n',coeff_t5);
            
            sum_Timm_Timm = sum(T(1:R-1,m,m) .* T(1:R-1,m,m));
            coeff_t3 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm -24*sum_Timn_Timn + 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn - 6*T(m,m,m)*T(m,n,n) + 6*T(n,n,n)*T(m,m,n) + 18*T(m,m,m)*T(n,n,n) - 36*T(m,m,n)*T(m,m,n) + 36*T(m,n,n)*T(m,n,n) - 18*T(m,n,n)*T(m,m,m) + 18*T(m,n,n)*T(m,m,m);
            %fprintf('%f\n',coeff_t3);
            
            coeff_t2 = 12*sum_Tijm_Tijn - 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*sum_Timm_Timn + 6*T(m,m,m)*T(n,n,n) - 12*T(m,m,m)*T(m,m,n) + 36*T(m,m,n)*T(m,n,n) - 18*T(m,n,n)*T(m,m,m) + 18*T(m,n,n)*T(m,m,n);
            %fprintf('%f\n',coeff_t2);
            
            coeff_t1 = 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn + 12*T(m,m,m)*T(m,n,n) - 6*T(m,m,m)*T(m,m,m) + 18*T(m,m,n)*T(m,m,n);
            %fprintf('%f\n',coeff_t1);
            
            coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T(m,m,m)*T(m,m,n);
            %fprintf('%f\n',coeff_t0);
            
            poly = [coeff_t6 coeff_t5 coeff_t4 coeff_t3 coeff_t2 coeff_t1 coeff_t0];
            rt = roots(poly);
            %disp(rt);
            rt = rt(imag(rt)==0);
            if isempty(rt) == 1
                fprintf('No real roots');
                flag = 1;
            end
       %end
    %end
            if flag == 1
                break;
            end
            
            [size_rt_x,size_rt_y] = size(rt);
            %num_roots = size(rt);
            obj_func_array = zeros(1,size_rt_x);
            for rti=1:size_rt_x
                c = cosd(rt(rti));
                s = sind(rt(rti));
                term1 = 3 .* sum(sum((c .* T(1:R-1,1:R-1,m) + s .* T(1:R-1,1:R-1,n)).^2));
                term2 = 3 .* sum(((c^2) .* T(1:R-1,m,m) + (s^2) .* T(1:R-1,n,n) + 2*c*s .* T(1:R-1,m,n)) .^ 2);
                term3 = ((c^3)*T(m,m,m) + (s^3)*T(n,n,n) + 3*(c^2)*s*T(m,m,n)+3*c*(s^2)*T(m,n,n)).^2;
                obj_func_array(rti) = term1 + term2 + term3;
            end
            
            [obj_max_val, obj_max_val_loc] = max(obj_func_array);
            theta = rt(obj_max_val_loc);
            %max_root(mi) = theta;
            
            
            %theta = max(rt);
            %theta = 8 * theta;
            %fprintf('%f\n',theta);
            %epsilon1=lowerBound+rand(1,1)*(upperBound-lowerBound);
            %theta = theta - epsilon1;
            %theta = 112;
            G{mi}(m,m) = cosd(theta);
            G{mi}(n,n) = cosd(theta);
            G{mi}(m,n) = -sind(theta);
            G{mi}(n,m) = sind(theta);
        %end
    %end   
    
end

    %ACap = ttm(tensor(A),{U{1}*U{1}',U{1}*U{1}',U{1}*U{1}'});
    
    %try this
    
    for mi=1:maxiter
        Q = Q * G{mi};
    end
    for mi=1:maxiter
        %TT = ttm(tensor(T),{U{mi}*U{mi}',U{mi}*U{mi}',U{mi}*U{mi}'});
        GT = G{mi}';
        TT = tensor(T);
        TT = ttm(TT,{GT,GT,GT});
        %TT = GT * 
        T = double(TT);
    end
    U = Q(:,1:R);
    AT = tensor(A);
    ACap = ttm(AT,{U*U',U*U',U*U'});
    
    %for mi=1:maxiter
        %Q = Q * G{mi};
    %end
    %for mi=1:maxiter        
        %TT = ttm(tensor(T),{G{mi}',G{mi}',G{mi}'});
        %T = tensor(TT);
    %end
    %U = Q(:,1:R);
    %ACap = ttm(T,{U*U',U*U',U*U'});
    
    FNormValue = norm(double(tenmat(A,1))-double(tenmat(ACap,1)),'fro')/(norm(double(tenmat(A,1)),'fro'));
    fnormarray(mmi) = FNormValue;
end
%U = Q(:,1:R);
AT = tensor(A);
%ACap = ttm(AT,{U',U',U'});
%ACap = ttm(ACap,{U,U,U});
for idx=1:299
    if(fnormarray(idx)<fnormarray(idx+1))
        fprintf('%d is greater than %d',idx+1,idx);
    %else
        %disp('no');
    end
end

end
