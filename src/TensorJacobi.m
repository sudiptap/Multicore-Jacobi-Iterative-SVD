function [ACap,fnormarray,flag,max_root,coeff_arr] = TensorJacobi(A,I,R)
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
%lowerBound = 0; upperBound = (2/I);
%epsilon=lowerBound+rand(1,1)*(upperBound-lowerBound);

%Step 2
maxiter = 3000;
theta = 117;
fnormarray = zeros(1,maxiter);
max_root = zeros(1,maxiter);
coeff_arr = cell(1,maxiter);
fprintf('iteration starts after this\n');

flag = 0;
%rng(5);
for mi=1:maxiter
    fprintf('i = %d\n',mi);
    %step 3
    G = eye(I);
    %for m=1:R
        %for n=(R+1):I
     m = randi([1 R]);
        n = randi([R+1 I]);
            %m = m_values(mi);
            %n = n_values(mi);
            %disp(m);disp(n);
            %m = 1; n = R+1;
            %algorithm 2 to maximize theta
            sum_Tijm_Tijn = sum(sum(T(1:R-1,1:R-1,m) .* T(1:R-1,1:R-1,n)));
            sum_Tinn_Timn = sum(T(1:R-1,n,n) .* T(1:R-1,m,n));
            coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T(n,n,n)*T(m,n,n);
            
            
            sum_Tijm_Tijn = sum(sum(T(1:R-1,1:R-1,m) .* T(1:R-1,1:R-1,n)));            
            sum_Timm_Timn = sum(T(1:R-1,m,m) .* T(1:R-1,m,n));            
            coeff_t4 = 6*sum_Tijm_Tijn - 12*sum_Tijm_Tijn - 12*sum_Tinn_Timn - 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*T(m,n,n)*T(n,n,n) - 6*T(n,n,n)*T(m,m,m) - 18*T(m,m,n)*T(m,n,n) - 36*T(m,n,n)*T(m,m,n) + 18*T(m,n,n)*T(n,n,n);
            
            sum_Tinn_Tinn = sum(T(1:R-1,n,n) .* T(1:R-1,n,n));
            sum_Tinn_Timm = sum(T(1:R-1,n,n) .* T(1:R-1,m,m));
            sum_Timn_Timn = sum(T(1:R-1,m,n) .* T(1:R-1,m,n));
            coeff_t5 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm - 24*sum_Timn_Timn + 6*T(n,n,n)*T(n,n,n) -12*T(n,n,n)*T(m,m,n) - 18*T(m,n,n)*T(m,n,n);
            
            sum_Timm_Timm = sum(T(1:R-1,m,m) .* T(1:R-1,m,m));
            coeff_t3 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm -24*sum_Timn_Timn + 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn - 6*T(m,m,m)*T(m,n,n) + 6*T(n,n,n)*T(m,m,n) + 18*T(m,m,m)*T(n,n,n) - 36*T(m,m,n)*T(m,m,n) + 36*T(m,n,n)*T(m,n,n) - 18*T(m,n,n)*T(m,m,m) + 18*T(m,n,n)*T(m,m,m);
            
            coeff_t2 = 12*sum_Tijm_Tijn - 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*sum_Timm_Timn + 6*T(m,m,m)*T(n,n,n) - 12*T(m,m,m)*T(m,m,n) + 36*T(m,m,n)*T(m,n,n) - 18*T(m,n,n)*T(m,m,m) + 18*T(m,n,n)*T(m,m,n);
            
            coeff_t1 = 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn + 12*T(m,m,m)*T(m,n,n) - 6*T(m,m,m)*T(m,m,m) + 18*T(m,m,n)*T(m,m,n);
            
            coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T(m,m,m)*T(m,m,n);
            
            poly = [coeff_t6 coeff_t5 coeff_t4 coeff_t3 coeff_t2 coeff_t1 coeff_t0];
            coeff_arr{mi} = poly;
            %disp(poly);
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
            %disp(rt);
            %compute theta that maximizes objective function
            [size_rt_x,size_rt_y] = size(rt);
            %num_roots = size(rt);
            obj_func_array = zeros(1,size_rt_x);
            for rti=1:size_rt_x
                c = cos(atan(rt(rti)));
                s = sin(atan(rt(rti)));
                term1 = 3 .* sum(sum((c .* T(1:R-1,1:R-1,m) + s .* T(1:R-1,1:R-1,n)).^2));
                term2 = 3 .* sum(((c^2) .* T(1:R-1,m,m) + (s^2) .* T(1:R-1,n,n) + 2*c*s .* T(1:R-1,m,n)) .^ 2);
                term3 = ((c^3)*T(m,m,m) + (s^3)*T(n,n,n) + 3*(c^2)*s*T(m,m,n)+3*c*(s^2)*T(m,n,n)).^2;
                obj_func_array(rti) = term1 + term2 + term3;
            end
            
            [obj_max_val, obj_max_val_loc] = max(obj_func_array);
            theta = atan(rt(obj_max_val_loc));
            max_root(mi) = theta;
            %theta = 45 * theta;
           
            %epsilon1=lowerBound+rand(1,1)*(upperBound-lowerBound);
            %theta = theta - epsilon1;
            G(m,m) = cos(theta);
            G(n,n) = cos(theta);
            G(m,n) = -sin(theta);
            G(n,m) = sin(theta);
        %end
    %end
    %all threads should stop here - in parallel implementation
    %step 5
    Q = Q * G;
    
    %step 6
    GT = G';
    TT = tensor(T);
    TT = ttm(TT,{GT,GT,GT});
    %TT = GT * 
    T = double(TT);
   
    %convergence criteria
    U = Q(:,1:R);
   
    AT = tensor(A);
    ACap = ttm(AT,{U*U',U*U',U*U'});
    
    
    fn = norm(double(tenmat(ACap,1)),'fro');
    %disp(fn);

    %disp(sptensor(ACap));
    clear AT;
    FNormValue = norm(double(tenmat(A,1))-double(tenmat(ACap,1)),'fro')/(norm(double(tenmat(A,1)),'fro'));
    fnormarray(mi) = FNormValue;
    %disp(FNormValue);
    clear ACap;
    %fprintf('fnorm value : %f',FNormValue);
    
end
U = Q(:,1:R);
AT = tensor(A);
%ACap = ttm(AT,{U',U',U'});
%ACap = ttm(ACap,{U,U,U});

ACap = ttm(AT,{U*U',U*U',U*U'});
end


