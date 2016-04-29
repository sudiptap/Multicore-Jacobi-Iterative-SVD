function [ACap,fnormarray,flag] = TensorJacobi_distributed_simulation1(A,I,R)
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
maxiter = 4;
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
for mmi = 1:2000
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
            
            poly = getPolynomialOne(T,m,n,R);
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
                c = cos(atan(rt(rti)));
                s = sin(atan(rt(rti)));
                term1 = 3 .* sum_matrix(((c .* T(1:R-1,1:R-1,m) + s .* T(1:R-1,1:R-1,n)).^2),m);
                term2 = 3 .* sum_vector((((c^2) .* T(1:R-1,m,m) + (s^2) .* T(1:R-1,n,n) + 2*c*s .* T(1:R-1,m,n)) .^ 2),m);
                term3 = ((c^3)*T(m,m,m) + (s^3)*T(n,n,n) + 3*(c^2)*s*T(m,m,n)+3*c*(s^2)*T(m,n,n)).^2;
                obj_func_array(rti) = term1 + term2 + term3;
            end
            
            [obj_max_val, obj_max_val_loc] = max(obj_func_array);
            theta = atan(rt(obj_max_val_loc));
            %max_root(mi) = theta;
            
            
            %theta = max(rt);
            theta = 1 * theta;
            %fprintf('%f\n',theta);
            %epsilon1=lowerBound+rand(1,1)*(upperBound-lowerBound);
            %theta = theta - epsilon1;
            %theta = 112;
            G{mi}(m,m) = cos(theta);
            G{mi}(n,n) = cos(theta);
            G{mi}(m,n) = -sin(theta);
            G{mi}(n,m) = sin(theta);
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
%for idx=1:299
    %if(fnormarray(idx)<fnormarray(idx+1))
        %fprintf('%d is greater than %d',idx+1,idx);
    %else
        %disp('no');
    %end
%end

end

function poly = getPolynomialOne(T,m,n,R)

sum_Tijn_Tijn = sum_matrix((T(1:R,1:R,n) .* T(1:R,1:R,n)),m);
sum_Tijm_Tijm = sum_matrix((T(1:R,1:R,m) .* T(1:R,1:R,m)),m);
sum_Tijm_Tijn = sum_matrix((T(1:R,1:R,m) .* T(1:R,1:R,n)),m);

%term1
p11 = [1 0 2 0 1];
p12 = [(-6*sum_Tijm_Tijn) (6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm) (6*sum_Tijm_Tijn)];
%p12 = [(-6*sum_Tijm*sum_Tijn) (6*sum_Tijn*sum_Tijn - 6*sum_Tijm*sum_Tijm) (6*sum_Tijm*sum_Tijn)];
p1 = conv(p11,p12);

%term2
p22 = zeros(1,5);
for i=1:R
    if (i ~= m)
        t1 = [6*T(i,n,n) 12*T(i,m,n) 6*T(i,m,m)];
        t2 = [-2*T(i,m,n) -2*T(i,m,m) + 2*T(i,n,n) 2*T(i,m,n)];
        %disp(conv(t1,t2));
        %disp(p22);
        p22 = p22 + conv(t1, t2);
    end
end

p21 = [1 0 1];
p2 = conv(p21,p22);

%term3
p31 = [2*T(n,n,n) 6*T(m,n,n) 6*T(m,m,n) 2*T(m,m,m)];
p32 = [(-3*T(m,n,n)) (3*T(n,n,n) - 6*T(m,m,n)) (-3*T(m,m,m) + 6*T(m,n,n)) (3*T(m,m,n))];
p3 = conv(p31, p32);


poly = p1 + p2 + p3;
end
