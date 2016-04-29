
function [ACap,fnormarray,flag,max_root,coeff_arr] = compute_poly_coefficients(A,I,R)
%disp(A(13,24,16));
G = rand(I,I);
rnd = 1;
for r=1:I
    for c=1:I        
        rnd = r/I * c/I;        
        G(r,c)=rnd;       				
    end
end
%disp(sparse(G));
T3 = ttm(tensor(A),{G,G,G});
for ti=1:I
    for tj=1:I
        for tk=1:I
            fprintf('%.29f\n', T3(ti,tj,tk));
        end
    end
end       
return;
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
maxiter = 1000;
theta = 117;
fnormarray = zeros(1,maxiter);
max_root = zeros(1,maxiter);
coeff_arr = cell(1,maxiter);
fprintf('iteration starts after this\n');
m_values = [10 3 10 4 2 10 2 5 9 6 4 2 7 3 7 6 5 7 7 4 5 3 5 5 4 8 7 9 4 5 3 7 10 7 5 6 5 8 8 8 3 2 7 6 5 1 2 8 2 10 8 8 7 7 10 9 3 4 1 9 1 7 9 7 2 10 5 2 4 5 5 8 4 8 10 3 8 6 5 9 10 6 9 4 9 7 4 4 7 5 9 10 8 5 1 1 3 5 6 5 10 3 8 6 9 3 10 7 1 2 6 2 9 1 5 3 9 3 5 3 1 3 10 9 4 2 4 1 10 10 10 4 1 7 5 1 7 7 6 10 8 9 10 7 3 7 4 2 10 2 10 1 6 8 5 1 3 7 1 3 3 6 3 1 9 9 5 10 10 3 5 10 4 1 1 10 4 2 5 2 7 5 3 5 3 9 3 9 7 4 4 4 1 8 9 1 9 10 4 4 4 7 3 6 8 7 5 1 9 1 7 5 10 10 9 1 8 10 6 10 6 5 10 6 4 8 9 10 8 3 4 10 3 2 7 2 1 4 2 1 7 8 1 5 5 6 3 1 7 7 9 7 8 2 2 8 3 5 3 3 1 10 6 1 8 9 1 7 7 10 6 8 6 4 4 10 8 8 2 1 9 6 5 8 4 1 8 10 3 4 2 3 3 8 2 5 8 2 8 5 9 2 7 2 7 9 9 1 3 8 7 7 8 8 10 8 7 9 4 5 6 2 6 10 4 6 3 8 4 7 7 4 5 10 3 9 1 5 7 8 4 4 6 1 8 4 1 1 2 4 10 5 6 9 6 6 10 8 4 7 6 7 1 2 3 10 1 3 5 4 9 4 1 4 10 8 3 3 5 9 1 10 3 2 4 3 5 2 6 2 10 2 4 8 9 8 5 5 2 9 3 10 7 7 10 1 10 2 9 7 8 8 3 2 1 1 1 4 5 2 1 3 8 7 5 3 8 5 7 8 6 3 4 5 10 3 2 4 3 6 6 1 5 7 3 9 6 7 9 8 3 1 9 6 8 9 4 8 8 10 2 1 10 9 4 1 10 2 8 8 3 2 9 5 7 7 5 9 9 6 5 1 8 3 3 4 10 2 6 5 3 2 3 3 10 5 6 2 1 2 8 10 2 8 1 1 6 10 2 2 1 9 5 3 5 10 3 10 1 5 10 6 7 4 10 3 4 10 2 5 9 8 4 10 6 9 1 4 2 8 6 2 4 1 6 3 10 10 10 2 4 4 5 2 7 8 3 3 2 5 9 4 8 4 3 4 7 2 7 1 6 6 10 9 3 10 2 1 7 10 9 9 4 1 6 4 9 2 10 1 6 5 5 10 10 4 4 8 5 10 10 3 7 10 7 2 4 3 4 10 5 5 10 9 3 6 4 5 6 8 10 8 8 10 6 5 8 4 3 3 4 2 9 1 3 10 10 4 9 7 8 8 2 1 5 4 4 8 2 10 7 10 6 2 10 2 3 1 4 2 8 9 1 5 4 10 5 6 3 8 9 10 4 9 5 7 9 6 2 7 9 7 6 7 2 4 6 7 5 7 8 4 10 1 3 10 4 6 8 8 7 5 4 3 7 10 6 4 5 2 2 10 6 3 10 8 5 2 2 9 5 4 4 8 4 9 1 9 9 4 6 6 3 9 3 4 8 8 7 3 8 4 3 6 8 10 2 5 6 9 4 6 2 6 1 9 10 10 7 6 6 1 3 10 3 7 6 9 8 7 3 10 1 8 6 5 1 9 5 5 9 3 7 3 8 5 7 5 5 6 7 4 8 3 1 10 2 5 6 3 1 4 2 6 5 1 4 10 5 4 3 6 9 2 2 9 4 10 6 5 7 3 1 4 8 4 2 5 2 7 4 8 1 5 4 8 10 4 3 10 6 1 4 10 6 4 3 8 8 1 7 6 9 10 8 1 2 4 8 3 2 4 9 9 9 9 10 4 5 8 4 7 3 3 6 5 5 2 4 9 4 10 5 2 1 4 1 1 3 5 10 6 8 7 1 5 8 9 1 4 3 5 9 8 9 10 2 5 6 2 10 1 3 10 1 8 9 8 10 10 9 9 6 1 5 8 1 2 9 1 5 6 2 3 2 1 5 3 1 5 9 1 5 6 3 6 6 3 7 1 10 3 10 9 9 8 7 8 10 6 3 4 7 8 4 5 1 1 3 5 2 2 1 1 6 7 2 9 7 2 5 6 7 10 5 7 7 8 1 1 7 9 8 7 4 7 8 7 1 6];
n_values = [13 29 18 27 13 24 30 28 15 21 17 11 24 21 12 16 18 16 20 28 26 16 18 15 21 19 29 15 22 30 11 19 23 17 30 21 29 12 13 23 27 11 12 20 30 20 28 22 26 18 27 24 16 14 15 22 20 30 19 26 20 14 16 22 26 19 19 21 11 15 25 27 22 16 27 22 19 18 12 26 28 24 19 22 20 15 14 29 11 29 19 28 27 14 14 30 16 21 30 17 13 15 18 15 22 19 14 19 13 21 12 21 16 27 17 16 27 29 28 15 27 30 21 22 22 21 25 14 12 27 24 19 26 17 11 15 13 18 27 29 23 23 20 21 18 12 23 26 30 15 22 18 19 21 29 15 20 12 15 23 21 16 30 13 24 21 21 12 27 26 15 20 27 26 13 15 16 18 24 25 30 13 17 12 19 20 19 29 19 19 24 19 25 17 20 17 28 21 24 28 23 16 27 19 20 11 12 15 18 19 13 28 14 13 14 12 29 22 25 23 18 20 15 20 26 11 22 20 19 26 15 21 11 20 23 13 28 12 20 16 18 25 17 28 19 29 17 17 25 23 20 21 21 11 24 18 17 26 21 30 20 25 12 24 19 17 15 28 11 28 30 29 14 29 28 14 29 25 20 30 29 19 24 12 29 30 30 14 25 18 13 21 22 26 28 22 22 12 11 21 16 27 13 30 15 11 12 13 30 16 15 14 30 15 22 21 11 16 20 22 25 26 25 19 29 13 13 21 25 24 24 25 24 28 26 30 23 28 26 22 14 19 22 29 27 26 29 15 12 16 24 22 21 14 17 30 11 17 28 12 11 26 29 13 27 12 15 18 11 23 21 18 21 25 13 15 26 15 14 17 29 13 14 19 15 13 26 28 18 27 29 14 24 30 17 11 19 29 27 26 22 30 21 15 30 29 24 25 20 18 14 22 13 26 30 25 22 28 12 14 24 12 26 16 13 22 23 12 26 19 30 28 27 28 18 13 17 11 22 17 21 14 27 11 16 30 21 11 14 19 12 26 16 21 24 13 21 28 17 22 20 19 29 21 28 27 23 18 30 17 28 26 29 26 16 19 19 24 28 23 16 16 14 18 12 25 25 30 28 12 23 27 21 13 28 24 17 15 14 23 28 19 11 27 12 23 22 28 14 23 17 18 20 25 26 12 29 27 27 11 25 25 17 30 26 12 21 20 24 26 26 24 29 26 24 27 24 17 29 11 12 21 21 25 15 15 30 15 29 27 30 14 27 26 22 23 20 19 13 18 27 28 23 28 16 28 29 24 25 24 19 21 30 11 15 11 13 11 14 11 14 28 22 11 21 22 21 18 22 14 30 19 17 23 21 21 14 30 20 22 15 21 26 26 12 25 18 18 14 13 15 30 14 24 20 15 30 30 26 24 21 13 13 21 11 22 20 27 12 16 30 19 18 16 18 22 17 17 25 17 23 15 11 17 25 28 15 23 28 16 19 16 13 26 30 27 22 16 23 20 17 23 27 24 17 21 19 12 30 14 21 13 22 28 22 14 16 11 23 25 16 23 23 30 14 12 13 30 15 21 25 28 11 23 27 21 26 29 12 21 11 12 11 20 11 14 15 17 14 20 21 15 17 16 22 18 22 30 24 18 30 12 19 30 29 25 11 17 23 29 21 30 17 16 29 17 19 23 17 30 11 28 14 22 29 25 17 14 15 28 24 11 25 28 28 20 24 25 25 26 17 17 27 20 20 26 24 27 24 11 17 14 19 27 13 11 17 16 30 17 24 14 13 25 21 20 26 27 13 22 12 28 12 12 19 11 24 19 19 13 27 17 11 18 23 21 24 11 24 25 15 14 22 24 16 12 14 28 25 20 18 22 22 18 27 30 11 22 29 25 15 28 18 26 21 23 20 21 18 20 13 19 26 17 27 11 11 15 14 12 28 18 19 21 29 14 11 27 14 23 26 14 13 26 11 27 24 18 24 12 20 15 11 20 18 17 30 13 12 22 29 13 26 21 22 14 11 13 23 18 29 25 29 21 17 24 16 12 20 11 17 24 27 21 22 19 17 29 21 16 18 17 29 16 13 23 12 18 28 29 30 21 23 13 24 14 21 16 11 11 18 14 25 21 14 11 29 23 13 16 25 19 13 15 18 28 13 19 17 11 22 27 25 23 13 19 28 11 25 25 13 20 26 15 30 25 27 23 29 14 12 20 11 13 18 20 30 28 26 22 29 28 27 12 24 18 30 27 30 18 17 18 17 26 28 30 23 25 28 25];
%size(m_values)
%size(m_values)

flag = 0;
%rng(5);
for mi=1:maxiter
    
    fprintf('i = %d\n',mi);
    %step 3
    G = eye(I);
    %for m=1:R
        %for n=(R+1):I
     %m = randi([1 R]);
        %n = randi([R+1 I]);
        m = m_values(mi);
        n = n_values(mi);
        %disp(m);
        %disp(n);
            %m = m_values(mi);
            %n = n_values(mi);
            %disp(m);disp(n);
            %m = 1; n = R+1;
            %algorithm 2 to maximize theta
            poly = getPolynomialTwo(T,m,n,R);
            
            
            
            %poly2 = getPolynomialTwo(T,m,n,R);
            %disp(poly1);
            %disp(poly2);
            %return;
          
            
            coeff_arr{mi} = poly;
            %disp(poly);
            rt = roots(poly);
            %disp(rt);
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
            %rt1 = zeros(size(rt));
            %disp(rt);
            for rti=1:size_rt_x
                rt(rti) = atan(rt(rti)); 
            end
            
            %num_roots = size(rt);
            obj_func_array = zeros(1,size_rt_x);
            for rti=1:size_rt_x
                %c = cos(atan(rt(rti)));
                %s = sin(atan(rt(rti)));
                c = cos(rt(rti));
                s = sin(rt(rti));
                term1 = 3 .* sum_matrix(((c .* T(1:R,1:R,m) + s .* T(1:R,1:R,n)).^2),m);
                term2 = 3 .* sum_vector((((c^2) .* T(1:R,m,m) + (s^2) .* T(1:R,n,n) + 2*c*s .* T(1:R,m,n)) .^ 2),m);
                term3 = ((c^3)*T(m,m,m) + (s^3)*T(n,n,n) + 3*(c^2)*s*T(m,m,n)+3*c*(s^2)*T(m,n,n)).^2;
                obj_func_array(rti) = term1 + term2 + term3;
            end
            
            [obj_max_val, obj_max_val_loc] = max(obj_func_array);
            %theta = atan(rt(obj_max_val_loc));
            theta = rt(obj_max_val_loc);
            %fprintf('best_root : %f',theta);
            %disp(theta);
            max_root(mi) = theta;
            %theta = 2 * theta;
           
            %epsilon1=lowerBound+rand(1,1)*(upperBound-lowerBound);
            %theta = theta - epsilon1;
            %disp('theta values');
            %disp(cos(theta));
            %disp(sin(theta))
            %disp('theta values');
            G(m,m) = cos(theta);
            G(n,n) = cos(theta);
            G(m,n) = -sin(theta);
            G(n,m) = sin(theta);
        %end
    %end
    %all threads should stop here - in parallel implementation
    %step 5
    %tic;
    Q = Q * G;
    %toc;
    %step 6
    GT = G';
    TT = tensor(T);
    TT = ttm(TT,{GT,GT,GT});
    %TT = GT * 
    T = double(TT);
    if mi==4
        for ti=1:I
            for tj=1:I
                for tk=1:I
                    fprintf('%.30f\n', T(ti,tj,tk));
                end
            end
        end
        return;
    end
   
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

function poly = getPolynomialTwo(T,m,n,R)

sum_Tijn_Tijn = sum_matrix((T(1:R,1:R,n) .* T(1:R,1:R,n)),m);
sum_Tijm_Tijm = sum_matrix((T(1:R,1:R,m) .* T(1:R,1:R,m)),m);
sum_Tijm_Tijn = sum_matrix((T(1:R,1:R,m) .* T(1:R,1:R,n)),m);
%temp = (T(1:R,1:R,m) .* T(1:R,1:R,n));disp(sparse(T(1:R,1:R,n)));%return;
sum_Tinn_Timn = sum_matrix((T(1:R,n,n) .* T(1:R,m,n)),m);
sum_Timn_Tinn = sum_Tinn_Timn;
sum_Tinn_Tinn = sum_vector((T(1:R,n,n) .* T(1:R,n,n)),m);
sum_Tinn_Timm = sum_vector((T(1:R,n,n) .* T(1:R,m,m)),m);
sum_Timm_Tinn = sum_Tinn_Timm;
sum_Timn_Timn = sum_vector((T(1:R,m,n) .* T(1:R,m,n)),m);
sum_Timm_Timm = sum_vector((T(1:R,m,m) .* T(1:R,m,m)),m);
sum_Timm_Timn = sum_vector((T(1:R,m,m) .* T(1:R,m,n)),m);

%disp('-----------coefficients makers-----------');
%disp(sum_Tijm_Tijn);
% disp(sum_Tinn_Timn);
% disp(sum_Tinn_Tinn);
% disp(sum_Tinn_Timm);
% disp(sum_Timn_Timn);
% disp(sum_Tijn_Tijn);
% disp(sum_Tijm_Tijm);
% disp(sum_Timm_Timn);
% disp(sum_Timm_Timm);



coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T(n,n,n)*T(m,n,n);

%coeff_t5 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm - 24*sum_Timn_Timn + 6*T(n,n,n)*T(n,n,n) -12*T(n,n,n)*T(m,m,n) - 18*T(m,n,n)*T(m,n,n) + 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm;
coeff_t5 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm -24*sum_Timn_Timn -12*sum_Timm_Tinn + 12*sum_Tinn_Tinn -18*T(m,n,n)*T(m,n,n) + 6*T(n,n,n)*T(n,n,n) - 12*T(m,m,n)*T(n,n,n);

%coeff_t4 = -6*sum_Tijm_Tijn - 12*sum_Tijm_Tijn - 12*sum_Tinn_Timn - 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*T(m,n,n)*T(n,n,n) - 6*T(n,n,n)*T(m,m,m) - 18*T(m,m,n)*T(m,n,n) - 36*T(m,n,n)*T(m,m,n) + 18*T(m,n,n)*T(n,n,n);
coeff_t4 = (-6) * sum_Tijm_Tijn + 24*sum_Tinn_Timn -36*sum_Timm_Timn -54*T(m,m,n)*T(m,n,n) + 30*T(m,n,n)*T(n,n,n) - 6*T(m,m,m)*T(n,n,n);

%coeff_t3 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm -24*sum_Timn_Timn + 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn - 6*T(m,m,m)*T(m,n,n) + 6*T(n,n,n)*T(m,m,n) + 18*T(m,m,m)*T(n,n,n) - 36*T(m,m,n)*T(m,m,n) + 36*T(m,n,n)*T(m,n,n) - 18*T(m,n,n)*T(m,m,m) + 18*T(m,n,n)*T(m,m,m);
coeff_t3 = 12*sum_Tijn_Tijn - 12*sum_Tijm_Tijm + 12*sum_Tinn_Tinn - 12*sum_Timm_Timm -24*T(m,n,n)*T(m,m,m) + 24*T(m,m,n)*T(n,n,n) - 36*T(m,m,n)*T(m,m,n) + 36*T(m,n,n)*T(m,n,n);

%coeff_t2 = 12*sum_Tijm_Tijn - 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*sum_Timm_Timn + 6*T(m,m,m)*T(n,n,n) - 12*T(m,m,m)*T(m,m,n) + 36*T(m,m,n)*T(m,n,n) - 18*T(m,n,n)*T(m,m,m) + 18*T(m,n,n)*T(m,m,n);
coeff_t2 = 6*sum_Tijm_Tijn - 24*sum_Timm_Timn + 36*sum_Timn_Tinn + 6*T(m,m,m)*T(n,n,n) -30*T(m,m,m)*T(m,m,n) + 54*T(m,m,n)*T(m,n,n);

%coeff_t1 = 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn + 12*T(m,m,m)*T(m,n,n) - 6*T(m,m,m)*T(m,m,m) + 18*T(m,m,n)*T(m,m,n);
coeff_t1 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm - 12*sum_Timm_Timm + 12*sum_Timm_Tinn + 24*sum_Timn_Timn -6*T(m,m,m)*T(m,m,m) + 12*T(m,m,m)*T(m,n,n) + 18*T(m,m,n)*T(m,m,n);

coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T(m,m,m)*T(m,m,n);

%disp(coeff_t6);
%disp(coeff_t5);
%disp(coeff_t4);
%disp(coeff_t3);
%disp(coeff_t2);
%disp(coeff_t1);
%disp(coeff_t0);

poly = [coeff_t6 coeff_t5 coeff_t4 coeff_t3 coeff_t2 coeff_t1 coeff_t0];
end


