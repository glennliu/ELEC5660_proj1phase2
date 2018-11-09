function s_des = trajectory_generator(ti, path)


if nargin > 1 % pre-process can be done here (given waypoints)
    % set parameters
    t0 = 0;
    T = 2;   %keyframe duration
    kr = 4; %derivative order
    st_1 = [-1;2];   %[p0 pt]
    st_2 = [2;0];   %[v0 vt]
    st_3 = [0;0];   %[a0 at]

    % induced params
    num_st = size(st_1,1);  %number of constrain
    N = 2*kr -1;

    % build matrix
    Q = zeros(N+1,N+1);
    A = zeros(num_st,N+1);
    f = zeros(N+1,1);
        
    for i=kr+1:N+1
        for j = kr+1:N+1
            Q(i,j) = ((i-1)*(i-2)*(i-3)*(i-4)*(j-1)*(j-2)*(j-3)*(j-4)/(i+j-9))*T^(i+j-9);
        end
    end
    
%     k1 = 0; % 1st derivation
    for k=0:1
        for i=0:N
            if i-k <0
            A(k*2+1,i+1) = 0;
            A(k*2+2,i+1) = 0;
            else
            A(k*2+1,i+1) = (factorial(i)/factorial(i-k))*t0^(i-k);
            A(k*2+2,i+1) = (factorial(i)/factorial(i-k))*T^(i-k);
            end
        end
    end
    
    [x,fval,exitflag,output,lamda] = quadprog(Q,f,[],[],A,[st_1;st_2]);
    
    % create output path
    frame = 100;
    s_des = zeros(frame+1,1);
    s_des(1,1) = st_1(1,1);
    for index=2:frame+1
        t = t0 +(index-1)*T/frame;
        s_des(index,1) = 0;
        for i=0:N
            s_des(index,1) = x(i+1)*t^(i)+s_des(index,1);
        end
    end
        
    disp('break');
    
else % output desired trajectory here (given time)
    warning('input errors for trajectory generator!')
    
end

end


