function s_des = trajectory_generator(ti, path)


if nargin > 1 % pre-process can be done here (given waypoints)
    % set parameters
    t0 = 0;
    T = 2;   % segment duration
    kr = 4; % derivative order
    frame = 100;    % number of frames
    
    % vias constrain
    st_1 = [-1;0;0;1];   % [x0 xt vx0 vxt ax0 axt]
    st_2 = [0;1;1;0];   % [y0 yt vy0 vyt ay0 ayt]
    st_3 = [1;0;0;2];   % [z0 zt vz0 vzt az0 azt]
    dim = 3;

    % induced params
%     dim = size(st_1,1)/2; % dimensions of path
    num_st = size(st_1,1);  %number of constraint
    N = 2*kr -1;

    % build matrix
    Q = zeros(N+1,N+1);
    A = zeros(num_st,N+1);
    f = zeros((N+1)*dim,1);
    s_des = zeros(frame+1,dim);

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
    
    % duplicate in dimension            
    A0 = zeros(num_st,N+1);
    Q0 = zeros(N+1,N+1);
    switch dim
        case 1
            Qtotal = Q;
            Atotal = A;
            s_des(1,1) = st_1(1,1);
        case 2
            Qtotal = [Q Q0
                      Q0 Q];
            Atotal = [A A0
                      A0 A];
            s_des(1,:) = [st_1(1,1),st_2(1,1)];
        case 3
            Qtotal = [Q Q0 Q0
                      Q0 Q Q0
                      Q0 Q0 Q];
            Atotal = [A A0 A0
                      A0 A A0
                      A0 A0 A];
            s_des(1,:) = [st_1(1,1),st_2(1,1),st_3(1,1)];
    end
    
    [x,fval,exitflag,output,lamda] = quadprog(Qtotal,f,[],[],Atotal,[st_1;st_2;st_3]);
    
    % create output path
    disp('debug');
    
    for index=2:frame+1
        t = t0 +(index-1)*T/frame;
        for i=0:N   % x path
            s_des(index,1) = x(i+1)*t^(i) + s_des(index,1);
            s_des(index,2) = x(i+N+2)*t^(i) + s_des(index,2);
            s_des(index,3) = x(i+2*N+3)*t^(i) + s_des(index,3);

        end
    end
        
    disp('break');
    
else % output desired trajectory here (given time)
    warning('input errors for trajectory generator!')
    
end

end


