function s_des = trajectory_generator(ti, path)


if nargin > 1 % pre-process can be done here (given waypoints)
    %% set parameters here
    t0 = 0;
    T = 2;   % segment duration
    kr = 4; % derivative order
    frame = 100;    % number of frames
    
    % derivative constraint
    % Notice: set acce st without vel st doesn't perform well
    st_x = [-1;0;0;0;0;0];   % [x0 xt vx0 vxt ax0 axt]
    st_y = [0;1;0;0;0;0];   % [y0 yt vy0 vyt ay0 ayt]
    st_z = [1;0;0;0;0;0];   % [z0 zt vz0 vzt az0 azt]
    
    % c_select = [p0 pt v0 vt a0 at]
    % 1 for choose, 0 for none
    % the same selection matrix for x,y,z dim
    c_select = [1 1 1 1 1 1]
    
    dim = 3;
%     st_order = size(st_x,1)/2;   % 2nd derivative st

    %% calculation
    % induced params
    num_st = size(st_x,1);  %number of constraint
    N = 2*kr -1;

    % build matrix
    Q = zeros(N+1,N+1);
    Ad = zeros(num_st,N+1);
    Ac = zeros((N+1)*num_st,1);
    f = zeros((N+1)*dim,1);
    s_des = zeros(frame+1,dim);

    for i=kr+1:N+1
        for j = kr+1:N+1
            Q(i,j) = ((i-1)*(i-2)*(i-3)*(i-4)*(j-1)*(j-2)*(j-3)*(j-4)/(i+j-9))*T^(i+j-9);
        end
    end

%   derivative constraint matrix
    for k=0:2
        for i=0:N
            if i-k <0
            Ad(k*2+1,i+1) = 0;
            Ad(k*2+2,i+1) = 0;
            else
            Ad(k*2+1,i+1) = (factorial(i)/factorial(i-k))*t0^(i-k);
            Ad(k*2+2,i+1) = (factorial(i)/factorial(i-k))*T^(i-k);
            end
        end
        
        if c_select(1,k*2+1) == 0
            Ad(k*2+1,:) = zeros(1,N+1);
        end
        
        if c_select(1,k*2+2) == 0
            Ad(k*2+2,:) = zeros(1,N+1);
        end
       
    end
    
%   continuity constraint

    
%   duplicate in dimension            
    A0 = zeros(num_st,N+1);
    Q0 = zeros(N+1,N+1);
    switch dim
        case 1
            Qtotal = Q;
            Atotal = Ad;
            s_des(1,1) = st_x(1,1);
        case 2
            Qtotal = [Q Q0
                      Q0 Q];
            Atotal = [Ad A0
                      A0 Ad];
            s_des(1,:) = [st_x(1,1),st_y(1,1)];
        case 3
            Qtotal = [Q Q0 Q0
                      Q0 Q Q0
                      Q0 Q0 Q];
            Atotal = [Ad A0 A0
                      A0 Ad A0
                      A0 A0 Ad];
            s_des(1,:) = [st_x(1,1),st_y(1,1),st_z(1,1)];
    end
    
    [x,fval,exitflag,output,lamda] = quadprog(Qtotal,f,[],[],Atotal,[st_x;st_y;st_z]);
    
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


