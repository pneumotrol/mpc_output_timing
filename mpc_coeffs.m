function coeff = mpc_coeffs(sysd,params)
    A = sysd.A;
    B = sysd.B;
    C = sysd.C;
    N = params.N;
    Q = params.Q;
    R = params.R;

    % dimensions
    l = size(C,1); % output
    m = size(B,2); % input
    n = size(A,2); % state

    % free output
    F_tmp = zeros(l*N,n);
    A_tmp = A;
    for i = 1:N
        F_tmp((1:l)+l*(i-1),:) = C*A_tmp;
        A_tmp = A_tmp*A;
    end
    coeff.F = F_tmp;

    % forced output
    G_tmp = zeros(l*N,m*N)';
    for i = 1:N
        A_tmp = eye(n);
        for j = i:N
            G_tmp((1:m)+m*(i-1),(1:l)+l*(j-1)) = (C*A_tmp*B)';
            A_tmp = A_tmp*A;
        end
    end
    coeff.G = G_tmp';

    % output weight (block diagonal)
    Q_tmp = zeros(l*N);
    for i = 1:N
        Q_tmp((1:l)+l*(i-1),(1:l)+l*(i-1)) = Q;
    end
    coeff.Q = Q_tmp;

    % input weight (block diagonal)
    R_tmp = zeros(m*N);
    for i = 1:N
        R_tmp((1:m)+m*(i-1),(1:m)+m*(i-1)) = R;
    end
    coeff.R = R_tmp;
end
