function [U] = MODSUP(X,m,k,Fw,omega,zeta)

%{
    > Assuming mass-normalised vectors, X'*M*X = I and X'*K*X = diag(lam0)
    > Assuming sorted omega values

    > omega is a vector of required solution frequencies [s^-1]
    > Fw is an input force in the same coordinates as X
    > m is the generalised mass (unit vector)
    > k is the generalised stiffness (k = diag(lam0), K*X = M*X*lam0)
    > zeta is an optional modal damping percentage constant/vector

    >> Velocity = U*im.*omega; Acceleration = -U.*omega^2
%}

% Initial variable definitions
    omega   = cvec(omega)';
	im      = complex(0,1);
    e       = exp(1.0);

    if ~exist('zeta','var')     b = zeros(size(m));
    else                        b = 2*(m.*omega).*zeta;
    end

% Modal force matrix
	f = X'*Fw;
    n = size(f,1);
    
% Truncated modal mass, stiffness and damping 
    if size(m,1) ~= n
        m = m(1:n);
        k = k(1:n);
        b = b(1:n,:);
    end
    
% Make expected f input match omega on calling MODSUP
    switch size(f,2)
        case 1
            f = repmat(f,1,length(omega));
        case numel(omega)
            % Do nothing. 
        otherwise
            error('size(Fw,2) must equal 1 or numel(omega)')                % Interpolation of input load before calling MODSUP can be used
    end

%{
% Modal solution, damped and undamped, naiive implementation
    for i = 1:length(omega)
        KSIw(:,i) = f(:,i)./(k - m*omega(i)^2);
        KSIc(:,i) = f(:,i)./(k + im*b(:,i).*omega(i) - m*omega(i)^2);       % Assuming b(:,i) is either scalar or vector of size n
    end
%}

% Modal solution, damped and undamped, KSI(mode,omega)
%    KSIw1 = f./(k - m.*omega.^2);                                           % Vectorised. Note: usually m == 1   
    KSIc1 = f./(k + im*b.*omega - m.*omega.^2);                             % For % crit modal c = 2*m*omega*zeta = ccrit*zeta = ccrit*G/2
   
% Modal superposition
    U = X*KSIc1;

end




