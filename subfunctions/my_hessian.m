function H = my_hessian(func,loc,h,varargin)
% ts 2015-09-21: 
% calculate Hessian matrix by central finite differences,
% any given function handle "func", which accepts parameters,
% any number of parameters (based on input "loc" and "h"),
% INPUT:
% func = function handle for function to evaluate
% loc = vector for parameter location where Hessian is evaluated
% h = step width vector (one element per parameter)
% varargin = additional data, if required by function
%
% input function must be of form: func(loc,varargin{:})
M = length(loc);
H = nan(M);
e = [1 1];
cen = func(loc,varargin{:}); % n=1
for a = 1:M % 2nd order diff. w.r.t. same parameter
    H(a,a) = (func(loc + sparse(1,a,1,1,M).*h,varargin{:}) -2*cen ...
        + func(loc - sparse(1,a,1,1,M).*h,varargin{:}))/(h(a)^2);
end % n=2*M
if M<2
    return
end
for a = 1:M-1 % diff. w.r.t. parameter a
    for b = a+1:M % diff. w.r.t. parameter b
        H(a,b) = (func(loc + sparse(e,[a b],e,1,M).*h,varargin{:})...
            -func(loc + sparse(e,[a b],[1 -1],1,M).*h,varargin{:})...
            -func(loc + sparse(e,[a b],[-1 1],1,M).*h,varargin{:})...
            +func(loc + sparse(e,[a b],-[1 1],1,M).*h,varargin{:})...
            )/(4*h(a)*h(b));
        H(b,a) = H(a,b);
    end
end % n=2*(M-1)*M
% n = 1 + 2*M^2
% M=parameters, n=calculations of func
% M   n
% 1   3
% 2   9
% 3  19
% 4  33
