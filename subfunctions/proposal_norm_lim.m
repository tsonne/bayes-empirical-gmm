function state = proposal_norm_lim(x0,sd,N,lb,ub)
state = x0 + sd.*randn(1,N);
while sum(state>=lb)+sum(state<=ub)~=2*N
    state = x0 + sd.*randn(1,N);
end