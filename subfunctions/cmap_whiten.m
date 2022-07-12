function map = cmap_whiten(cmap,wpower)
% fade input colormap to white (such that cmap(1,:)-->white).
% INPUT:
%        wpower = exponent of whitening percentage, (DEF: 1)
%                    if 0 < wpower < 1 --> more whitening
%                    if 1 < wpower     --> less whitening
%
% 2016-05-18 tsonne
assert(nargin>0,'need at least colormap input')
if nargin<2, wpower=1; end
assert(wpower>=0,'wpower must be >= 0.');

SZ = size(cmap);
w  = ones(SZ);
g  = linspace(1,0,SZ(1))';
g  = repmat(g,1,3);
g  = g.^wpower; % nonlinear scaling
map = (1-g).*cmap + g.*w;