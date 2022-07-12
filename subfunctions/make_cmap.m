function map = make_cmap(Cols,ColPos,N)
% make colormap according to given input,
% INPUT:
%    Cols   = (n,3)-matrix of discrete colortriplets
%             they will be anchorpoints, interpolate between them
%    ColPos = relative position vector, (lenght=n)
%             determines order and space between given colors
%    N      = (optional) number of colors to produce (DEF: 256)
% 2016-05-18 tsonne
if nargin<3, N=256; end
assert(nargin>1,'need Cols and ColPos!')
nr = size(Cols,1);
assert(nr==numel(ColPos),'Given colors and position vector mismatch!')

cp = (ColPos - min(ColPos))/(max(ColPos)-min(ColPos));
cp(cp<0)=0;
cp(cp>1)=1;
map = interp1(cp(:),Cols,linspace(0,1,N));
map(map<0)=0;
map(map>1)=1;