function imopt = opt_process(im0,imw,w)

CSparam.Itnlim = 50;	% default number of iterations
CSparam.gradToll = 1e-30;	% step size tollerance stopping criterea (not used)

CSparam.l1Smooth = 1e-15;	% smoothing parameter of L1 norm
CSparam.pNorm = 1;  % type of norm to use (i.e. L1 L2 etc)

% line search parameters
CSparam.lineSearchItnlim = 100;
CSparam.lineSearchAlpha = 0.01;
CSparam.lineSearchBeta = 0.6;
CSparam.lineSearchT0 = 1;

CSparam.TV=TVOP;
CSparam.data = im0;
CSparam.TVWeight= w;

% CSparam.BW = power(edge(:,:,slice),1);
% gaussianFilter = fspecial('gaussian', [5, 5], 20);
% CSparam.BW = imfilter(abs(CSparam.BW), gaussianFilter, 'circular', 'conv');
% CSparam.BW = CSparam.BW/max(CSparam.BW(:));
% CSparam.BW = 0.9*power(edge(:,:,slice),14);
CSparam.BW = imw;
x0 = im0; 
for n=1:1:5
    x1=CGx(x0,CSparam);
    x0=x1;
end

imopt = x1;

end