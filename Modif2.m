function fun=Modif2(fun,IndBckgrnd)
%    global a2 x BH1 esg
%    myfunc= @(a1,a2) a1*x.^3+a2*x
%    plot(x,myfunc(a1,a2))
% set background to NaN 
    fun(IndBckgrnd)= NaN;
%normalize to [0,1]
    fun = (fun -min(fun(:)))./(max(fun(:))-min(fun(:)));
    imagesc(fun)
    set(gca,'xtick',[])
    set(gca,'ytick',[])

end