function fun=Modif5(fun,IndBckgrnd,ImPhysio,alpha1,alpha2)
% set background to NaN
fun(IndBckgrnd)= NaN;
%normalize to [0,1]
fun = (fun -min(fun(:)))./(max(fun(:))-min(fun(:)));
fun = uint8(fun*255);
%Im=(fun*alpha1+ImPhysio*alpha2*(1-alpha1))/(alpha1+alpha2*(1-alpha1));
%imagesc(Im)
%imagesc(cat(3,zeros(size(fun)),fun*alpha1,zeros(size(fun))));

imagesc(cat(3,ImPhysio*alpha2,fun*alpha1,zeros(size(fun))));
%imagesc(fun*alpha1);

hold on
%iim2=imagesc(cat(3,ImPhysio,zeros(size(fun)),zeros(size(fun))));
%set(iim2,'AlphaData',alpha2);
set(gca,'xtick',[])
set(gca,'ytick',[])
hold off

end