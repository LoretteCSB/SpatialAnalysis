function fun=Modif6(fun_Promo,fun_Inhib,IndBckgrnd,ImPhysio,alpha1,alpha2,alpha3)
fun = fun_Promo+fun_Inhib;
fun(IndBckgrnd)= NaN;% set background to NaN

delta_Promo=max(fun_Promo(:)) -min(fun_Promo(:));
delta_Inhib=max(fun_Inhib(:)) -min(fun_Inhib(:));


%{
if (delta_Promo )
    if ( max(-fun_Inhib(:))-min(-fun_Inhib(:)) )
        imagesc(cat(3,ImPhysio*alpha1,...
            ( fun_Promo-min(fun_Promo(:) ))./delta_Promo*255*alpha2,...
            (-fun_Inhib-min(-fun_Inhib(:)))./delta_Inhib*255*alpha3));
        %normalize to [0,1]
        fun = (fun -min(fun(:)))./(max(fun(:))-min(fun(:)));
    else
        imagesc(cat(3,ImPhysio*alpha1,...
            ( fun_Promo-min(fun_Promo(:) ))./delta_Promo*255*alpha2,...
            zeros(size(fun_Promo))));
        %normalize to [0,1]
        fun = (fun -min(fun(:)))./(max(fun(:))-min(fun(:)));
    end
else
    if ( max(-fun_Inhib(:))-min(-fun_Inhib(:)) )
        imagesc(cat(3,ImPhysio*alpha1,...
             zeros(size(fun_Promo)),...
            (-fun_Inhib-min(-fun_Inhib(:)))./delta_Inhib*255*alpha3));
        %normalize to [0,1]
        fun = (fun -min(fun(:)))./(max(fun(:))-min(fun(:)));
    else
        imagesc(cat(3,ImPhysio*alpha1,...
            zeros(size(fun_Promo)),...
            zeros(size(fun_Promo))));
        
    end
    
end
%}
fun(fun<0)=0;
fun = 255*(fun -min(fun(:)))./(max(fun(:))-min(fun(:)));

%imshow(cat(3,uint8(ImPhysio*alpha1),uint8(fun*alpha1,zeros(size(fun_Promo)))))
imshow(cat(3,uint8(ImPhysio*alpha1),uint8(fun*alpha1),uint8(zeros(size(fun_Promo)))))

colorbar

hold on
set(gca,'xtick',[])
set(gca,'ytick',[])
hold off


fun = uint8(fun*255);


end