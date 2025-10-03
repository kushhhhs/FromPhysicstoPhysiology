function TA=MakeTissueAreas(S)
% this function creates rectangular tissue areas in 2D, based on the
% settings of:
% S.Ytop S.Ybottom  vertical limits at the arterial side and at the venous side,  
% S.Xleft and S.Xright left and right limits of the whole brain
% S.N_Tissue_Areas the number of tissue areas next to each other
% this function assumes flow from top to bottom, such that the rectangles
% indicate parallel tissue areas

% % % awidth=(S.Xright-S.Xleft)/S.N_Tissue_Areas;
% % % for ta=1:S.N_Tissue_Areas
% % % 	TA(ta).Ytop=S.Ytop;
% % % 	TA(ta).Ybottom=S.Ybottom;
% % % 	TA(ta).Xleft=S.Xleft+(ta-1)*awidth;
% % % 	TA(ta).Xright=TA(ta).Xleft+awidth; 
% % % 	% any other derived geometric parameters here
% % % 	TA(ta).Area=(TA(ta).Ytop-TA(ta).Ybottom)*(TA(ta).Xright-TA(ta).Xleft);
% % % end
% % % end

% october 2022: have this more general als Vertices (corners) n-by-2 points
% [x1 y1, ..., xn yn]. Here, make these squares

awidth=(S.Xright-S.Xleft)/S.N_Tissue_Areas;
for ta=1:S.N_Tissue_Areas
	Ytop=S.Ytop;
	Ybottom=S.Ybottom;
	Xleft=S.Xleft+(ta-1)*awidth;
	Xright=Xleft+awidth;
	TA(ta).Vertices=[Xleft,Ybottom; Xleft,Ytop;Xright,Ytop;Xright,Ybottom];
	% any other derived geometric parameters here
	TA(ta).AreaSize=polyarea(TA(ta).Vertices(:,1),TA(ta).Vertices(:,2));
end
end