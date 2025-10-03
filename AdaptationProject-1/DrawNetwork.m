function hfig = DrawNetwork(S,fignum)

% FP2P: this is very basic drawing of the topology, you could extend this
% to color-coded or line width-coded plots of radius, flow, WSS, ...(if you
% do the adapt assigment)
% you could mark the nodes that have a clot (if you do the clot assignment)
% you could also make multiple subplots in a single plot
% or make a movie of your plots, showing the dynamics


hfig=figure(fignum);clf; hold on 

for i=1:S.nin
	style='ob'; 
    plot(S.IN(i).x,S.IN(i).y,style)
	
end
for j=1:S.nie
	 xx=[S.IN(S.IE(j).nodes).x];yy=[S.IN(S.IE(j).nodes).y]; %% xx are x values of left and right node
	 plot(xx,yy,'-b')
	 
end
v=find([S.IN.nsources]); % these indices in IN are connected to SE
for k=1:length(v)
	 plot(S.IN(v(k)).x,S.IN(v(k)).y,'or')
	 
end


%% get the labeling of the nodes and elements
% remove this or adapt if you want
for i=1:S.nin
		S.IN(i).x=S.IN(i).pos(1); % couldn't get the indexing right below, so the hard way here
		S.IN(i).y=S.IN(i).pos(2);
		text(S.IN(i).x,S.IN(i).y,num2str(i),'Color','black','FontSize',14)
end
for j=1:S.nie
		 xx=[S.IN(S.IE(j).nodes).x]; % pos(1) instead of x doesn't work
		 yy=[S.IN(S.IE(j).nodes).y]; %% xx are x values of left and right node
		 xm=(xx(1)+xx(2))/2; ym=(yy(1)+yy(2))/2;
		 text(xm,ym,num2str(j),'Color','blue','FontSize',12)
end

%% Draw the tissue areas if these exist, this is for Transport
if isfield(S,'TA')
    for ta=1:S.N_Tissue_Areas
        hlarea=line(S.TA(ta).Vertices([1:end,1],1),S.TA(ta).Vertices([1:end,1],2));
        set(hlarea,'color',[0.5 0.5 0],'LineStyle','-.');

    end
end


daspect([ 1 1 1]) % data aspect ratio same scale in all directions
hold off

