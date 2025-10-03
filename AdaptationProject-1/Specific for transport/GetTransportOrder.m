function [S,TransportOrder]=GetTransportOrder(S)
% this function sorts out the order of internal elements for calculation of contrast transport
% basic principle: we can calculate exit function from entrance function if
% all exit functions of upstream ('mother') connected vessel are known.
% S.IE gets an extra field 'mothers', a vector of all other internal
% elements that are immediately upstream. The inflow concentration is the
% flow-weighted average of the outflow concentrations of all mothers
% this function assumes that C(t) is known for all INFLOW source elements.
% this function is based on input and output C(t) defined for
% elements and not nodes. 

TransportOrder=zeros(1,S.nie); % predefine this vector
TransportOrderDone=zeros(1,S.nie); % indices into S.IE whose order in transport calculations has been determined
iTO=0; % index into TransportOrder
ie=0;
while iTO<S.nie
	if ie<S.nie
		ie=ie+1;
		if ~TransportOrderDone(ie) 
			% figure out which is the upstream node
			j=1+(S.IE(ie).Q<0); % j=1 for forward flow, 2 for backflow
			%figure out which elements are connected to the upstream node
			in=S.IE(ie).nodes(j); % this is the node we are dealing with
			coniev=S.IN(in).ie; % all elements connected to this node 
			coniev=coniev(find(coniev~=ie)); % remove the current element from this vector
			% coniev has all the connected ie to the relevant node, but not
			% all these connected ie bring flow to the node, we need to figure
			% out which ie are really upstream ('mother'), rather than being
			% 'sister'. This is the case if the 'remote' node P is higher than
			% the P in our in . 
			mothers=[]; % we will store the mothers for later easier use
			TransportReady=1; % switch to zero if one of the mothers has not been done
			for k=1:length(coniev)
				nodes=S.IE(coniev(k)).nodes;	
				j=find(nodes~=in);% one or two, index to the remote node
				ismother=S.IN(nodes(j)).P>S.IN(in).P;
				if ismother
					mothers=[mothers,coniev(k)];
					if ~TransportOrderDone(coniev(k)) 
						TransportReady=0; break
					end
				end
			end
			if TransportReady
				S.IE(ie).mothers=mothers; % store the indices into IE of mother (upstream) vessels
				TransportOrderDone(ie)=1; 
				iTO=iTO+1;
				TransportOrder(iTO)=ie; % add this element to the calculation order for transport
			end
		end
	else
		if iTO<S.nie
			ie=1; % redo the sequence
		end
	end
		
end
