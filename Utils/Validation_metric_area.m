function AREA_diff=Validation_metric_area(y1,y2,disp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation metric (AREA)
%% INPUTS
% y1 : simulated data (n-by-1 vector, where n : # simulated data)
% y2  : experimental data (m-by-1 vector, where m : # experimental data)
%% OUTPUTS
% AREA : simulated data
% 2016.05.07: Coded by SS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical cumulative distribution function (ECDF)
[f1,x1] = ecdf(y1); [f2,x2] = ecdf(y2);
[x1 f1]=stairs(x1,f1); [x2 f2]=stairs(x2,f2);

%% PRE-PROCESSING FOR COMPUTING AREA METRICS (LIKE INTERPOLATION)
x1=x1(2:end); x2=x2(2:end);
f1=f1(2:end); f2=f2(2:end);

X=[x1 f1 ones(size(x1)); x2 f2 2*ones(size(x2))];
[Y,I] = sort(X(:,1));
X=X(I,:);

A=[]; B=[];
for i=1:2:size(X,1)
    IND1=find(x1==X(i,1) & X(i,3)==1);
    IND2=find(x2==X(i,1) & X(i,3)==2);
    
    if ~isempty(IND1), A(i:i+1,1)=f1(IND1);
    else
        if i==1, A(i:i+1,1)=0;
        else A(i:i+1,1)=A(i-1,1); end
    end
    
    if ~isempty(IND2), B(i:i+1,1)=f2(IND2);
    else
        if i==1, B(i:i+1,1)=0;
        else B(i:i+1,1)=B(i-1,1); end
    end
end

%% COMPUTE THE INDICE ACCORDING TO EACH SECTIONS
Index_n_area=A<=B;
Index_n_area(1)=Index_n_area(2); Index_n_area(end)=Index_n_area(end-1);

INDEX=[]; n_area=1; n_sign=0;
for i=1:size(A,1)-1
    INDEX(i,1)=n_area; INDEX(i,2)=n_sign;
    
    if ~(Index_n_area(i)==Index_n_area(i+1))
        n_area=n_area+1;
        if n_sign==0 % A>B
            n_sign=1;
        else % A<B
            n_sign=0;
        end
        
    end
    
end

%% COMPUTE AREA
AREA=[];
for i=1:n_area
    IND=find(INDEX(:,1)==i);
    if size(IND,1)==1, Area_A=0; Area_B=0;
    else Area_A=trapz(X(IND,1),A(IND)); Area_B=trapz(X(IND,1),B(IND)); end
    AREA(i)=Area_A-Area_B;
    
end
AREA_diff = sum(abs(AREA));

%% LOCALIZE THE ARES SECTION
if disp
%     figure;
%     plot(x1,f1,'b-','linewidth',2); hold on
%     plot(x2,f2,'r:','linewidth',2);
    
    for i=1:n_area
        IND=find(INDEX(:,1)==i);
        [a b]=stairs(X(IND,1),B(IND));
        [c d]=stairs(X(IND,1),A(IND));
        
        x=[a' fliplr(a')];yy=[b' fliplr(d')];
        fill_obj=fill(x,yy,'b');
        set(fill_obj,'facecolor',[0.867 0.949 1],'edgecolor','none','facealpha',0.3);
    end
    
%     set(gca,'fontsize',20,'fontweight','bold');
end
end